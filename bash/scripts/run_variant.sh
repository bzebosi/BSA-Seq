#!/bin/bash

# Load configuration
source "$(dirname "$0")/variant_config.sh"

# Fail a pipeline if any command in it fails
set -o pipefail

# -------------------------------------------------------------------------------------------------------
# log functions for updates and status
# -------------------------------------------------------------------------------------------------------
logmsg() { echo "$(date '+%Y-%m-%d %H:%M:%S'): $*" ; }


# -------------------------------------------------------------------------------------------------------
# Validate project and its directory
# -------------------------------------------------------------------------------------------------------
get_project_dir(){
    local project=$1
    local dir=${sample_loc[$project]:-}

    # check if the project was provided
    if [[ -z ${project} ]]; then
        logmsg "Error : No project specified."
        return 1
    fi 
    
    # check if the project directory exits
    if [[ -d ${dir} ]]; then
        echo ${dir} 
    else
        logmsg "Error: Project directory does not exit ${dir}"
        return 1
    fi
}

# -------------------------------------------------------------------------------------------------------
# Create one or more directories
# -------------------------------------------------------------------------------------------------------
create_dir() {
    local dir

    # Loop through each directory provided
    for dir in "$@"; do
        
        # check if the directory already exists
        if [[ -d ${dir} ]]; then
            logmsg "Directory already exists: ${dir}"

            # create the directory
            elif mkdir -p ${dir}; then
                logmsg "Directory created: ${dir}"
            
            # Report an error if the directory failed to create
            else
                logmsg "Error: Failed to create directory: ${dir}"
                return 1
        fi
    done
}

# -------------------------------------------------------------------------------------------------------
# download genome functions
# . looks up URL in genome_urls[]
# -------------------------------------------------------------------------------------------------------

download_genome(){
    local gbase=$1
    local url=${genome_urls[$gbase]:-}
    local genome=${ref_dir}/${gbase}.fa.gz

    # Create reference directory
    create_dir ${ref_dir}

    # Check if a URL was found for the genome
    if [[ -z $url ]]; then
        logmsg "Error: No URL found for genome: $gbase"
        return 1
    fi

    logmsg "Checking genome $gbase : $genome"

    # Skip download if genome already exists and is not empty
    if [[ -s ${genome} ]]; then
        logmsg "Genome ${genome} already exits - skipping download."
    else
        # Download genome
        logmsg "Downloading ${gbase} genome.."
        if wget -O ${genome} ${url}; then
            logmsg "${gbase} sucessfully downloaded"
        else 
            logmsg "Error: download failed or URL missing for ${gbase}"
            return 1
        fi
    fi
    
}

# -------------------------------------------------------------------------------------------------------
# index genome
# -------------------------------------------------------------------------------------------------------

index_genome(){
    local gbase=$1
    local genome_gz=${ref_dir}/${gbase}.fa.gz
    local genome=${ref_dir}/${gbase}.fa
    local idx_dir=${ref_dir}/indexes
    local mmi=${idx_dir}/${gbase}.mmi
    local fai=${idx_dir}/${gbase}.fa.fai
    local link_fa=${idx_dir}/${gbase}.fa

    create_dir ${idx_dir}

    # check genome 
    if [[ ! -s ${genome} && ! -s ${genome_gz} ]]; then 
        logmsg "Error: ${genome_gz} not found or empty - cannot index"
        return 1
    fi

    # unzip genome if needed
    if [[ ! -s ${genome} &&  -s ${genome_gz} ]]; then 
        logmsg "Unzipping ${genome_gz}"
        if gunzip -c ${genome_gz} > ${genome}; then
            logmsg "Unzipped genome created: ${genome}"
        else
            logmsg "Error: Failed to unzip ${genome_gz}" && exit 1
        fi
    fi

    # index the genome using minimap2
    logmsg "Starting Minimap2 indexing for ${genome}.."
    if [[ -s ${mmi} ]]; then
        logmsg "Minimap2 index already exists for ${gbase} - skipping indexing"
    else
        logmsg "Building minimap2 index for ${gbase}..."
        if minimap2 -t ${threads} -d ${mmi} ${genome} ; then
            logmsg "Minimap2 index completed ${mmi}"
        else
            logmsg "Error: minimap2 index failed for ${genome}" && exit 1
        fi
    fi

    # index the genome using samtools faidx .fai and .gzi
    if [[ -s ${fai} ]]; then
        logmsg "samtools index .fai for ${gbase} already exists - skip indexing"
    else
        if [[ -s ${genome}.fai ]]; then
            mv -u ${genome}.fai ${fai} && logmsg "moved ${genome}.fai to ${fai}"
        else
            logmsg "Running samtools faidx on ${genome}.."
            if samtools faidx ${genome} ; then
                logmsg "samtools faidx complete for ${gbase}"
                if [[ -f ${genome}.fai ]]; then
                    mv -u ${genome}.fai ${fai} && logmsg "moved ${genome}.fai to ${fai}"
                fi
            else
                logmsg "samtools index failed for ${gbase}." && exit 1 
            fi
        fi
    fi

    # create fasta symlink besides the indexs
    if [[ -e ${link_fa} ]]; then
        logmsg "symlink already exists"
    else
        if ln -s ${genome} ${link_fa}; then
            logmsg "Symlink created ${link_fa} for ${genome}"
        else
            logmsg "Error: failed to created symlink for genome: ${gbase}." && exit 1 
        fi
    fi
}


# -------------------------------------------------------------------------------------------------------
# trim read
# -------------------------------------------------------------------------------------------------------

trim_reads(){
    local project=$1
    local project_dir=$(get_project_dir ${project}) || return 1
    local reads_dir=${project_dir}/reads
    local trim_dir=${project_dir}/trim_reads
    local stats_dir=${project_dir}/stats
    local reports_dir=${stats_dir}/reports
    local summary=${reports_dir}/fastp_trim_summary.tsv

    create_dir ${stats_dir} ${trim_dir} ${reports_dir} || return 1

    # Initialize summary file 
    if [[ ! -s ${summary} ]]; then
        # Add a header to the summary file
        echo -e "Sample\tRead_Type\tTotal_Reads_Before\tGC_Content_Before\tTotal_Reads_After\t\
        GC_Content_After\tPassed_Filter_Reads\tLow_Quality_Reads\tToo_Many_N_Reads\t\
        Too_Short_Reads\tDuplication_Rate\tInsert_Size_Peak" > ${summary}

        logmsg "Summary file created: ${summary}"
    fi

    local R1
    for R1 in ${reads_dir}/*_R1.fq.gz; do
        local fname=$(basename ${R1})

        local sbase=${fname%_R1.fq.gz}
        local R2=${reads_dir}/${sbase}_R2.fq.gz
        local O1=${trim_dir}/${sbase}_trim_R1.fq.gz
        local O2=${trim_dir}/${sbase}_trim_R2.fq.gz
        local ht=${reports_dir}/${sbase}_fastp_report.html
        local jt=${reports_dir}/${sbase}_fastp_report.json
        local metrics=""
        local read_type

        # Check matching R2
        if [[ -s "${R2}" ]]; then
            read_type="PE"
        else
            read_type="SE"
        fi

        # Skip if already summarized
        if awk -F '\t' -v sample="${sbase}" -v type="${read_type}" '
            NR > 1 && $1 == sample && $2 == type && NF == 12 {found = 1}
            END {exit !found}
        ' "${summary}"; then
            logmsg "${sbase} already completed. Skipping."
            continue
        fi

        # Run fastp

        logmsg "Starting fastp for ${sbase} (${read_type})"
        if [[ "${read_type}" == "PE" ]]; then
            if [[ ! -s "${O1}" || ! -s "${O2}" || ! -s "${jt}" ]]; then
                if fastp -i ${R1} -I ${R2} -o ${O1} -O ${O2} --detect_adapter_for_pe \
                    --thread ${threads} -h ${ht} -j ${jt}; then
                    logmsg "fastp complete for ${sbase}"
                else
                    logmsg "ERROR : fastp failed ${sbase}"
                    continue
                fi
            else
                logmsg "$(basename "${O1}") and $(basename "${O2}") already exist and trimmed. skip fastp..."
            fi

        else
            if [[ ! -s "${O1}" || ! -s "${jt}" ]]; then
                if fastp -i ${R1} -o ${O1} --thread ${threads} -h ${ht} -j ${jt} ; then
                    logmsg "fastp complete for ${sbase}"
                else
                    logmsg "ERROR : fastp failed ${sbase}"
                    continue
                fi
            else
                logmsg "$(basename "${O1}") already exists and is trimmed. Skipping fastp..."
            fi
        fi

        # Extract fastp metrics
        if [[ -s "${jt}" ]]; then
            metrics="$(
                jq -r '
                    [.summary.before_filtering.total_reads, .summary.before_filtering.gc_content,
                    .summary.after_filtering.total_reads, .summary.after_filtering.gc_content,
                    .filtering_result.passed_filter_reads, .filtering_result.low_quality_reads,
                    .filtering_result.too_many_N_reads, .filtering_result.too_short_reads,
                    .duplication.rate, .insert_size.peak ] | @tsv' "${jt}"
                )" || metrics=""
        fi

        # Check metrics
        if [[ -z "${metrics}" ]]; then
            logmsg "ERROR: Could not extract fastp metrics for ${sbase}"
            continue
        fi

        # Add summary row
        echo -e "${sbase}\t${read_type}\t${metrics}" >> "${summary}"

        logmsg "fastp complete for ${sbase} (${read_type})"
    done

    logmsg "Read trimming completed for project: ${project}"
}

# -------------------------------------------------------------------------------------------------------
# Map Reads 
# -------------------------------------------------------------------------------------------------------

map_reads(){
    local project=$1
    local gbase=$2
    local project_dir=$(get_project_dir "${project}") || return 1
    local trim_dir=${project_dir}/trim_reads
    local idx_dir=${ref_dir}/indexes
    local bam_dir=${project_dir}/bam

    local stats_dir=${project_dir}/stats
    local reports_dir=${stats_dir}/reports
    local plots_dir=${stats_dir}/plots

    local snps_dir=${project_dir}/snps_variants
    local snps_vcf=${snps_dir}/snps_vcf
    local snps_tsv=${snps_dir}/snps_tsv

    local svs_dir=${project_dir}/svs_variants
    local svs_vcf=${svs_dir}/svs_vcf
    local svs_tsv=${svs_dir}/svs_tsv


    create_dir "${bam_dir}" "${stats_dir}" "${reports_dir}" "${plots_dir}" || return 1
    create_dir "${snps_dir}" "${snps_vcf}" "${snps_tsv}" "${svs_dir}" "${svs_vcf}" "${svs_tsv}" || return 1
    # create reference paths
    local idx_mmi=${idx_dir}/${gbase}.mmi
    local genome_fa=${idx_dir}/${gbase}.fa

    if [[ -s "${idx_mmi}" && -s "${genome_fa}" ]]; then
        logmsg "Using index: ${idx_mmi} and FASTA: ${genome_fa}"
    else
        logmsg "Error: index or FASTA missing for ${gbase}." && exit 1 
    fi

    # Check if overall coverage file exists, create if missing
    local overall_coverage="${reports_dir}/overall_coverage.tsv"

    # Initialize coverage summary
    if [[ ! -s "${overall_coverage}" ]]; then
        if echo -e "Genome\tSample\tRead_Type\tCoverage\tTotal_Reads\tMapped_Reads\tProperly_Paired\tMapped_%" > "${overall_coverage}"; then

            logmsg "Successfully created ${overall_coverage}"
        else
            logmsg "Error creating ${overall_coverage}." && exit 1
        fi
    fi

    local O1

    for O1 in "${trim_dir}"/*_trim_R1.fq.gz; do
        local sbase=$(basename "${O1}" "_trim_R1.fq.gz")
        local O2=${trim_dir}/${sbase}_trim_R2.fq.gz
        local tag=${gbase}_${sbase}
        local bam_out=${bam_dir}/${tag}.bam
        local read_type
        local vcf_out=${snps_vcf}/${tag}.vcf.gz

        # check if paired or single end
        if [[ -s "${O2}" ]]; then
            read_type="PE"
        else
            read_type="SE"
        fi

        # Run Minimap2 alignment

        # Skip if BAM and BAM index already exists
        if [[ -s "${bam_out}" ]]; then
            logmsg "${sbase} already mapped to ${gbase}. Skipping mapping"

            if [[ ! -s "${bam_out}.bai" ]]; then
                logmsg "BAM index missing. Indexing ${bam_out}"

                if ! samtools index "${bam_out}"; then
                    logmsg "Indexing of ${bam_out} failed" && exit 1
                fi
            fi

        else
            # Map paired-end reads
            if [[ "${read_type}" == "PE" ]]; then
                logmsg "PE mapping ${sbase} to ${gbase} started"

                if minimap2 -ax sr -t "${threads}" "${idx_mmi}" "${O1}" "${O2}" | samtools sort -@ "${threads}" -o "${bam_out}"; then
                    logmsg "PE alignment and sorting of ${bam_out} completed"
                else
                    logmsg "ERROR: PE mapping failed for ${bam_out}" && exit 1
                fi
            # Map single-end reads
            else
                logmsg "SE mapping ${sbase} to ${gbase} started"

                if minimap2 -ax sr -t "${threads}" "${idx_mmi}" "${O1}" | samtools sort -@ "${threads}" -o "${bam_out}"; then

                    logmsg "SE alignment and sorting for ${bam_out} completed"
                else
                    logmsg "ERROR: SE mapping failed for ${bam_out}" && exit 1
                fi
            fi

            # Index BAM file
            logmsg "samtools indexing ${bam_out} started"

            if samtools index "${bam_out}"; then
                logmsg "${bam_out} successfully indexed"
            else
                logmsg "Indexing of ${bam_out} failed" && exit 1
            fi
        fi

        logmsg "Read mapping completed for ${sbase}"

        if awk -F '\t' -v genome="${gbase}" -v sample="${sbase}" \
            '$1 == genome && $2 == sample {found=1} END {exit !found}' "${overall_coverage}"; then

            logmsg "Coverage for ${tag} already present. Skipping stats.."
        else
            # Run samtools flagstat and extract relevant values
            local stats total_reads mapped_reads properly_paired depth

            logmsg "Flagstat for ${bam_out}"
            stats=$(samtools flagstat "${bam_out}")

            # Extract read statistics
            total_reads=$(echo "${stats}" | awk '/in total/ {print $1; exit}')

            mapped_reads=$(echo "${stats}" | awk '$4 == "mapped" {print $1; exit}')

            if [[ ${read_type} == "PE" ]]; then
                properly_paired=$(echo "${stats}" | awk '/properly paired/ {print $1; exit}')
            else
                # For single-end data, paired reads are not applicable
                properly_paired="NA"
            fi

            # Check read statistics
            if [[ -z "${total_reads}" || -z "${mapped_reads}" || "${total_reads}" -eq 0 ]]; then
                logmsg "Error: Could not retrieve valid read statistics for ${bam_out}"
                continue
            fi

            local alignment_percentage
            if alignment_percentage=$(awk -v total="${total_reads}" -v mapped="${mapped_reads}" \
                'BEGIN {printf "%.2f", (mapped/total)*100}'); then
                logmsg "Alignment percentage calculated"
            fi

            # compute coverage depth using samtoools - depth
            logmsg "computing coverage depth ${bam_out} started"
            if depth=$(samtools depth -a "${bam_out}" |
                awk '{sum += $3} END {if (NR > 0) printf "%.2f", sum / NR; else print 0}'); then

                logmsg "Coverage calculated for ${bam_out}: ${depth}x"
            else
                logmsg "Coverage calculation failed for ${bam_out}"
                continue
            fi

            # Append structured tab-separated results
            echo -e "${gbase}\t${sbase}\t${read_type}\t${depth}\t${total_reads}\t${mapped_reads}\t${properly_paired}\t${alignment_percentage}" >> "${overall_coverage}"
        fi


        # Variant calling with bcftools mpileup + call
        logmsg  "mpileup and variant calling ${bam_out} against ${gbase} started."
        if [[ -s "${vcf_out}" ]]; then
            logmsg "${vcf_out} already exists. Skipping."
        else
            if ! bcftools mpileup --ignore-RG -f "${genome_fa}" "${bam_out}" --threads "${threads}" \
                | bcftools call -m -v -Oz --threads "${threads}" -o "${vcf_out}"; then
                logmsg "bcftools mpileup for ${tag} failed" && exit 1
            fi
            logmsg "mpileup and variant calling for ${tag} completed."
        fi

        if [[ -s "${vcf_out}.csi" || -s "${vcf_out}.tbi" ]]; then
            logmsg "VCF index for ${vcf_out} already exists. Skipping indexing."
        else
            # index vcf files
            logmsg "indexing ${vcf_out} started."
            if ! bcftools index "${vcf_out}"; then
                logmsg "bcftools index for ${vcf_out} failed" && exit 1
            fi
            logmsg "indexing ${vcf_out} completed."
        fi

        # bcftool stats
        local bstats="${reports_dir}/${tag}_bcfstats.tsv"
        logmsg "Generating stats for VCF: ${vcf_out}"
        if [[ -s ${bstats} ]]; then
            logmsg "Stats file already exists: ${bstats}."
        else
            if ! bcftools stats "${vcf_out}" > "${bstats}"; then
                logmsg "ERROR: bcftools stats failed on ${vcf_out}" && exit 1
            else 
                logmsg "bcftools stats written to ${bstats}."
            fi     
        fi

        # decompress .vcf.gz and pipe to bcftools
        local snp_table="${snps_tsv}/${tag}_snps.tsv"
        logmsg "Creating the Final SNP table for ${snp_table} started."

        if [[ -s ${snp_table} ]]; then
            logmsg "${snp_table} already exists. Skipping."
        else
            echo -e "CHROM\tPOS\tREF\tALT\tQUAL\tDP\tFref\tRref\tFalt\tRalt" > "${snp_table}"

            if bgzip -d -c "${vcf_out}" | grep -E '^#|^chr[0-9]+\b' | 
                bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%DP\t[%DP4{0}]\t[%DP4{1}]\t[%DP4{2}]\t[%DP4{3}]\n' >> "${snp_table}"; then
                logmsg "SNP table created for ${vcf_out}."
            else
                logmsg "bcf filter for ${tag} failed" && exit 1
            fi
        fi

        if [[ ${sv_call} == "true" ]]; then
            local manta_dir=${svs_dir}/manta_svs
            local manta_run="${manta_dir}/${tag}_svs"

            # make sure all dirs exist
            create_dir "${manta_dir}" "${svs_vcf}" "${svs_tsv}" "${manta_run}" || return 1

            # Structural variant calling with Manta
            logmsg "Running Manta for: ${tag} started."

            local vcfs=( "$manta_run"/results/variants/*.vcf.gz )

            if [[ -f "${vcfs[0]}" ]]; then
                logmsg "Manta SVs already exist for ${tag}, skipping Manta."
            
            # Configure Manta and run manta
            elif [[ -f "${manta_run}/runWorkflow.py" ]]; then
                logmsg "Manta already configured for ${tag}. Running workflow."

                if ! "${manta_run}/runWorkflow.py" -m local -j "${threads}" \
                    > "${manta_run}/mantaWorkflow_${tag}.log" 2>&1; then
                    logmsg "Manta workflow failed. Check log: ${manta_run}/mantaWorkflow_${tag}.log"
                    exit 1
                fi
            else
                if ! configManta.py --bam "${bam_out}" --referenceFasta "${genome_fa}" \
                    --runDir "${manta_run}" > "${manta_run}/configManta_${tag}.log" 2>&1; then
                    logmsg "Manta configuration failed for ${tag}."
                    exit 1
                fi

                    if ! "${manta_run}/runWorkflow.py" -m local -j "${threads}" \
                        > "${manta_run}/mantaWorkflow_${tag}.log" 2>&1; then
                        logmsg "Manta workflow failed. Check log: ${manta_run}/mantaWorkflow_${tag}.log"
                        exit 1
                    fi
            fi

            # copy VCFs
            local vrt_dir="${manta_run}/results/variants"
            for svf in candidateSmallIndels.vcf.gz candidateSmallIndels.vcf.gz.tbi \
                candidateSV.vcf.gz candidateSV.vcf.gz.tbi diploidSV.vcf.gz diploidSV.vcf.gz.tbi; do  

                local source="${vrt_dir}/${svf}"
                local destination="${svs_vcf}/${tag}_${svf}"

                if [[ -e "${destination}" ]] ; then  
                    logmsg "${destination} already exists—skipping copy."
                else
                    if ! cp "${source}" "${destination}"; then 
                        logmsg "Failed to copy ${source}"
                        continue
                    fi
                fi
            done

            # Extract important SV fields
            local vf
            for vf in "${svs_vcf}"/*.vcf.gz; do
                local sv_tsv=${svs_tsv}/$(basename ${vf} .vcf.gz).tsv
                
                logmsg "Extracting important SV info from: $(basename ${sv_tsv})"
                # Add headers and extract fields
                if [[ -s ${sv_tsv} ]] ; then
                    logmsg "${sv_tsv} already exists—skipping."
                else
                    echo -e "CHROM\tPOS\tREF\tALT\tQUAL\tFILTER\tSVTYPE\tSVLEN\tEND" > ${sv_tsv}
                    # decompress .vcf.gz and Extract specific fields from the diplodVCF file
                    if bgzip -d -c "${vf}" \
                        | grep -E '^#|^chr[0-9]+\b' \
                        | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/END\n' >> "${sv_tsv}"; then
                        
                        logmsg "$(basename "${sv_tsv}") successfully created."
                    else
                        logmsg "$(basename "${sv_tsv}") failed"
                    fi
                fi
            done
        else
            logmsg "Structural Variant calling for ${tag} not needed. ...skipping."
        fi
    done
}

# before any samples…
for gx in "${goi[@]}"; do
    download_genome "${gx}"
    index_genome   "${gx}"
done

# trim once per sample
for sx in "${samples[@]}"; do
    if project_dir="$(get_project_dir "${sx}")"; then
        logmsg "${sx} ready (dir: ${project_dir})"
        trim_reads "${sx}"
    else
        logmsg "Skipping '${sx}' — missing/empty in sample_loc"
    fi
done


# genome mapping
for gx in "${goi[@]}"; do
    logmsg "Mapping all samples to genome: ${gx}"

    for sx in "${samples[@]}"; do
        if project_dir="$(get_project_dir "${sx}")"; then
            map_reads "${sx}" "${gx}"
        else
            logmsg "Skipping '${sx}' — missing/empty in sample_loc"
        fi
    done
done
