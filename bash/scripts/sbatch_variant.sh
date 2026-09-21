#!/bin/bash

#SBATCH --nodes=1  # Number of tasks
#SBATCH --cpus-per-task=8   # Number of cores
#SBATCH --output="variant_%j.out"
#SBATCH --error="variant_%j.error"
#SBATCH -A leiboff_lab           # Lab account (priority)
#SBATCH -p leiboff_lab           # Lab partition
#SBATCH -w cerebro               # Target cerebro node (optional)

# Define variables
env_name="bsa_seq_pkgs"
script=/nfs5/BPP/Leiboff_Lab/Brian/scripts/bsa-seq-pipeline/scripts/config/run_variant.sh

# Make allocated CPUs available to pipeline
export threads="${SLURM_CPUS_PER_TASK:-8}"

# Load Conda environment management
source /local/cqls/software/x86_64/miniforge3/etc/profile.d/conda.sh

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
logmsg() { echo "$(timestamp): $*"; }

# Create Conda environment if it doesn't exist
if ! conda env list | awk '{print $1}' | grep -Fxq "${env_name}"; then
    logmsg "Creating Conda environment '${env_name}'..."

    conda create \
        --name "${env_name}" \
        --channel conda-forge --channel bioconda --channel defaults --strict-channel-priority \
        libboost=1.85 minimap2 fastp jq manta=1.6.0 vcftools samtools bcftools python tectonic \
        htslib seqkit gatk4 snpeff sift4g -y || {
            logmsg "ERROR: Failed to create Conda environment: ${env_name}"
            exit 1
        }
else
    logmsg "Conda environment '${env_name}' already exists."
fi


# Activate the environment
conda activate "${env_name}" || {
    logmsg "Error: Failed to activate Conda environment '${env_name}'."
    exit 1
}

# Job information
logmsg "Job ID: ${SLURM_JOB_ID}"
logmsg "Node: ${SLURM_NODELIST}"
logmsg "Threads: ${threads}"

# Export environment details for reproducibility
conda env export -n "${env_name}" > "${env_name}.yaml"

# Run processing script
if [ -x "${script}" ]; then
    logmsg "Running script: ${script}"
    "${script}" && logmsg "Successfully ran ${script}."
else
    logmsg "Error: Script not found or not executable: ${script}"
    conda deactivate
    exit 1
fi

# Deactivate and clean up environment
conda deactivate

logmsg "SLURM job completed."