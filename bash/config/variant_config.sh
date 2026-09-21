#!/bin/bash
# ============================================================================================================================
# VARIANT CALLING CONFIG
# ============================================================================================================================

# Uses Slurm allocation if available, otherwise defaults to 8
threads=${SLURM_CPUS_PER_TASK:-8}

# Run structural variant calling with Manta
sv_call="true" # set to "false" to skip structural variant calling

# ---------------------------------------------------------------------------------------------------------------------------------
# Reference genome location
# ---------------------------------------------------------------------------------------------------------------------------------
ref_dir="/nfs5/BPP/Leiboff_Lab/Brian/genomes_tmp/genomes"


# ---------------------------------------------------------------------------------------------------------------------------------
# Available reference genomes
# Genomes of interest - edit keys & URLs
# ---------------------------------------------------------------------------------------------------------------------------------
declare -A genome_urls=(
    [b73]="https://download.maizegdb.org/Genomes/B73/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz"
    [mo17]="https://download.maizegdb.org/Genomes/Mo17/Zm-Mo17-REFERENCE-CAU-1.0/Zm-Mo17-REFERENCE-CAU-1.0.fa.gz"
    [w22]="https://download.maizegdb.org/Genomes/W22/Zm-W22-REFERENCE-NRGENE-2.0/Zm-W22-REFERENCE-NRGENE-2.0.fa.gz"
    [p39]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-P39-REFERENCE-NAM-1.0/Zm-P39-REFERENCE-NAM-1.0.fa.gz"
    [ki11]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-Ki11-REFERENCE-NAM-1.0/Zm-Ki11-REFERENCE-NAM-1.0.fa.gz"
    [a188]="https://download.maizegdb.org/Genomes/A188/Zm-A188-REFERENCE-KSU-1.0/Zm-A188-REFERENCE-KSU-1.0.fa.gz"
    [a632]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-A632-REFERENCE-CAAS_FIL-1.0/Zm-A632-REFERENCE-CAAS_FIL-1.0.fa.gz"
    [oh43]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-Oh43-REFERENCE-NAM-1.0/Zm-Oh43-REFERENCE-NAM-1.0.fa.gz"
    [cml247]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-CML247-REFERENCE-NAM-1.0/Zm-CML247-REFERENCE-NAM-1.0.fa.gz"
    [tx303]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-Tx303-REFERENCE-NAM-1.0/Zm-Tx303-REFERENCE-NAM-1.0.fa.gz"
    [ky21]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-Ky21-REFERENCE-NAM-1.0/Zm-Ky21-REFERENCE-NAM-1.0.fa.gz"
    [oh7b]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-Oh7B-REFERENCE-NAM-1.0/Zm-Oh7B-REFERENCE-NAM-1.0.fa.gz"
    [tzi8]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-Tzi8-REFERENCE-NAM-1.0/Zm-Tzi8-REFERENCE-NAM-1.0.fa.gz"
    [nc350]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-NC350-REFERENCE-NAM-1.0/Zm-NC350-REFERENCE-NAM-1.0.fa.gz"
    [cml277]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-CML277-REFERENCE-NAM-1.0/Zm-CML277-REFERENCE-NAM-1.0.fa.gz"
    [b97]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-B97-REFERENCE-NAM-1.0/Zm-B97-REFERENCE-NAM-1.0.fa.gz"
    [cml322]="https://download.maizegdb.org/Genomes/NAM_Founders/Zm-CML322-REFERENCE-NAM-1.0/Zm-CML322-REFERENCE-NAM-1.0.fa.gz"
)

# ---------------------------------------------------------------------------------------------------------------------------------
# Available sample locations
# ---------------------------------------------------------------------------------------------------------------------------------
declare -A sample_loc=(
    [q1]="/nfs5/BPP/Leiboff_Lab/Brian/bsa_seq/bds1"
    [q2]="/nfs5/BPP/Leiboff_Lab/Brian/bsa_seq/S348"
    [q3]="/nfs5/BPP/Leiboff_Lab/Brian/bsa_seq/S344"
)


# ---------------------------------------------------------------------------------------------------------------------------------
# Genomes to process
# ---------------------------------------------------------------------------------------------------------------------------------

declare -a goi=( b73 )

# ---------------------------------------------------------------------------------------------------------------------------------
# samples to process
# ---------------------------------------------------------------------------------------------------------------------------------

declare -a samples=( q2 q3 )

# *************** End of user setting **************************************