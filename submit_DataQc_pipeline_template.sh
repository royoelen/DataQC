#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=6G
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --job-name="DataQc"

# These are needed modules in UT HPC to get singularity and Nextflow running. Replace with appropriate ones for your HPC.
module load java-1.8.0_40
module load singularity/3.5.3
module load squashfs/4.4

# If you follow the eQTLGen phase II cookbook and analysis folder structure,
# some of the following paths are pre-filled.
# https://github.com/eQTLGen/eQTLGen-phase-2-cookbook/wiki/eQTLGen-phase-II-cookbook

# We set the following variables for nextflow to prevent writing to your home directory (and potentially filling it completely)
# Feel free to change these as you wish.
export SINGULARITY_CACHEDIR=../../singularitycache
export NXF_HOME=../../nextflowcache

# Disable pathname expansion. Nextflow handles pathname expansion by itself.
set -f

# Define paths
nextflow_path=../../tools # folder where Nextflow executable is

# Genotype data
bfile_path=[full path to your input genotype files without .bed/.bim/.fam extension]

# Other data
exp_path=[full path to your gene expression matrix]
gte_path=[full path to your genotype-to-expression file]
exp_platform=[expression platform name: HT12v3/HT12v4/HuRef8/RNAseq/AffyU219/AffyHumanExon]
cohort_name=[name of the cohort]
genome_build="GRCh37"
output_path=../output # Output path

# Additional settings and optional arguments for the command

# --GenOutThresh [numeric threshold]
# --GenSdThresh [numeric threshold]
# --ExpSdThresh [numeric threshold]
# --ContaminationArea [number between 0 and 90, default 30]
# --ContaminationSlope [number between 0 and 90, default 45]
# --InclusionList [file with the list of samples to restrict the analysis]
# --ExclusionList [file with the list of samples to remove from the analysis]
# --preselected_sex_check_vars "data/Affy6_pruned_chrX_variant_positions.txt"
# --AdditionalCovariates [file with additional covariates. First column should be `SampleID`]
# --gen_qc_steps 'WGS'
# --fam [PLINK .fam file. Takes precedence over .fam file supplied with `--bfile`]
# --plink_executable [path to plink executable (PLINK v1.90b6.26 64-bit)]
# --plink2_executable [path to plink2 executable (PLINK v2.00a3.7LM 64-bit Intel)]
# --reference_1000g_folder [path to folder with 1000G reference data]
# --chain_path [folder with hg19->hg38 and hg38->hg19 chain files]

# Command:
NXF_VER=21.10.6 ${nextflow_path}/nextflow run DataQC.nf \
--bfile ${bfile_path} \
--expfile ${exp_path} \
--gte ${gte_path} \
--exp_platform ${exp_platform} \
--cohort_name ${cohort_name} \
--genome_build ${genome_build} \
--outdir ${output_path}  \
-profile slurm,singularity \
-resume
