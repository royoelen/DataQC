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

# Define paths
# If you follow eQTLGen phase II cookbook, you can use some provided default paths

nextflow_path=../../tools # folder where Nextflow executable is

geno_path=[full path to your input genotype folder]
exp_path=[full path to your gene expression matrix]
gte_path=[full path to your genotype-to-expression file]
exp_platform=[expression platform name: HT12v3/HT12v4/HuRef8/RNAseq/AffyU219/AffyHumanExon]
cohort_name=[name of the cohort]
output_path=../output # Output path, can be kept as is

# Optional arguments for the command
# --GenOutThresh [numeric threshold]
# --GenSdThresh [numeric threshold]
# --ExpSdThresh [numeric threshold]
# --ContaminationArea [number between 0 and 90, default 30]
# --ExclusionList [file with the list of samples to remove from the analysis]

# Command:
NXF_VER=21.10.6 ${nextflow_path}/nextflow run DataQC.nf \
--bfile ${geno_path} \
--expfile ${exp_path} \
--gte ${gte_path} \
--exp_platform ${exp_platform} \
--cohort_name ${cohort_name} \
--outdir ${output_path}  \
-profile slurm,singularity \
-resume
