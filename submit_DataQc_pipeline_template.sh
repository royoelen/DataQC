#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=6G
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=[your e-mail@email.com]
#SBATCH --job-name="DataQc"

# These are needed modules in UT HPC to get singularity and Nextflow running. Replace with appropriate ones for your HPC.
module load java-1.8.0_40
module load singularity/3.5.3
module load squashfs/4.4

# Define paths
nextflow_path=[full path to your Nextflow executable]

geno_path=[full path to your input genotype folder]
exp_path=[full path to your raw gene expression matrix]
gte_path=[full path to your genotype-to-expression file]
exp_platform=[expression platform name e.g. HT12v3 or RNAseq]
cohort_name=[name of the cohort]
output_path=[name of the output path]

# Command
${nextflow_path}/nextflow run DataQC.nf \
--bfile ${geno_path} \
--expfile ${exp_path} \
--gte ${gte_path} \
--exp_platform ${exp_platform} \
--cohort_name ${cohort_name} \
--outdir ${output_path}  \
-profile slurm \
-resume
