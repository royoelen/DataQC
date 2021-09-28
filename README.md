# eQTLGen data QC pipeline

Automatic data quality check and processing for unimputed genotype data and unprocessed gene expression data.

Performs the following main steps:

- 

## Usage information

### Input files

- Unimputed genotype file in plink .bed format (https://www.cog-genomics.org/plink/1.9/input#bed). Has to be in **hg19**.
- Raw, unprocessed gene expression matrix. Tab-delimited file, genes/probes in the rows, samples in the columns.
    - For Illumina arrays, probe ID has to be Illumina ArrayAddress.
    - First column has header - and contains either ENSEMBL gene IDs (RNA-seq), Illumina ArrayAddress IDs for probes or probe names for corresponding Affymetrix array.
    - For RNA-seq, gene ID has to be stable ENSEMBL gene ID (ENSEMBL v75).
    - For Affymetrix arrays: information TBA.
- Genotype-to-expression linking file (gte). Tab-delimited file, no header, 2 columns: sample ID in genotype data, corresponding sample ID in gene expression data.

### Running the data QC command

This is example using SLURM scheduler:

```
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
--cohort_name [MyCohortName] \
--outdir ${output_path}  \
-profile slurm \
-resume

```

### Output

Pipeline makes the following output (most relevant files outlined):

```
|--output
    |--outputfolder_exp
        |--exp_data_QCd
            |--exp_data_preprocessed.txt
        |--exp_PCs
            |--exp_PCs.txt
        |--exp_data_summary
            |--raw_gene_summary.txt
            |--processed_gene_summary.txt
            |--...
        |--exp_plots
            |--...
    |--outputfolder_gen
        |--gen_data_QCd
            |--*_ToImputation.bed
            |--*_ToImputation.bim
            |--*_ToImputation.fam
            |--*_ToImputation.log
        |--gen_PCs
            |--GenotypePCs.txt
        |--gen_data_summary
            |--1000G_PC_projections.txt
            |--...
        |--gen_plots
            |--...
    |--pipeline_info
        |--DataQc_report.html
        |--...
    |--Report_DataQc.html
    |--CovariatesPCs.txt
```

#### Steps to take

1. Investigate the file `Report_DataQc.html`, fix any issues with the genotype/gene expression data according to the plots and instructions. 

2. If data had any quality issues, re-run the pipeline when issues are removed, check the `Report_DataQc.html` again.

When all issues are solved:

3. Files: `output/outputfolder_gen/gen_data_QCd/*_ToImputation.bed`, `output/outputfolder_gen/gen_data_QCd/*_ToImputation.bim`, `output/outputfolder_gen/gen_data_QCd/*_ToImputation.fam` are the filtered and QCd genotype files which need to be the input for imputation pipeline `https://gitlab.com/eqtlgen-group/eqtlgen-imputation-pipeline`

4. File: `output/outputfolder_exp/exp_data_QCd/exp_data_preprocessed.txt` is the filtered, QCd and preprocessed file which needs to be the input for per-cohort preparations and encoding pipeline `https://gitlab.com/eqtlgen-group/PerCohortPreparations`.

5. File: `output/CovariatePCs.txt` is the covariate file which contains 10 first genotype PCs and 50 first gene expression PCs. This should be used as an input for per-cohort preparations and encoding pipeline `https://gitlab.com/eqtlgen-group/PerCohortPreparations`.

**TODO!** Combine with cell metric PGS calculation pipeline/output.

5. Share with central site the following files:

    - `output/outputfolder_exp/exp_data_summary/raw_gene_summary.txt`: gene expression summary statistics (mean, median, sd, min, max, Shapiro test P) before processing, used in central site to filter out lowly expressed genes, containing outliers, etc.
    - `output/outputfolder_exp/exp_data_summary/processed_gene_summary.txt`: gene expression summary statistics (mean, median, sd, min, max, Shapiro test P) after processing, used in central site to filter out lowly expressed genes, containing outliers, etc.
    - `output/Report_DataQc.html`: Data QC report, used in central site to check the per-cohort QC plots. 
    - `output/pipeline_info/DataQc_report.html`: Technical pipline runtime report, used in central site for debugging, if such need arises. 
