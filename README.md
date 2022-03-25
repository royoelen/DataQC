# eQTLGen data QC pipeline

Automatic data quality check and processing for unimputed genotype data and unprocessed gene expression data.

Performs the following main steps:

- Genotypes
    - Standard SNP QC filtering (call-rate>0.95, Hardy-Weinberg P>1e-6, MAF>0.01).
    - Individual-level missingness filter <0.05.
    - Comparison of reported and genetic check, removal of mismatched samples.
    - Removal of samples with unclear genetic sex.
    - Removal of samples with excess heterozygosity (+/-3SD from the mean).
    - Removal of related samples (3d degree relatives). From each pair of related samples, one is kept in the data.
    - Visualisation of samples in the genetic reference space, instructions how to remove or split the data in case of multi-ancestry samples.
    - Removal of in-sample genetic outliers.
    - Calculates 10 first genetic principal components (PCs), used in analyses as covariates to correct for population stratification.
- Gene expression:
    - Aligns sample IDs between genotype samples and gene expression samples.
    - Filters in blood-expressed genes.
    - Replaces array probe IDs with gene IDs.
    - Iteratively checks for expression outliers and removes samples which fail this check.
    - Appropriately normalises the data according to the expression platform used and applies additional inverse normal transformation.
    - Calculates 100 first expression PCs, used in analyses as covariates to correct for unknown variation.
    - Calculates the expression summary statistics for each gene, used to do *post-hoc* QC and gene filtering in the meta-analysis.
- Additional steps:
    - Removes samples whose genetic sex does not match with the expression of sex-specific genes.
    - Reorders the genotype samples into random order.
    - Organises all the QCd data into the standard folder format.
    - Provides commented `html` QC report which should be used to get an overview of the quality of the data.

## Usage information

### Input files

- Unimputed genotype file in plink .bed format (https://www.cog-genomics.org/plink/1.9/input#bed). Genome build to be in **hg19**. It is advisable that .fam file also includes observed sex for all samples (format: males=1, females=2), so that pipeline does extra check on that. However, if this information is not available for all samples, pipeline just skips this check.
- Raw, unprocessed gene expression matrix. Tab-delimited file, genes/probes in the rows, samples in the columns.
    - First column has header "-".
    - For Illumina arrays, probe ID has to be Illumina ArrayAddress.
    - For RNA-seq, gene ID has to be stable ENSEMBL gene ID (ENSEMBL v75).
    - For Affymetrix arrays we expect that gene **expression matrix has already gone through the standard preprocessing** and is in the same format as was used in eQTLGen phase 1 analyses (incl. array probe names).
- Genotype-to-expression linking file (gte). Tab-delimited file, no header, 2 columns: sample ID in genotype data, corresponding sample ID in gene expression data.

### Additional settings

There are three arguments which can be used to adjust certain outlier detection thresholds. These should be adjusted after initial run with the default settings and after investigating the diagnostic plots in the `Report_DataQc_[cohort name].html`. Then the pipeline should be re-run with adjusted settings.

`--GenOutThresh` Threshold for declaring genotype sample genetic outlier, based on LOF "outlierness" metric. Default is 0.4.

`--GenSdThresh` Threshold for declaring genotype sample genetic outlier, based on the deviation from the means of first two genetic PCs. Defaults to 3 SD from the mean.

`--ExpSdThresh` Threshold for declaring expression sample outlier, based on the deviation from the means of first two expression PCs. Defaults to 4 SD from the mean.

`--ContaminationArea` Threshold for declaring samples as contaminated on XIST vs Y-chr gene expression plot. Defaults to 30 degrees, meaning that samples which have high expression of both, X-chr and Y-chr genes are likely contaminated.


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

# Optional arguments for the command
# --GenOutThresh [numeric threshold]
# --GenSdThresh [numeric threshold]
# --ExpSdThresh [numeric threshold]
# --ContaminationArea [numeric threshold]

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
    |--Report_DataQc_[cohort name].html
    |--CovariatesPCs.txt
```

#### Steps to take

1. Investigate the file `Report_DataQc*.html`, fix any issues with the genotype/gene expression data according to the plots and instructions. You can adjust three arguments of the pipeline to adjust certain outlier detection thresholds according to your data.

2. If data had any quality issues, re-run the pipeline when issues are removed, check the `Report_DataQc.html` again.

When all issues are solved:

3. Files: `output/outputfolder_gen/gen_data_QCd/*_ToImputation.bed`, `output/outputfolder_gen/gen_data_QCd/*_ToImputation.bim`, `output/outputfolder_gen/gen_data_QCd/*_ToImputation.fam` are the filtered and QCd genotype files which need to be the input for imputation pipeline `https://gitlab.com/eqtlgen-group/eqtlgen-imputation-pipeline`

4. The whole folder `output` should be specified as an input for per-cohort preparations and encoding pipeline `https://gitlab.com/eqtlgen-group/PerCohortPreparations`. This pipeline automatically uses the processed, QCd expression data and covariate file to run data encoding and partial derivative calculation. It then organises the encoded matrices for sharig with central site. It also extracts some QC files for sharing with the central site:

- `output/Report_DataQc_[cohort name].html`: most important data QC report, used in central site to check the per-cohort QC information.
- `output/outputfolder_exp/exp_data_summary/raw_gene_summary.txt`: gene expression summary statistics (mean, median, sd, min, max, nr of unique values, Shapiro test P) before normalisation, used in central site to filter out lowly expressed genes, genes having outliers, etc.
- `output/outputfolder_exp/exp_data_summary/processed_gene_summary.txt`: gene expression summary statistics (mean, median, sd, min, max, nr of unique values, Shapiro test P) after normalisation, used in central site to filter out lowly expressed genes, genes having outliers, etc.
- `output/outputfolder_gen/plots/*`, `output/outputfolder_exp/plots/*`: Separate diagnostic plots which can be later used in the manuscript materials.
- `output/pipeline_info/DataQc_report.html`: Technical pipeline runtime report, it can be used in central site for debugging.

