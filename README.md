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

### Requirements for the system

- Have access to HPC with multiple cores.
- Have Bash >=3.2 installed.
- Have Java >=8 installed.
- Have Slurm scheduler managing the jobs in the HPC.
- HPC has Singularity installed and running.

### Setup of the pipeline

You can either clone it by using git (if available in HPC):

`git clone https://github.com/eQTLGen/DataQC.git`

Or just download this from the gitlab/github download link and unzip.

### Input files

- Unimputed genotype file in plink .bed format (https://www.cog-genomics.org/plink/1.9/input#bed). Genome build has to be in **hg19/GRCh37 (default)** or **hg38/GRCh38**. It is advisable that .fam file also includes observed sex for all samples (format: males=1, females=2), so that pipeline does extra check on that. However, if this information is not available for all samples, pipeline just skips this check. Input path has to be without `.bed/.bim/.fam` extension.
- Raw, unprocessed gene expression matrix. Tab-delimited file, genes/probes in the rows, samples in the columns.
    - First column has header "-".
    - For Illumina arrays, probe ID has to be Illumina ArrayAddress.
    - For RNA-seq, gene ID has to be stable ENSEMBL gene ID (ENSEMBL v75).
    - For Affymetrix arrays we expect that gene **expression matrix has already gone through the standard preprocessing** and is in the same format as was used in eQTLGen phase 1 analyses (incl. array probe names).
- Genotype-to-expression linking file (gte). Tab-delimited file, no header, 2 columns: sample ID in genotype data, corresponding sample ID in gene expression data.

### Required inputs

`--cohort_name`                 Name of the cohort.

`--bfile`                       Path to the unimputed genotype files in plink bed/bim/fam format (without extensions bed/bim/fam).

`--genome_build`                Genome build of the cohort. Either hg19, GRCh37, hg38 or GRCh38. Defaults to hg19.

`--expfile`                     Path to the un-preprocessed gene expression matrix (genes/probes in the rows, samples in the columns). Can be from RNA-seq experiment or from array. NB! For Affymetrix arrays (AffyU219, AffyExon) we assume that standard preprocessing and normalisation is already done.

`--gte`                         Genotype-to-expression linking file. Tab-delimited, no header. First column: sample ID for genotype data. Second column: corresponding sample ID for gene expression data. Can be used to filter samples from the analysis.

`--exp_platform`                Gene expression platform. HT12v3, HT12v4, HuRef8, RNAseq, AffyU219, AffyHumanExon.

`--outdir`                      Path to the output directory.
    
### Additional settings

There are five arguments which can be used to adjust certain outlier detection thresholds. These should be adjusted after initial run with the default settings and after investigating the diagnostic plots in the `Report_DataQc_[cohort name].html`. Then the pipeline should be re-run with adjusted settings.


`--GenOutThresh` Threshold for declaring genotype sample genetic outlier, based on LOF "outlierness" metric. Default is 0.4.

`--GenSdThresh` Threshold for declaring genotype sample genetic outlier, based on the deviation from the means of first two genetic PCs. Defaults to 3 SD from the mean.

`--ExpSdThresh` Threshold for declaring expression sample outlier, based on the deviation from the means of first two expression PCs. Defaults to 4 SD from the mean.

`--ContaminationArea` Threshold for declaring samples as contaminated on XIST vs Y-chr gene expression plot. Defaults to 30 degrees, meaning that samples which have high expression of both, X-chr and Y-chr genes are likely contaminated.

Optional arguments:

`--AdditionalCovariates` Tab-separated file with additional external covariates deemed to be relevant for eQTL mapping in given dataset. First column must have header "SampleID" and following columns must include corresponding covariates with informative headers (E.g. "GenotypeBatch", etc.). Categorical covariates must be specified in the text format (E.g. "Batch1", "Batch2", "Batch3"), not encoded as numbers. Pipeline does one-hot encoding for you. Numerical covariates are allowed as well. If specified, this file should include covariate information for each eQTL sample which passes QC and NAs are not allowed. 

`--InclusionList` File with the genotype IDs to keep in the analysis (one per row). Useful for e.g. keeping in only the samples which have part of the biobank, etc. By default, pipeline keeps all samples in. No header needed.

`--ExclusionList` File with the genotype IDs to remove from the analysis (one per row). Useful for removing part of the samples from the analysis if these are from different ancestry. In case of the overlap between inclusion list and exclusion list, intersect is kept in the analysis. No header needed.


### Running the data QC command

Go to folder `dataqc` and modify the Slurm script template `submit_DataQc_pipeline_template.sh` with your input paths. Below is an example template for Slurm scheduler. Some of the paths are pre-filled, assuming that you follow [eQTLGen phase II cookbook](https://github.com/eQTLGen/eQTLGen-phase-2-cookbook/wiki/eQTLGen-phase-II-cookbook) and its recommended folder structure, however you can also use custom paths.

```bash
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

# Define paths 
nextflow_path=../../tools # folder where Nextflow executable is

geno_path=[full path to your input genotype files without .bed/.bim/.fam extension]
exp_path=[full path to your gene expression matrix]
gte_path=[full path to your genotype-to-expression file]
exp_platform=[expression platform name: HT12v3/HT12v4/HuRef8/RNAseq/AffyU219/AffyHumanExon]
cohort_name=[name of the cohort]
output_path=../output # Output path

# Additional settings and optional arguments for the command

# --GenOutThresh [numeric threshold]
# --GenSdThresh [numeric threshold]
# --ExpSdThresh [numeric threshold]
# --ContaminationArea [number between 0 and 90, default 30]
# --AdditionalCovariates [file with additional covariates]
# --InclusionList [file with the list of samples to restrict the analysis]
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
```

You can save the modified script version to informative name, e.g. `submit_DataQc_[**CohortName_PlatformName**].sh`.

Then submit the job `sbatch submit_DataQc_[**CohortName_PlatformName**].sh`. This initiates pipeline, makes analysis environment (using singularity or conda) and automatically submits the steps in correct order and parallel way. Separate `work` directory is made to the folder and contains all interim files.


### Monitoring and debugging

- Monitoring:
  - Monitor the `slurm-***.out` log file and check if all the steps finish without error. Trick: command `watch tail -n 20 slurm-***.out` helps you to interactively monitor the status of the jobs.
  - Use `squeue -u [YourUserName]` to see if individual tasks are in the queue.
- If the pipeline crashes (e.g. due to walltime), you can just resubmit the same script after the fixes. Nextflow does not rerun completed steps and continues only from the steps which had not completed.
- When the work has finished, download and check the job report. This file  is automatically written to your output folder `pipeline_info` subfolder, for potential errors or warnings. E.g. `output/pipeline_info/DataQcReport.html`.
- When you need to do some debugging, then you can use the last section of aforementioned report to figure out in which subfolder from `work` folder the actual step was run. You can then navigate to this folder and investigate the following hidden files:
  - `.command.sh`: script which was submitted
  - `.command.log`: log file for seeing the analysis outputs/errors.
  - `.command.err`: file which lists the errors, if any.


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
- `output/pipeline_info/DataQc_report.html`: Technical pipeline runtime report. It is not automatically added to shared folder, however it can be used in the central site for debugging, if need arises.

## Acknowledgements

Genotype QC and covariate preparations make extensive use of the [bigsnpr package](https://privefl.github.io/bigsnpr/) and [plink 2](https://www.cog-genomics.org/plink/2.0/).

Gene expression processing makes use of [preprocesscore](https://bioconductor.org/packages/release/bioc/html/preprocessCore.html) and [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) R packages.

### Citation

[Privé, F., Luu, K., Blum, M. G. B., McGrath, J. J., &#38; Vilhjálmsson, B. J. (2020). Efficient toolkit implementing best practices for principal component analysis of population genetic data. <i>Bioinformatics</i>, <i>36</i>(16), 4449–4457. https://doi.org/10.1093/BIOINFORMATICS/BTAA520](https://academic.oup.com/bioinformatics/article/36/16/4449/5838185)

[Chang, C. C., Chow, C. C., Tellier, L. C. A. M., Vattikuti, S., Purcell, S. M., &#38; Lee, J. J. (2015). Second-generation PLINK: Rising to the challenge of larger and richer datasets. GigaScience, 4(1). https://doi.org/10.1186/s13742-015-0047-8](https://academic.oup.com/gigascience/article/4/1/s13742-015-0047-8/2707533)

[Bolstad B (2021). preprocessCore: A collection of pre-processing functions. R package version 1.56.0, https://github.com/bmbolstad/preprocessCore.](https://bioconductor.org/packages/release/bioc/html/preprocessCore.html)

[Robinson MD, McCarthy DJ, Smyth GK (2010). “edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.” Bioinformatics, 26(1), 139-140. doi: 10.1093/bioinformatics/btp616.](10.1093/bioinformatics/btp616)

### Contacts

For this Nextflow pipeline: urmo.vosa at gmail.com.