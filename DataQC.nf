def helpMessage() {
    log.info"""
    =======================================================
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\'
        |\\ | |__  __ /  ` /  \\ |__) |__         }  {
        | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                              `._,._,\'
     DataQC v${workflow.manifest.version}
    =======================================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run DataQC.nf \
        --bfile EstBB_HT12v3\
        --expfile EstBB_HT12v3_exp.txt\
        --gte gte_EstBB_HT12v3.txt\
        --exp_platform HT12v3\
        --cohort_name EstBB_HT12v3\
        --outdir EstBB_HT12v3_PreImputationQCd\
        -profile slurm\
        -resume

    Mandatory arguments:
      --bfile                       Path to the unimputed plink bgen files (without extensions bed/bim/fam).
      --expfile                     Path to the unpreprocessed gene expression matrix (genes/probes in the rows, samples in the columns). Can be from RNA-seq experiment or from array.
      --gte                         Genotype-to-expression linking file. Tab-delimited, no header. First column: sample ID for genotype data. Second column: corresponding sample ID for gene expression data. 
      --exp_platform                Indicator indicating the gene expression platform. HT12v3, HT12v4, RNAseq, AffyU219, AffyExon.
      --outdir                      Path to the output directory.
      --Sthresh                     "Outlierness" score threshold for excluding ethnic outliers. Defaults to 4 but should be adjusted according to visual inspection.
      --ExpSdThreshold              Standard deviation threhshold for excluding gene expression outliers. By default, samples away by 3 SDs from the median of PC1 are removed.

    """.stripIndent()
}

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


// Define input channels
Channel
    .from(params.bfile)
    .map { study -> [file("${study}.bed"), file("${study}.bim"), file("${study}.fam")]}
    .set { bfile_ch }

Channel
    .from(params.expfile)
    .map { study -> [file("${study}.txt")]}
    .set { expfile_ch }

Channel
    .from(params.gte)
    .map { study -> [file("${study}.txt")]}
    .set { gte_ch }

Channel
    .fromPath('bin/Report_template.Rmd')
    .set { report_ch }

params.sthresh = 0.4
params.exp_platform = ''
params.cohort_name = ''
params.outdir = ''

// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'
GenotypeGC v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']            = 'DataQC'
summary['Pipeline Version']         = workflow.manifest.version
summary['Run Name']                 = custom_runName ?: workflow.runName
summary['PLINK bfile']              = params.bfile
summary['S threshold']              = params.sthresh
summary['Expression matrix']        = params.expfile
summary['GTE file']                 = params.gte
summary['Max Memory']               = params.max_memory
summary['Max CPUs']                 = params.max_cpus
summary['Max Time']                 = params.max_time
summary['Cohort name']              = params.cohort_name
summary['Expression platform']      = params.exp_platform
summary['Output dir']               = params.outdir
summary['Working dir']              = workflow.workDir
summary['Container Engine']         = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']             = "$HOME"
summary['Current user']             = "$USER"
summary['Current path']             = "$PWD"
summary['Working dir']              = workflow.workDir
summary['Script dir']               = workflow.projectDir
summary['Config Profile']           = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']            = params.awsregion
   summary['AWS Queue']             = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "========================================="

process GenotypeQC {

    tag {GenotypeQC}

    input:
      set file(bfile), file(bim), file(fam) from bfile_ch
      val s_stat from params.sthresh

    output:
      path ('outputfolder_gen') into output_ch_genotypes
      file 'outputfolder_gen/gen_data_QCd/*_ToImputation.fam' into gen_samples

      """
      Rscript --vanilla $baseDir/bin/GenQcAndPosAssign.R  \
      --target_bed ${bfile} \
      --sample_list $baseDir/data/unrelated_reference_samples_ids.txt \
      --pops $baseDir/data/1000G_pops.txt \
      --S_threshold ${s_stat} \
      --output outputfolder_gen
      """
}

process GeneExpressionQC {

    tag {GeneExpressionQC}

    input:
      file exp_mat from expfile_ch
      file gte from gte_ch
      file gen_samples from gen_samples
      val exp_platform from params.exp_platform

    output:
      path ('outputfolder_exp') into output_ch_geneexpression

      script:
      if (exp_platform == 'HT12v3')
      """
      Rscript --vanilla $baseDir/bin/ProcessExpression.R  \
      --expression_matrix ${exp_mat} \
      --genotype_to_expression_linking ${gte} \
      --genotype_samples ${gen_samples} \
      --platform ${exp_platform} \
      --emp_probe_mapping $baseDir/data/EmpiricalProbeMatching_Illumina_HT12v3_20170808.txt \
      --output outputfolder_exp
      """
      else if (exp_platform == 'HT12v4')
      """
      Rscript --vanilla $baseDir/bin/ProcessExpression.R  \
      --expression_matrix ${exp_mat} \
      --genotype_to_expression_linking ${gte} \
      --genotype_samples ${gen_samples} \
      --platform ${exp_platform} \
      --emp_probe_mapping $baseDir/data/EmpiricalProbeMatching_Illumina_HT12v4_20170808.txt \
      --output outputfolder_exp
      """
      else if (exp_platform == 'RNAseq')
      """
      Rscript --vanilla $baseDir/bin/ProcessExpression.R  \
      --expression_matrix ${exp_mat} \
      --genotype_to_expression_linking ${gte} \
      --genotype_samples ${gen_samples} \
      --platform ${exp_platform} \
      --emp_probe_mapping NULL \
      --output outputfolder_exp
      """
}

process RenderReport {

    tag {RenderReport}

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
      path output_gen from output_ch_genotypes
      path output_exp from output_ch_geneexpression
      path report from report_ch
      val exp_platform from params.exp_platform

    output:
      path ('outputfolder_gen/*') into output_ch2
      path ('outputfolder_exp/*') into output_ch3
      path ('Report_DataQc.html') into report_ch2
      path ('CovariatePCs.txt') into combined_covariates

      """
      # Make combined covariate file
      Rscript --vanilla $baseDir/bin/MakeCovariateFile.R

      # Make report
      cp -L ${report} notebook.Rmd

      R -e 'library(rmarkdown);rmarkdown::render("notebook.Rmd", "html_document", 
      output_file = "Report_DataQc.html", params = list(dataset_name = "${params.cohort_name}", platform = "${exp_platform}"))'
      """
}
