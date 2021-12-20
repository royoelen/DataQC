library(data.table)
library(preprocessCore)
library(edgeR)
library(ggplot2)
library(optparse)
library(patchwork)

setDTthreads(8)

# TODO: for gene QC report on RNA-seq, consider including the following metrics from GTEx:
# Genes were selected based on expression thresholds of >0.1 TPM in at least 20% of samples and ≥6 reads in at least 20% of samples.
# For calculation of TPM: for 19,942 genes extract the total length of gene CDs, based on latest ENSEMBL
# https://gist.github.com/slowkow/c6ab0348747f86e2748b


# Argument parser
option_list <- list(
    make_option(c("-e", "--expression_matrix"), type = "character",
    help = "Unprocessed gene expression matrix from array or RNA-seq experiment. Samples in columns, genes/probes in the rows."),
    make_option(c("-l", "--genotype_to_expression_linking"), type = "character",
    help = "Genotype-to-expression linking file. No header. First column- sample IDs in the genotype file, second column- corresponding sample IDs in the gene expression matrix."),
    make_option(c("-g", "--genotype_samples"), type = "character",
    help = ".fam file with the samples for which there is genotype data available. Used for filtering those in."),
    make_option(c("-p", "--platform"), type = "character",
    help = "Gene expression platform. This determines the normalization method and replaces probes with best-matching genes based on empirical probe mapping. One of: HT12v3, HT12v4, RNAseq, AffyU291, AffyHuEx."),
    make_option(c("-m", "--emp_probe_mapping"), type = "character",
    help = "Empirical probe matching file. Used to link the best array probe to each blood-expressed gene."),
    make_option(c("-s", "--sd"), type = "double", default = 4,
    help = "Standard deviation threshold for removing expression samples. By default, samples away 3 SDs from the median of PC1 are removed."),
    make_option(c("-o", "--output"), type = "character",
    help = "Output folder where to put preprocessed data matrix, expression PCs, etc.")
    )

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

# Debug
print(args$expression_matrix)
print(args$genotype_to_expression_linking)
print(args$genotype_samples)
print(args$platform)
print(args$emp_probe_mapping)
print(args$sd)
print(args$output)

# Make output folder structure
dir.create(args$output)
dir.create(paste0(args$output, "/exp_plots"))
dir.create(paste0(args$output, "/exp_data_QCd"))
dir.create(paste0(args$output, "/exp_PCs"))
dir.create(paste0(args$output, "/exp_data_summary"))

#############
# functions #
#############
Z_transform <- function(x){
    z <- (x - mean(x))/sd(x)
    return(z)
    }

center_data <- function(x){
    z <- (x - mean(x))
    return(z)
    }

INT_transform <- function(x){
    int <- qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
    return(int)
    }

comp_cv <- function(x){sd(x) / mean(x)}
shap_test <- function(x){
    if (length(unique(x)) > 1) {
        return(shapiro.test(x)$p.value)
    } else {
        return(NA)
    }
}

illumina_array_preprocess <- function(exp, gte, gen){
    # Leave in only probes for which there is empirical probe mapping info and convert data into matrix
    emp <- fread(args$emp_probe_mapping)
    emp <- emp[, c(1, 2), with = FALSE]
    colnames(exp)[1] <- "Probe"

    exp$Probe <- as.character(exp$Probe)
    emp$Probe <- as.character(emp$Probe)

    exp <- merge(exp, emp, by = "Probe")
    exp <- as.data.frame(exp)
    rownames(exp) <- exp[, ncol(exp)]
    exp <- exp[, -ncol(exp)]
    exp <- exp[, -1]
    exp <- as.matrix(exp)

    # Remove samples which are not in the gte or in genotype data
    gte <- fread(args$genotype_to_expression_linking, header = FALSE)
    gte$V1 <- as.character(gte$V1)
    gte$V2 <- as.character(gte$V2)

    geno_fam <- fread(args$genotype_samples, header = FALSE)
    gte <- gte[gte$V1 %in% geno_fam$V2,]
    gte <- gte[gte$V2 %in% as.character(colnames(exp)),]

    message(paste(nrow(gte), "overlapping samples in the gte file AND genotype data."))

    exp <- exp[, colnames(exp) %in% gte$V2]
    exp <- exp[, base::order(colnames(exp))]
    gte <- gte[base::order(gte$V2), ]

    if(!all(colnames(exp) == gte$V2)){stop("Something went wrong in matching genotype and expression IDs. Please debug!")}
    colnames(exp) <- gte$V1

    # quantile normalization
    and_n <- normalize.quantiles(exp, copy = FALSE)
    colnames(and_n) <- colnames(exp)
    rownames(and_n) <- rownames(exp)

    # log2 transformation (not needed because INT is applied)
    # and_n <- log2(and_n)
    message(paste(ncol(and_n), "samples in normalised expression matrix."))

    return(and_n)
}

RNAseq_preprocess <- function(exp, gte, gen){
    # Leave in only probes for which there is empirical probe mapping info and convert data into matrix
    colnames(exp)[1] <- "Probe"

    exp$Probe <- as.character(exp$Probe)

    exp <- as.data.frame(exp)
    rownames(exp) <- exp$Probe
    exp <- exp[, -1]
    exp <- as.matrix(exp)

    # Remove samples which are not in the gte or in genotype data
    gte <- fread(args$genotype_to_expression_linking, header = FALSE)
    gte$V1 <- as.character(gte$V1)
    gte$V2 <- as.character(gte$V2)

    geno_fam <- fread(args$genotype_samples, header = FALSE)
    gte <- gte[gte$V1 %in% geno_fam$V2,]
    gte <- gte[gte$V2 %in% as.character(colnames(exp)),]

    message(paste(nrow(gte), "overlapping samples in the gte file AND genotype data."))

    exp <- exp[, colnames(exp) %in% gte$V2]
    exp <- exp[, base::order(colnames(exp))]
    gte <- gte[base::order(gte$V2), ]

    if(!all(colnames(exp) == gte$V2)){stop("Something went wrong in matching genotype and expression IDs. Please debug!")}

    colnames(exp) <- gte$V1

    exp_obj <- DGEList(counts = exp)

    # TMM-normalized counts
    exp_obj <- calcNormFactors(exp_obj, method = "TMM")

    exp_n <- cpm(exp_obj, log = FALSE)

    print("RNAseq_preprocess")
    print(str(exp_n))

    print(min(exp_n))
    print(max(exp_n))
    print(any(is.na(exp_n)))

    # log2 transformation (+ add 0.25 for solving issues with log2(0)) (not needed because INT is applied)
    # and_n <- log2(and_n + 0.25)
    # message(paste(ncol(and_n), "samples in normalised expression matrix."))

    return(exp_n)
}

exp_summary <- function(x){

    per_gene_mean <- apply(x, 1, mean)
    per_gene_median <- apply(x, 1, median)
    per_gene_min <- apply(x, 1, min)
    per_gene_max <- apply(x, 1, max)
    per_gene_sd <- apply(x, 1, sd)
    unique_values <- apply(x, 1, function(x) length(unique(x)))
    per_gene_shapiro <- apply(x, 1, shap_test)

    gene_summary <- data.table(gene = rownames(x), 
    mean = per_gene_mean,
    median = per_gene_median,
    min = per_gene_min,
    max = per_gene_max,
    sd = per_gene_sd,
    nr_unique_values = unique_values,
    shapiro_P = per_gene_shapiro)

    return(gene_summary)
}

IterativeOutlierDetection <- function(input_exp, sd_threshold = 1, platform = c("HT12v3", "HT12v4", "RNAseq", "AffyU291", "AffyHuEx")) {
  and <- input_exp
  list_ggplots <- list()
  summary_pcs <- list()
  platform <- match.arg(platform)
  message(paste0("Expression platform is ", platform, ", preprocessing the data..."))

  it_round <- 0
  plot_it_round <- 0
  nr_outliers <- 1

  message(paste("Removing samples which deviate more than", sd_threshold, "SDs from the mean values of PC1 or PC2."))
  message("Starting iterative outlier detection...")

  while (nr_outliers > 0) {

    if (platform %in% c("HT12v3", "HT12v4")){
      and_p <- illumina_array_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
      and_p <- log2(and_p + .Machine$double.eps)
      #and_p <- apply(and_p, 1, INT_transform)
      #and_p <- t(and_p)
      #and_p <- apply(and_p, 1, center_data)
      #and_p <- t(and_p)

    } else if(platform %in% c("RNAseq")){
      and_p <- RNAseq_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
      and_p <- log2(and_p + .Machine$double.eps)
      #and_p <- apply(and_p, 1, INT_transform)
      #and_p <- t(and_p)
      #and_p <- apply(and_p, 1, center_data)
      #and_p <- t(and_p)
    } ## TODO: if needed. Add methods for Affymetrix arrays.
    message("Data preprocessed!")

    it_round <- it_round + 1
    pcs <- prcomp(t(and_p), center = FALSE, scale. = FALSE)
    message("PCA calculation finished!")

    print(str(pcs))

    PCs <- as.data.table(pcs$x)
    PCs$sample <- rownames(pcs$x)

    print(PCs)

    importance <- pcs$sdev^2 / sum(pcs$sdev^2)
    summary_pcs[[it_round]] <- data.table(PC = paste0("PC", 1:50), explained_variance = importance[1:50])

    PCs$outlier <- "no"
    PCs1_mean <- mean(PCs$PC1)
    PCs1_sd <- sd(PCs$PC1)
    PCs[(PCs$PC1 > PCs1_mean + sd_threshold * PCs1_sd) | (PCs$PC1 < PCs1_mean - sd_threshold * PCs1_sd)]$outlier <- "yes"

    PCs2_mean <- mean(PCs$PC2)
    PCs2_sd <- sd(PCs$PC2)
    PCs[(PCs$PC2 > PCs2_mean + sd_threshold * PCs2_sd) | (PCs$PC2 < PCs2_mean - sd_threshold * PCs2_sd)]$outlier <- "yes"

    plot_it_round <- plot_it_round + 1
    list_ggplots[[plot_it_round]] <- ggplot(PCs, aes(x = PC1, y = PC2, colour = outlier)) +
      geom_point(alpha = 0.3) +
      theme_bw() +
      scale_colour_manual(values = c("no" = "black", "yes" = "red")) +
      ggtitle(paste0(it_round, ". iteration round"))
    
    plot_it_round <- plot_it_round + 1
    list_ggplots[[plot_it_round]] <- ggplot(PCs, aes(x = PC3, y = PC4, colour = outlier)) +
      geom_point(alpha = 0.3) +
      theme_bw() +
      scale_colour_manual(values = c("no" = "black", "yes" = "red")) +
      ggtitle(paste0(it_round, ". iteration round"))
    
    plot_it_round <- plot_it_round + 1
    list_ggplots[[plot_it_round]] <- ggplot(PCs, aes(x = PC5, y = PC6, colour = outlier)) +
      geom_point(alpha = 0.3) +
      theme_bw() +
      scale_colour_manual(values = c("no" = "black", "yes" = "red")) +
      ggtitle(paste0(it_round, ". iteration round"))
    
    plot_it_round <- plot_it_round + 1
    list_ggplots[[plot_it_round]] <- ggplot(PCs, aes(x = PC7, y = PC8, colour = outlier)) +
      geom_point(alpha = 0.3) +
      theme_bw() +
      scale_colour_manual(values = c("no" = "black", "yes" = "red")) +
      ggtitle(paste0(it_round, ". iteration round"))

    print(str(PCs))

    non_outliers <- PCs[!PCs$outlier %in% c("yes"), ]$sample

    print(str(non_outliers))

    non_outliers <- gte[as.character(gte$V1) %in% non_outliers, ]$V2
  
    print(str(non_outliers))
    # Remove outlier samples from the unprocessed data
    and <- and[, colnames(and) %in% c(non_outliers, "Feature"), with = FALSE]

    nr_outliers <- length(PCs[PCs$outlier %in% c("yes"), ]$sample)

    if (nr_outliers > 0) {
      message(paste0("Iteration round ", it_round, ". Removed ", nr_outliers, " outlier(s). Re-processing the data and running another round of PCA."))
      summary_table_temp <- data.table(Stage = paste("After removal of expression outliers in ", it_round, " iteration."), Nr_of_features = nrow(and), Nr_of_samples = ncol(and))
      summary_table <- rbind(summary_table, summary_table_temp)
    } else if (nr_outliers == 0) {
      message(paste0("Iteration round ", it_round, ". No outliers detected. Finalizing interactive outlier detection."))
    }
  }
  return(list(exp_mat = and, plots = list_ggplots, summary_pcs = summary_pcs))
}

############
# Analysis #
############
# Read in raw expression matrix
and <- fread(args$expression_matrix)
colnames(and)[1] <- "Feature"
message(paste("Initially:", nrow(and), "genes/probes and ", ncol(and), "samples"))

summary_table <- data.table(Stage = "Unprocessed matrix", Nr_of_features = nrow(and), Nr_of_samples = ncol(and))

# Remove samples which are not in the gte or in genotype data
gte <- fread(args$genotype_to_expression_linking, header = FALSE)
geno_fam <- fread(args$genotype_samples, header = FALSE)
gte <- gte[gte$V1 %in% geno_fam$V2, ]

and <- and[, colnames(and) %in% c("Feature", gte$V2), with = FALSE]

summary_table_temp <- data.table(Stage = "Samples with available genotype info", Nr_of_features = nrow(and), Nr_of_samples = ncol(and))
summary_table <- rbind(summary_table, summary_table_temp)

if (!args$platform %in% c("HT12v3", "HT12v4", "RNAseq", "AffyU291", "AffyHuEx")){stop("Platform has to be one of HT12v3, HT12v4, RNAseq, AffyU291, AffyHuEx")}

iterative_outliers <- IterativeOutlierDetection(and, sd_threshold = args$sd, platform = args$platform) 

print(str(and))

# Keep in the original data only non-outlier samples

exp_non_outliers <- colnames(iterative_outliers$exp_mat)
print("!!!!!!!!!!!!!!!!!!!!!!!!!")
print(head(exp_non_outliers))
#exp_non_outliers <- gte[gte$V1 %in% exp_non_outliers, ]$V2
and <- and[, colnames(and) %in% c("Feature", exp_non_outliers), with = FALSE]

print(str(and))

# Final re-process, re-calculate PCs, re-visualise and write out
if (args$platform %in% c("HT12v3", "HT12v4")){
and_p <- illumina_array_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
}
if (args$platform %in% c("RNAseq")){
and_p <- RNAseq_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
}

# Apply inverse normal transformation and scaling to normalised data.
#and_p <- apply(and_p, 1, Z_transform) # No Z-transform as data will be forced to normal distribution anyway
and_p <- apply(and_p, 1, INT_transform)
and_p <- t(and_p)
#and_p <- apply(and_p, 1, center_data)
#and_p <- t(and_p)

summary_table_temp <- data.table(Stage = "After removal of all expression outliers", Nr_of_features = nrow(and_p), Nr_of_samples = ncol(and_p))
summary_table <- rbind(summary_table, summary_table_temp)

pcs <- prcomp(t(and_p), center = FALSE, scale. = FALSE)
PCs <- data.table(Sample = rownames(pcs$x), pcs$x)

p5 <- ggplot(PCs, aes(x = PC1, y = PC2)) + geom_point(alpha = 0.3) + theme_bw()
p6 <- ggplot(PCs, aes(x = PC3, y = PC4)) + geom_point(alpha = 0.3) + theme_bw()
p7 <- ggplot(PCs, aes(x = PC5, y = PC6)) + geom_point(alpha = 0.3) + theme_bw()
p8 <- ggplot(PCs, aes(x = PC7, y = PC8)) + geom_point(alpha = 0.3) + theme_bw()

p <- (iterative_outliers$plots[[1]] | iterative_outliers$plots[[2]]) / (iterative_outliers$plots[[3]] | iterative_outliers$plots[[4]])
ggsave(paste0(args$output, "/exp_plots/PCA_before.png"), height = 10, width = 11, units = "in", dpi = 300, type = "cairo")

p <- (p5 | p6) / (p7 | p8)
ggsave(paste0(args$output, "/exp_plots/PCA_after.png"), height = 10, width = 10, units = "in", dpi = 300, type = "cairo")

and_p2 <- as.data.table(and_p)
and_p2 <- data.table(`-` = rownames(and_p), and_p2)
fwrite(and_p2, paste0(args$output, "/exp_data_QCd/exp_data_preprocessed.txt"), sep = "\t", quote = FALSE)

# Write out importance of final PCs
importance <- pcs$sdev^2 / sum(pcs$sdev^2)
summary_pcs <- data.table(PC = paste0("PC", 1:100), explained_variance = importance[1:100])
summary_pcs$PC <- factor(summary_pcs$PC, levels = paste0("PC", 1:100))

fwrite(PCs[, c(1:101)], paste0(args$output, "/exp_PCs/exp_PCs.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
fwrite(summary_pcs, paste0(args$output, "/exp_data_summary/", "summary_pcs.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

p <- ggplot(summary_pcs, aes(x = PC, y = explained_variance)) + geom_bar(stat = "identity") + theme_bw() +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(args$output, "/exp_plots/PCA_final_scree_plot.png"), height = 6, width = 17, units = "in", dpi = 300, type = "cairo")

# Summary statistics ----
message("Calculating descriptive summary statistiscs for every gene...")
# Per gene, calculate mean, median, min, max, sd, and shapiro test P-value
# before preprocessing
message("Before normalisation.")
# Remove sample outliers
sample_non_outliers <- colnames(and_p)
sample_non_outliers <- gte[gte$V1 %in% sample_non_outliers, ]$V2

and <- and[, colnames(and) %in% c("Feature", sample_non_outliers), with = FALSE]

and <- as.data.frame(and)
colnames(and)[1] <- "gene"
and$gene <- as.character(and$gene)
rownames(and) <- and$gene
and <- and[, -1]

print(str(and))

gene_summary <- exp_summary(and)

emp <- fread(args$emp_probe_mapping)
emp <- emp[, c(1, 2), with = FALSE]

gene_summary <- merge(as.data.table(gene_summary), emp, by.x = "gene", by.y = "Probe")

gene_summary <- gene_summary[, c(ncol(gene_summary), 1:(ncol(gene_summary) - 1)), with = FALSE]
colnames(gene_summary)[1:2] <- c("gene", "probe")

fwrite(gene_summary, paste0(args$output, "/exp_data_summary/", "raw_gene_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# on fully processed expression matrix
message("After normalization.")

gene_summary <- exp_summary(as.data.table(and_p))
fwrite(gene_summary, paste0(args$output, "/exp_data_summary/", "processed_gene_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# Write out summary table about features and samples
fwrite(summary_table, paste0(args$output, "/exp_data_summary/", "summary_table.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
