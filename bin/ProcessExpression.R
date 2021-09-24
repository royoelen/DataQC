library(data.table)
library(preprocessCore)
library(edgeR)
library(ggplot2)
library(optparse)
library(patchwork)

setDTthreads(8)

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
    make_option(c("-o", "--output"), type = "character",
    help = "Output folder where to put preprocessed data matrix, expression PCs, etc.")
    )

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

# Debug
print(args$expression_matrix)
print(args$genotype_to_expression_linking)
print(args$genotype_samples)
print(args$emp_probe_mapping)
print(args$output)

# Make output folder structure
dir.create(args$output)
dir.create(paste0(args$output, "/exp_plots"))
dir.create(paste0(args$output, "/exp_data_QCd"))
dir.create(paste0(args$output, "/exp_PCs"))
dir.create(paste0(args$output, "/exp_data_summary"))

# functions
Z_transform <- function(x){
    z <- (x - mean(x))/sd(x)
    return(z)}

INT_transform <- function(x){
    int <- qnorm((rank(x, na.last = "keep")-0.5)/sum(!is.na(x)))
    return(int)
}
comp_cv <- function(x){sd(x) / mean(x)}
shap_test <- function(x){shapiro.test(x)$p.value}

illumina_HT12v3_preprocess <- function(exp, gte, gen){
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
    geno_fam <- fread(args$genotype_samples, header = FALSE)
    gte <- gte[gte$V1 %in% geno_fam$V2,]
    
    message(paste(nrow(gte), "overlapping samples in the gte file AND genotype data"))

    exp <- exp[, colnames(exp) %in% gte$V2]

    # quantile normalization
    and_n <- normalize.quantiles(exp, copy = FALSE)
    colnames(and_n) <- colnames(exp)
    rownames(and_n) <- rownames(exp)

    # log2 transformation
    and_n <- log2(and_n)

    return(and_n)
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
    unique_values = unique_values,
    shapiro_P = per_gene_shapiro)

    return(gene_summary)
}

# Read in raw expression matrix
and <- fread(args$expression_matrix)
colnames(and)[1] <- "Feature"
message(paste("Initially:", nrow(and), "genes/probes and ", ncol(and), "samples"))

summary_table <- data.table(Stage = "Unprocessed matrix", Nr_of_features = nrow(and), Nr_of_samples = ncol(and))


if (!args$platform %in% c("HT12v3", "HT12v4", "RNAseq", "AffyU291", "AffyHuEx")){stop("Platform has to be one of HT12v3, HT12v4, RNAseq, AffyU291, AffyHuExs")}

# 1.1 data is already appropriately processed ----

# 1.2 Raw Illumina HT12v3 expression array -----
if (args$platform == "HT12v3"){
and_p <- illumina_HT12v3_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
}

summary_table_temp <- data.table(Stage = "After inital preprocessing", Nr_of_features = nrow(and_p), Nr_of_samples = ncol(and_p))
summary_table <- rbind(summary_table, summary_table_temp)

# 1.3 Raw RNA-seq count matrix ----


# QC ----
## Calculate PCs and visualise

pcs <- prcomp(t(and_p), center = FALSE, scale. = FALSE)
PCs <- as.data.table(pcs$x)
PCs$sample <- rownames(pcs$x)

importance <- pcs$sdev^2/sum(pcs$sdev^2)
summary_pcs <- data.table(PC = paste0("PC", 1:50), explained_variance = importance[1:50])

summary_pcs$PC <- factor(summary_pcs$PC, levels = paste0("PC", 1:50))
fwrite(summary_pcs, paste0(args$output, "/exp_data_summary/", "summary_raw_pcs.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

p <- ggplot(summary_pcs, aes(x = PC, y = explained_variance)) + geom_bar(stat = "identity") + theme_bw()
ggsave(paste0(args$output, "/exp_plots/PCA_raw_scree_plot.png"), height = 6, width = 17, units = "in", dpi = 400, type = "cairo")

## Remove outliers (PC 1 and 2, remove what lays out 3 SDs)
PCs$outlier <- "no"
PCs1_med <- median(PCs$PC1)
PCs1_sd <- sd(PCs$PC1)
PCs[(PCs$PC1 > PCs1_med + 3 * PCs1_sd) | (PCs$PC1 < PCs1_med - 3 * PCs1_sd)]$outlier <- "yes"

PCs2_med <- median(PCs$PC2)
PCs2_sd <- sd(PCs$PC2)
PCs[(PCs$PC2 > PCs2_med + 3 * PCs2_sd) | (PCs$PC2 < PCs2_med - 3 * PCs2_sd)]$outlier <- "yes"

p1 <- ggplot(PCs, aes(x = PC1, y = PC2, colour = outlier)) + geom_point(alpha = 0.3) + theme_bw() + scale_colour_manual(values = c("no" = "black", "yes" = "red")) + ggtitle("Before outlier removal\nnormalised\nlog-transformed")

# Remove outliers from primary data
non_outliers <- PCs[!PCs$outlier %in% c("yes"), ]$sample
and <- and[, colnames(and) %in% c(non_outliers, "Feature"), with = FALSE]

## Re-process, re-calculate PCs, re-visualise and write out
and_p <- illumina_HT12v3_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
and_p <- apply(and_p, 1, Z_transform)
and_p <- apply(and_p, 2, INT_transform)
and_p <- t(abi_p)

summary_table_temp <- data.table(Stage = "After removal of expression outliers", Nr_of_features = nrow(and_p), Nr_of_samples = ncol(and_p))
summary_table <- rbind(summary_table, summary_table_temp)

pcs <- prcomp(t(and_p), center = FALSE, scale. = FALSE)
PCs <- data.table(Sample = rownames(pcs$x), pcs$x)

p2 <- ggplot(PCs, aes(x = PC1, y = PC2)) + geom_point(alpha = 0.3) + theme_bw() + ggtitle("After outlier removal\nnormalised\nlog-transformed\nZ-transformed\ninverse normal transformed")

p <- p1 + p2
ggsave(paste0(args$output, "/exp_plots/PCA_before_and_after.png"), height = 8, width = 15, units = "in", dpi = 400, type = "cairo")

fwrite(and_p, paste0(args$output, "/exp_data_QCd/exp_data_preprocessed.txt"), sep = "\t", quote = FALSE)

# Write out importance of PCs
importance <- pcs$sdev^2/sum(pcs$sdev^2)
summary_pcs <- data.table(PC = paste0("PC", 1:50), explained_variance = importance[1:50])
summary_pcs$PC <- factor(summary_pcs$PC, levels = paste0("PC", 1:50))

fwrite(PCs[, c(1:51)], paste0(args$output, "/exp_PCs/exp_PCs.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
fwrite(summary_pcs, paste0(args$output, "/exp_data_summary/", "summary_pcs.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

p <- ggplot(summary_pcs, aes(x = PC, y = explained_variance)) + geom_bar(stat = "identity") + theme_bw()
ggsave(paste0(args$output, "/exp_plots/PCA_final_scree_plot.png"), height = 6, width = 17, units = "in", dpi = 400, type = "cairo")

# Summary statistics ----
# Per gene, calculate mean, median, min, max, sd, skewness, centrality and shapiro test P-value
# before preporcessing
and <- as.data.frame(and)
colnames(and)[1] <- "gene"
and$gene <- as.character(and$gene)
and <- and[, -1]

gene_summary <- exp_summary(and)

emp <- fread(args$emp_probe_mapping)
emp <- emp[, c(1, 2), with = FALSE]

gene_summary <- merge(as.data.table(gene_summary), emp, by.x = "gene", by.y = "Probe")

gene_summary <- gene_summary[, c(ncol(gene_summary), 1:(ncol(gene_summary) - 1)), with = FALSE]
colnames(gene_summary)[1:2] <- c("gene", "probe")

fwrite(gene_summary, paste0(args$output, "/exp_data_summary/", "raw_gene_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# on fully processed expression matrix
gene_summary <- exp_summary(and_p)
fwrite(gene_summary, paste0(args$output, "/exp_data_summary/", "processed_gene_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# Write out summary table about features and samples
fwrite(summary_table, paste0(args$output, "/exp_data_summary/", "summary_table.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
