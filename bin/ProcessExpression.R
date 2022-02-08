library(data.table)
library(preprocessCore)
library(edgeR)
library(ggplot2)
library(optparse)
library(patchwork)
library(MASS)
library(dplyr)

setDTthreads(8)

# Argument parser
option_list <- list(
    make_option(c("-e", "--expression_matrix"), type = "character",
    help = "Unprocessed gene expression matrix from array or RNA-seq experiment. Samples in columns, genes/probes in the rows."),
    make_option(c("-l", "--genotype_to_expression_linking"), type = "character",
    help = "Genotype-to-expression linking file. No header. First column- sample IDs in the genotype file, second column- corresponding sample IDs in the gene expression matrix."),
    make_option(c("-p", "--platform"), type = "character",
    help = "Gene expression platform. This determines the normalization method and replaces probes with best-matching genes based on empirical probe mapping. One of: HT12v3, HT12v4, HuRef8, RNAseq, AffyU219, AffyHumanExon."),
    make_option(c("-m", "--emp_probe_mapping"), type = "character",
    help = "Empirical probe matching file. Used to link the best array probe to each blood-expressed gene."),
    make_option(c("-s", "--sd"), type = "double", default = 4,
    help = "Standard deviation threshold for removing expression samples. By default, samples away 4 SDs from the median of PC1 are removed."),
    make_option(c("-c", "--contamination_area"), type = "double", default = 0.3,
                help = "Area that marks likely contaminated samples based on sex-chromosome gene expression. Must be an angle between 0 and 90. The angle represents the total area around the y = x function."),
    make_option(c("-i", "--sex_info"), type = "character",
    help = "File with sex information. Plink2 --check-sex filtered output."),
    make_option(c("-f", "--geno_filter"), type = "character",
    help = "File with filtered genotype samples. Plink fam file."),
    make_option(c("-o", "--output"), type = "character",
    help = "Output folder where to put preprocessed data matrix, expression PCs, etc.")
    )

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

# Debug
print(args$expression_matrix)
print(args$genotype_to_expression_linking)
print(args$geno_filter)
print(args$sex_info)
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

illumina_array_preprocess <- function(exp, gte, gen, normalize = TRUE){
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

    geno_fam <- fread(args$sex_info, header = TRUE)
    gte <- gte[gte$V1 %in% geno_fam$IID,]
    gte <- gte[gte$V2 %in% as.character(colnames(exp)),]

    message(paste(nrow(gte), "overlapping samples in the gte file AND genotype data."))

    exp <- exp[, colnames(exp) %in% gte$V2]
    exp <- exp[, base::order(colnames(exp))]
    gte <- gte[base::order(gte$V2), ]

    if(!all(colnames(exp) == gte$V2)){stop("Something went wrong in matching genotype and expression IDs. Please debug!")}
    colnames(exp) <- gte$V1

    if (normalize == TRUE){
    # quantile normalization
    exp_n <- normalize.quantiles(exp, copy = FALSE)
    colnames(exp_n) <- colnames(exp)
    rownames(exp_n) <- rownames(exp)
    }else{exp_n <- exp}

    # log2 transformation (not needed because INT is applied)
    # and_n <- log2(and_n)
    message(paste(ncol(exp_n), "samples in normalised expression matrix."))

    return(exp_n)
}

RNAseq_preprocess <- function(exp, gte, gen, normalize = TRUE){
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
    exp <- abs(as.matrix(exp))

    # Remove samples which are not in the gte or in genotype data
    gte <- fread(args$genotype_to_expression_linking, header = FALSE)
    gte$V1 <- as.character(gte$V1)
    gte$V2 <- as.character(gte$V2)

    geno_fam <- fread(args$sex_info, header = TRUE)
    gte <- gte[gte$V1 %in% geno_fam$IID,]
    gte <- gte[gte$V2 %in% as.character(colnames(exp)),]

    message(paste(nrow(gte), "overlapping samples in the gte file AND genotype data."))

    exp <- exp[, colnames(exp) %in% gte$V2]
    exp <- exp[, base::order(colnames(exp))]
    gte <- gte[base::order(gte$V2), ]

    if(!all(colnames(exp) == gte$V2)){stop("Something went wrong in matching genotype and expression IDs. Please debug!")}

    colnames(exp) <- gte$V1

    # Remove genes with no variance
    gene_variance <- data.frame(gene = rownames(exp), gene_variance = apply(exp, 1, var))
    exp <- exp[!rownames(exp) %in% gene_variance[gene_variance$gene_variance == 0, ]$gene, ]

    if (normalize == TRUE){
    # TMM-normalized counts
    exp_n <- DGEList(counts = exp)
    exp_n <- calcNormFactors(exp_n, method = "TMM")
    exp_n <- cpm(exp_n, log = FALSE)
    }else{exp_n <- exp}

    # log2 transformation (+ add 0.25 for solving issues with log2(0)) (not needed because INT will be applied)
    # and_n <- log2(and_n + 0.25)
    message(paste(ncol(exp_n), "samples in normalised expression matrix."))

    return(exp_n)
}

Affy_preprocess <- function(exp, gte, gen){
    # Leave in only probes for which there is empirical probe mapping info and convert data into matrix
    message("For Affymetrix arrays we assume that input expression matrix is already appropriately normalised and transformed.")
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

    geno_fam <- fread(args$sex_info, header = TRUE)
    gte <- gte[gte$V1 %in% geno_fam$IID,]
    gte <- gte[gte$V2 %in% as.character(colnames(exp)),]

    message(paste(nrow(gte), "overlapping samples in the gte file AND genotype data."))

    exp <- exp[, colnames(exp) %in% gte$V2]
    exp <- exp[, base::order(colnames(exp))]
    gte <- gte[base::order(gte$V2), ]

    if(!all(colnames(exp) == gte$V2)){stop("Something went wrong in matching genotype and expression IDs. Please debug!")}

    colnames(exp) <- gte$V1

    # Remove genes with no variance
    gene_variance <- data.frame(gene = rownames(exp), gene_variance = apply(exp, 1, var))
    exp <- exp[!rownames(exp) %in% gene_variance[gene_variance$gene_variance == 0, ]$gene, ]

    # No normalisation done, it assumes that you have already done this.
    message(paste(ncol(exp_n), "samples in normalised expression matrix."))

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

IterativeOutlierDetection <- function(input_exp, sd_threshold = 1, platform = c("HT12v3", "HT12v4", "HuRef8", "RNAseq", "AffyU291", "AffyHuEx")) {
  and <- input_exp
  list_ggplots <- list()
  summary_pcs <- list()
  outliers <- c()
  platform <- match.arg(platform)
  message(paste0("Expression platform is ", platform, ", preprocessing the data..."))

  it_round <- 0
  plot_it_round <- 0
  nr_outliers <- 1

  message(paste("Removing samples which deviate more than", sd_threshold, "SDs from the mean values of PC1 or PC2."))
  message("Starting iterative outlier detection...")

  while (nr_outliers > 0) {

    if (platform %in% c("HT12v3", "HT12v4", "HuRef8")){
      and_p <- illumina_array_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
      and_p <- log2(and_p)
      #and_p <- apply(and_p, 1, INT_transform)
      #and_p <- t(and_p)
      #and_p <- apply(and_p, 1, center_data)
      #and_p <- t(and_p)
    } else if(platform %in% c("RNAseq")){
      and_p <- RNAseq_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
      and_p <- log2(and_p + 0.25)
      #and_p <- apply(and_p, 1, INT_transform)
      #and_p <- t(and_p)
      #and_p <- apply(and_p, 1, center_data)
      #and_p <- t(and_p)
    } else if(platform %in% c("AffyU219", "AffyHumanExon")){
      and_p <- Affy_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
      #and_p <- apply(and_p, 1, INT_transform)
      #and_p <- t(and_p)
      #and_p <- apply(and_p, 1, center_data)
      #and_p <- t(and_p)
    }
    message("Data preprocessed!")

    it_round <- it_round + 1
    pcs <- prcomp(t(and_p), center = FALSE, scale. = FALSE)
    message("PCA calculation finished!")

    PCs <- as.data.table(pcs$x)
    PCs$sample <- rownames(pcs$x)

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

    non_outliers <- PCs[!PCs$outlier %in% c("yes"), ]$sample

    non_outliers <- gte[as.character(gte$V1) %in% non_outliers, ]$V2
  
    # Remove outlier samples from the unprocessed data
    and <- and[, colnames(and) %in% c(non_outliers, "Feature"), with = FALSE]

    nr_outliers <- length(PCs[PCs$outlier %in% c("yes"), ]$sample)

    if (nr_outliers > 0) {
      message(paste0("Iteration round ", it_round, ". Removed ", nr_outliers, " outlier(s). Re-processing the data and running another round of PCA."))
      summary_table_temp <- data.table(Stage = paste("After removal of expression outliers in ", it_round, " iteration."), Nr_of_features = nrow(and), Nr_of_samples = ncol(and))
      summary_table <- rbind(summary_table, summary_table_temp)

      outliers <- c(outliers, PCs[PCs$outlier %in% c("yes"), ]$sample)
    } else if (nr_outliers == 0) {
      message(paste0("Iteration round ", it_round, ". No outliers detected. Finalizing interactive outlier detection."))
    }
    if (it_round == 1){
      PCs_it1 <- PCs
    } else {PCs_it1[PCs_it1$sample %in% outliers, ]$outlier <- "yes"}
  }

  return(list(exp_mat = and_p, plots = list_ggplots, summary_pcs = summary_pcs, PCs_first = PCs_it1))
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
geno_fam <- fread(args$sex_info, header = TRUE)
gen_filter <- fread(args$geno_filter, header = FALSE)
gte <- gte[gte$V1 %in% geno_fam$IID, ]
gte <- gte[gte$V1 %in% gen_filter$V2, ]

and <- and[, colnames(and) %in% c("Feature", gte$V2), with = FALSE]

summary_table_temp <- data.table(Stage = "Samples which overlap with QCd genotype info", Nr_of_features = nrow(and), Nr_of_samples = ncol(and) - 1)
summary_table <- rbind(summary_table, summary_table_temp)

if (!args$platform %in% c("HT12v3", "HT12v4", "HuRef8", "RNAseq", "AffyU219", "AffyHumanExon")){stop("Platform has to be one of HT12v3, HT12v4, HuRef8, RNAseq, AffyU291, AffyHuEx")}

iterative_outliers <- IterativeOutlierDetection(and, sd_threshold = args$sd, platform = args$platform) 

# Keep in the original data only non-outlier samples
exp_non_outliers <- colnames(iterative_outliers$exp_mat)
exp_non_outliers <- gte[gte$V1 %in% exp_non_outliers, ]$V2
# Remove outlier samples from MDS
and <- and[, colnames(and) %in% c("Feature", exp_non_outliers), with = FALSE]

summary_table_temp <- data.table(Stage = "After removal of all expression outliers (PCA)", Nr_of_features = nrow(and), Nr_of_samples = ncol(and))
summary_table <- rbind(summary_table, summary_table_temp)

message("Running MDS")
if (args$platform %in% c("HT12v3", "HT12v4", "HuRef8")){
  and_p <- illumina_array_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
  and_p <- log2(and_p)
}
if (args$platform %in% c("RNAseq")){
  and_p <- illumina_array_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
  and_p <- log2(and_p + 0.25)
}
if (args$platform %in% c("AffyU219", "AffyHumanExon")){
  and_p <- Affy_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
}

dist <- cor(and_p, method = "pearson")
mds <- isoMDS(1 - dist, k = 2)
mds <- as.data.frame(mds$points)
colnames(mds) <- paste("MDS coordinate", 1:2)
mds$Sample <- colnames(and_p)

# Find samples deviating from the mean MDS1 and MDS2
mds$outlier <- "no"
mean_mds1 <- mean(mds$`MDS coordinate 1`)
sd_mds1 <- sd(mds$`MDS coordinate 1`)
mean_mds2 <- mean(mds$`MDS coordinate 2`)
sd_mds2 <- sd(mds$`MDS coordinate 2`)
if (nrow(mds[mds$`MDS coordinate 1` > mean_mds1 + args$sd * sd_mds1 |
mds$`MDS coordinate 1` < mean_mds1 - args$sd * sd_mds1 |
mds$`MDS coordinate 2` > mean_mds2 + args$sd * sd_mds2 |
mds$`MDS coordinate 2` < mean_mds2 - args$sd * sd_mds2, ]) > 0){

mds[mds$`MDS coordinate 1` > mean_mds1 + args$sd * sd_mds1 |
mds$`MDS coordinate 1` < mean_mds1 - args$sd * sd_mds1 |
mds$`MDS coordinate 2` > mean_mds2 + args$sd * sd_mds2 |
mds$`MDS coordinate 2` < mean_mds2 - args$sd * sd_mds2, ]$outlier <- "yes"

}

# Add sex info
sex <- fread(args$sex_info, header = FALSE)
sex <- sex[, c(2, 4), with = FALSE]
colnames(sex) <- c("Sample", "Sex")
mds <- merge(mds, sex, by = "Sample")

p <- ggplot(mds, aes(x = `MDS coordinate 1`, `MDS coordinate 2`, colour = outlier, shape = Sex)) +
geom_point(alpha = 0.3) +
theme_bw() +
scale_colour_manual(values = c("no" = "black", "yes" = "red"))

ggsave(paste0(args$output, "/exp_plots/MDS_before.png"), height = 5, width = 5.5, units = "in", dpi = 300, type = "cairo")
ggsave(paste0(args$output, "/exp_plots/MDS_before.pdf"), height = 5, width = 5.5, units = "in", dpi = 300)

exp_non_outliers <- gte[gte$V1 %in% mds[mds$outlier == "no", ]$Sample, ]$V2
and <- and[, colnames(and) %in% c("Feature", exp_non_outliers), with = FALSE]

summary_table_temp <- data.table(Stage = "After removal of all expression outliers (MDS)", Nr_of_features = nrow(and), Nr_of_samples = ncol(and))
summary_table <- rbind(summary_table, summary_table_temp)

# Final re-process, re-calculate PCs, re-visualise and write out
if (args$platform %in% c("HT12v3", "HT12v4", "HuRef8")){
and_pp <- illumina_array_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
} else if (args$platform %in% c("RNAseq")){
and_pp <- RNAseq_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
} else if (args$platform %in% c("AffyU219", "AffyHumanExon")){
and_pp <- Affy_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
}

# Visualise the expression of X-specific and Y-specific genes
xist <- and_pp[rownames(and_pp) == "ENSG00000229807", ]
y_genes <- c("ENSG00000234795", "ENSG00000237048", "ENSG00000275866", "ENSG00000239225", "ENSG00000169789", "ENSG00000129862", "ENSG00000273693", "ENSG00000243040", "ENSG00000198692",
"ENSG00000233699", "ENSG00000254488", "ENSG00000236424", "ENSG00000234414", "ENSG00000099715", "ENSG00000180910", "ENSG00000280969", "ENSG00000215560", "ENSG00000229236",
"ENSG00000223637", "ENSG00000131007", "ENSG00000176728", "ENSG00000012817", "ENSG00000182415", "ENSG00000099725", "ENSG00000184895", "ENSG00000129824", "ENSG00000067646",
"ENSG00000183878", "ENSG00000154620", "ENSG00000241859", "ENSG00000099721", "ENSG00000092377", "ENSG00000205944", "ENSG00000168757", "ENSG00000197038", "ENSG00000232808",
"ENSG00000233803", "ENSG00000229549", "ENSG00000242389", "ENSG00000165246", "ENSG00000215580", "ENSG00000131002")
y_genes <- and_pp[rownames(and_pp) %in% y_genes, ]

y_mean <- apply(y_genes, 2, mean)
y_genes <- data.frame(sample = colnames(y_genes), xist = xist, y_genes = y_mean)

geno_fam_f <- geno_fam[, c(2, 4), with = FALSE]
colnames(geno_fam_f) <- c("sample", "Sex")
geno_fam_f$Sex <- as.character(geno_fam_f$Sex)

y_genes <- merge(y_genes, geno_fam_f, by = "sample")
max_exp <- max(y_genes$y_genes, y_genes$xist)

y_genes$expressionSex <- case_when(
  y_genes$y_genes > y_genes$xist ~ 1,
  y_genes$y_genes < y_genes$xist ~ 2
)

y_genes$mismatch <- case_when(
  y_genes$Sex == 0 ~ "unknown",
  y_genes$expressionSex == y_genes$Sex ~ "no",
  y_genes$expressionSex != y_genes$Sex ~ "yes"
)

x_expression_median <- median(y_genes[y_genes$Sex == 1 & y_genes$expressionSex == 1, "y_genes"])
y_expression_median <- median(y_genes[y_genes$Sex == 2 & y_genes$expressionSex == 2, "xist"])

lower_slope <- tan((45 - args$contamination_area / 2) / 180*pi)
upper_slope <- tan((45 + args$contamination_area / 2) / 180*pi)

y_genes$contaminated <- case_when(
  (y_genes$y_genes > ((y_genes$xist - x_expression_median) * lower_slope + y_expression_median)
    & y_genes$y_genes < ((y_genes$xist - x_expression_median) * upper_slope + y_expression_median)) ~ "yes",
  TRUE ~ "no"
)

#
# y_genes$mismatch <- "no"
#
# y_genes$mismatch[y_genes$Sex == 0] <- "unknown"
#
# if (nrow(y_genes[(y_genes$y_genes > y_genes$xist & y_genes$Sex == 2) | (y_genes$y_genes < y_genes$xist & y_genes$Sex == 1), ]) > 0){
# y_genes[(y_genes$y_genes > y_genes$xist & y_genes$Sex == 2) | (y_genes$y_genes < y_genes$xist & y_genes$Sex == 1), ]$mismatch <- "yes"
# }

exclusion_zone <- tibble(x = c(x_expression_median, max_exp)) %>%
  mutate(lower_bound = (x - x_expression_median) * lower_slope + y_expression_median,
         upper_bound = (x - x_expression_median) * upper_slope + y_expression_median)

base_plot <- ggplot(data=exclusion_zone, aes(x = x, ymin = lower_bound, ymax = upper_bound)) +
  geom_ribbon(fill = "blue", alpha = 0.2) +
  geom_segment(aes(x = 0, y = 0, xend = max_exp, yend = max_exp), linetype = 2, colour = "blue") +
  geom_point(data = y_genes, inherit.aes = F, alpha = 0.3, aes(col = mismatch, shape = Sex, x = xist, y = y_genes)) +
  scale_colour_manual(values = c("no" = "black", "unknown" = "orange", "yes" = "red")) +
  theme_bw() + ylab("mean of Y genes") + xlab("XIST")

ggsave(paste0(args$output, "/exp_plots/SexSpecificGenes.png"), height = 5, width = 6, units = "in", dpi = 300, type = "cairo")
ggsave(paste0(args$output, "/exp_plots/SexSpecificGenes.pdf"), height = 5, width = 6, units = "in", dpi = 300)

# Filter out potential sex mismatches
and_pp <- and_pp[, colnames(and_pp) %in% y_genes[y_genes$mismatch == "no", ]$sample]

summary_table_temp <- data.table(Stage = "Samples after removal of sex errors", Nr_of_features = nrow(and_pp), Nr_of_samples = ncol(and_pp))
summary_table <- rbind(summary_table, summary_table_temp)

# Apply inverse normal transformation to normalised data.
#and_p <- apply(and_p, 1, Z_transform) # No Z-transform as data will be forced to normal distribution anyway
and_p <- apply(and_pp, 1, INT_transform)
and_p <- t(and_p)

pcs <- prcomp(t(and_p), center = FALSE, scale. = FALSE)
PCs <- data.table(Sample = rownames(pcs$x), pcs$x)

# Add info about sex
message("Outline sex on the plots")

PCs_int <- PCs
PCs_norm <- iterative_outliers$PCs_first

PCs_int <- merge(PCs_int, sex, by = "Sample")
PCs_norm <- merge(PCs_norm, sex, by.x = "sample", by.y = "Sample")

p5 <- ggplot(PCs_int, aes(x = PC1, y = PC2, shape = Sex)) + geom_point(alpha = 0.3) + theme_bw()
p6 <- ggplot(PCs_int, aes(x = PC3, y = PC4, shape = Sex)) + geom_point(alpha = 0.3) + theme_bw()
p7 <- ggplot(PCs_int, aes(x = PC5, y = PC6, shape = Sex)) + geom_point(alpha = 0.3) + theme_bw()
p8 <- ggplot(PCs_int, aes(x = PC7, y = PC8, shape = Sex)) + geom_point(alpha = 0.3) + theme_bw()

p1 <- ggplot(PCs_norm, aes(x = PC1, y = PC2, colour = outlier, shape = Sex)) +
      geom_point(alpha = 0.3) +
      theme_bw() +
      scale_colour_manual(values = c("no" = "black", "yes" = "red")) +
      ggtitle(paste0(1, ". iteration round"))
p2 <- ggplot(PCs_norm, aes(x = PC3, y = PC4, colour = outlier, shape = Sex)) +
      geom_point(alpha = 0.3) +
      theme_bw() +
      scale_colour_manual(values = c("no" = "black", "yes" = "red")) +
      ggtitle(paste0(1, ". iteration round"))
p3 <- ggplot(PCs_norm, aes(x = PC5, y = PC6, colour = outlier, shape = Sex)) +
      geom_point(alpha = 0.3) +
      theme_bw() +
      scale_colour_manual(values = c("no" = "black", "yes" = "red")) +
      ggtitle(paste0(1, ". iteration round"))
p4 <- ggplot(PCs_norm, aes(x = PC7, y = PC8, colour = outlier, shape = Sex)) +
      geom_point(alpha = 0.3) +
      theme_bw() +
      scale_colour_manual(values = c("no" = "black", "yes" = "red")) +
      ggtitle(paste0(1, ". iteration round"))

p <- (p1 | p2) / (p3 | p4)

ggsave(paste0(args$output, "/exp_plots/PCA_before.png"), height = 9, width = 10.5, units = "in", dpi = 300, type = "cairo")
ggsave(paste0(args$output, "/exp_plots/PCA_before.pdf"), height = 9, width = 10.5, units = "in", dpi = 300)

p <- (p5 | p6) / (p7 | p8)
ggsave(paste0(args$output, "/exp_plots/PCA_after.png"), height = 7.5, width = 9, units = "in", dpi = 300, type = "cairo")
ggsave(paste0(args$output, "/exp_plots/PCA_after.pdf"), height = 7.5, width = 9, units = "in", dpi = 300)

# Convert to HASE format
and_p2 <- as.data.table(t(and_p))
and_p2 <- data.table(`ID` = colnames(and_p), and_p2)

fwrite(and_p2, paste0(args$output, "/exp_data_QCd/exp_data_preprocessed.txt"), sep = "\t", quote = FALSE)

# Write out importance of final PCs
importance <- pcs$sdev^2 / sum(pcs$sdev^2)
summary_pcs <- data.table(PC = paste0("PC", 1:100), explained_variance = importance[1:100])
summary_pcs$PC <- factor(summary_pcs$PC, levels = paste0("PC", 1:100))

fwrite(PCs[, c(1:min(101, nrow(PCs))), with = F], paste0(args$output, "/exp_PCs/exp_PCs.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
fwrite(summary_pcs, paste0(args$output, "/exp_data_summary/", "summary_pcs.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

p <- ggplot(summary_pcs, aes(x = PC, y = explained_variance)) + geom_bar(stat = "identity") + theme_bw() +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(args$output, "/exp_plots/PCA_final_scree_plot.png"), height = 6, width = 17, units = "in", dpi = 300, type = "cairo")
ggsave(paste0(args$output, "/exp_plots/PCA_final_scree_plot.pdf"), height = 6, width = 17, units = "in", dpi = 300)

# Summary statistics ----
message("Calculating descriptive summary statistics for every gene...")
# Per gene, calculate mean, median, min, max, sd, and shapiro test P-value
# before preprocessing
message("Before normalisation.")
# Remove sample outliers
sample_non_outliers <- colnames(and_p)
sample_non_outliers <- gte[gte$V1 %in% sample_non_outliers, ]$V2

and <- and[, colnames(and) %in% c("Feature", sample_non_outliers), with = FALSE]

if (args$platform %in% c("HT12v3", "HT12v4", "HuRef8")){
and_pp <- illumina_array_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples, normalize = FALSE)
} else if (args$platform %in% c("RNAseq")){
and_pp <- RNAseq_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples, normalize = FALSE)
} else if (args$platform %in% c("AffyU219", "AffyHumanExon")){
and_pp <- Affy_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
}

gene_summary <- exp_summary(and_pp)

fwrite(gene_summary, paste0(args$output, "/exp_data_summary/", "raw_gene_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# on fully processed expression matrix
message("After normalization.")
if (args$platform %in% c("HT12v3", "HT12v4", "HuRef8")){
and_pp <- illumina_array_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples, normalize = TRUE)
} else if (args$platform %in% c("RNAseq")){
and_pp <- RNAseq_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples, normalize = TRUE)
} else if (args$platform %in% c("AffyU291", "AffyHuEx")){
and_pp <- Affy_preprocess(and, args$genotype_to_expression_linking, args$genotype_samples)
}

gene_summary <- exp_summary(and_pp)
fwrite(gene_summary, paste0(args$output, "/exp_data_summary/", "processed_gene_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# Write out summary table about features and samples
fwrite(summary_table, paste0(args$output, "/exp_data_summary/", "summary_table.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
