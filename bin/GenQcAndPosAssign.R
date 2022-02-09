#!/usr/bin/env Rscript

library(bigsnpr)
library(dplyr)
library(ggplot2)
library(data.table)
library(optparse)
library(patchwork)
library(stringr)
library(rmarkdown)
library(Cairo)

# Argument parser
option_list <- list(
    make_option(c("-t", "--target_bed"), type = "character",
    help = "Name of the target genotype file (bed/bim/fam format). Required file extension: .bed."),
    make_option(c("-g", "--gen_exp"), type = "character",
    help = "Tab-delimited genotype-to-expression sample ID linking file."),
    make_option(c("-s", "--sample_list"), type = "character",
    help = "Path to the file listing unrelated samples for reference data (tab-delimited .txt)."),
    make_option(c("-p", "--pops"), type = "character",
    help = "Path to the file indicating the population for each sample in reference data."),
    make_option(c("-a", "--pruned_variants_sex_check"), type = "character",
    help = "Path to a file with pruned X-chromosome variants to use in the sex check"),
    make_option(c("-o", "--output"), type = "character", help = "Folder with all the output files."),
    make_option(c("-S", "--S_threshold"), default = 0.4,
    help = "Numeric threshold to declare samples outliers, based on the genotype PCs. Defaults to 0.4 but should always be visually checked and changed, if needed."),
    make_option(c("-d", "--SD_threshold"), default = 0.4,
    help = "Numeric threshold to declare samples outliers, based on the genotype PCs. Defaults to 0.4 but should always be visually checked and changed, if needed.")
    )

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

# Remove the check of parallel blas
options(bigstatsr.check.parallel.blas = FALSE)

# Report settings
print(args$target_bed)
print(args$gen_exp)
print(args$sample_list)
print(args$pops)
print(args$pruned_variants_sex_check)
print(args$output)
print(args$S_threshold)
print(args$SD_threshold)

bed_simplepath <- stringr::str_replace(args$target_bed, ".bed", "")

# Make output folder structure
dir.create(args$output)
dir.create(paste0(args$output, "/gen_plots"))
dir.create(paste0(args$output, "/gen_data_QCd"))
dir.create(paste0(args$output, "/gen_PCs"))
dir.create(paste0(args$output, "/gen_data_summary"))

# Download plink2 executable
download_plink2("plink", AVX2 = FALSE)
# Download plink 1.9 executable
download_plink("plink")

# Download subsetted 1000G reference
bedfile <- download_1000G("data")

# Target data
## Original file
message("Read in target data.")
target_bed <- bed(args$target_bed)
summary_table <- data.frame(stage = "Raw file", Nr_of_SNPs = target_bed$ncol, Nr_of_samples = target_bed$nrow)

## Keep in only samples which are present in genotype-to-expression file
gte <- fread(args$gen_exp, sep = "\t", header = FALSE)
samples_to_include_gte <- data.frame(FID = target_bed$.fam$family.ID, IID = target_bed$.fam$sample.ID)
samples_to_include_gte <- samples_to_include_gte[samples_to_include_gte$IID %in% gte$V1, ]

fwrite(samples_to_include_gte, "SamplesToIncludeGte.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
# Remove samples not in GTE
system(paste0("plink/plink2 --bfile ", bed_simplepath, " --output-chr 26 --keep SamplesToIncludeGte.txt --make-bed --threads 4 --out ", bed_simplepath, "_filtered"))

temp_QC <- data.frame(stage = "Samples in genotype-to-expression file", Nr_of_SNPs = target_bed$ncol, Nr_of_samples = nrow(samples_to_include_gte))
summary_table <- rbind(summary_table, temp_QC)

# Do SNP and sample missingness QC on raw genotype bed
message("Do SNP and genotype QC.")
snp_plinkQC(
  plink.path = "plink/plink2",
  prefix.in = paste0(bed_simplepath, "_filtered"),
  prefix.out = paste0(bed_simplepath, "_QC"),
  file.type = "--bfile",
  maf = 0.01,
  geno = 0.05,
  mind = 0.05,
  hwe = 1e-6,
  autosome.only = FALSE,
  extra.options = paste0("--output-chr 26 --not-chr 0 25-26 --set-all-var-ids ", r"(@:#[b37]\$r,\$a)"),
  verbose = TRUE
)

# Read in reference and target genotype data
ref_bed <- bed("data/1000G_phase3_common_norel.bed")
# Read in QCd target genotype data
target_bed <- bed(paste0(bed_simplepath, "_QC.bed"))
temp_QC <- data.frame(stage = "SNP CR>0.95; HWE P>1e-6; MAF>0.01; GENO<0.05; MIND<0.05", Nr_of_SNPs = target_bed$ncol, Nr_of_samples = target_bed$nrow)
summary_table <- rbind(summary_table, temp_QC)

## Assert that all IIDs are unique
if (any(duplicated(target_bed$fam$sample.ID))) {
  stop("Individual sample IDs should be unique. Exiting...")
}

sex_check_data_set_chromosomes <- unique(target_bed$map$chromosome)

sex_check_out_path <- paste0(args$output, "/gen_data_QCd/SexCheck.txt")
sex_check_removed_out_path <- paste0(args$output, "/gen_data_QCd/SexCheckRemoved.txt")
sex_check_samples <- target_bed$fam

if (23 %in% sex_check_data_set_chromosomes) {

  # Do sex check
  message("Do sex check.")

  # Split x if needed

  pruned_variants_sex_check <- args$pruned_variants_sex_check

  if (!is.null(pruned_variants_sex_check)
    && pruned_variants_sex_check != ""
    && file.exists(pruned_variants_sex_check)
    && nrow(fread(pruned_variants_sex_check, header = F)) > 0) {

    message("Using predefined pruned variants for sex-check:")
    message(pruned_variants_sex_check)

    system(paste0(
      "plink/plink --bfile ", bed_simplepath, "_QC", " --extract range ", pruned_variants_sex_check,
      " --maf 0.05 --make-bed --out ", bed_simplepath, "_split"))

  } else {

    message("Not using predefined pruned variants for sex-check")

    system(paste0(
      "plink/plink --bfile ", bed_simplepath, "_QC",
      " --chr X --maf 0.05 --split-x hg19 no-fail --make-bed --out ", bed_simplepath, "_split"))

  }

  ## Pruning
  system(paste0("plink/plink2 --bfile ", bed_simplepath, "_split",
                " --rm-dup 'exclude-mismatch' --indep-pairwise 20000 200 0.2 --out check_sex_x"))
  ## Sex check
  system(paste0("plink/plink --bfile ", bed_simplepath, "_split --extract check_sex_x.prune.in --check-sex"))

  #stop("check output")

  ## If there is sex info in the fam file for all samples then remove samples which fail the sex check or genotype-based F is >0.2 & < 0.8
  sexcheck <- fread("plink.sexcheck")
  ## Annotate samples who have clear sex
  sexcheck$F_PASS <- !(sexcheck$F > 0.2 & sexcheck$F < 0.8)
  temp_QC <- data.frame(stage = "Sex check (0.2<F<0.8)",
                        Nr_of_SNPs = target_bed$ncol,
                        Nr_of_samples = sum(sexcheck$F_PASS))
  summary_table <- rbind(summary_table, temp_QC)

  sexcheck$MATCH_PASS <- case_when(sexcheck$PEDSEX == 0 ~ T,
                                   sexcheck$STATUS == "PROBLEM" ~ F,
                                   TRUE ~ T)

  sexcheck$PASS <- sexcheck$MATCH_PASS & sexcheck$F_PASS

  if (any(sexcheck$PEDSEX %in% c(1, 2))) {

    temp_QC <- data.frame(stage = "Sex check (reported and genetic sex mismatch)",
                          Nr_of_SNPs = target_bed$ncol,
                          Nr_of_samples = sum(sexcheck$PASS))
    summary_table <- rbind(summary_table, temp_QC)

  } else {
    message("No sex info in the .fam file.")
  }

  print(sexcheck$PEDSEX)

  sex_cols <- c("0" = "black", "1" = "orange", "2" = "blue")

  p <- ggplot(sexcheck, aes(x = F, fill = factor(PEDSEX))) +
    geom_histogram(position="stack", color = "black", alpha = 0.5) +
    scale_fill_manual(values = sex_cols, breaks = c("0", "1", "2"), labels = c("Unknown", "Male", "Female"), name = "Reported sex") +
    geom_vline(xintercept = c(0.2, 0.8), colour = "red", linetype = 2) + theme_bw()

  ggsave(paste0(args$output, "/gen_plots/SexCheck.png"), p, type = "cairo", height = 7 / 2, width = 9, units = "in", dpi = 300)
  ggsave(paste0(args$output, "/gen_plots/SexCheck.pdf"), p, height = 7 / 2, width = 9, units = "in", dpi = 300)

  fwrite(sexcheck, sex_check_out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  fwrite(sexcheck[!sexcheck$PASS,], sex_check_removed_out_path, sep = "\t", quote = FALSE, row.names = FALSE)

} else {
  warning("No X chromosome present. Skipping sex-check...")

  sexcheck <- sex_check_samples[,c(1,2,5)]
  colnames(sexcheck) <- c("FID", "IID", "PEDSEX")
  sexcheck$PEDSEX_COPY <- sexcheck$PEDSEX
  sexcheck$STATUS <- NA_character_
  sexcheck$F <- NA_real_

  fwrite(sexcheck, sex_check_out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  file.create(sex_check_removed_out_path)
}

# Remove sex chromosomes
snp_plinkQC(
  plink.path = "plink/plink2",
  prefix.in = paste0(bed_simplepath, "_QC"),
  prefix.out = paste0(bed_simplepath, "_QC", "_QC"),
  file.type = "--bfile",
  maf = 0.01,
  geno = 0.05,
  mind = 0.05,
  hwe = 1e-6,
  autosome.only = TRUE,
  extra.options = paste0("--output-chr 26 --remove ", sex_check_removed_out_path),
  verbose = TRUE
)

# replace the previous QC version
system(paste0("rm ", bed_simplepath, "_QC.*"))
system(paste0("mv ", bed_simplepath, "_QC_QC.bed ", bed_simplepath, "_QC.bed"))
system(paste0("mv ", bed_simplepath, "_QC_QC.bim ", bed_simplepath, "_QC.bim"))
system(paste0("mv ", bed_simplepath, "_QC_QC.fam ", bed_simplepath, "_QC.fam"))

# Read in again QCd target genotype data
target_bed <- bed(paste0(bed_simplepath, "_QC.bed"))
temp_QC <- data.frame(stage = "Removed X/Y", Nr_of_SNPs = target_bed$ncol, Nr_of_samples = nrow(target_bed$fam))
summary_table <- rbind(summary_table, temp_QC)

# Do heterozygosity check
message("Do heterozygosity check.")

# Prune variants
system(paste0("plink/plink2 --bfile ", bed_simplepath, "_QC --rm-dup 'exclude-mismatch' --indep-pairwise 50 1 0.2"))

system(paste0("plink/plink2 --bfile ", bed_simplepath, "_QC --extract plink2.prune.in --het"))
het <- fread("plink2.het")
het$het_rate <- (het$OBS_CT - het$`O(HOM)`) / het$OBS_CT

het_fail_samples <- het[het$het_rate < mean(het$het_rate) - 3 * sd(het$het_rate) | het$het_rate > mean(het$het_rate) + 3 * sd(het$het_rate), ]

# Get the indices of those samples that passed heterozygozity check
indices_of_het_failed_samples <- match(het_fail_samples, target_bed$fam$sample.ID)
indices_of_het_passed_samples <- rows_along(target_bed)[-indices_of_het_failed_samples]

temp_QC <- data.frame(stage = "Excess heterozygosity (mean+/-3SD)", Nr_of_SNPs = target_bed$ncol, Nr_of_samples = length(indices_of_het_passed_samples))
summary_table <- rbind(summary_table, temp_QC)
temp_QC <- data.frame(stage = "Excess heterozygosity (mean+/-3SD)", Nr_of_SNPs = target_bed$ncol, Nr_of_samples = nrow(het) - nrow(het_fail_samples))
summary_table <- rbind(summary_table, temp_QC)

p <- ggplot(het, aes(x = het_rate)) + geom_histogram(color = "#000000", fill = "#000000", alpha = 0.5) + 
xlab("Heterozygosity rate") + 
geom_vline(xintercept = c(mean(het$het_rate), mean(het$het_rate) + 3 * sd(het$het_rate), mean(het$het_rate) - 3 * sd(het$het_rate)), linetype = 2, colour = "red") + 
theme_bw()

ggsave(paste0(args$output, "/gen_plots/HetCheck.png"), type = "cairo", height = 7 / 2, width = 9, units = "in", dpi = 300)
ggsave(paste0(args$output, "/gen_plots/HetCheck.pdf"), height = 7 / 2, width = 9, units = "in", dpi = 300)

# Project the data on QCd 1000G reference
message("Projecting samples to 1000G reference.")
unrelated_ref_samples <- fread(args$sample_list)
unrelated_ref_samples <- as.numeric(unrelated_ref_samples$ind.row)

proj_PCA <- bed_projectPCA(
  obj.bed.ref = ref_bed,
  ind.row.ref = unrelated_ref_samples,
  obj.bed.new = target_bed,
  k = 10,
  strand_flip = TRUE,
  join_by_pos = TRUE,
  match.min.prop = 0.5,
  build.new = "hg19",
  build.ref = "hg19",
  liftOver = NULL,
  verbose = TRUE,
  ncores = 4
)

## Visualise PCs
abi <- as.data.frame(proj_PCA$OADP_proj)
colnames(abi) <- paste0("PC", 1:10)

PCs_ref <- predict(proj_PCA$obj.svd.ref)
abi2 <- as.data.frame(PCs_ref)
colnames(abi2) <- paste0("PC", 1:10)

abi2$sample <- ref_bed$.fam$sample.ID[unrelated_ref_samples]
abi2 <- abi2[, c(11, 1:10)]

pops <- fread(args$pops)
pops <- pops[, c(2, 6, 7)]
abi2 <- merge(abi2, pops, by.x = "sample", by.y = "SampleID")
abi2 <- abi2[, c(1, 12, 13, 2:11)]

abi <- data.frame(sample = target_bed$fam$sample.ID, Population = "Target", Superpopulation = "Target", abi)

abi$type <- "Target"
abi2$type <- "1000G"

combined <- rbind(abi, abi2)

combined$Superpopulation <- factor(combined$Superpopulation, levels = c("Target", "EUR", "EAS", "AMR", "SAS", "AFR"))

p00 <- ggplot(combined, aes(x = PC1, y = PC2, alpha = type)) + 
geom_point() + 
theme_bw() + 
scale_alpha_manual(values = c("Target" = 1, "1000G" = 0)) + 
ggtitle("Target sample projections\nin 1000G PC space")

combined_h <- combined[combined$Superpopulation == "Target", ]
#combined_h$Superpopulation <- factor(combined_h$Superpopulation, levels = "Target")
#combined_h$Type <- factor(combined_h$Type, levels = "Target")

p0 <- ggplot(combined_h, aes(x = PC1, y = PC2)) + 
geom_point() + 
theme_bw() + 
ggtitle("Target sample projections\nzoomed in")

p1 <- ggplot(combined, aes(x = PC1, y = PC2, colour = Superpopulation, alpha = type)) + 
geom_point() + 
theme_bw() + 
scale_color_manual(values = c("Target" = "black", "EUR" = "blue", 
"EAS" = "goldenrod", "AMR" = "lightgrey", "SAS" = "orange", "AFR" = "red")) +
scale_alpha_manual(values = c("Target" = 1, "1000G" = 0.2))

p2 <- ggplot(combined, aes(x = PC3, y = PC4, colour = Superpopulation, alpha = type)) + 
geom_point() + theme_bw() + 
scale_color_manual(values = c("Target" = "black", "EUR" = "blue", 
"EAS" = "goldenrod", "AMR" = "lightgrey", "SAS" = "orange", "AFR" = "red")) +
scale_alpha_manual(values = c("Target" = 1, "1000G" = 0.2))

p3 <- ggplot(combined, aes(x = PC5, y = PC6, colour = Superpopulation, alpha = type)) + 
geom_point() + theme_bw() + 
scale_color_manual(values = c("Target" = "black", "EUR" = "blue", 
"EAS" = "goldenrod", "AMR" = "lightgrey", "SAS" = "orange", "AFR" = "red")) +
scale_alpha_manual(values = c("Target" = 1, "1000G" = 0.2))

p4 <- ggplot(combined, aes(x = PC7, y = PC8, colour = Superpopulation, alpha = type)) + 
geom_point() + theme_bw() + 
scale_color_manual(values = c("Target" = "black", "EUR" = "blue", 
"EAS" = "goldenrod", "AMR" = "lightgrey", "SAS" = "orange", "AFR" = "red")) +
scale_alpha_manual(values = c("Target" = 1, "1000G" = 0.2))

p5 <- ggplot(combined, aes(x = PC9, y = PC10, colour = Superpopulation, alpha = type)) + 
geom_point() + theme_bw() + 
scale_color_manual(values = c("Target" = "black", "EUR" = "blue", 
"EAS" = "goldenrod", "AMR" = "lightgrey", "SAS" = "orange", "AFR" = "red")) +
scale_alpha_manual(values = c("Target" = 1, "1000G" = 0.2))

p <- p00 + p0 + p1 + p2 + p3 + p4 + p5 + plot_layout(nrow = 4)

ggsave(paste0(args$output, "/gen_plots/SamplesPCsProjectedTo1000G.png"), type = "cairo", height = 20, width = 9.5 * 1.6, units = "in", dpi = 300)
ggsave(paste0(args$output, "/gen_plots/SamplesPCsProjectedTo1000G.pdf"), height = 20, width = 9.5 * 1.6, units = "in", dpi = 300)
fwrite(abi[, -c(2, 3, ncol(abi))], paste0(args$output, "/gen_data_summary/1000G_PC_projections.txt"), sep = "\t", quote = FALSE )

## Assign each sample to the superpopulation
message("Assign each sample to 1000G superpopulation.")
### Calculate distance of each sample to all samples per each population
target_samples <- abi[, -c(2, 3, ncol(abi))]

#### Use 3 PCs
target_samples <- target_samples[, c(1:4)]
rownames(target_samples) <- target_samples$sample
target_samples <- target_samples[, -1]

population_assign_res <- data.frame(sample = rownames(target_samples), abi = rownames(target_samples))

#### EUR
for(population in c("EUR", "EAS", "AMR", "SAS", "AFR")){
abi_e <- abi2[abi2$Superpopulation == population, ]
head(abi_e)

sup_pop_samples <- abi_e[, -c(2, 3, ncol(abi_e))]

sup_pop_samples <- sup_pop_samples[, c(1:4)]
rownames(sup_pop_samples) <- sup_pop_samples$sample
sup_pop_samples <- sup_pop_samples[, -1]

head(sup_pop_samples)

comb <- rbind(target_samples, sup_pop_samples)
head(comb)
distance <- as.matrix(dist(comb, method = "euclidean"))
head(distance)
distance <- distance[c(1:nrow(target_samples)), -c(1:nrow(target_samples))]
head(distance)

head(rowMeans(distance))

distance <- data.frame(sample = rownames(target_samples), MeanDistance = rowMeans(distance))
colnames(distance)[2] <- population

population_assign_res <- cbind(population_assign_res, distance[, -1])

print(paste("distance:", population))

}

colnames(population_assign_res)[3:ncol(population_assign_res)] <- c("EUR", "EAS", "AMR", "SAS", "AFR")
fwrite(population_assign_res[, -1], paste0(args$output, "/gen_data_summary/PopAssignResults.txt"), sep = "\t", quote = FALSE )

# Find related samples
message("Find related samples.")
related <- snp_plinkKINGQC(
  plink2.path = "plink/plink2",
  bedfile.in = paste0(bed_simplepath, "_QC.bed"),
  thr.king = 2^-4.5,
  make.bed = FALSE,
  ncores = 4
)

### Do PCA on target data
message("Find genetic outliers.")

print(related)
### First remove one related sample from each pair and those failing heterozygosity check
# TODO, check if this is correct that only one out of two is removed
ind.rel <- match(unique(c(related$IID2, het_fail_samples$IID)), target_bed$fam$sample.ID)

print(unique(c(related$IID2, het_fail_samples$IID)))

ind.norel <- rows_along(target_bed)[-ind.rel]

temp_QC <- data.frame(stage = "Relatedness: thr. KING>2^-4.5", Nr_of_SNPs = target_bed$ncol, Nr_of_samples = length(ind.norel))
summary_table <- rbind(summary_table, temp_QC)

message("Find genetic outliers: do PCA on QCd target data.")
### PCA
target_pca <- bed_autoSVD(target_bed, ind.row = ind.norel, k = 10, ncores = 4)

### Find outlier samples
prob <- bigutilsr::prob_dist(target_pca$u, ncores = 4)
S <- prob$dist.self / sqrt(prob$dist.nn)

# Put threshold for outlier samples, this is by default 0.4!
Sthresh <- args$S_threshold

p <- ggplot() +
  geom_histogram(aes(S), color = "#000000", fill = "#000000", alpha = 0.5) +
  scale_x_continuous(breaks = 0:5 / 5, limits = c(0, NA)) +
  scale_y_sqrt(breaks = c(10, 100, 500)) +
  theme_bigstatsr() +
  labs(x = "Statistic of outlierness", y = "Frequency (sqrt-scale)") +
  geom_vline(aes(xintercept = Sthresh), colour = "red", linetype = 2)

ggsave(paste0(args$output, "/gen_plots/PC_dist_outliers_S.png"), type = "cairo", height = 7 / 2, width = 9, units = "in", dpi = 300)
ggsave(paste0(args$output, "/gen_plots/PC_dist_outliers_S.pdf"), height = 7 / 2, width = 9, units = "in", dpi = 300)

# Visualise PCs, outline individual outlier samples
PCs <- predict(target_pca)

PCs <- as.data.frame(PCs)
colnames(PCs) <- paste0("PC", 1:10)
PCs$S <- S

PCs$outlier_ind <- "no"

if (any(PCs$S > Sthresh)) {
  PCs[PCs$S > Sthresh, ]$outlier_ind <- "yes"
}
PCs$sd_outlier <- "no"
sd_outlier_selection <- (PCs$PC1 > mean(PCs$PC1) + args$SD_threshold * sd(PCs$PC1) | PCs$PC1 < mean(PCs$PC1) - args$SD_threshold * sd(PCs$PC1)) | (PCs$PC2 > mean(PCs$PC2) + args$SD_threshold * sd(PCs$PC2) | PCs$PC2 < mean(PCs$PC2) - args$SD_threshold * sd(PCs$PC2))
if (any(sd_outlier_selection)) {
  PCs[sd_outlier_selection, ]$sd_outlier <- "yes"
}

PCs$outlier <- "no"
if (nrow(PCs[PCs$outlier_ind == "yes" & PCs$sd_outlier == "no", ]) > 0){
PCs[PCs$outlier_ind == "yes" & PCs$sd_outlier == "no", ]$outlier <- "S outlier"}
if(nrow(PCs[PCs$outlier_ind == "no" & PCs$sd_outlier == "yes", ]) > 0){
PCs[PCs$outlier_ind == "no" & PCs$sd_outlier == "yes", ]$outlier <- "SD outlier"}
if(nrow(PCs[PCs$outlier_ind == "yes" & PCs$sd_outlier == "yes", ]) > 0){
PCs[PCs$outlier_ind == "yes" & PCs$sd_outlier == "yes", ]$outlier <- "S and SD outlier"
}
# For first 2 PCs also remove samples which deviate from the mean

p1 <- ggplot(PCs, aes(x = PC1, y = PC2, colour = outlier)) + theme_bw() + geom_point(alpha = 0.5) + scale_color_manual(values = c("no" = "black", "SD outlier" = "#d79393", "S outlier" = "red", "S and SD outlier" = "firebrick")) +
geom_vline(xintercept = c(mean(PCs$PC1) + 3 * sd(PCs$PC1), mean(PCs$PC1) - 3 * sd(PCs$PC1)), colour = "firebrick", linetype = 2) + 
geom_hline(yintercept = c(mean(PCs$PC2) + 3 * sd(PCs$PC2), mean(PCs$PC2) - 3 * sd(PCs$PC2)), colour = "firebrick", linetype = 2)
p2 <- ggplot(PCs, aes(x = PC3, y = PC4, colour = outlier)) + theme_bw() + geom_point(alpha = 0.5) + scale_color_manual(values = c("no" = "black", "SD outlier" = "#d79393", "S outlier" = "red", "S and SD outlier" = "firebrick"))
p3 <- ggplot(PCs, aes(x = PC5, y = PC6, colour = outlier)) + theme_bw() + geom_point(alpha = 0.5) + scale_color_manual(values = c("no" = "black", "SD outlier" = "#d79393", "S outlier" = "red", "S and SD outlier" = "firebrick"))
p4 <- ggplot(PCs, aes(x = PC7, y = PC8, colour = outlier)) + theme_bw() + geom_point(alpha = 0.5) + scale_color_manual(values = c("no" = "black", "SD outlier" = "#d79393", "S outlier" = "red", "S and SD outlier" = "firebrick"))
p5 <- ggplot(PCs, aes(x = PC9, y = PC10, colour = outlier)) + theme_bw() + geom_point(alpha = 0.5) + scale_color_manual(values = c("no" = "black", "SD outlier" = "#d79393", "S outlier" = "red", "S and SD outlier" = "firebrick"))

p <- p1 + p2 + p3 + p4 + p5 + plot_layout(nrow = 3)

ggsave(paste0(args$output, "/gen_plots/PCA_outliers.png"), type = "cairo", height = 10 * 1.5, width = 9 * 1.5, units = "in", dpi = 300)
ggsave(paste0(args$output, "/gen_plots/PCA_outliers.pdf"), height = 10 * 1.5, width = 9 * 1.5, units = "in", dpi = 300)

# Filter out related samples and outlier samples, write out QCd data
message("Filter out related samples and outlier samples, write out QCd data.")
ind.row <- ind.norel[PCs$outlier == "no"]
samples_to_include <- data.frame(family.ID = target_bed$.fam$family.ID[ind.row], sample.IDD2 = target_bed$.fam$sample.ID[ind.row])

temp_QC <- data.frame(stage = paste0("Outlier samples: thr. S>", Sthresh, " PC1/PC2 SD deviation thresh ", args$SD_threshold), Nr_of_SNPs = target_bed$ncol, Nr_of_samples = nrow(samples_to_include))
summary_table <- rbind(summary_table, temp_QC)

fwrite(data.table::data.table(samples_to_include), "SamplesToInclude.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
# Remove samples
system(paste0("plink/plink2 --bfile ", bed_simplepath, "_QC --output-chr 26 --keep SamplesToInclude.txt --make-bed --threads 4 --out ", args$output, "/gen_data_QCd/", bed_simplepath, "_ToImputation"))

# Rerun PCA on QCd data
bed_qc <- bed(paste0(args$output, "/gen_data_QCd/", bed_simplepath, "_ToImputation.bed"))
target_pca_qcd <- bed_autoSVD(bed_qc, k = 10, ncores = 4)

PCsQ <- predict(target_pca_qcd)

PCsQ <- as.data.frame(PCsQ)

colnames(PCsQ) <- paste0("PC", 1:10)
rownames(PCsQ) <- bed_qc$.fam$sample.ID

# Visualise
p1 <- ggplot(PCsQ, aes(x = PC1, y = PC2)) + theme_bw() + geom_point(alpha = 0.5)
p2 <- ggplot(PCsQ, aes(x = PC3, y = PC4)) + theme_bw() + geom_point(alpha = 0.5)
p3 <- ggplot(PCsQ, aes(x = PC5, y = PC6)) + theme_bw() + geom_point(alpha = 0.5)
p4 <- ggplot(PCsQ, aes(x = PC7, y = PC8)) + theme_bw() + geom_point(alpha = 0.5)
p5 <- ggplot(PCsQ, aes(x = PC9, y = PC10)) + theme_bw() + geom_point(alpha = 0.5)
p <- p1 + p2 + p3 + p4 + p5 + plot_layout(nrow = 3)

ggsave(paste0(args$output, "/gen_plots/Target_PCs_postQC.png"), type = "cairo", height = 10 * 1.5, width = 9 * 1.3, units = "in", dpi = 300)
ggsave(paste0(args$output, "/gen_plots/Target_PCs_postQC.pdf"), height = 10 * 1.5, width = 9 * 1.3, units = "in", dpi = 300)

# Write out
fwrite(PCsQ, paste0(args$output, "/gen_PCs/GenotypePCs.txt"), row.names = TRUE, sep = "\t", quote = FALSE)

# Visualise loadings
plot(target_pca_qcd, type = "loadings", loadings = 1:10, coeff = 0.6)
ggsave(paste0(args$output, "/gen_plots/Target_PCs_postQC_Loadings.png"), type = "cairo", height = (5 * 7) * 0.7, width = (5 * 7) * 0.7, units = "in", dpi = 300)

# Write out final summary
fwrite(summary_table, paste0(args$output, "/gen_data_summary/summary_table.txt"), sep = "\t", quote = FALSE)

# Reorder the samples and write out the sample file
message("Shuffle sample order.")
rows <- sample(nrow(samples_to_include))
samples_to_include2 <- samples_to_include[rows, ]
fwrite(samples_to_include2, "ShuffledSampleOrder.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

system(paste0("plink/plink2 -bfile ", args$output, "/gen_data_QCd/", bed_simplepath, "_ToImputation ",
"--indiv-sort f ShuffledSampleOrder.txt ",
"--make-bed ",
"--out ", args$output, "/gen_data_QCd/", bed_simplepath, "_ToImputation_temp"))

system(paste0("mv ", args$output, "/gen_data_QCd/", bed_simplepath, "_ToImputation_temp.bed ",
args$output, "/gen_data_QCd/", bed_simplepath, "_ToImputation.bed"))
system(paste0("mv ", args$output, "/gen_data_QCd/", bed_simplepath, "_ToImputation_temp.bim ",
args$output, "/gen_data_QCd/", bed_simplepath, "_ToImputation.bim"))
system(paste0("mv ", args$output, "/gen_data_QCd/", bed_simplepath, "_ToImputation_temp.fam ",
args$output, "/gen_data_QCd/", bed_simplepath, "_ToImputation.fam"))
