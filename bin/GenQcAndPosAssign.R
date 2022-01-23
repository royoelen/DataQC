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
    make_option(c("-s", "--sample_list"), type = "character", 
    help = "Path to the file listing unrelated samples for reference data (tab-delimited .txt)."),
    make_option(c("-p", "--pops"), type = "character", 
    help = "Path to the file inidicating the population for each sample in reference data."),
    make_option(c("-o", "--output"), type = "character", help = "Folder with all the output files."),
    make_option(c("-S", "--S_threshold"), default = 0.4, 
    help = "Numeric threshold to declare samples outliers, based on the genotype PCs. Defaults to 0.4 but should always be visually checked and changed, if needed.")
    )

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)


# Debug
print(args$target_bed)
print(args$sample_list)
print(args$pops)
print(args$output)
print(args$S_threshold)

bed_simplepath <- stringr::str_replace(args$target_bed, ".bed", "")

# Make output folder structure
dir.create(args$output)
dir.create(paste0(args$output, "/gen_plots"))
dir.create(paste0(args$output, "/gen_data_QCd"))
dir.create(paste0(args$output, "/gen_PCs"))
dir.create(paste0(args$output, "/gen_data_summary"))

# Remove the check of parallel blas
options(bigstatsr.check.parallel.blas = FALSE)

# Download plink2 executable
plink2 <- download_plink2("plink", AVX2 = FALSE)

# Download subsetted 1000G reference
bedfile <- download_1000G("data")

# Target data
## Original file
target_bed <- bed(args$target_bed)
summary_table <- data.frame(stage = "Raw file", Nr_of_SNPs = target_bed$ncol, Nr_of_samples = target_bed$nrow)

## Do SNP and sample QC on raw genotype bed
### SNP QC
#snp_plinkQC(
#  plink.path = "plink/plink2",
#  prefix.in = bed_simplepath,
#  file.type = "--bfile",
#  maf = 0.01,
#  geno = 0.05,
#  mind = 0.05,
#  hwe = 1e-6,
#  autosome.only = TRUE,
#  extra.options = "",
#  verbose = TRUE
#)

### Find related samples
related <- snp_plinkKINGQC(
  plink2.path = "plink/plink2",
  bedfile.in = paste0(bed_simplepath, "_QC.bed"),
  thr.king = 2^-4.5,
  make.bed = FALSE,
  ncores = 4
)

# Read in reference and target genotype data
ref_bed <- bed("data/1000G_phase3_common_norel.bed")
target_bed <- bed(paste0(bed_simplepath, "_QC.bed"))
temp_QC <- data.frame(stage = "SNP CR>0.95; HWE P>1e-6; MAF>0.01; MIND<0.05", Nr_of_SNPs = target_bed$ncol, Nr_of_samples = target_bed$nrow)
summary_table <- rbind(summary_table, temp_QC)

### Do PCA on target data
### First remove one related sample from each pair
ind.rel <- match(c(related$IID2), target_bed$fam$sample.ID)
ind.norel <- rows_along(target_bed)[-ind.rel]

temp_QC <- data.frame(stage = "Relatedness: thr. KING>2^-4.5", Nr_of_SNPs = target_bed$ncol, Nr_of_samples = length(ind.norel))
summary_table <- rbind(summary_table, temp_QC)

### PCA
target_pca <- bed_autoSVD(target_bed, ind.row = ind.norel, k = 10, ncores = 4)

### Find outlier samples
prob <- bigutilsr::prob_dist(target_pca$u, ncores = 4)
S <- prob$dist.self / sqrt(prob$dist.nn)

print("!!!!!!!!!!S:")
print(S)

p1 <- ggplot() +
  geom_histogram(aes(S), color = "#000000", fill = "#000000", alpha = 0.5) +
  scale_x_continuous(breaks = 0:5 / 5, limits = c(0, NA)) +
  scale_y_sqrt(breaks = c(10, 100, 500)) +
  theme_bigstatsr() +
  labs(x = "Statistic of outlierness", y = "Frequency (sqrt-scale)") 

# Put threshold for outlier samples, this is by default 0.4!
Sthresh <- args$S_threshold

print("test123")

p2 <- ggplot() +
  geom_histogram(aes(S), color = "#000000", fill = "#000000", alpha = 0.5) +
  scale_x_continuous(breaks = 0:5 / 5, limits = c(0, NA)) +
  scale_y_sqrt(breaks = c(10, 100, 500)) +
  theme_bigstatsr() +
  labs(x = "Statistic of outlierness", y = "Frequency (sqrt-scale)") +
  geom_vline(aes(xintercept = Sthresh), colour = "red", linetype = 2)

p <- p2
ggsave(paste0(args$output, "/gen_plots/PC_dist_outliers_S.png"), type = "cairo", height = 7 / 2, width = 9, units = "in", dpi = 300)

print("test123")

# Visualise PCs, outline outlier samples
PCs <- predict(target_pca)

print(str(PCs))

PCs <- as.data.frame(PCs)
colnames(PCs) <- paste0("PC", 1:10)
PCs$S <- S

print("test345")

print(str(PCs))

PCs$outlier <- "no"
if (any(PCs$S > Sthresh)) {
  PCs[PCs$S > Sthresh, ]$outlier <- "yes"
}

print(str(PCs))

plotLabels <- F
uniqueShapes <- intToUtf8(c(97:122, 65:90), multiple = T)
PCs$thisLabel <- ""

if (sum(PCs$outlier == "yes", na.rm = T) <= length(uniqueShapes) && plotLabels) {
  PCs[PCs$outlier == "yes", ]$thisLabel <- as.character(uniqueShapes[1:sum(PCs$outlier == "yes", na.rm = T)])
}

print(PCs[PCs$outlier == "yes", ])

p1 <- ggplot(PCs, aes(x = PC1, y = PC2, colour = outlier, label = thisLabel)) + theme_bw() + geom_point() + geom_text(hjust = 0, nudge_x = 1, size = 4, show.legend = F) + scale_color_manual(values = c("no" = "black", "yes" = "red"))
p2 <- ggplot(PCs, aes(x = PC3, y = PC4, colour = outlier, label = thisLabel)) + theme_bw() + geom_point() + geom_text(hjust = 0, nudge_x = 1, size = 4, show.legend = F) + scale_color_manual(values = c("no" = "black", "yes" = "red"))
p3 <- ggplot(PCs, aes(x = PC5, y = PC6, colour = outlier, label = thisLabel)) + theme_bw() + geom_point() + geom_text(hjust = 0, nudge_x = 1, size = 4, show.legend = F) + scale_color_manual(values = c("no" = "black", "yes" = "red"))
p4 <- ggplot(PCs, aes(x = PC7, y = PC8, colour = outlier, label = thisLabel)) + theme_bw() + geom_point() + geom_text(hjust = 0, nudge_x = 1, size = 4, show.legend = F) + scale_color_manual(values = c("no" = "black", "yes" = "red"))
p5 <- ggplot(PCs, aes(x = PC9, y = PC10, colour = outlier, label = thisLabel)) + theme_bw() + geom_point() + geom_text(hjust = 0, nudge_x = 1, size = 4, show.legend = F) + scale_color_manual(values = c("no" = "black", "yes" = "red"))

p <- p1 + p2 + p3 + p4 + p5 + plot_layout(nrow = 3)

ggsave(paste0(args$output, "/gen_plots/PCA_outliers.png"), type = "cairo", height = 10 * 1.5, width = 9 * 1.5, units = "in", dpi = 300)

# Filter out related samples and outlier samples, write out QCd data
ind.row <- ind.norel[S < Sthresh]
samples_to_include <- data.frame(family.ID = target_bed$.fam$family.ID[-ind.row], sample.IDD2 = target_bed$.fam$sample.ID[-ind.row])

temp_QC <- data.frame(stage = paste0("Outlier samples: thr. S>2^", Sthresh), Nr_of_SNPs = target_bed$ncol, Nr_of_samples = target_bed$nrow - nrow(samples_to_include))
summary_table <- rbind(summary_table, temp_QC)

snp_plinkRmSamples(plink.path = "plink/plink2",
bedfile.in = paste0(bed_simplepath, "_QC.bed"),
bedfile.out = paste0(args$output, "/gen_data_QCd/", bed_simplepath, "_ToImputation.bed"),
df.or.files = samples_to_include)

# Rerun PCA on QCd data
bed_qc <- bed(paste0(args$output, "/gen_data_QCd/", bed_simplepath, "_ToImputation.bed"))
target_pca_qcd <- bed_autoSVD(bed_qc, k = 10, ncores = 4)

PCsQ <- predict(target_pca_qcd)

PCsQ <- as.data.frame(PCsQ)

colnames(PCsQ) <- paste0("PC", 1:10)
rownames(PCsQ) <- bed_qc$.fam$sample.ID

# Visualise
p1 <- ggplot(PCsQ, aes(x = PC1, y = PC2)) + theme_bw() + geom_point()
p2 <- ggplot(PCsQ, aes(x = PC3, y = PC4)) + theme_bw() + geom_point()
p3 <- ggplot(PCsQ, aes(x = PC5, y = PC6)) + theme_bw() + geom_point()
p4 <- ggplot(PCsQ, aes(x = PC7, y = PC8)) + theme_bw() + geom_point()
p5 <- ggplot(PCsQ, aes(x = PC9, y = PC10)) + theme_bw() + geom_point()
p <- p1 + p2 + p3 + p4 + p5 + plot_layout(nrow = 3)

ggsave(paste0(args$output, "/gen_plots/Target_PCs_postQC.png"), type = "cairo", height = 10 * 1.5, width = 9 * 1.3, units = "in", dpi = 300)

# Write out
fwrite(PCsQ, paste0(args$output, "/gen_PCs/GenotypePCs.txt"), row.names = TRUE, sep = "\t", quote = FALSE)

# Visualise loadings
plot(target_pca_qcd, type = "loadings", loadings = 1:10, coeff = 0.6)

ggsave(paste0(args$output, "/gen_plots/Target_PCs_postQC_Loadings.png"), type = "cairo", height = (5 * 7) * 0.7, width = (5 * 7) * 0.7, units = "in", dpi = 300)

# Project the data on QCd 1000G reference
unrelated_ref_samples <- fread(args$sample_list)
unrelated_ref_samples <- as.numeric(unrelated_ref_samples$ind.row)

proj_PCA <- bed_projectPCA(
  obj.bed.ref = ref_bed,
  ind.row.ref = unrelated_ref_samples,
  obj.bed.new = bed_qc,
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

abi <- data.frame(sample = bed_qc$fam$sample.ID, Population = "Target", Superpopulation = "Target", abi)

abi$type <- "Target"
abi2$type <- "1000G"

combined <- rbind(abi, abi2)

combined$Superpopulation <- factor(combined$Superpopulation, levels = c("Target", "EUR", "EAS", "AMR", "SAS", "AFR"))

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

p <- p1 + p2 + p3 + p4 + p5 + plot_layout(nrow = 3)

ggsave(paste0(args$output, "/gen_plots/SamplesPCsProjectedTo1000G.png"), type = "cairo", height = 11 * 1.5, width = 9 * 1.5, units = "in", dpi = 300)

fwrite(summary_table, paste0(args$output, "/gen_data_summary/summary_table.txt"), sep = "\t", quote = FALSE)
fwrite(abi[, -c(2, 3, ncol(abi))], paste0(args$output, "/gen_data_summary/1000G_PC_projections.txt"), sep = "\t", quote = FALSE )

## Assign each sample to the superpopulation
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
