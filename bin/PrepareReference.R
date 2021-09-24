library(bigsnpr)
library(data.table)

# Remove the check of parallel blas
options(bigstatsr.check.parallel.blas = FALSE)

# This part follows the QC steps as specified in the bigsnpr manual: https://privefl.github.io/bigsnpr/articles/bedpca.html
# Download plink2 executable
plink2 <- download_plink2("data")

# Download subsetted 1000G reference
bedfile <- download_1000G("data")

# Find relatives
rel <- snp_plinkKINGQC(
  plink2.path = "tools/plink2",
  bedfile.in = bedfile,
  thr.king = 2^-4.5,
  make.bed = FALSE,
  ncores = 8
)

ref_bed <- bed("data/1000G_phase3_common_norel.bed")

# Filter out related sample pairs
ind.rel <- match(c(rel$IID1, rel$IID2), ref_bed$fam$sample.ID)
ind.norel <- rows_along(ref_bed)[-ind.rel]

# Calculate PCA
obj.svd <- bed_autoSVD(ref_bed, ind.row = ind.norel, k = 20, ncores = 8)

# Find outlier samples
prob <- bigutilsr::prob_dist(obj.svd$u, ncores = 8)
S <- prob$dist.self / sqrt(prob$dist.nn)

# Remove sample (S>0.5 as shown in manual) and redo PCA
ind.row <- ind.norel[S < 0.5]

fwrite(as.data.frame(ind.row), "data/unrelated_reference_samples_ids.txt", sep = "\t", quote = FALSE, row.names = FALSE)