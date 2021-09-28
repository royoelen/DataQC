library(data.table)

gen_cov <- fread("outputfolder_gen/gen_PCs/GenotypePCs.txt")
exp_cov <- fread("outputfolder_exp/exp_PCs/exp_PCs.txt")

colnames(gen_cov) <- c("SampleID", paste0("GenPC", 1:10))
colnames(exp_cov) <- c("SampleID", paste0("ExpPC", 1:50))

cov <- merge(gen_cov, exp_cov, by = "SampleID")

fwrite(cov, "CovariatesPCs.txt", sep = "\t", quote = FALSE, row.names = FALSE)