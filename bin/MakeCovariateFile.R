library(data.table)

args <- commandArgs(trailingOnly = TRUE)

setDTthreads(1)

gen_cov <- fread("outputfolder_gen/gen_PCs/GenotypePCs.txt")
exp_cov <- fread("outputfolder_exp/exp_PCs/exp_PCs.txt")

colnames(gen_cov) <- c("SampleID", paste0("GenPC", 1:10))
colnames(exp_cov) <- c("SampleID", paste0("ExpPC", 1:min(100, ncol(exp_cov) - 1)))

cov <- merge(gen_cov, exp_cov, by = "SampleID")

if (length(args) > 0){


    add_cov <- fread(args[1])
    add_cov <- add_cov[complete.cases(add_cov), ]

    if (!all(apply(add_cov[, -1, with = FALSE], 2, is.numeric))){
        stop("Error: all covariates have to be numeric and not include NA's. Binary covariates have to be encoded as 0/1.")
    }

    if (all(cov$SampleID %in% add_cov$SampleID)){
        cov <- merge(cov, add_cov, by = "SampleID")
        message("Additional covariates included to the CovariatePCs.txt.")

    } else {
        stop("Error: you have less covariates in the additional covariate file than in the PC file! Each sample has to be in the covariate file!")
    }

}

fwrite(cov, "CovariatePCs.txt", sep = "\t", quote = FALSE, row.names = FALSE)
