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

    add_cov_sup <- add_cov[, sapply(add_cov, is.character), with = FALSE]

    if (ncol(add_cov_sup) > 1){
    add_cov_sup <- dcast(data = melt(add_cov_sup, id.vars = "SampleID"), SampleID ~ variable + value, length)

    indik <- sapply(add_cov, is.numeric)
    indik[1] <- TRUE

    add_cov_sup2 <- add_cov[, indik, with = FALSE]

    add_cov <- merge(add_cov_sup, add_cov_sup2, by = "SampleID")

    }

    if (all(cov$SampleID %in% add_cov$SampleID)){
        cov <- merge(cov, add_cov, by = "SampleID")
        message("Additional covariates included to the CovariatePCs.txt.")

    } else {
        stop("Error: you have less covariates in the additional covariate file than in the PC file! Each sample has to be in the covariate file!")
    }

}

fwrite(cov, "CovariatePCs.txt", sep = "\t", quote = FALSE, row.names = FALSE)
