#!/usr/bin/env Rscript

### ANCESTRY: PHENOTYPE & COVARIATE FILES ### 
# Usage: Rscript 10b_prep_files.R AFR
# ANC - EUR, AFR

ARGS <- commandArgs(trailingOnly = TRUE)
ANC <- ARGS[1]

library(data.table)

HDIR <- "/project/knathans_tecac/jenny/breast/analysis"
WDIR <- file.path(HDIR, ANC, "cov")
ID_DIR <- file.path(HDIR, "ids")
dir.create(WDIR, showWarnings = FALSE)
setwd(WDIR)

# Input files
covar_file <- "/static/PMBB/PMBB-Release-2024-3.0/Phenotypes/3.0/PMBB-Release-2024-3.0_covariates.txt"
id_file <- file.path(HDIR, ANC, "step1", paste0(ANC, "_common_snps_qc.id"))
eigenvec_file <- file.path(WDIR, paste0(ANC, "_common_snps_qc.eigenvec"))
case_file <- file.path(ID_DIR, paste0(ANC, "_case_ids.txt"))
control_file <- file.path(ID_DIR, paste0(ANC, "_control_ids.txt"))
exclude_file <- file.path(WDIR, "fid_iid_remove.txt")

### PHENOTYPE ### 
ids <- read.table(id_file, header = FALSE, col.names = c("FID", "IID"))
exclude_ids <- read.table(exclude_file, header = FALSE, col.names = c("FID", "IID"))

ids$combined_id <- paste(ids$FID, ids$IID, sep = "_")
exclude_ids$combined_id <- paste(exclude_ids$FID, exclude_ids$IID, sep = "_")
ids <- ids[!ids$combined_id %in% exclude_ids$combined_id, ]
ids$combined_id <- NULL

cases <- readLines(case_file)
controls <- readLines(control_file)

pheno <- data.frame(
  FID = ids$FID,
  IID = ids$IID,
  Y1 = ifelse(ids$IID %in% cases, 1, 
             ifelse(ids$IID %in% controls, 0, NA))
)

write.table(pheno, file = "phenotype.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)


### COVARIATES ###
covar <- read.table(covar_file, header = TRUE)
pca <- read.table(eigenvec_file, header = FALSE, 
                 col.names = c("FID", "IID", paste0("PC", 1:10)))

setDT(ids); setDT(covar); setDT(pca)

result <- merge(ids, pca, by = c("FID", "IID"), all.x = TRUE)
result <- merge(result, 
               covar[, .(IID = person_id, Batch = Batch, Age = Sample_age)],
               by = "IID", all.x = TRUE)

# Set NA for missing values
cols <- c(paste0("PC", 1:10), "Batch", "Age")
result[, (cols) := lapply(.SD, function(x) ifelse(is.na(x), "NA", x)), .SDcols = cols]

fwrite(result[, .(FID, IID, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, Batch, Age)],
       file = "covariates.txt", sep = "\t", quote = FALSE)
