#!/usr/bin/env Rscript
library(data.table)
library(qqman)
library(ggplot2)

# Script to create genome-wide Manhattan and QQ plots from merged REGENIE output
# This script greatly increases the speed of plotting by subsampling a small fraction
# of the non-significant SNPs

# INPUT: INFILE - merged .regenie file (all chromosomes)
# OUTPUT: OUTFILE - output prefix for plots

# ------------------------------ get.GIF ------------------------------ #
get.GIF <- function(p.vec)
    # compute genomic inflation factor
    # input: p.vec (numeric), vector of p-values
    # output: GIF
{
    median(qchisq(1 - p.vec, 1)) / qchisq(0.5, 1)
}
# --------------------------------------------------------------------- #

# -------------------------- gg.qq ------------------------------------ #
gg.qq <- function(p.vec, x.pos, y.pos)
    # recreation of qq from qqman, but using ggplot
    # input:
    #   p.vec (numeric), vector of p-values
    #   x.pos, y.pos (numeric), x/y coordinates for the position of the GIF text
    # output:
    #   p1, a ggplot object
{
    if ((!is.numeric(p.vec)))
        stop("input must be numeric")

    p.vec <- p.vec[!is.na(p.vec) & !is.nan(p.vec) & !is.null(p.vec) & is.finite(p.vec) &
                       p.vec < 1 & p.vec > 0]
    o <- -log10(sort(p.vec, decreasing = FALSE))
    e <- -log10(ppoints(length(p.vec)))

    txt <- paste("lambda", round(get.GIF(p.vec), 2), sep = " == ")
    p1 <- ggplot(data = data.frame(o = o, e = e), aes(x = e, y = o)) +
        geom_point() +
        theme_minimal() +
        xlab("Expected") +
        ylab("Observed") +
        geom_abline(slope = 1, intercept = 0, color = "red") +
        annotate("text", x = x.pos, y = y.pos,
                 label = txt,
                 parse = T,
                 size = unit(10, "pt"))
    return(p1)
}
# --------------------------------------------------------------------- #

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    print("Usage: Rscript genome_wide_regenie_plot.R INFILE OUTFILE")
    print("INFILE should be a merged REGENIE format text file (all chromosomes)")
    print("OUTFILE should be the output prefix (without extension)")
    stop()
}

INFILE <- args[1]
OUTFILE <- args[2]

print(paste("Input file:", INFILE))
print(paste("Output prefix:", OUTFILE))

# Check if input file exists
if (!file.exists(INFILE)) {
    stop(paste("ERROR: Input file not found:", INFILE))
}

# PARAMS
# the threshold for p-values; sub-sample any values with p > p.thresh
p.thresh <- 0.01

# proportion of snps to keep
k.thresh <- 0.2

# Output files
MANOUTFILE_PNG <- paste(OUTFILE, "_manhattan.png", sep = "")
QQOUTFILE_PNG <- paste(OUTFILE, "_qq.png", sep = "")
MANOUTFILE_PDF <- paste(OUTFILE, "_manhattan.pdf", sep = "")
QQOUTFILE_PDF <- paste(OUTFILE, "_qq.pdf", sep = "")

# ---

print("Reading association testing data...")

dat <- fread(INFILE, header = TRUE)
dat <- as.data.frame(dat)

print(paste("Total SNPs loaded:", nrow(dat)))

# Handle X chromosome conversion
if (any(dat$CHROM == "X")) {
    dat$CHROM[which(dat$CHROM == "X")] <- 23
}

# Data type conversions
dat$ID <- as.character(dat$ID)
dat$CHROM <- as.integer(dat$CHROM)
dat$GENPOS <- as.integer(dat$GENPOS)
dat$CHISQ <- as.numeric(dat$CHISQ)
dat$LOG10P <- as.numeric(dat$LOG10P)
dat$INFO <- as.numeric(dat$INFO)
dat$PVAL <- as.numeric(10^(-dat$LOG10P))

# Remove any rows with missing essential data
dat <- dat[!is.na(dat$CHROM) & !is.na(dat$GENPOS) & !is.na(dat$PVAL), ]

print(paste("SNPs after filtering:", nrow(dat)))

# Create chromosome labels
chrlabs <- as.character(seq(1:length(unique(dat$CHROM))))

print("Data loading complete!")

# ---

print("Subsampling association data...")

# Subsample non-significant SNPs
nonsig_indices <- which(dat$PVAL > p.thresh)
x <- length(nonsig_indices)

if (x > 0) {
    s <- sample(dat$ID[nonsig_indices], round(k.thresh * x), replace = FALSE)
    print(paste("Subsampled", length(s), "non-significant SNPs out of", x))
} else {
    s <- character(0)
    print("No non-significant SNPs to subsample")
}

# Create subset for plotting (significant SNPs + subsampled non-significant)
plot_dat <- dat[(dat$ID %in% s) | (dat$PVAL <= p.thresh), ]

print(paste("Total SNPs for plotting:", nrow(plot_dat)))
print("Subsampling complete!")

# ---

print("Creating Manhattan plot...")

# PNG version
png(MANOUTFILE_PNG, height = 600, width = 1200)
manhattan(plot_dat,
          chr = "CHROM", bp = "GENPOS", p = "PVAL", snp = "ID", 
          col = c("gray", "black"),
          logp = TRUE, chrlabs = chrlabs,
          main = paste("Manhattan Plot - Genome-wide Association"))
dev.off()

# PDF version
pdf(MANOUTFILE_PDF, height = 6, width = 12)
manhattan(plot_dat,
          chr = "CHROM", bp = "GENPOS", p = "PVAL", snp = "ID", 
          col = c("gray", "black"),
          logp = TRUE, chrlabs = chrlabs,
          main = paste("Manhattan Plot - Genome-wide Association"))
dev.off()

print("Manhattan plot complete!")

# ---

print("Creating QQ plot...")

# PNG version
png(QQOUTFILE_PNG, height = 600, width = 600)
print(gg.qq(dat$PVAL, 1, 3))
dev.off()

# PDF version
pdf(QQOUTFILE_PDF, height = 6, width = 6)
print(gg.qq(dat$PVAL, 1, 3))
dev.off()

print("QQ plot complete!")

# ---

print("====== PLOTTING SUMMARY ======")
print(paste("Manhattan plot (PNG):", MANOUTFILE_PNG))
print(paste("Manhattan plot (PDF):", MANOUTFILE_PDF))
print(paste("QQ plot (PNG):", QQOUTFILE_PNG))
print(paste("QQ plot (PDF):", QQOUTFILE_PDF))

# Print some summary statistics
print(paste("Total SNPs analyzed:", nrow(dat)))
print(paste("SNPs with p < 0.05:", sum(dat$PVAL < 0.05)))
print(paste("SNPs with p < 0.01:", sum(dat$PVAL < 0.01)))
print(paste("SNPs with p < 1e-5:", sum(dat$PVAL < 1e-5)))
print(paste("SNPs with p < 5e-8:", sum(dat$PVAL < 5e-8)))
print(paste("Genomic inflation factor (lambda):", round(get.GIF(dat$PVAL), 3)))

print("Script completed successfully!")
