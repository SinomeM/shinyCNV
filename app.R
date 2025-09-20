library(data.table)
library(ggplot2)
library(shiny)

cargs <- commandArgs(trailingOnly = T)
# working dir, cnvs, samples, snps


# quickly check inputs
if (!dir.exists(cargs[1])) stop('Provided folder not found!')
if (!file.exists(cargs[2])) stop('CNVs table not found!')
if (!file.exists(cargs[3])) stop('Samples table not found!')
if (!file.exists(cargs[4])) stop('SNPs table not found!')

# load data (CNVs, samples and SNPs) and initialise 'vo' column if necessary
cnvs <- fread(cargs[2])
if (!'vo' %in% colnames(cnvs)) cnvs[, vo := -9]
if (!'CN' %in% colnames(cnvs)) cnvs[, CN := NA]
if (!'length' %in% colnames(cnvs)) cnvs[, length := -9]
if (!'numsnp' %in% colnames(cnvs)) cnvs[, numsnp := -9]
cnvs <- cnvs[, .(sample_ID, chr, start, end, numsnp, length, GT, CN, vo)]
samples <- fread(cargs[3])[, .(sample_ID, file_path_tabix)]
snps <- fread(cargs[4])
setorder(cnvs, chr, start)