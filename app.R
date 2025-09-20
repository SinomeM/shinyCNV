# CNV visualization tool

# This app takes as input:
#  - a folder to save results
#  - a table of CNVs (sample_ID, chr, start, end, GT, vo, CN, length, numsnp)
#  - a table linking each sample_ID to a tabix-indexed file (sample_ID, file_path_tabix)
#  - a table of filtered SNP (Name, Chr, Position)

# The tabix-indexed files should contain is expected to be a tab-separated file
# with the following columns (no header): chr, position, position, LRR, BAF, LRR_adj
# (if LRR_adj is not present, the app will use LRR instead).

# The purpose of the app is to visualize each CNV in a plot with LRR/BAF values
# in order to visually validate the call. Each CNV can be marked as either true
# false or unknown.
# The input CNV table can also be filtered based on different metrics (GT, minimum
# length or number of SNPs, previous visual inspection outcome).
# Also the app can focus on a specific locus ad show all CNVs in the region with
# minimum overlap with the fixed locus.

# This version is a major update of the previous app_old.R, with added features:
#  - new UI layout based on bslib (possibly built around page_sidebar())
#  - Interactive plots instead of static (can be zoomed in/out)
#  - Ability to update CNV coordinates (start/end)
#  - Ability to toggle the SNPs filtering on and off
#  - Other minor improvements

library(data.table)
library(ggplot2)
library(shiny)
library(bslib)

# Set this to true when running from command line
if (F) {
  cargs <- commandArgs(trailingOnly = T)
  # working dir, cnvs, samples, snps

  # quickly check inputs
  if (!dir.exists(cargs[1])) stop('Provided folder not found!')
  if (!file.exists(cargs[2])) stop('CNVs table not found!')
  if (!file.exists(cargs[3])) stop('Samples table not found!')
  if (!file.exists(cargs[4])) stop('SNPs table not found!')
  # Load data (CNVs, samples and SNPs)
  cnvs <- fread(cargs[2])
  samples <- fread(cargs[3])[, .(sample_ID, file_path_tabix)]
  snps <- fread(cargs[4])
}

# For testing purposes only
if (T) {
  cnvs <- fread('./data/cnvs.txt')
  samples <- fread('data/samples_list.txt')
  snps <- fread('data/hd_1kG_hg19.snppos.filtered.gz')
}

# Initialise 'vo' column if necessary and add empty columns if not present
if (!'vo' %in% colnames(cnvs)) cnvs[, vo := -9]
if (!'CN' %in% colnames(cnvs)) cnvs[, CN := NA]
if (!'length' %in% colnames(cnvs)) cnvs[, length := -9]
if (!'numsnp' %in% colnames(cnvs)) cnvs[, numsnp := -9]
cnvs <- cnvs[, .(sample_ID, chr, start, end, numsnp, length, GT, CN, vo)]
setorder(cnvs, chr, start)


# UI function ----

# Define the app main page using bslib::page_sidebar()

# The layout is a central page with CNV table one the top, CNV plot in the middle
# and the buttons to validate CNVs and move around at the bottom

# In the sidebar, there are options to filter CNVs, change plot height etc

ui <- page_sidebar(
  title = 'CNVs Visual Validation',
  sidebar = sidebar(
    title = 'Filters CNVs, Select locus and Settings',
    position = 'left'
  )
)



# Server function ----

# The main function of the server is to filter the CNV table if needed, load
# one CNV line and create the plot. The plot needs to be interactive, meaning it
# must be possible to zoom in and out. The plot is made of two panels, showing the
# LRR and BAF values for each SNP in the region respectively. The CNVs coordinates
# are marked with horizontal dashed lines. If a fixed locus in selected,
# the locus boundaries are also marked with vertical dashed lines.

# The boundaries updating functionalist will be added later

server <- function(input, output, session) {
  
}