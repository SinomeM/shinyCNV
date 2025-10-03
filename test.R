library(ggplot2)
library(data.table)
library(plotly)

cnvs <- fread('./data/cnvs.txt')
samples <- fread('data/samples_list.txt')
snps <- fread('data/hd_1kG_hg19.snppos.filtered.test.gz')
snps[, ':=' (Name = as.character(Name),
            Chr = as.character(Chr),
            Position = as.integer(Position))]

load_sample_snps <- function(tabix_path, chr) {
  # Use system tabix to extract SNPs for the chromosome
  cmd <- paste0("tabix ", tabix_path, " ", chr, ":", 0,
                        "-", 275000000) # 275 Mbp is larger than chromosome 1
  snp_dt <- tryCatch({
    fread(cmd = cmd, header = FALSE)
  }, error = function(e) {
    data.table() # return empty table on error
  })
  # Assign column names if data is present
  if (ncol(snp_dt) >= 6) {
    setnames(snp_dt, c("chr", "start", "end", "LRR", "BAF", "LRR_adj"))
    snp_dt[, ':=' (chr = as.character(chr),
                   start = as.numeric(start),
                   end = as.numeric(end),
                   LRR = as.numeric(LRR),
                   BAF = as.numeric(BAF),
                   LRR_adj = as.numeric(LRR_adj))]
  }
  if (ncol(snp_dt) == 5) {
    setnames(snp_dt, c("chr", "start", "end", "LRR", "BAF"))
    snp_dt[, ':=' (chr = as.character(chr),
                   start = as.numeric(start),
                   end = as.numeric(end),
                   LRR = as.numeric(LRR),
                   BAF = as.numeric(BAF))]
  }
  return(unique(snp_dt))
}

i <- 1
cnv <- cnvs[i]

tabix_path <- samples[sample_ID == cnv$sample_ID, file_path_tabix]
chr <- cnv$chr

# Load SNPs for the sample and chromosome
snp_dt <- load_sample_snps(tabix_path, chr)
snps_chr <- snps[Chr == chr, ]

pl <- ggplot(snp_dt, aes(x = start, y = LRR_adj)) +
        geom_point() + theme_bw()
ggplotly(pl)

ply <- plot_ly(snp_dt, x = ~start, y = ~LRR_adj,
               type = "scatter", mode = "markers",
               marker = list(color = "black", opacity = 0.5),
               text = ~paste("Position:", start, "<br>LRR:", LRR_adj),
               hoverinfo = "text") |>
  layout(
    xaxis = list(range = c(24e6, 25e6)),
    yaxis = list(range = c(-2, 2), fixedrange = TRUE)
  )
ply
