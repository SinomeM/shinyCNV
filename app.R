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

library(DT)
library(bslib)
library(shiny)
library(ggplot2)
library(data.table)


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
if (!'length' %in% colnames(cnvs)) cnvs[, length := end-start+1]
if (!'numsnp' %in% colnames(cnvs)) cnvs[, numsnp := NA]
cnvs <- cnvs[, .(sample_ID, chr, start, end, numsnp, length, GT, CN, vo)]
setorder(cnvs, chr, start)


# UI function ----

# Define the app main page using bslib::page_sidebar()

# The layout is a central page with CNV table one the top, CNV plot in the middle
# and the buttons to validate CNVs and move around at the bottom

# In the sidebar, there are options to filter CNVs, change plot height etc

ui <- fluidPage(
  # Sidebar layout, contains filtering and settings
  layout_sidebar(
    sidebar = sidebar(
      position = 'left',
      # Various settings
      fluidRow(
        h4('Settings'),
        textInput('project_name', 'Project Name', ''),
        selectInput('snp_filtering', 'Filter SNPs based on input SNP table?',
                    c('yes', 'no'), 'yes'),
        sliderInput('plot_height', 'Change plot height',
                    min = 400, max = 1000, value = 750, step = 50)
      ),
      # CNVs filtering
      fluidRow(
        h4('Filter the CNV table'),
        selectInput('vo_filter', 'Filter CNV previous VI',
                    c('new', 'true', 'false', 'unkown', 'all'), 'all'),
        selectInput('gt_filter', 'Filter CNV GT', c('dels', 'dups', 'both'), 'both'),
        textInput("min_len_filter", "Minimum CNV length", '0'),
        textInput("max_len_filter", "Maximum CNV length", '10000000'),
        textInput("min_snp_filter", "Minimum number of SNPs", '0'),
        actionButton('run_filtering', 'Run Filtering', class = 'btn-primary')
      ),
      # Fixed locus
      fluidRow(
        h4('Select fixed locus'),
        textInput('locus_chr', 'Locus Chromosome', ''),
        textInput('locus_start', 'Locus Start Position', ''),
        textInput('locus_end', 'Locus End Position', ''),
        textInput('min_overlap', 'Minimum overlap with locus (bp)', '0')
      )
    ),
    # CNV table at the top of the main page
    div(
      DTOutput('cnv_table')
    ),
    # CNV plot
    div(
       plotOutput('cnv_plot')
    ),
    # Buttons to validate CNVs and move around
    div(
      fluidRow(
        actionButton("true", "True", class = "btn-success"),
        actionButton("false", "False", class = "btn-danger"),
        actionButton("unk", "Unkown", class = 'btn-warning'),
        actionButton("err", "Error"),
        actionButton("prv", "Previous", class = "btn-default"),
        actionButton("nxt", "Next", class = "btn-default"),
        textOutput('progress')
      )
    )
  )
)



# Server function ----

# The main function of the server is to filter the CNV table if needed, load
# one CNV line and create the plot. The plot needs to be interactive, meaning it
# must be possible to zoom in and out and each dot can be selected clicking on it.
# The plot is made of two panels, showing the
# LRR and BAF values for each SNP in the region respectively. The CNVs coordinates
# are marked with horizontal dashed lines. If a fixed locus in selected,
# the locus boundaries are also marked with vertical dashed lines.

# The boundaries updating functionalities will be added later

server <- function(input, output, session) {
  # 1. Initialize reactive values
  r_state <- reactiveValues(
    cnvs = cnvs,                    # original CNV table
    filtered_cnvs = cnvs,           # filtered version
    current_idx = 1                 # current CNV index
  )

  # 2. Filtering logic
  observeEvent(input$run_filtering, {
    filtered <- r_state$cnvs

    # Convert inputs to numeric and check if they are valid
    min_len <- as.numeric(input$min_len_filter)
    max_len <- as.numeric(input$max_len_filter)
    min_snp <- as.numeric(input$min_snp_filter)

    # Apply filters based on input values
    if (!is.null(input$min_len_filter)) {
      filtered <- filtered[length >= min_len]
    }
    if (!is.null(input$max_len_filter)) {
      filtered <- filtered[length <= max_len]
    }
    if (!is.null(input$min_snp_filter)) {
      filtered <- filtered[numsnp >= min_snp]
    }
    if (!is.null(input$gt_filter)) {
      if (input$gt_filter == 'dels') {
        filtered <- filtered[GT == 1, ]
      } else if (input$gt_filter == 'dups') {
        filtered <- filtered[GT == 2, ]
      }
    }
    if (!is.null(input$vo_filter)) {
      if (input$vo_filter == 'new') {
        filtered <- filtered[vo == -9, ]
      } else if (input$vo_filter == 'true') {
        filtered <- filtered[vo == 1, ]
      } else if (input$vo_filter == 'false') {
        filtered <- filtered[vo == 2, ]
      } else if (input$vo_filter == 'unkown') {
        filtered <- filtered[vo == 3, ]
      }
    }

    # Update filtered table and reset index
    r_state$filtered_cnvs <- filtered
    r_state$current_idx <- 1
  })

  # 3. Navigation and validation logic
  observeEvent(input$nxt, {
    if (r_state$current_idx < nrow(r_state$filtered_cnvs)) {
      r_state$current_idx <- r_state$current_idx + 1
    }
  })
  observeEvent(input$prv, {
    if (r_state$current_idx > 1) {
      r_state$current_idx <- r_state$current_idx - 1
    }
  })

  # Validation buttons
  observeEvent(input$true, {
    idx <- r_state$current_idx
    row <- r_state$filtered_cnvs[idx]
    r_state$cnvs[sample_ID == row$sample_ID & 
                 chr == row$chr & 
                 start == row$start & 
                 end == row$end, 
                 vo := 1]
    if (r_state$current_idx < nrow(r_state$filtered_cnvs)) {
      r_state$current_idx <- r_state$current_idx + 1
    }
  })
  observeEvent(input$false, {
    idx <- r_state$current_idx
    row <- r_state$filtered_cnvs[idx]
    r_state$cnvs[sample_ID == row$sample_ID & 
                 chr == row$chr & 
                 start == row$start & 
                 end == row$end, 
                 vo := 2]
    if (r_state$current_idx < nrow(r_state$filtered_cnvs)) {
      r_state$current_idx <- r_state$current_idx + 1
    }
  })
  observeEvent(input$unk, {
    idx <- r_state$current_idx
    row <- r_state$filtered_cnvs[idx]
    r_state$cnvs[sample_ID == row$sample_ID & 
                 chr == row$chr & 
                 start == row$start & 
                 end == row$end, 
                 vo := 3]
    if (r_state$current_idx < nrow(r_state$filtered_cnvs)) {
      r_state$current_idx <- r_state$current_idx + 1
    }
  })
  observeEvent(input$err, {
    idx <- r_state$current_idx
    row <- r_state$filtered_cnvs[idx]
    r_state$cnvs[sample_ID == row$sample_ID & 
                 chr == row$chr & 
                 start == row$start & 
                 end == row$end, 
                 vo := -7]
    if (r_state$current_idx < nrow(r_state$filtered_cnvs)) {
      r_state$current_idx <- r_state$current_idx + 1
    }
  })

  # 4. Render current CNV table row
  output$cnv_table <- renderDT({
    current_row <- r_state$filtered_cnvs[r_state$current_idx]
    datatable(
      current_row,
      options = list(dom = 't', 
                    paging = FALSE,
                    searching = FALSE,
                    info = FALSE),
      rownames = FALSE
    )
  })

  # 5. Progress text
  output$progress <- renderText({
    total <- nrow(r_state$filtered_cnvs)
    idx <- r_state$current_idx
    percent <- if (total > 0) round(100 * (idx-1) / total, 1) else 0
    sprintf("CNV %d out of %d (%.1f%%)", idx, total, percent)
  })

  # 6. Helper function to load SNPs for a sample and chromosome using tabix
  load_sample_snps <- function(tabix_path, chr) {
    # Use system tabix to extract SNPs for the chromosome (UPDATE COMMAND)
    cmd <- sprintf("tabix %s %s", shQuote(tabix_path), shQuote(chr))
    snp_dt <- tryCatch({
      fread(cmd = cmd, header = FALSE)
    }, error = function(e) {
      data.table() # return empty table on error
    })
    # Assign column names if data is present
    if (ncol(snp_dt) >= 6) {
      setnames(snp_dt, c("chr", "start", "end", "LRR", "BAF", "LRR_adj"))
    }
    snp_dt
  }

  # 7. CNV plot
  output$cnv_plot <- renderPlot({
    # Get current CNV
    cnv <- r_state$filtered_cnvs[r_state$current_idx]
    if (is.null(cnv) || nrow(cnv) == 0) {
      plot.new()
      title("No CNV selected")
      return()
    }

    # Get sample tabix file path
    sample_row <- samples[sample_ID == cnv$sample_ID]
    tabix_path <- sample_row$file_path_tabix
    chr <- cnv$chr

    # Load SNPs for the sample and chromosome
    snp_dt <- load_sample_snps(tabix_path, chr)

    # Placeholder plot for LRR and BAF
    par(mfrow = c(2, 1))
    plot(snp_dt$start, snp_dt$LRR, pch = 20, col = "blue",
         main = "LRR values", xlab = "Position", ylab = "LRR")
    abline(v = c(cnv$start, cnv$end), lty = 2, col = "red")
    plot(snp_dt$start, snp_dt$BAF, pch = 20, col = "green",
         main = "BAF values", xlab = "Position", ylab = "BAF")
    abline(v = c(cnv$start, cnv$end), lty = 2, col = "red")
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)