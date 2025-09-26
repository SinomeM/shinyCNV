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
      title = 'Filters CNVs, Select locus and Settings',
      position = 'left',
      # Various settings
      fluidRow(
        h6('Settings'),
        textInput('project_name', 'Project Name', ''),
        selectInput('snp_filtering', 'Filter SNPs based on input SNP table?',
                    c('yes', 'no'), 'yes'),
        sliderInput('plot_height', 'Change plot height',
                    min = 400, max = 1000, value = 750, step = 50)
      ),
      # CNVs filtering
      fluidRow(
        h6('Filter the CNV table'),
        selectInput('vo_filter', 'Filter CNV previous VI',
                    c('new', 'true', 'false', 'unk', 'other', 'all'), 'all'),
        selectInput('gt_filter', 'Filter CNV GT', c('dels', 'dups', 'both'), 'both'),
        textInput("min_len_filter", "Minimum CNV length", '0'),
        textInput("max_len_filter", "Maximum CNV length", '10000000'),
        textInput("min_snp_filter", "Minimum number of SNPs", '0')
      ),
      # Fixed locus
      fluidRow(
        h6('Select fixed locus'),
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
        actionButton("prev", "Previous", class = "btn-default"),
        actionButton("next", "Next", class = "btn-default"),
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

# The boundaries updating functionalist will be added later

server <- function(input, output, session) {
  # 1. Reactive values to store state
  r_state <- reactiveValues(
    cnvs = cnvs,
    filtered_cnvs = cnvs,
    current_idx = 1
  )

  # 2. Filtering logic
  observeEvent(
    {
      input$vo_filter
      input$gt_filter
      input$min_len_filter
      input$max_len_filter
      input$min_snp_filter
      input$locus_chr
      input$locus_start
      input$locus_end
      input$min_overlap
    },
    {
      # TODO: Implement filtering logic here
      # Update r_state$filtered_cnvs based on filters
      # Reset r_state$current_idx to 1 if needed
    }
  )

  # 3. Navigation logic
  observeEvent(input$next, {
    if (r_state$current_idx < nrow(r_state$filtered_cnvs)) {
      r_state$current_idx <- r_state$current_idx + 1
    }
  })
  observeEvent(input$prev, {
    if (r_state$current_idx > 1) {
      r_state$current_idx <- r_state$current_idx - 1
    }
  })
  observeEvent(input$true, {
    idx <- r_state$current_idx
    row <- r_state$filtered_cnvs[idx, ]
    r_state$cnvs[sample_ID == row$sample_ID & chr == row$chr & start == row$start & end == row$end, vo := 1]
    if (r_state$current_idx < nrow(r_state$filtered_cnvs)) {
      r_state$current_idx <- r_state$current_idx + 1
    }
  })
  observeEvent(input$false, {
    idx <- r_state$current_idx
    row <- r_state$filtered_cnvs[idx, ]
    r_state$cnvs[sample_ID == row$sample_ID & chr == row$chr & start == row$start & end == row$end, vo := 2]
    if (r_state$current_idx < nrow(r_state$filtered_cnvs)) {
      r_state$current_idx <- r_state$current_idx + 1
    }
  })
  observeEvent(input$unk, {
    idx <- r_state$current_idx
    row <- r_state$filtered_cnvs[idx, ]
    r_state$cnvs[sample_ID == row$sample_ID & chr == row$chr & start == row$start & end == row$end, vo := 3]
    if (r_state$current_idx < nrow(r_state$filtered_cnvs)) {
      r_state$current_idx <- r_state$current_idx + 1
    }
  })
  observeEvent(input$err, {
    idx <- r_state$current_idx
    row <- r_state$filtered_cnvs[idx, ]
    r_state$cnvs[sample_ID == row$sample_ID & chr == row$chr & start == row$start & end == row$end, vo := -7]
    if (r_state$current_idx < nrow(r_state$filtered_cnvs)) {
      r_state$current_idx <- r_state$current_idx + 1
    }
  })

# 4. CNV table
  output$cnv_table <- renderDT({
    datatable(
      r_state$filtered_cnvs[r_state$current_idx],
      class = 'cell-border stripe',
      options = list(dom = "t", paging = FALSE, info = FALSE, searching = FALSE),
      rownames = F
    )
  })
}


# Run the app ----
shinyApp(ui = ui, server = server)