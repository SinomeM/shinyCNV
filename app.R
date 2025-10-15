library(DT)
library(bslib)
library(shiny)
library(ggplot2)
library(data.table)
library(plotly)

# Controls whether the inputs are taken from command line or from test data
testing <- TRUE

if (testing == FALSE) {
  # working dir, cnvs, samples, snps
  cargs <- commandArgs(trailingOnly = T)

  # Check inputs
  if (!dir.exists(cargs[1])) warning('Provided folder not found!')
  if (!file.exists(cargs[2])) stop('CNVs table not found!')
  if (!file.exists(cargs[3])) stop('Samples table not found!')
  if (!file.exists(cargs[4])) stop('SNPs table not found!')
  # Load data (CNVs, samples and SNPs)
  wkdir <- cargs[1]
  cnvs <- fread(cargs[2])
  samples <- fread(cargs[3])[, .(sample_ID, file_path_tabix)]
  snps <- fread(cargs[4])
}

if (testing) {
  wkdir <- './tmp'
  cnvs <- fread('./data/cnvs.txt')
  samples <- fread('./data/samples_list.txt')
  snps <- fread('./data/hd_1kG_hg19.snppos.filtered.test.gz')
  snps[, ':=' (Name = as.character(Name),
              Chr = as.character(Chr),
              Position = as.integer(Position))]
}

# Initialise 'vo' column if necessary and add empty columns if not present
if (!'vo' %in% colnames(cnvs)) cnvs[, vo := -9]
if (!'CN' %in% colnames(cnvs)) cnvs[, CN := NA]
if (!'length' %in% colnames(cnvs)) cnvs[, length := end-start+1]
if (!'numsnp' %in% colnames(cnvs)) cnvs[, numsnp := NA]
cnvs <- cnvs[, .(sample_ID, chr, start, end, numsnp, length, GT, CN, vo)]
setorder(cnvs, chr, start)
cnvs[, ':=' (original_start = start,
             original_end = end)]

# Check if the provided wkdir exists, if not create it (recursive set to
# false, if it's needed it's likely the user made a mistake)
if (!dir.exists(wkdir)) dir.create(wkdir, recursive = F)


# UI function ----

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .shiny-notification {
        position: fixed;
        top: 1rem;
        left: 1rem;
      }
      .shiny-notification-warning {
        background-color: #c62828;
        color: #ffffff;
        font-size: 2rem;
        font-weight: 600;
        box-shadow: 0 0 20px rgba(198, 40, 40, 0.7);
      }
    ")),
    tags$script(HTML("
      $(document).on('keydown', function(e) {
        if ($(e.target).is('input, textarea, select, button')) return;
        switch (e.key) {
          case 't':
          case 'T':
            $('#btn_true').click();
            break;
          case 'f':
          case 'F':
            $('#btn_false').click();
            break;
          case 'u':
          case 'U':
            $('#btn_unk').click();
            break;
          case 'n':
          case 'N':
            $('#btn_nxt').click();
            break;
          case 'p':
          case 'P':
            $('#btn_prv').click();
            break;
          case 'e':
          case 'E':
            $('#btn_err').click();
            break;
          case 'r':
          case 'R':
            $('#btn_ref').click();
            break;
        }
      });
    "))
  ),
  layout_sidebar(
    sidebar = sidebar(
      position = 'left',
      # Various settings
      fluidRow(
        h4('Settings'),
        textInput('project_name', 'Project Name', ''),
        selectInput('snp_filtering', 'Show only filtered SNPs?',
                    c('yes', 'no'), 'yes'),
        sliderInput('zoom_lvl', 'Base zoom level', min = 2, max = 20, value = 10, step = 2)
      ),
      # CNVs filtering
      fluidRow(
        h4('Filter the CNV table'),
        selectInput('vo_filter', 'Filter CNV previous VI',
                    c('all', 'new', 'true', 'false', 'unkown', 'error'), 'all'),
        selectInput('gt_filter', 'Filter CNV GT', c('both', 'dels', 'dups'), 'both'),
        textInput("min_len_filter", "Minimum CNV length (bp)", '25000'),
        textInput("max_len_filter", "Maximum CNV length (bp)", '10000000'),
        textInput("min_snp_filter", "Minimum number of SNPs", '5'),
        actionButton('run_filtering', 'Simple Filtering', class = 'btn-primary')
      ),
      # Fixed locus
      fluidRow(
        h4('Select fixed locus'),
        textInput('locus_name', 'Locus Name', ''),
        textInput('locus_chr', 'Locus Chromosome', ''),
        textInput('locus_start', 'Locus Start Position (bp)', ''),
        textInput('locus_end', 'Locus End Position (bp)', ''),
        textInput('min_overlap', 'Minimum overlap with locus (proportion)', '0'),
        actionButton('run_fixed_locus', 'Fixed Locus Mode',
                     class = 'btn-primary')
      ),
      # Select CNVs overlapping region
      fluidRow(
        h4('Select CNVs overlapping a region'),
        textInput('region_chr', 'Region Chromosome', ''),
        textInput('region_start', 'Region Start Position (bp)', ''),
        textInput('region_end', 'Region End Position (bp)', ''),
        textInput('min_iou', 'Minimum IOU with region', '0'),
        actionButton('run_select_region', 'Select CNVs in Region',
                     class = 'btn-primary')
      )
    ),
    # CNV table at the top of the main page
    div(
      DTOutput('cnv_table')
    ),
    # CNV plot
    div(
       #plotOutput('cnv_plot')
       plotlyOutput('cnv_plotly')
    ),
    # Buttons to validate CNVs and move around
    div(
      fluidRow(
        actionButton("btn_true", "True", class = "btn-success"),
        actionButton("btn_false", "False", class = "btn-danger"),
        actionButton("btn_unk", "Unkown", class = 'btn-warning'),
        actionButton("btn_err", "Error"),
        actionButton("btn_prv", "Previous", class = "btn-default"),
        actionButton("btn_nxt", "Next", class = "btn-default"),
        actionButton("btn_ref", "Update CNV coordinates", class = "btn-info"),
        textOutput('progress')
      ),
      fluidRow(
        verbatimTextOutput("selected")
      )
    )
  )
)


# Server function ----

server <- function(input, output, session) {
  # 1. Initialize reactive values
  r_state <- reactiveValues(
    cnvs = cnvs,                       # original CNV table
    filtered_cnvs = cnvs,              # filtered version
    current_idx = 1,                   # current CNV index
    fixed_locus = FALSE,               # whether in fixed locus mode
    select_region = FALSE,             # whether in select region mode
    snp_pos = integer(),               # positions of filtered SNPs
    filtered_snps_chr = data.table(),  # filtered SNPs for the current chromosome
    plot_retrigger = 0                 # to force plot re-rendering
  )

  # 2. Filtering logic
  observeEvent(input$run_filtering, {
    # always start from the original CNV table
    filtered <- r_state$cnvs

    # Convert inputs to integer and check if they are valid
    min_len <- as.integer(input$min_len_filter)
    max_len <- as.integer(input$max_len_filter)
    min_snp <- as.integer(input$min_snp_filter)

    # Apply filters based on input values
    if (!is.na(input$min_len_filter)) {
      filtered <- filtered[length >= min_len]
    }
    if (!is.na(input$max_len_filter)) {
      filtered <- filtered[length <= max_len]
    }
    if (!is.na(input$min_snp_filter)) {
      filtered <- filtered[numsnp >= min_snp]
    }
    if (!is.null(input$gt_filter) && input$gt_filter != 'both' &&
        !is.na(input$gt_filter)) {
      if (input$gt_filter == 'dels') {
        filtered <- filtered[GT == 1, ]
      } else if (input$gt_filter == 'dups') {
        filtered <- filtered[GT == 2, ]
      }
    }
    if (!is.null(input$vo_filter) && input$vo_filter != 'all' &&
        !is.na(input$vo_filter)) {
      if (input$vo_filter == 'new') {
        filtered <- filtered[vo == -9, ]
      } else if (input$vo_filter == 'true') {
        filtered <- filtered[vo == 1, ]
      } else if (input$vo_filter == 'false') {
        filtered <- filtered[vo == 2, ]
      } else if (input$vo_filter == 'unkown') {
        filtered <- filtered[vo == 3, ]
      } else if (input$vo_filter == 'error') {
        filtered <- filtered[vo == -7, ]
      }
    }

    # Update filtered table and reset index
    r_state$filtered_cnvs <- filtered
    r_state$current_idx <- 1
  })

  # 3. Navigation and validation logic
  observeEvent(input$btn_nxt, {
    r_state$current_idx <- r_state$current_idx + 1
  })
  observeEvent(input$btn_prv, {
    r_state$current_idx <- r_state$current_idx - 1
  })

  # Validation buttons
  observeEvent(input$btn_true, {
    idx <- r_state$current_idx
    r_state$filtered_cnvs[idx, vo := 1]
    r_state$current_idx <- r_state$current_idx + 1
  })
  observeEvent(input$btn_false, {
    idx <- r_state$current_idx
    r_state$filtered_cnvs[idx, vo := 2]
    r_state$current_idx <- r_state$current_idx + 1
  })
  observeEvent(input$btn_unk, {
    idx <- r_state$current_idx
    r_state$filtered_cnvs[idx, vo := 3]
    r_state$current_idx <- r_state$current_idx + 1
  })
  observeEvent(input$btn_err, {
    idx <- r_state$current_idx
    r_state$filtered_cnvs[idx, vo := -7]
    r_state$current_idx <- r_state$current_idx + 1
  })

  # 4. Render current CNV table row
  output$cnv_table <- renderDT({
    r_state$plot_retrigger  # to force re-rendering when needed
    current_row <- r_state$filtered_cnvs[r_state$current_idx]
    datatable(
      current_row[, .(sample_ID, chr, "start (Mbp)" = round(start/1e6, 2),
                      "end (Mbp)" = round(end/1e6, 2), numsnp,
                      "length (Mbp)" = round(length/1e6, 2), GT, CN, vo)],
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
  # Cache SNPs to avoid reloading them multiple times
  snp_cache <- new.env(parent = emptyenv())

  load_sample_snps <- function(tabix_path, chr) {
    # Check cache first
    key <- paste(tabix_path, chr, sep = "::")
    if (exists(key, envir = snp_cache, inherits = FALSE))
      return(get(key, envir = snp_cache))

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
                     start = as.integer(start),
                     end =   as.integer(end),
                     LRR = as.numeric(LRR),
                     BAF = as.numeric(BAF),
                     LRR_adj = as.numeric(LRR_adj))]
    }
    if (ncol(snp_dt) == 5) {
      setnames(snp_dt, c("chr", "start", "end", "LRR", "BAF"))
      snp_dt[, ':=' (chr = as.character(chr),
                     start = as.integer(start),
                     end =   as.integer(end),
                     LRR = as.numeric(LRR),
                     BAF = as.numeric(BAF))]
    }

    # Store in cache
    assign(key, snp_dt, envir = snp_cache)

    return(unique(snp_dt))
  }

  # 7. CNV plot
  # plotly version
  output$cnv_plotly <- renderPlotly({
    r_state$plot_retrigger  # to force re-rendering when needed

    # if no CNVs or index out of bounds, return empty plot
    if (nrow(r_state$filtered_cnvs) == 0 ||
        r_state$current_idx > nrow(r_state$filtered_cnvs)) {
      return(plotly_empty())
    }

    # Get current CNV, if null or empty return empty plot
    cnv <- r_state$filtered_cnvs[r_state$current_idx]
    if (is.null(cnv) || nrow(cnv) == 0) {
      return(plotly_empty())
    }

    # Load data for the sample in the entire chromosome
    tabix_path <- samples[sample_ID == cnv$sample_ID, file_path_tabix]
    chr <- cnv$chr
    snp_dt <- load_sample_snps(tabix_path, chr)

    # Clean/filtered SNPs in the chromosome
    if (r_state$current_idx == 1) {
      r_state$filtered_snps_chr <- snps[Chr == chr,]
      r_state$snp_pos <- r_state$filtered_snps_chr[, unique(Position)]
    }
    # chromosome changed, reload filtered SNPs and clear cache
    else if (r_state$filtered_snps_chr[, Chr][1] != chr) {
      if (length(ls(envir = snp_cache)) > 0)
        rm(list = ls(envir = snp_cache), envir = snp_cache)

      r_state$filtered_snps_chr <- snps[Chr == chr,]
      r_state$snp_pos <- r_state$filtered_snps_chr[, unique(Position)]
    }

    # Select filtered SNPs if requested
    if (!is.null(input$snp_filtering) && input$snp_filtering == "yes") {
      snp_dt <- snp_dt[start %in% r_state$snp_pos, ]
    }
    # Empty plot if no snps are left
    if (nrow(snp_dt) == 0) return(plotly_empty())

    # Clamp LRR values into [-1.5, 1.5]
    if ("LRR" %in% names(snp_dt)) {
      snp_dt[, LRR := pmax(pmin(LRR, 1.5), -1.5)]
    }
    if ("LRR_adj" %in% names(snp_dt)) {
      snp_dt[, LRR_adj := pmax(pmin(LRR_adj, 1.5), -1.5)]
    }

    # Use LRR_adj if present, else LRR
    lrr_col <- if ("LRR_adj" %in% names(snp_dt)) "LRR_adj" else "LRR"

    # Default zoom window (ADD, center on locus in fixed locus and region modes)
    cnv_len <- if (!is.na(cnv$length)) cnv$length else (cnv$end - cnv$start + 1)
    flank <- input$zoom_lvl * cnv_len
    window_start <- max(cnv$start - flank, 0)
    window_end <- cnv$end + flank

    # Get all CNVs on the  chromosome for this sample (before filtering)
    all_cnvs <- r_state$cnvs[sample_ID == cnv$sample_ID & chr == cnv$chr, ]

    # Create the CNVs outlines. Note that the current CNV has a thin border
    # while the other CNVs have no border
    lrr_outlines <- list()
    baf_outlines <- list()
    for (i in 1:all_cnvs[, .N]) {
      line <- all_cnvs[i]
      # Find if this CNV matches the current one
      is_current <- (line$sample_ID == cnv$sample_ID && 
                     line$start == cnv$start && 
                     line$end == cnv$end)
      lrr_outlines[[i]] <- list(
        type = "rect",
        xref = "x", x0 = line$start, x1 = line$end,
        yref = "y", y0 = -1.5, y1 = 1.5,
        fillcolor = "rgba(255, 0, 0, 0.1)",
        line = list(width = ifelse(is_current, 0.15, 0))
      )
      baf_outlines[[i]] <- list(
        type = "rect",
        xref = "x", x0 = line$start, x1 = line$end,
        yref = "y", y0 = 0, y1 = 1,
        fillcolor = "rgba(0, 0, 255, 0.1)",
        line = list(width = ifelse(is_current, 0.15, 0))
      )
    }

    # Add an additional outline (in light gray) for the locus/region if relevant
    if (r_state$fixed_locus) {
      ix <- length(lrr_outlines) + 1
      lrr_outlines[[ix]] <- list(
        type = "rect",
        xref = "x", x0 = as.integer(input$locus_start),
        x1 = as.integer(input$locus_end),
        yref = "y", y0 = -1.5, y1 = 1.5,
        fillcolor = "rgba(200, 200, 200, 0.3)",
        line = list(width = 0.5)
      )
      baf_outlines[[ix]] <- list(
        type = "rect",
        xref = "x", x0 = as.integer(input$locus_start),
        x1 = as.integer(input$locus_end),
        yref = "y", y0 = 0, y1 = 1,
        fillcolor = "rgba(200, 200, 200, 0.3)",
        line = list(width = 0.4)
      )
    }
    if (r_state$select_region) {
      ix <- length(lrr_outlines) + 1
      lrr_outlines[[ix]] <- list(
        type = "rect",
        xref = "x", x0 = as.integer(input$region_start),
        x1 = as.integer(input$region_end),
        yref = "y", y0 = -1.5, y1 = 1.5,
        fillcolor = "rgba(200, 200, 200, 0.3)",
        line = list(width = 0.3)
      )
      baf_outlines[[ix]] <- list(
        type = "rect",
        xref = "x", x0 = as.integer(input$region_start),
        x1 = as.integer(input$region_end),
        yref = "y", y0 = 0, y1 = 1,
        fillcolor = "rgba(200, 200, 200, 0.3)",
        line = list(width = 0.3)
      )
    }

    # LRR plot
    p1 <- plot_ly(
      snp_dt, x = ~start, y = as.formula(paste0("~", lrr_col)),
      type = "scatter", mode = "markers",
      name = "LRR",
      marker = list(color = "red", opacity = 0.8),
      hoverinfo = "text"
    ) %>%
      layout(
        xaxis = list(range = c(window_start, window_end),
                     autorange = FALSE, title = "Position (Mbp)"),
        yaxis = list(range = c(-1.5, 1.5), fixedrange = TRUE, title = "LRR"),
        shapes = lrr_outlines
      )

    # BAF plot
    p2 <- plot_ly(
      snp_dt, x = ~start, y = ~BAF,
      type = "scatter", mode = "markers",
      name = "BAF",
      marker = list(color = "blue", opacity = 0.8),
      hoverinfo = "text"
    ) %>%
      layout(
        xaxis = list(range = c(window_start, window_end), autorange = FALSE),
        yaxis = list(range = c(0, 1), fixedrange = TRUE, title = "BAF"),
        shapes = baf_outlines
      )

    # Combine LRR and BAF plots into a single figure with shared x-axis
    fig <- subplot(p2, p1, nrows = 2, shareX = TRUE, titleY = TRUE) %>%
      layout(dragmode = "select")  %>%
      config(
        modeBarButtonsToRemove = c(
          'lasso2d', 'autoScale2d', 'hoverClosestCartesian', 'zoomin2d',
          'hoverCompareCartesian', 'toggleSpikelines', 'toImage', 'zoomout2d'),
        displaylogo = FALSE
      ) %>%
      event_register("plotly_selected")
    fig
  })

  # 8. Save table
  observeEvent(r_state$current_idx, {
    req(nrow(r_state$filtered_cnvs) > 0)

    pname <- trimws(input$project_name %||% "")

    # Throw a warning if project name is not provided
    if (!nzchar(pname)) {
      showNotification("Please provide a project name before continuing.", type = "warning")
      return()
    }

    out_path <- file.path(wkdir, paste0(pname, ".tsv"))
    fwrite(r_state$filtered_cnvs, out_path, sep = "\t")
  }, ignoreNULL = TRUE)

  # 9. Fixed locus
  observeEvent(input$run_fixed_locus, {
    loc_name <- as.character(input$locus_name)
    loc_chr <- as.integer(input$locus_chr)
    loc_start <- as.integer(input$locus_start)
    loc_end <- as.integer(input$locus_end)
    min_overlap <- as.numeric(input$min_overlap)

    min_len <- as.integer(input$min_len_filter)
    max_len <- as.integer(input$max_len_filter)
    min_snp <- as.integer(input$min_snp_filter)

    # if something is not provided, raise a warning
    if (is.null(loc_name) || loc_name == '' ||
        is.null(loc_chr) || is.na(loc_chr) ||
        is.null(loc_start) || is.na(loc_start) ||
        is.null(loc_end) || is.na(loc_end) ||
        is.null(min_overlap) || is.na(min_overlap) ||
        min_overlap < 0 || min_overlap > 1) {
      showNotification("Locus info missing or wrong format.",
                       type = "warning")
      return()
    }

    # If all inputs are valid, select the putative carriers from the original CNV table
    cnvs <- r_state$cnvs

    # first apply filters if provided
    if (!is.na(min_len)) {
      cnvs <- cnvs[length >= min_len]
    }
    if (!is.na(max_len)) {
      cnvs <- cnvs[length <= max_len]
    }
    if (!is.na(min_snp)) {
      cnvs <- cnvs[numsnp >= min_snp]
    }
    if (!is.null(input$gt_filter) && input$gt_filter != 'both' &&
        !is.na(input$gt_filter)) {
      if (input$gt_filter == 'dels') {
        cnvs <- cnvs[GT == 1, ]
      } else if (input$gt_filter == 'dups') {
        cnvs <- cnvs[GT == 2, ]
      }
    }

    # then select the CNVs overlapping the locus and compute the overlap
    cnvs <- cnvs[chr == loc_chr & 
                 start <= loc_end & 
                 end >= loc_start, ]
    cnvs[, overlap := (pmin(end, loc_end) - pmax(start, loc_start) + 1) /
             (loc_end - loc_start + 1)]
    cnvs <- cnvs[overlap >= min_overlap, ]

    # assign the locus boundaries to the CNVs
    cnvs[, ':=' (locus_name = loc_name,
                 chr = loc_chr,
                 start = loc_start,
                 end = loc_end)]

    # Update filtered table and reset index
    cnvs <- unique(cnvs[, .(sample_ID, chr, start, end, locus_name, GT)])
    cnvs[, ':=' (numsnp = NA, CN = NA, length = NA, vo = -9)]
    r_state$filtered_cnvs <- cnvs
    r_state$current_idx <- 1
    r_state$fixed_locus <- TRUE
  })

  # 10. Select region
  observeEvent(input$run_select_region, {
    reg_chr <- as.integer(input$region_chr)
    reg_start <- as.integer(input$region_start)
    reg_end <- as.integer(input$region_end)
    min_iou <- as.numeric(input$min_iou)

    min_len <- as.integer(input$min_len_filter)
    max_len <- as.integer(input$max_len_filter)
    min_snp <- as.integer(input$min_snp_filter)

    # if something is not provided, raise a warning
    if (is.null(reg_chr) || is.na(reg_chr) ||
        is.null(reg_start) || is.na(reg_start) ||
        is.null(reg_end) || is.na(reg_end) ||
        is.null(min_iou) || is.na(min_iou) ||
        min_iou < 0 || min_iou > 1) {
      showNotification("Region info missing or wrong format.",
                       type = "warning")
      return()
    }

    # If all inputs are valid, select the putative carriers from the original CNV table
    cnvs <- r_state$cnvs

    # first apply filters if provided
    if (!is.na(min_len)) {
      cnvs <- cnvs[length >= min_len]
    }
    if (!is.na(max_len)) {
      cnvs <- cnvs[length <= max_len]
    }
    if (!is.na(min_snp)) {
      cnvs <- cnvs[numsnp >= min_snp]
    }
    if (!is.null(input$gt_filter) && input$gt_filter != 'both' &&
        !is.na(input$gt_filter)) {
      if (input$gt_filter == 'dels') {
        cnvs <- cnvs[GT == 1, ]
      } else if (input$gt_filter == 'dups') {
        cnvs <- cnvs[GT == 2, ]
      }
    }

    # then select the CNVs overlapping the locus and compute the overlap
    cnvs <- cnvs[chr == reg_chr & 
                 start <= reg_end & 
                 end >= reg_start, ]
    cnvs[, iou := (pmin(end, reg_end) - pmax(start, reg_start) + 1) /
                    (pmax(end, reg_end) - pmin(start, reg_start) + 1)]
    cnvs <- cnvs[iou >= min_iou, ]

    # Update filtered table and reset index
    r_state$filtered_cnvs <- cnvs
    r_state$current_idx <- 1
    r_state$select_region <- TRUE
  })

  # 11. Boundaries update
  output$selected <- renderPrint({
    d <- event_data("plotly_selected")
    if (is.null(d)) 'No points selected'
    else paste0('New start and stop positions: ',
                min(d$x), ' - ', max(d$x), '. ',
                'New number of SNPs: ', nrow(d))
  })

  observeEvent(input$btn_ref, {
    idx <- r_state$current_idx
    d <- event_data("plotly_selected")
    if (!is.null(d)) {
      # update start and end positions
      r_state$filtered_cnvs[idx, ':=' (start = min(d$x), end = max(d$x))]
      # update length and numsnp
      r_state$filtered_cnvs[idx, ':=' (length = end - start + 1, numsnp = nrow(d))]
      # retrigger plot and table rendering
      r_state$plot_retrigger <- r_state$plot_retrigger + 1
    }
  })

}

# Run the app ----
shinyApp(ui = ui, server = server)