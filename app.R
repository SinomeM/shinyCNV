library(data.table)
library(ggplot2)
library(shiny)
# library(CNValidatron)

cargs <- commandArgs(trailingOnly = T)
# working dir
# cnvs
# samples
# snps

source('./R/plotting_function.R')

# quickly check inputs
if (!dir.exists(cargs[1])) stop('Provided folder not found!')
if (!file.exists(cargs[2])) stop('CNVs table not found!')
if (!file.exists(cargs[3])) stop('Samples table not found!')
if (!file.exists(cargs[4])) stop('SNPs table not found!')

# load data (CNVs, samples and SNPs) and initialise 'vo' column if necessary
cnvs <- fread(cargs[2])
if (!'vo' %in% colnames(cnvs)) cnvs[, vo := -9]
cnvs <- cnvs[, .(sample_ID, chr, start, end, numsnp, length, GT, CN, vo)]
samples <- fread(cargs[3])[, .(sample_ID, file_path_tabix)]
snps <- fread(cargs[4])
setorder(cnvs, chr, start)


# UI function ----
ui <- fluidPage(
  titlePanel("CNVs Visual Validation"),

  sidebarLayout(
    sidebarPanel(
      textOutput('prog'),
      br(),
      fluidRow(
        textAreaInput('project', 'Project Name', 'visual_inspection'),
        sliderInput('pl_h', 'Change plot height',
                    min = 400, max = 1000, value = 750, step = 50),
        br(),
        actionButton("true", "True", class = "btn-success"),
        actionButton("false", "False", class = "btn-danger"),
        actionButton("unk", "Unkown", class = 'btn-warning'),
        actionButton("err", "Error"),
        actionButton('ref', 'Needs Boundary Ref')
        ),
      fluidRow(
        actionButton("prev", "Previous", class = "btn-default"),
        actionButton("nxt", "Next", class = "btn-default"),
        actionButton('save', 'Save results')
        ),
      br(),
      fluidRow(
        h5('Use the following fields to filter CNVs'),
        selectInput('vo_f', 'Filter CNV previous VI',
                    c('new', 'true', 'false', 'unk', 'other', 'all'), 'all'),
        selectInput('gt_f', 'Filter CNV GT', c('dels', 'dups', 'both'), 'both'),
        textInput("min_len", "Minimum CNV length", '0'),
        textInput("max_len", "Maximum CNV length", '10000000'),
        textInput("min_snp", "Minimum number of SNPs", '0'),
        h5('Fixed locus? Select input the details in the following fields'),
        textInput("locus", "Locus name", '0'),
        selectInput("loc_chr", "Locus chr", 0:24, '0'),
        textInput("loc_st", "Locus start", '0'),
        textInput("loc_en", "Locus end", '0'),
        sliderInput("loc_min_overlap", "Minimum overlap CNV over locus", min = 0,
                    max = 0.75, value = 0, step = 0.05),
        h5('After setting all required filters, click the "Run Filtering!" and wait for the updated CNV count'),
        actionButton('run', 'Run Filtering!', class = 'btn-success')
      )
    ),

    mainPanel(
      tableOutput('cnv_line'),
      plotOutput("pl"),
      textOutput('empty_dt')
    )
  )
)

# Server function ----
server <- function(input, output, session) {

  # initialise dt, I guess it could be done better
  r_dt <- reactiveValues()

  r_dt$i <- 1
  r_dt$cnvs <- cnvs
  r_dt$line <- cnvs[1]

  r_dt$vo <- c(1:4, -7, -9)
  r_dt$gt <- 1:2
  r_dt$min_len <- 0
  r_dt$min_snp <- 0

  r_dt$loc_st <- 0
  r_dt$loc_en <- 0
  r_dt$loc_chr <- 0
  r_dt$loc_min_overlap <- 0

  r_dt$empty_dt <- F

  # Observers, filter the CNV table
  observeEvent(input$vo_f, {
    if (input$vo_f == 'true')  r_dt$vo <- 1
    if (input$vo_f == 'false') r_dt$vo <- 2
    if (input$vo_f == 'unk')   r_dt$vo <- 3
    if (input$vo_f == 'new')   r_dt$vo <- -9
    if (input$vo_f == 'other') r_dt$vo <- c(4, -7)
    if (input$vo_f == 'all')   r_dt$vo <- c(1:4, -7, -9)
  })


  observeEvent(input$gt_f, {
    if (input$gt_f == 'dels') r_dt$gt <- 1
    if (input$gt_f == 'dups') r_dt$gt <- 2
    if (input$gt_f == 'both') r_dt$gt <- 1:2
  })

  observeEvent(input$min_len, {
    r_dt$min_len <- as.integer(input$min_len)
  })

  observeEvent(input$min_snp, {
    r_dt$min_snp <- as.integer(input$min_snp)
  })

  observeEvent(input$loc_chr, {
    r_dt$loc_chr <- as.integer(input$loc_chr)
  })

  observeEvent(input$loc_st, {
    r_dt$loc_st <- as.integer(input$loc_st)
  })

  observeEvent(input$loc_en, {
    r_dt$loc_en <- as.integer(input$loc_en)
  })


  observeEvent(input$run, {
    r_dt$cnvs <- cnvs
    #
    # update cnv table if a locus is selected
    if(r_dt$loc_chr != 0 & r_dt$loc_st != 0 & r_dt$loc_en != 0) {
      lloc <- data.table(locus = input$locus, chr = r_dt$loc_chr,
                         start = r_dt$loc_st, end = r_dt$loc_en)
      r_dt$cnvs <- QCtreeCNV::select_stitch_calls(r_dt$cnvs, lloc,
                                                  minoverlap = 0)
      if (nrow(r_dt$cnvs) != 0) {
        r_dt$cnvs[, ':=' (loc_st = lloc$start, loc_en = lloc$end)]
        r_dt$cnvs[, c('gap', 'stitch', 'densnp') := NULL]
        r_dt$cnvs <- r_dt$cnvs[overlap >= input$loc_min_overlap, ]
        #print(r_dt$cnvs)
      }
    }
    #
    # update cnv table and line
    if (nrow(r_dt$cnvs) == 0) {
      r_dt$cnvs <- cnvs[1]
      r_dt$empty_dt <- T
    }
    else {
      r_dt$cnvs <- r_dt$cnvs[GT %in% r_dt$gt & vo %in% r_dt$vo &
                             between(length, r_dt$min_len, r_dt$max_len)
                             & numsnp >= r_dt$min_snp, ]

      r_dt$line <- r_dt$cnvs[r_dt$i]
    }
  })


  # Observers, update the table and move index based on buttons press

  ## true, increase i and update 'vo' to 1
  observeEvent(input$true, {
    r_dt$cnvs[r_dt$i, vo := 1]
    r_dt$i <- r_dt$i + 1
    r_dt$line <- r_dt$cnvs[r_dt$i]
  })

  ## false, increase i and update 'vo' to 2
  observeEvent(input$false, {
    r_dt$cnvs[r_dt$i, vo := 2]
    r_dt$i <- r_dt$i + 1
    r_dt$line <- r_dt$cnvs[r_dt$i]
  })

  ## uknown, increase i and update 'vo' to 3
  observeEvent(input$unk, {
    r_dt$cnvs[r_dt$i, vo := 3]
    r_dt$i <- r_dt$i + 1
    r_dt$line <- r_dt$cnvs[r_dt$i]
  })

  ## refine, increase i and update 'vo' to 4
  observeEvent(input$ref, {
    r_dt$cnvs[r_dt$i, vo := 4]
    r_dt$i <- r_dt$i + 1
    r_dt$line <- r_dt$cnvs[r_dt$i]
  })

  ## error, increase i and update 'vo' to -7
  observeEvent(input$err, {
    r_dt$cnvs[r_dt$i, vo := -7]
    r_dt$i <- r_dt$i + 1
    r_dt$line <- r_dt$cnvs[r_dt$i]
  })

  ## previous, decrease i
  observeEvent(input$prev, {
    r_dt$i <- r_dt$i - 1
    r_dt$line <- r_dt$cnvs[r_dt$i]
  })

  ## next, increase i
  observeEvent(input$nxt, {
    r_dt$i <- r_dt$i + 1
    r_dt$line <- r_dt$cnvs[r_dt$i]
  })


  # Plot and table line for the current CNV

  output$prog <- renderText({
    paste0('CNV ', r_dt$i, ' out of ', r_dt$cnvs[, .N], '. ',
            round((r_dt$i - 1) / r_dt$cnvs[, .N] * 100), '% completed.')
  })

  output$cnv_line <- renderTable({
    r_dt$line[, ]
  })

  output$pl <- renderPlot({
    # copied from CNValidatron nut now in ./R/
    plot_cnv(r_dt$line, samples[sample_ID == r_dt$line[, sample_ID], ],
             snps)
  }, width = function() input$pl_h * 1.2, height = function() input$pl_h)

  output$empty_dt <- renderText({
    if (r_dt$empty_dt) paste0('NO CNVs match the selected filters!!')
    else ''
  })

  # Save results when asked

  observeEvent(input$save, {
    fwrite(r_dt$cnvs, paste0(cargs[1], '/', input$project, '_vi_res.txt'), sep = '\t')
  })

  # automatically save every 25 inspections
  observe({
    if (r_dt$i %% 25 == 0)
      fwrite(r_dt$cnvs, paste0(cargs[1], '/', input$project, '_vi_res.txt'), sep = '\t')
  })
}


# Run the app ----
shinyApp(ui = ui, server = server)
