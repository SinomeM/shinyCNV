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

# load data (CNVs, samples and SNPs) and initialise 'vo' column if necessary
cnvs <- fread(cargs[2])
if (!'vo' %in% colnames(cnvs)) cnvs[, vo := -9]
cnvs <- cnvs[, .(sample_ID, chr, start, end, numsnp, length, GT, vo)]
samples <- fread(cargs[3])[, .(sample_ID, file_path_tabix)]
snps <- fread(cargs[4])
setorder(cnvs, chr, start)


# UI function ----
ui <- fluidPage(
  titlePanel("CNV Visual Validation"),

  sidebarLayout(
    sidebarPanel(
      textOutput('prog'),
      br(),
      tableOutput('cnv_line'),
      br(),
      fluidRow(
        selectInput('vo_f', 'Filter CNV previous VI (select this first)',
                    c('new', 'true', 'false', 'unk', 'other', 'all'), 'all'),
        selectInput('gt_f', 'Filter CNV GT', c('dels', 'dups', 'both'), 'both')
      ),
      br(),
      br(),
      fluidRow(
        actionButton("true", "True", class = "btn-success"),
        actionButton("false", "False", class = "btn-danger"),
        actionButton("unk", "Unkown", class = 'btn-warning'),
        actionButton("err", "Error"),
        actionButton('ref', 'Boundary Refinement Needed')
      ),
      fluidRow(
        actionButton("prev", "Previous", class = "btn-default"),
        actionButton("nxt", "Next", class = "btn-default")
      ),
      br(),
      br(),
      br(),
      fluidRow(
        textAreaInput('project', 'Project Name', 'visual_inspection'),
        actionButton('save', 'Save results')
      ),
      br(),
      br(),
      fluidRow(
        sliderInput('pl_h', 'Change plot height',
                    min = 400, max = 1000, value = 750)
      )
    ),

    mainPanel(
      plotOutput("pl")
    )
  )
)

# Server function ----
server <- function(input, output, session) {

  # initialise dt, I guess it could be done better
  r_dt <- reactiveValues()
  r_dt$i <- 1
  r_dt$all_cnvs <- cnvs
  r_dt$cnvs <- cnvs
  r_dt$line <- cnvs[1]

  # Observers, filter the CNV table
  observeEvent(input$vo_f, {
    if (input$vo_f == 'true') r_dt$cnvs <- r_dt$all_cnvs[vo == 1, ]
    if (input$vo_f == 'false') r_dt$cnvs <- r_dt$all_cnvs[vo == 2, ]
    if (input$vo_f == 'unk') r_dt$cnvs <- r_dt$all_cnvs[vo == 3, ]
    if (input$vo_f == 'new') r_dt$cnvs <- r_dt$all_cnvs[vo == -9, ]
    if (input$vo_f == 'other') r_dt$cnvs <- r_dt$all_cnvs[!vo %in% c(1:3, -9), ]
    if (input$vo_f == 'all') r_dt$cnvs <- r_dt$all_cnvs
    r_dt$line <- r_dt$cnvs[r_dt$i]
  })

  observeEvent(input$gt_f, {
    gt <- ifelse(input$gt_f == 'dels', 1, 2)
    r_dt$cnvs <- r_dt$cnvs[GT == gt, ]
    r_dt$line <- r_dt$cnvs[r_dt$i]
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
            round(r_dt$i / r_dt$cnvs[, .N] * 100), '% completed.')
  })

  output$cnv_line <- renderTable({
    r_dt$line[, !c('sample_ID')]
  })

  output$pl <- renderPlot({
    # copied from CNValidatron nut now in ./R/
    plot_cnv(r_dt$line, samples[sample_ID == r_dt$line[, sample_ID], ],
             snps, tmp_plot = 1)
  }, width = function() input$pl_h * 1.1, height = function() input$pl_h)


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
