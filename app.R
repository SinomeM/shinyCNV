library(data.table)
library(ggplot2)
library(shiny)

# all this could be passed with command line arguments
setwd('~/Documents/CNValidatron_fl/')
devtools::load_all()
setwd('~/Documents/shinyCNV/simpler_app/')

cnvs <- fread('~/Documents/CNValidatron_trained_models/third_iteration/selected_cnvs.txt')[,
          .(sample_ID, chr, start, end, numsnp, GT)]
cnvs[, vo := -9]
samples <- fread('~/Documents/CNValidatron_trained_models/third_iteration/samples.txt')[,
             .(sample_ID, file_path_tabix)]
snps <- fread('~/Documents/UKB_data/data/snppos_filtered.txt')
setorder(cnvs, chr, start)
# until here


# if using the full table I can add back the selector for GT and vo in the sidebar

# add a message to keep track of the progress

#
ui <- fluidPage(
  titlePanel("CNV Visual Validation"),
  
  sidebarLayout(
    sidebarPanel(
      # fluidRow(
      #   # for the moment only cnv file is provided and plots are already made
      #   # beforehand outside of the app
      #   #fileInput("fi_samples", "Select samples list file"),
      #   #fileInput("fi_loci", "Select loci list file"),
      #   fileInput("fi_cnvs", "Select putative CNVs file"),
      #   #fileInput("fi_snps", "Select filtered SNPs file"),
      # ),
      # br(),
      # br(),
      fluidRow(
        actionButton("true", "True", class = "btn-success"),
        actionButton("false", "False", class = "btn-danger"),
        actionButton("unk", "Unkown", class = 'btn-warning'),
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
                    min = 400, max = 1000, value = 500)
      )
    ),
    
    mainPanel(
      # plotOutput("imag", width = '30%', height = '20%'),
      tableOutput('cnv_line'),
      br(),
      plotOutput("pl")
    )
  )
)

#
server <- function(input, output, session) {
  
  # initialise dt
  r_dt <- reactiveValues()
  r_dt$i <- 1
  r_dt$cnvs <- cnvs
  r_dt$line <- cnvs[1]
  
  # cnvs_r <- reactive({
  #   if (!is.null(input$fi_cnvs)) {
  #     tmp <- input$fi_cnvs
  #     cnvs <- fread(tmp$datapath, header = T)
  #     if (!"vi" %in% colnames(cnvs)) cnvs[, vi := -9]
  #     if (!"ix" %in% colnames(cnvs)) cnvs[, ix := 1:nrow(cnvs)]
  #     cnvs
  #   }
  # })
  # at the moment CNV table is loaded "manually" (outside the proper app)
  
  # Observers, update the table and move index based on buttons press
  
  ## true, increase i and update 'vo' to 1
  observeEvent(input$true, {
    cnvs[r_dt$i, vo := 1]
    r_dt$i <- r_dt$i + 1
    r_dt$line <- cnvs[r_dt$i]
  })
  
  ## false, increase i and update 'vo' to 2
  observeEvent(input$false, {
    cnvs[r_dt$i, vo := 2]
    r_dt$i <- r_dt$i + 1
    r_dt$line <- cnvs[r_dt$i]
  })
  
  ## uknown, increase i and update 'vo' to 3
  observeEvent(input$unk, {
    cnvs[r_dt$i, vo := 3]
    r_dt$i <- r_dt$i + 1
    r_dt$line <- cnvs[r_dt$i]
  })
  
  ## refine, increase i and update 'vo' to 4
  observeEvent(input$ref, {
    cnvs[r_dt$i, vo := 4]
    r_dt$i <- r_dt$i + 1
    r_dt$line <- cnvs[r_dt$i]
  })
  
  ## previous, decrease i
  observeEvent(input$prev, {
    r_dt$i <- r_dt$i - 1
    r_dt$line <- cnvs[r_dt$i]
  })
  
  ## next, increase i
  observeEvent(input$nxt, {
    r_dt$i <- r_dt$i + 1
    r_dt$line <- cnvs[r_dt$i]
  })
  

  
  # Plot and table line for the current CNV
  
  output$imag <- renderImage({
    pt <- r_dt$line[, plot_id]
    list(src = pt)
  }, deleteFile = F)
  
  output$cnv_line <- renderTable({
    r_dt$line
  })
  
  output$pl <- renderPlot({
    # from CNValidatron
    plot_cnv(r_dt$line, samples[sample_ID == r_dt$line[, sample_ID], ],
             snps, tmp_plot = 1, shrink_lrr = 0.1)
  }, width = function() input$pl_h * 1.1, height = function() input$pl_h)
  
  # Save results when asked
  observeEvent(input$save, {
    fwrite(r_dt$cnvs, paste0('./', input$project, 'vi_res.txt'), sep = '\t')
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)
