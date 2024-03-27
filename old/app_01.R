#

library(shiny)

ui <- fluidPage(

    # Application title
    titlePanel("Shiny CNV Visual Inspection tool TESTING"),

    # multiple tabs, one for each major step:
    # data load, putative CNVs filtering, and Visual Inspection
    mainPanel(
      tabsetPanel(
        id = "sidetab",
        
        tabPanel("Import data",
                 sidebarPanel(
                   fileInput("samples", "Select samples list file"),
                   fileInput("loci", "Select loci list file"),
                   fileInput("cnvs", "Select putative CNVs file"),
                   fileInput("snps", "Select filtered SNPs file"),
                   
                   actionButton("dttest", "Check files", class = "btn-warning"),
                   actionButton("dtload", "Continue", class = "btn-success")
                 ),
                 mainPanel(
                   textOutput("dataload")
                 )
          ),
        
        tabPanel("Select Putative CNVs",
                 sidebarPanel(),
                 mainPanel()
          ),
        
        tabPanel("Visual Inspection",
                 sidebarPanel(),
                 mainPanel()
          )
        ),
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  
  # step 1, load and check main files
  msgsnps <- eventReactive(input$dttest, {

    if (is.null(input$snps)) snpmsg <- "No SNP file selected!"
    # else if (!all(c("Name", "Chr", "Position") 
    #               %in% colnames(input$snps))) 
    #        snpmsg <- "SNPs file in correct format"
    #   else snpmsg <- "prova"

  })
  
  msgother <- eventReactive(input$dttest, {

    # check all files formats
    lociok <- T; cnvsok <- T; sampleok <- T
    if (!all(c("locus", "chr", "start", "end") 
             %in% colnames(input$loci))) lociok <- F
    # stop("Some required columns are missing from loci object")
    if (!all(c("sample_ID", "chr", "start", "end", "CN") 
             %in% colnames(input$cnvs))) cnvsok <- F
    # stop("Some required columns are missing from CNV calls object")
    if (!all(c("sample_ID", "file_path", "file_path_tabix") 
             %in% colnames(input$samples))) sampleok <- F
    # stop("Some required columns are missing from samples list object")
    
    if (sum(lociok, cnvsok, sampleok) == 3)
      checkmsg <- paste0("All good! Click continue and proceed to ",
                         "CNVs selection.")
    else checkmsg <- "Something is wrong"

  })
  
  output$dataload <- renderText({
    # newline is being ignored
    c(msgsnps(), "ptova", msgother())
  }, sep = "\n")
  
  
  # step 2, select calls
  
  
  # step 3, visual inspection 


}

# Run the application 
shinyApp(ui = ui, server = server)
