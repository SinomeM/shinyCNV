library(shiny)
library(data.table)

if (file.exists("tmp.txt")) file.remove("tmp.txt")
if (file.exists("tmp1.txt")) file.remove("tmp1.txt")


ui <- fluidPage(
  
  # Application title
  titlePanel("Shiny CNV Visual Inspection tool TESTING"),
  
  # multiple tabs, one for each major step:
  # data load, putative CNVs filtering, and Visual Inspection
  mainPanel(
    tabsetPanel(
      id = "sidetab",
      
      # First panel
      tabPanel("Import data",
               sidebarPanel(
                 fileInput("fi_samples", "Select samples list file"),
                 fileInput("fi_loci", "Select loci list file"),
                 fileInput("fi_cnvs", "Select putative CNVs file"),
                 fileInput("fi_snps", "Select filtered SNPs file"),
                 
                 actionButton("loadOK", "Proceed", class = "btn-success")
               ),
               mainPanel(
                 verbatimTextOutput("sampload", placeholder = T),
                 verbatimTextOutput("lociload", placeholder = T),
                 verbatimTextOutput("cnvsload", placeholder = T),
                 verbatimTextOutput("snpsload", placeholder = T)
               )
      ),
      
      # Second panel
      tabPanel("Select Putative CNVs",
               sidebarPanel(
                 selectInput("si_locus", "Select Locus", "None"),
                 selectInput("si_gt", "Select putative CNVs type", 
                             c("Duplications", "Deletions", "Both"),
                             selected = "Both"),
                 selectInput("si_eval", "Select putative CNVs state",
                             c("New", "True", "False", "Unkknown", "Error", "All"),
                             selected = "New"),
                 
                 actionButton("filterOK", "Proceed", class = "btn-success")
               ),
               mainPanel(
                 tableOutput("cnvs")
               )
      ),
      
      # Third panel
      tabPanel("Visual Inspection",
               # would be nice to have the buttons centered, use table instead of
               # row?
               sidebarPanel(
                 fluidRow(
                   actionButton("true", "True", class = "btn-success")
                 ),
                 fluidRow(
                   actionButton("false", "False", class = "btn-danger")
                 ),
                 fluidRow(
                   actionButton("unk", "Unkwown")
                 ),
                 fluidRow(
                   actionButton("err", "Error", class = "btn-warning")
                 ),
                 fluidRow(
                   # empty line???
                 ),
                 fluidRow(
                   actionButton("prev", "Previous")
                 ),
                 fluidRow(
                   actionButton("next", "Next")
                 ),
                 fluidRow(
                   actionButton("save", "Save Results")
                 )
               ),
               mainPanel()
      )
    ),
  )
)

server <- function(input, output) {
  
  # data load and checks
   
  # samples list
  output$sampload <- renderPrint({
    
    "No loci file selected"
    
    if (!is.null(input$fi_samples)) {
      tmp <- input$fi_samples
      if (!is.null((tmp))) {
        samples <- fread(tmp$datapath, header = T)
        if (!all(c("sample_ID", "file_path_tabix")
                 %in% colnames(samples)))
          "Samples list loaded but NOT in the correct format!!!"
        else "Samples list successfully loaded"
      }
    }
    
  })

  # loci list
  output$lociload <- renderPrint({
    
    "No loci file selected"
    
    if (!is.null(input$fi_loci)) {
      tmp <- input$fi_loci
      if (!is.null((tmp))) {
        loci <- fread(tmp$datapath, header = T)
        if (!all(c("locus", "chr", "start", "end")
                 %in% colnames(loci)) | 
            "GT" %in% colnames(loci))
          "Loci list loaded but NOT in the correct format!!!"
        else "Loci list successfully loaded"
      }
    }
    
  })
  
  # putative CNVs
  output$cnvsload <- renderPrint({
    
    "No CNVs file selected"
    
    if (!is.null(input$fi_cnvs)) {
      msg <- ""
      tmp <- input$fi_cnvs
      if (!is.null((tmp))) {
        cnvs <- fread(tmp$datapath, header = T)
        colnames(cnvs)
        if (!all(c("sample_ID", "chr", "start", "end", "GT")
                 %in% colnames(cnvs)))
          msg <- "Putative CNVs loaded but NOT in the correct format!!!"
        else {
          if (!"vi" %in% colnames(cnvs)) cnvs[, vi := -9]
          if (!"ix" %in% colnames(cnvs)) cnvs[, ix := 1:nrow(cnvs)]
          msg <- "Putative CNVs successfully loaded"
        
        # save CNV table as tmp.txt to be loaded in the second tab
        fwrite(cnvs, "tmp.txt")
        # return message to be printed
        }
        msg
      }
    }
    
  })
  
  # # SNPs list
  output$snpsload <- renderPrint({
    
    "No SNPs file selected"
    
    if (!is.null(input$fi_snps)) {
      tmp <- input$fi_snps
      if (!is.null((tmp))) {
        snps <- fread(tmp$datapath, header = T)
        colnames(snps)
        if (!all(c("chr", "postion") %in% colnames(snps)))
          "SNPs loaded but NOT in the correct format!!!"
        else "SNPs successfully loaded"
      }
    }
    
  })
  
  
  # filter putative CNVs
  
  # update loci list
  observeEvent(input$loadOK, {
    
    if (file.exists("tmp.txt")) {
      cnvs <- fread("tmp.txt")
      updateSelectInput(inputId = "si_locus", 
                        choices = unique(cnvs$locus))
    }
    
  })
  
  # update putative cnv list
  output$cnvs <- renderTable({
    
    cnvs <- fread("tmp.txt")
    s_loc <- input$si_locus
    
    if (input$si_gt == "Deletions") s_gt <- 1
    if (input$si_gt == "Duplications") s_gt <- 2
    if (input$si_gt == "Both") s_gt <- c(1,2)

    c("New", "True", "False", "Unkknown", "Error", "All")
    if (input$si_eval == "New") s_eval <- -9
    if (input$si_eval == "True") s_eval <- 1
    if (input$si_eval == "False") s_eval <- 2
    if (input$si_eval == "Unkknown") s_eval <- 3
    if (input$si_eval == "Error") s_eval <- 4
    if (input$si_eval == "All") s_eval <- 5
    
    cnvs <- cnvs[GT %in% s_gt & locus %in% s_loc & vi %in% s_eval,]
    
    # save CNV table as tmp.txt to be loaded in the third tab
    observeEvent(input$filterOK, {
      fwrite(cnvs, "tmp1.txt")
    })
    
    cnvs[, .(sample_ID, chr, start, end, GT, CN, locus, vi, ix)]
    
  }, striped = T, hover = T, boarded = T)


}

# Run the application 
shinyApp(ui = ui, server = server)
