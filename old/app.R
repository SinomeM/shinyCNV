library(shiny)
library(data.table)
library(ggplot2)

plotcnv <- function(cnvs, loci, samples_file, snps = NA) {

  reg_len = 2000000
  i <- cnvs[1, ix]

  cnvi <- cnvs[ix == i, ]

  s <- cnvi[, sample_ID]
  loc <- cnvi[, locus]
  ll <- QCtreeCNV:::getline_locus(loci[locus == loc,])
  r_st <- as.integer(ll[3]) - reg_len
  r_en <- as.integer(ll[3]) + reg_len

  fp <- samples_file[sample_ID == s, file_path_tabix][1]
  if (length(fp) == 0) stop("Can't find sample ", s,
                            " in samples list file")
  raw <- QCtreeCNV:::get_region_tabix(ll[2], r_st, r_en, fp)
  colnames(raw) <- c("Chr", "Position","Position2", "LRR", "BAF", "LRRadj")
  if(!is.na(snps)) {
    raw <- merge(raw, snps, by = c("Chr", "Position"))
  }
  raw[between(as.integer(Position), as.integer(ll[3]), as.integer(ll[4])),
      core := T][is.na(core), core := F]

  m <- 1000000
  cc <- ifelse(cnvi[, GT] == 1, "red", "green")
  vv <- ifelse(cnvi[, GT] == 1, -1.5, 1.5)
  pl_lrr <-
    ggplot() +
    geom_point(data = raw, aes(Position/m, LRR, colour=core)) +
    scale_color_manual(values = c("black", "blue")) +
    geom_segment(aes(x=r_st/m, xend=r_en/m, y=0, yend=0),
                 linetype = "dashed", size=0.15, alpha=0.75) +
    xlab("Position (Mbp)") +
    guides(colour = "none") +
    geom_segment(data=cnvs[ix == i, ], aes(y=vv, yend=vv, x=start/m, xend=end/m),
                 colour=cc, arrow=arrow(angle=90, ends="both",
                                        length=unit(0.05, "inches"))) +
    ylim(-2, 2) +
    xlim(r_st/m, r_en/m) +
    theme_classic()
  pl_baf <-
    ggplot() +
    geom_point(data = raw, aes(Position/m, BAF, colour=core)) +
    scale_color_manual(values = c("black", "blue")) +
    geom_segment(aes(x=r_st/m, xend=r_en/m, y=0.5, yend=0.5),
                 linetype = "dashed", size=0.15, alpha=0.75) +
    xlab("Position (Mbp)") +
    guides(colour = "none") +
    ylim(0, 1) +
    xlim(r_st/m, r_en/m) +
    theme_classic()

  pl <- cowplot::plot_grid(pl_baf, pl_lrr, nrow = 2)

  pl

}

plotcnv(fread("data/put_cnv.txt")[, ix := 1:2], fread("data/loci.txt"),
        fread("data/samples_list.txt"))


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
               sidebarPanel(width = 6,
                        fileInput("fi_samples", "Select samples list file"),
                        fileInput("fi_loci", "Select loci list file"),
                        fileInput("fi_cnvs", "Select putative CNVs file"),
                        fileInput("fi_snps", "Select filtered SNPs file"),

                        actionButton("loadOK", "Proceed", class = "btn-success")),
               mainPanel(
                 verbatimTextOutput("sampload", placeholder = T),
                 verbatimTextOutput("lociload", placeholder = T),
                 verbatimTextOutput("cnvsload", placeholder = T),
                 verbatimTextOutput("snpsload", placeholder = T)
               )
      ),

      # Second panel
      tabPanel("Select Putative CNVs",
               sidebarPanel(width = 6,
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
                 tableOutput("cnvstb")
               )
      ),

      # Third panel
      tabPanel("Visual Inspection",
               # would be nice to have the buttons centered, use table instead of
               # row?
               sidebarPanel(width = 6,
                 fluidRow(
                   column(5,
                     actionButton("true", "True", class = "btn-success"),
                     actionButton("unk", "Unkwown"),
                     actionButton("prev", "Previous"),
                     downloadButton("save", "Save Results")),
                   column(5,
                     actionButton("false", "False", class = "btn-danger"),
                     actionButton("err", "Error", class = "btn-warning"),
                     actionButton("next", "Next"))
                 )),
               mainPanel(plotOutput("activeplot"))
      )))
)

server <- function(input, output) {

  # data load and checks
  samples <- reactive({
    if (!is.null(input$fi_samples)) {
      tmp <- input$fi_samples
      fread(tmp$datapath, header = T)
    }
  })

  loci <- reactive({
    if (!is.null(input$fi_loci)) {
      tmp <- input$fi_loci
      fread(tmp$datapath, header = T)
    }
  })

  cnvs <- reactive({
    if (!is.null(input$fi_cnvs)) {
      tmp <- input$fi_cnvs
      cnvs <- fread(tmp$datapath, header = T)
      if (!"vi" %in% colnames(cnvs)) cnvs[, vi := -9]
      if (!"ix" %in% colnames(cnvs)) cnvs[, ix := 1:nrow(cnvs)]
      cnvs
    }
  })

  snps <- reactive({
    if (!is.null(input$fi_snps)) {
      tmp <- input$fi_snps
      fread(tmp$datapath, header = T)
    }
  })

  output$sampload <- renderPrint({

    samples <- samples()

    if (is.null(samples)) "No samples file selected"
    else
      if (!all(c("sample_ID", "file_path_tabix")
               %in% colnames(samples)))
        "Samples list loaded but NOT in the correct format!!!"
      else "Samples list successfully loaded"
  })

  output$lociload <- renderPrint({

    loci <- loci()

    if (is.null(loci)) "No loci file selected"
    else
      if (!all(c("locus", "chr", "start", "end")
                %in% colnames(loci)) | "GT" %in% colnames(loci))
        "Loci list loaded but NOT in the correct format!!!"
      else "Loci list successfully loaded"
  })

  output$cnvsload <- renderPrint({
    cnvs <- cnvs()

    if (is.null(cnvs))
      "No CNVs file selected"
    else
      if (!all(c("sample_ID", "chr", "start", "end", "GT")
               %in% colnames(cnvs)))
        "Putative CNVs loaded but NOT in the correct format!!!"
      else {
        # update loci list
        updateSelectInput(inputId = "si_locus",
                          choices = unique(cnvs()$locus))
        "Putative CNVs successfully loaded"
      }

  })

  output$snpsload <- renderPrint({

    snps <- snps()

    if (is.null(snps)) "No SNPs file selected"
    else
      if (!all(c("chr", "postion") %in% colnames(snps)))
        "SNPs loaded but NOT in the correct format!!!"
      else "SNPs successfully loaded"
  })


  # filter putative CNVs

  cnvstb <- reactive({

    if (input$loadOK != 0) {

      cnvs <- cnvs()
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

      cnvs[GT %in% s_gt & locus %in% s_loc & vi %in% s_eval, 
             .(sample_ID, chr, start, end, GT, CN, locus, vi, ix)]
    }

  })

  output$cnvstb <- renderTable(cnvstb(), striped = T, hover = T, boarded = T)


  # produce plot and save inspection results

  activeplot <- reactive({

    if (input$filterOK != 0) {
      if (!is.null(cnvstb()) & !is.null(loci()) & !is.null(samples()))
        if (!is.null(snps())) plotcnv(cnvs(), loci(), samples(), snps())
        else plotcnv(cnvs(), loci(), samples())
    }

  })

  output$activeplot <- renderPlot(activeplot())

}

# Run the application
shinyApp(ui = ui, server = server)
