library(shiny)
library(plotly)

ui <- fluidPage(
  selectizeInput(
    inputId = "cities", 
    label = "Select a city", 
    choices = unique(txhousing$city), 
    selected = "Abilene",
    multiple = TRUE
  ),
  plotlyOutput(outputId = "p")
)

server <- function(input, output, session, ...) {
  output$p <- renderPlotly({
    req(input$cities)
    if (identical(input$cities, "")) return(NULL)
    p <- ggplot(data = filter(txhousing, city %in% input$cities)) + 
      geom_line(aes(date, median, group = city))
    height <- session$clientData$output_p_height
    width <- session$clientData$output_p_width
    ggplotly(p, height = height, width = width)
  })
}

shinyApp(ui, server)


# Relevant examples from https://plotly-r.com/linking-views-with-shiny.html
plotly_example("shiny", "drag_lines")
plotly_example("shiny", "event_data_parcoords")