library(shiny)
library(DT)

ui <- fluidPage(
  actionButton("update_button", "Update Column"),
  DTOutput("data_table")
)

server <- function(input, output, session) {
  # Example dataset
  data <- reactiveValues(df = data.frame(
    ID = 1:5,
    Value = c(10, 20, 30, 40, 50)
  ))
  
  # Update column value on button press
  observeEvent(input$update_button, {
    # Update logic here, for example, doubling the 'Value' column
    data$df$Value <- data$df$Value * 2
  })
  
  # Render the updated dataset
  output$data_table <- renderDataTable({
    datatable(data$df)
  })
}

shinyApp(ui, server)