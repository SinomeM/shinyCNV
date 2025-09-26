library(shiny)
library(bslib)

ui <- fluidPage(
  theme = bs_theme(version = 5),
  layout_sidebar(
    sidebar = sidebar(
      title = "Controls",
      width = 280,
      open = "open",
      position = "left",
      selectInput("x", "X", choices = names(mtcars))
    ),
    div("Main content goes here")
  )
)

server <- function(input, output, session) {}
shinyApp(ui, server)
