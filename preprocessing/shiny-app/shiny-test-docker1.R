library(shiny)

# Define UI
ui <- fluidPage(
  titlePanel("Simple Shiny App"),
  sidebarLayout(
    sidebarPanel(
      textInput("name", "Enter your name:", ""),
      actionButton("greet", "Greet!")
    ),
    mainPanel(
      verbatimTextOutput("greeting")
    )
  )
)

# Define server
server <- function(input, output) {
  observeEvent(input$greet, {
    name <- input$name
    output$greeting <- renderText({
      if (name != "") {
        paste0("Hello, ", name, "!")
      } else {
        "Please enter your name."
      }
    })
  })
}

# Run the app
shinyApp(ui = ui, server = server)
