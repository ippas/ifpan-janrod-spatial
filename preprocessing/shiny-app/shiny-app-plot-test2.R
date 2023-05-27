# create plots
library(ggplot2)
library(shiny)

p1 <- ggplot(mtcars, aes(mpg, cyl)) + geom_point()
p2 <- ggplot(mtcars, aes(hp, cyl)) + geom_point()
p3 <- ggplot(mtcars, aes(disp, cyl)) + geom_point()
p4 <- ggplot(mtcars, aes(mpg, disp)) + geom_point()

# Define the UI
ui <- fluidPage(
  titlePanel("Individual Plots Display"),
  sidebarLayout(
    sidebarPanel(
      checkboxInput("display_plot1", "Display Plot 1", value = TRUE),
      checkboxInput("display_plot2", "Display Plot 2", value = TRUE),
      checkboxInput("display_plot3", "Display Plot 3", value = TRUE),
      checkboxInput("display_plot4", "Display Plot 4", value = TRUE)
    ),
    mainPanel(
      uiOutput("plot1_ui"),
      uiOutput("plot2_ui"),
      uiOutput("plot3_ui"),
      uiOutput("plot4_ui")
    )
  )
)

# Define the server
server <- function(input, output) {
  
  output$plot1 <- renderPlot({ p1 })
  output$plot2 <- renderPlot({ p2 })
  output$plot3 <- renderPlot({ p3 })
  output$plot4 <- renderPlot({ p4 })
  
  output$plot1_ui <- renderUI({
    if (input$display_plot1) {
      fluidRow(column(6, plotOutput("plot1")))
    }
  })
  
  output$plot2_ui <- renderUI({
    if (input$display_plot2) {
      fluidRow(column(6, plotOutput("plot2")))
    }
  })
  
  output$plot3_ui <- renderUI({
    if (input$display_plot3) {
      fluidRow(column(6, plotOutput("plot3")))
    }
  })
  
  output$plot4_ui <- renderUI({
    if (input$display_plot4) {
      fluidRow(column(6, plotOutput("plot4")))
    }
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
