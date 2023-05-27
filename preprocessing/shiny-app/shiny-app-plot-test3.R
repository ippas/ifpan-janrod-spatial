# create plots
library(ggplot2)
library(shiny)

p1 <- plot_list$S6269Nr1 
p2 <- plot_list$S6269Nr3
p3 <- plot_list$S7788Nr1
p4 <- plot_list$S7788Nr2
# Define the UI
ui <- fluidPage(
  titlePanel("Individual Plots Display"),
  sidebarLayout(
    sidebarPanel(
      selectInput("columns", "Number of Columns:", choices = 1:4),
      checkboxInput("display_plot1", "Display Plot 1", value = TRUE),
      checkboxInput("display_plot2", "Display Plot 2", value = TRUE),
      checkboxInput("display_plot3", "Display Plot 3", value = TRUE),
      checkboxInput("display_plot4", "Display Plot 4", value = TRUE)
    ),
    mainPanel(
      uiOutput("plots_ui", width = "1500px")
    )
  )
)
# Define the server
server <- function(input, output) {
  
  output$plot1 <- renderPlot({ p1 })
  output$plot2 <- renderPlot({ p2 })
  output$plot3 <- renderPlot({ p3 })
  output$plot4 <- renderPlot({ p4 })
  
  output$plots_ui <- renderUI({
    columns_per_plot <- 12 / as.numeric(input$columns) # calculate the number of Bootstrap columns per plot
    plots <- lapply(1:4, function(i) {
      if (input[[paste0('display_plot', i)]]) {
        column(columns_per_plot, plotOutput(paste0('plot', i), height = "350px", width = "380px"))
        # column(columns_per_plot, plotOutput(paste0('plot', i)))
      }
    })
    do.call(fluidRow, plots) # create a fluidRow with the desired number of columns per plot
  })
}


# Run the application 
shinyApp(ui = ui, server = server)

