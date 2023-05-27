library(ggplot2)
library(shiny)
library(patchwork)


p1 <- plot_list$S6269Nr1 
p2 <- plot_list$S6269Nr3
p3 <- plot_list$S7788Nr1
p4 <- plot_list$S7788Nr2
p5 <- plot_list$S7788Nr3
p6 <- plot_list$S7788Nr11
p7 <- plot_list$S6230Nr3
p8 <- plot_list$S6230Nr4
p9 <- plot_list$S6269Nr2
p10 <- plot_list$S6269Nr4
p11 <- plot_list$S7788Nr15
p12 <- plot_list$S7788Nr16

# assuming p1, p2, p3, p4 are predefined ggplot objects

# Define the UI
ui <- fluidPage(
  titlePanel("Individual Plots Display"),
  sidebarLayout(
    sidebarPanel(
      numericInput("width_plot", "Width:", 600, min = 100, step = 10),
      numericInput("height_plot", "Height:", 600, min = 100, step = 10),
      selectInput("ncol", "Number of columns", choices = 1:4),
      checkboxInput("display_plot1", "Display Plot 1", value = TRUE),
      checkboxInput("display_plot2", "Display Plot 2", value = TRUE),
      checkboxInput("display_plot3", "Display Plot 3", value = TRUE),
      checkboxInput("display_plot4", "Display Plot 4", value = TRUE),
      checkboxInput("display_plot5", "Display Plot 5", value = TRUE),
      checkboxInput("display_plot6", "Display Plot 6", value = TRUE),
      checkboxInput("display_plot7", "Display Plot 7", value = TRUE),
      checkboxInput("display_plot8", "Display Plot 8", value = TRUE),
      checkboxInput("display_plot9", "Display Plot 9", value = TRUE),
      checkboxInput("display_plot10", "Display Plot 10", value = TRUE),
      checkboxInput("display_plot11", "Display Plot 11", value = TRUE),
      checkboxInput("display_plot12", "Display Plot 12", value = TRUE)
    ),
    mainPanel(
      plotOutput("combinedPlot", width = "100%", height = "100%")
    )
  )
)

# Define the server
server <- function(input, output) {
  
  output$combinedPlot <- renderPlot({
    plots_to_combine <- list()
    if(input$display_plot1) {plots_to_combine[[1]] <- p1}
    if(input$display_plot2) {plots_to_combine[[2]] <- p2}
    if(input$display_plot3) {plots_to_combine[[3]] <- p3}
    if(input$display_plot4) {plots_to_combine[[4]] <- p4}
    if(input$display_plot5) {plots_to_combine[[5]] <- p5}
    if(input$display_plot6) {plots_to_combine[[6]] <- p6}
    if(input$display_plot7) {plots_to_combine[[7]] <- p7}
    if(input$display_plot8) {plots_to_combine[[8]] <- p8}
    if(input$display_plot9) {plots_to_combine[[9]] <- p9}
    if(input$display_plot10) {plots_to_combine[[10]] <- p10}
    if(input$display_plot11) {plots_to_combine[[11]] <- p11}
    if(input$display_plot12) {plots_to_combine[[12]] <- p12}
    
    # Filter NULL values from the list
    plots_to_combine <- Filter(Negate(is.null), plots_to_combine)
    
    # If there are no plots to display, return NULL
    if(length(plots_to_combine) == 0) return(NULL)
    
    combined_plot <- wrap_plots(plots_to_combine, ncol = as.numeric(input$ncol))
    print(combined_plot)
  }, width = function() { as.numeric(input$width) }, height = function() { as.numeric(input$height) })
}

# Run the application 
shinyApp(ui = ui, server = server)
