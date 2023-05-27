# Define the UI
ui <- fluidPage(
  titlePanel("Interactive mtcars plot"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("x", "Choose x variable:", choices = colnames(mtcars)),
      selectInput("y", "Choose y variable:", choices = colnames(mtcars)),
      selectInput("num_columns",
                  "Number of Columns to Display:",
                  choices = c(4, 8, 2))
    ),
    

    
    mainPanel(
      plotOutput("scatterPlot")
    )
  )
)

# Define the server
server <- function(input, output) {
  output$scatterPlot <- renderPlot({
    spatial_gene_plot(spatial_data = risperidone_st_data_half,
                      type_data = "raw_data",
                      gene = "Sgk1",
                      samples =  c(samples_saline, samples_risperidone),
                      min_percentile = 0.00,
                      max_percentile = 1,
                      size = 0.8,
                      ncol = as.numeric(input$num_columns),
                      normalization = T)
    
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

spatial_gene_plot(spatial_data = risperidone_st_data_half,
                  type_data = "raw_data",
                  gene = "Sgk1",
                  samples =  c(samples_saline, samples_risperidone),
                  min_percentile = 0.00,
                  max_percentile = 1,
                  size = 0.8,
                  ncol = 4,
                  normalization = T) %>% .[[1]]

p1


spatial_feature_plot(spatial_data = risperidone_st_data_half,
                     type_data = "quantile_normalize",
                     peak_id = "risperidone-peak-19962",
                     # samples =  c(samples_saline[-1], samples_risperidone[-3]),
                     samples = c(samples_saline, samples_risperidone),
                     min_percentile = 0,
                     max_percentile = 1,
                     normalization = T,
                     return_list = T) -> plot_list
  plot_layout(ncol = 4) 

p1[1]


p1 %>% as.list() -> p_list
print(p_list)
p1[[1,2,3]]

p_list <- list(sample1 = p1[[1]],
               sample2 = p1[[2]])

p_list
p1 %>% class

# Create the individual plots
plot1 <- ggplot(data = mtcars) + geom_point(aes(x = mpg, y = disp))
plot2 <- ggplot(data = iris) + geom_bar(aes(x = Species, fill = Species))

# Set names for each plot

names(plot1) <- "Plot 1"
names(plot2) <- "Plot 2"

my_patchwork <- wrap_plots(plot1, plot2, names = c("Plot 1", "Plot 2"))

# Create the patchwork object
my_patchwork <- plot1 + plot2

my_patchwork[["Plot 1"]]

my_patchwork[[which(names(my_patchwork) == "Plot 1")]]
