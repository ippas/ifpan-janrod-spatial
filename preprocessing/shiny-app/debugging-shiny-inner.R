library(shiny)
library(DT)
library(dplyr)


# Define UI
ui <- fluidPage(
  titlePanel("Data Table Example"),
  mainPanel(
    dataTableOutput("filter_statistics")
  )
)

# Define server logic
server <- function(input, output) {
  
  # Function to provide the filtered data (for demonstration, we are using the mtcars dataset here)
  filtered_data <- reactive({
    mtcars
  })
  
  # Rendering the DataTable
  output$filter_statistics <- renderDataTable({
    datatable_data <- filtered_data() %>%
      mutate_if(is.numeric, round, 2)
    
    datatable(datatable_data,
              options = list(paging = TRUE,
                             scrollX = TRUE,
                             scrollY = TRUE,
                             server = FALSE,
                             pageLength = 20,
                             buttons = c('csv', 'excel'),
                             columnDefs = list(list(targets = '_all', className = 'dt-center'))
              ),
              extensions = 'Buttons',
              filter = 'top',
              selection = 'single',
              rownames = FALSE
    ) 
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
