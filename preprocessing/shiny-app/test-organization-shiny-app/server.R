# Define server logic required to generate and plot a random distribution 
server <- shinyServer(function(input, output) {
    output$hist <- renderPlot({
        # Generate a histogram with the requested number of observations
        dist <- rnorm(input$slider)
        hist(dist)
    })
})

