# Define UI for your application
shinyUI(
ui <- fluidPage(
        titlePanel("HELLO REAL WORLD!!!!!"),
        sidebarLayout(
            sidebarPanel(
                # Add input widgets here
                sliderInput("slider", "Number of observations:", 1, 100, 50)
            ),
            mainPanel(
                # Add output elements here
                plotOutput("hist")
            )
        )
    )
)
