# Reactive expression to create a vector of selected clusters
selected_clusters <- reactive({
  # Initialize an empty vector
  selected <- c()
  
  # Loop through each checkbox and check if it is selected
  for (i in 0:19) {
    if (input[[paste("cluster", i, sep = "")]]) {
      selected <- c(selected, i)
    }
  }
  
  # Return the vector of selected clusters
  selected
})
