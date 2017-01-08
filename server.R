library(shiny)
source('./read_ddpcr_data.R')

function(input, output) {
  output$plotFacet <- reactivePlot(function(){
    inFile <- input$input_files
    print(as.numeric(input$k))
    if (is.null(inFile))
      return("Please enter data.") # inFile will be NULL initially
    
    dat = parse_files_from_list(inFile$datapath, k = as.numeric(input$k), fileNames = inFile$name)
    plotFacet(dat)
  })
}

