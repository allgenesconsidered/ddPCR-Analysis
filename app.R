library(shiny)
library(ggplot2)
source('read_ddpcr_data.R')


ui <- fluidPage(
  titlePanel("ddPCR Analysis"),
  fluidRow(
    column(2, wellPanel(
      fileInput('plate_map', 'Input Plate Map (Optional)',
                accept=c('text/csv', 
                         'text/comma-separated-values,text/plain', 
                         '.csv'), multiple = T),
      fileInput('input_files', 'Choose CSV Files to Analyze',
                accept=c('text/csv', 
                         'text/comma-separated-values,text/plain', 
                         '.csv'), multiple = T),
      tags$hr(),
      h4('CSV Options'),
      checkboxInput('header', 'Header', TRUE),
      radioButtons('sep', 'Separator',
                   c(Comma=',',
                     Semicolon=';',
                     Tab='\t'),
                   ','),
      selectizeInput('k', 'K-means', choices = c(
        `4` = 4, `5` = 5, `6` = 6
      ), multiple = F), width = 2
    )),
    column(10,
           tabsetPanel(
           tabPanel( "Facet Graph",
                     downloadButton('download_graph', 'Download .PNG'),
                     plotOutput('plotFacet')
           ),
           tabPanel( "Data Table",
                     downloadButton('download_csv', 'Download .CSV'),
                     tableOutput('clean_csv')
           ))
    )
  )
)


server <- function(input, output) {
  
  dat <- reactiveValues(
  )
  
  
  #' Generic function for loading in the dataset. 
  loaded_data <- function(){
    inFile <- input$input_files
    if (is.null(inFile))
      return("Please enter data.") # inFile will be NULL initially
    dat$loaded_data = parse_files_from_list(inFile$datapath, k = as.numeric(input$k), fileNames = inFile$name)
    return(dat$loaded_data)
  }
  
  #' Function to output the facet plot. Data will be reloaded when
  #' new data is input into the two fileInput() UI functions. 
  #' TODO : Make the facet diplay height dynamic based on the number of
  #' graphs loaded (more col = smaller height, more rows = larger height).
  output$plotFacet <- renderPlot({
    loaded_data <- as.data.frame(loaded_data())
    if (!is.null(input$plate_map)){
      MAP = processPlatemap(input$plate_map$datapath)
      print(plotFacet(loaded_data, MAP))
    }
    print(plotFacet(loaded_data))
  }, height = 1000)
  
  clean_ddPCR_data <- function(){
    loaded_data <- as.data.frame(dat$loaded_data)
    MAP = processPlatemap(input$plate_map$datapath)
    outputPlate <- cleanPlate(loaded_data, MAP)
  }
  
  output$download_csv <- downloadHandler(
    filename = function() { paste('analyzed_ddPCR_data', '.csv', sep='') },
    content = function(file) {
      write.csv(clean_ddPCR_data(), file, row.names = F)}
  )
  
  output$clean_csv <- renderTable({
    clean_ddPCR_data()
  })
}

shinyApp(ui,server)




