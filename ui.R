library(shiny)

fluidPage(
  titlePanel("ddPCR Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput('input_files', 'Choose CSV File(s)',
                accept=c('text/csv', 
                         'text/comma-separated-values,text/plain', 
                         '.csv'), multiple = T),
      tags$hr(),
      checkboxInput('header', 'Header', TRUE),
      radioButtons('sep', 'Separator',
                   c(Comma=',',
                     Semicolon=';',
                     Tab='\t'),
                   ','),
      selectizeInput('k', 'K-means', choices = c(
        `2` = 2, `3` = 3, `4` = 4, `5` = 5, `6` = 6, `7` = 7
        ), multiple = F)
      ),
    mainPanel(
      plotOutput('plotFacet')
    )
  )
)


