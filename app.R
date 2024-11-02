# Load necessary packages
library(shiny)
library(shinydashboard)
library(shinythemes)
library(XML)
library(openxlsx)
library(DT)
library(shinycssloaders)

# Define UI
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "HPA XML Parser"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Upload & Parse", tabName = "upload", icon = icon("upload")),
      menuItem("humanBrain Data", tabName = "brain", icon = icon("brain")),
      menuItem("Tissue Data", tabName = "tissue", icon = icon("dna")),
      menuItem("Download", tabName = "download", icon = icon("download")),
      menuItem("Select Theme", tabName = "theme", icon = icon("paint-brush"))
    )
  ),
  dashboardBody(
    tabItems(
      # Upload & Parse Tab
      tabItem(tabName = "upload",
              fluidRow(
                box(title = "Upload XML File", status = "primary", solidHeader = TRUE, width = 12,
                    fileInput("xmlfile", "Choose XML File",
                              accept = c(".xml")),
                    actionButton("parse", "Parse XML", icon = icon("play"), class = "btn-success")
                )
              ),
              fluidRow(
                box(title = "Instructions", status = "info", solidHeader = TRUE, width = 12,
                    "1. Upload your XML file by clicking the 'Choose XML File' button.",
                    br(),
                    "2. Click the 'Parse XML' button to process the file.",
                    br(),
                    "3. Navigate to the 'humanBrain Data' and 'Tissue Data' tabs to preview the data.",
                    br(),
                    "4. Go to the 'Download' tab to download the results as an Excel file."
                )
              )
      ),
      # humanBrain Data Tab
      tabItem(tabName = "brain",
              fluidRow(
                box(title = "humanBrain Data Preview", status = "primary", solidHeader = TRUE, width = 12,
                    DTOutput("brainData") %>% withSpinner(color = "#3c8dbc")
                )
              )
      ),
      # Tissue Data Tab
      tabItem(tabName = "tissue",
              fluidRow(
                box(title = "Tissue Data Preview", status = "primary", solidHeader = TRUE, width = 12,
                    DTOutput("tissueData") %>% withSpinner(color = "#3c8dbc")
                )
              )
      ),
      # Download Tab
      tabItem(tabName = "download",
              fluidRow(
                box(title = "Download Data", status = "primary", solidHeader = TRUE, width = 12,
                    downloadButton("downloadData", "Download Excel File", class = "btn-primary")
                )
              )
      ),
      # Theme Selector Tab
      tabItem(tabName = "theme",
              fluidRow(
                box(title = "Select Theme", status = "primary", solidHeader = TRUE, width = 12,
                    shinythemes::themeSelector()
                )
              )
      )
    )
  )
)

# Define Server logic
server <- function(input, output, session) {
  
  rv <- reactiveValues()
  
  # Observe when 'Parse XML' button is clicked
  observeEvent(input$parse, {
    req(input$xmlfile)
    xmlfile <- input$xmlfile$datapath
    document <- xmlTreeParse(xmlfile, useInternalNodes = TRUE)
    rootnode <- xmlRoot(document)
    
    showModal(modalDialog(
      title = "Parsing XML File",
      "Please wait while the file is being processed...",
      footer = NULL
    ))
    
    # Parse data for 'humanBrain' assayType
    brainData <- parse_assayType(rootnode, "humanBrain")
    
    # Parse data for 'tissue' assayType
    tissueData <- parse_assayType(rootnode, "tissue")
    
    # Save data to reactive values
    rv$brainData <- brainData
    rv$tissueData <- tissueData
    
    removeModal()
    
    # Render data tables for preview
    output$brainData <- renderDT({
      if (!is.null(rv$brainData)) {
        datatable(rv$brainData, options = list(pageLength = 10, scrollX = TRUE))
      } else {
        datatable(data.frame(Message = "No humanBrain data found."), options = list(dom = 't'))
      }
    })
    output$tissueData <- renderDT({
      if (!is.null(rv$tissueData)) {
        datatable(rv$tissueData, options = list(pageLength = 10, scrollX = TRUE))
      } else {
        datatable(data.frame(Message = "No tissue data found."), options = list(dom = 't'))
      }
    })
    
    # Download handler for Excel file
    output$downloadData <- downloadHandler(
      filename = function() {
        paste0("parsed_data_", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        wb <- createWorkbook()
        if (!is.null(rv$brainData)) {
          addWorksheet(wb, "humanBrain")
          writeData(wb, "humanBrain", rv$brainData)
        }
        if (!is.null(rv$tissueData)) {
          addWorksheet(wb, "tissue")
          writeData(wb, "tissue", rv$tissueData)
        }
        saveWorkbook(wb, file, overwrite = TRUE)
      }
    )
  })
  
  # Function to parse XML based on assayType
  parse_assayType <- function(rootnode, assayType) {
    # Get all 'rnaExpression' nodes with the specified assayType
    xpath_expr <- paste0("//rnaExpression[@assayType='", assayType, "']")
    rnaExpressions <- getNodeSet(rootnode, xpath_expr)
    
    data_list <- list()
    
    for (rnaExpr in rnaExpressions) {
      data_nodes <- getNodeSet(rnaExpr, ".//data")
      for (data_node in data_nodes) {
        tissue_node <- getNodeSet(data_node, ".//tissue")[[1]]
        region_text <- xmlValue(tissue_node)
        organ_attr <- xmlGetAttr(tissue_node, "organ")
        
        RNASamples <- getNodeSet(data_node, ".//RNASample")
        for (rnaSample in RNASamples) {
          attrs <- xmlAttrs(rnaSample)
          ls <- c(attrs, Region = region_text, Organ = organ_attr)
          data_list[[length(data_list) + 1]] <- ls
        }
      }
    }
    
    if (length(data_list) == 0) {
      return(NULL)
    }
    
    # Get all unique attribute names
    all_names <- unique(unlist(lapply(data_list, names)))
    
    # Ensure all data entries have the same columns
    data_list <- lapply(data_list, function(x) {
      x[setdiff(all_names, names(x))] <- NA
      x[all_names]
    })
    
    df <- do.call(rbind, lapply(data_list, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    # Set column names
    colnames(df) <- all_names
    return(df)
  }
  
}

# Run the application
shinyApp(ui, server)

