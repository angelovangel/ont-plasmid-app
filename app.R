library(shiny)
library(shinyWidgets)
library(shinyjs)
library(bslib)
library(bsicons)
library(shinyjs)
library(tibble)
library(stringr)
library(dplyr)
library(processx)
library(shinyFiles)
library(shinybusy)
library(digest)
library(readxl)
library(digest)
library(reactable)


# Shiny app to run the ONT plasmid pipeline using pueue (BCL)

bin_on_path = function(bin) {
  exit_code = suppressWarnings(system2("command", args = c("-v", bin), stdout = FALSE))
  return(exit_code == 0)
}

sidebar <- sidebar(
  title = 'Inputs',
  selectizeInput('pipeline', 'Assembly pipeline', choices = c('bacteria', 'plasmid'), selected = 'plasmid'),
  fileInput('upload', 'Upload sample sheet', multiple = F, accept = c('.xlsx', '.csv'), placeholder = 'xlsx or csv file'),
  shinyDirButton("fastq_folder", "Select fastq_pass folder", title ='Please select a fastq_pass folder from a run', multiple = F),
  checkboxInput('report', 'Faster html report', value = T),
  tags$hr(),
  actionButton('start', 'Start pipeline'),
  tags$hr(),
  actionButton('attach', 'Attach session')
)

ui <- page_navbar(
  useShinyjs(),
  fillable = T,
  title = 'ONT plasmid assembly app',
  theme = bs_theme(bootswatch = 'yeti', primary = '#7B241C'),
  sidebar = sidebar,
  nav_panel(
    title = '',
    card(
      max_height = '250px',
      reactableOutput('table')
    ),
    card(
      max_height = '500px',
      verbatimTextOutput('stdout')
    )
  )
)

server <- function(input, output, session) {
  # check if pueue and process-ont.sh are on path
  if (!bin_on_path('pueue')) {
    notify_failure('pueue not found', position = 'center-bottom')
  } else if (!bin_on_path('process-ontseq.sh')) {
    notify_failure('process-ontseq.sh not found', position = 'center-bottom')
  } else {
    notify_success('pipeline is ready', position = 'center-bottom')
  }
  
  empty_df <- data.frame(
    id = NA,
    status = NA,
    command = NA,
    label = NA,
    start = NA
  )

  # reactives
  samplesheet <- reactive({
    file <- input$upload
  })
  
  selected <- reactive({
    getReactableState('table', 'selected')
  })
  
  pu_status <- reactive({
    invalidateLater(2000, session = session)
    j <- system2('pueue', args = c('status', '-j'), stdout = T)
    l <- jsonlite::fromJSON(j)
    df <- purrr::map_df(l$tasks, dplyr::bind_rows)
    if(nrow(df) == 0) {
      empty_df
    } else {
      df %>% 
        dplyr::select(id, status, command, label, start) %>% 
        unique() %>% 
        mutate(
          status = unlist(status),
          command = str_trunc(command, 70)
        )
    }
    #purrr::map_df(l$tasks, dplyr::bind_rows) %>% select(id, status, command, label, start, end) %>% unique() %>% mutate(status = unlist(status))
  })
  
  # dir choose management --------------------------------------
  default_path <- Sys.getenv('DEFAULT_PATH')
  volumes <- c(ont_data = default_path, getVolumes()())
  shinyDirChoose(input, "fastq_folder", 
                 roots = volumes,
                 session = session,
                 restrictions = system.file(package = "base")) 
  
  # build arguments for main call and display them on stdout at the same time
  output$stdout <- renderPrint({
    if (is.integer(input$fastq_folder)) {
      cat("No fastq folder selected\n")
    } else if (is.null(samplesheet()$datapath)) {
      cat("No samplesheet uploaded")
    } else {
      # hard set fastq folder and build arguments
      selectedFolder <<- parseDirPath(volumes, input$fastq_folder)
      nbarcodes <<- length(list.files(path = selectedFolder, pattern = "barcode*", recursive = F))
      htmlreport <- if_else(input$report, '-r', '')
      
      arguments <<- c('-p', selectedFolder, '-c', samplesheet()$datapath, htmlreport)  
      
      #:) remove empty strings
      #arguments <- arguments[arguments != ""] 
      cat(
        'Selected folder:\n', selectedFolder, '\n', '-------\n\n',
        'Number of barcodes:\n', nbarcodes, '\n',  '-------\n\n',
        'Command:\n',
        'process-ontseq.sh', arguments)
      
    }
  })
  
  # outputs
  output$table <- renderReactable({
    reactable(
      empty_df, pagination = FALSE, highlight = TRUE, height = 200, compact = T, 
      fullWidth = T, selection = 'single', onClick = 'select', defaultSelected = 1,
      theme = reactableTheme(
        rowSelectedStyle = list(backgroundColor = "#eee", boxShadow = " inset 2px 0 0 0 #7B241C")
      ),
      columns = list(
        id = colDef(minWidth = 25),
        status = colDef(minWidth = 40),
        command = colDef(minWidth = 200),
        label = colDef(minWidth = 40),
        start = colDef(minWidth = 60, format = colFormat(datetime = T, locales = 'en-GB'))
      )
    )
  })
  
  observe({
    updateReactable('table', data = pu_status(), selected = selected())
  })
  
}

shinyApp(ui, server)

