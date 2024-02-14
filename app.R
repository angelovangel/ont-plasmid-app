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
library(shinyalert)
library(digest)
library(readxl)
library(digest)
library(reactable)
library(lubridate)


# Shiny app to run the ONT plasmid pipeline using pueue (BCL)

bin_on_path = function(bin) {
  exit_code = suppressWarnings(system2("command", args = c("-v", bin), stdout = FALSE))
  return(exit_code == 0)
}

sidebar <- sidebar(
  title = '',
  selectizeInput('pipeline', 'Assembly pipeline', choices = c('bacteria', 'plasmid'), selected = 'plasmid'),
  fileInput('upload', 'Upload sample sheet', multiple = F, accept = c('.xlsx', '.csv'), placeholder = 'xlsx or csv file'),
  shinyDirButton("fastq_folder", "Select fastq_pass folder", title ='Please select a fastq_pass folder from a run', multiple = F),
  checkboxInput('report', 'Faster html report', value = T),
  #tags$hr(),
  textInput('session_name', 'Name for new session', value = 'assembly', placeholder = 'to identify later'),
  actionButton('start', 'Start pipeline'),
  #tags$hr(),
  actionButton('show_session', 'Show session pane'),
  actionButton('ctrlc', 'Send ctrl-c to session'),
  actionButton('kill', 'Kill session'),
  downloadButton('log', 'Save log to file')
)

ui <- page_navbar(
  useShinyjs(),
  useShinyalert(),
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

server_fallback <- function(input, output, session) {
  shinybusy::report(
    title = 'App not configured', 
    text = 'Please start pueued service and make sure process-ontseq.sh is in your path', 
    type = 'failure')
}
server <- function(input, output, session) {
  # check if pueue and process-ont.sh are on path
  if (!bin_on_path('tmux')) {
    notify_failure('tmux not found', position = 'center-bottom')
  } else if (!bin_on_path('ont-plasmid.sh')) {
    notify_failure('ont-plasmid.sh not found', position = 'center-bottom')
  } else if (!bin_on_path('nextflow')) {
    notify_failure('nextflow not found', position = 'center-bottom')
  } else {
    notify_success('Server is ready', position = 'center-bottom')
  }
  
  
  empty_df <- data.frame(
    session_id = NA,
    started = NA,
    runtime = NA,
    command = NA,
    #active = NA,
    attached = NA,
    session_path = NA
  )
  
  # reactives
  tmux_sessions <- reactive({
    invalidateLater(2000, session)
    oldw <- getOption("warn")
    options(warn = -1)
    tmuxinfo <- system2("bin/helper.sh", stdout = TRUE, stderr = TRUE)
    options(warn = oldw)
    
    if (any(str_detect(tmuxinfo, 'no server|error'))) {
      empty_df
    } else {
      data.frame(
        session_id = str_split_i(tmuxinfo, " ", 2),
        started = str_split_i(tmuxinfo, " ", 1) %>% as.numeric() %>% as.POSIXct(),
        runtime = NA,
        command = str_split_i(tmuxinfo, " ", 5),
        #active = str_split_i(tmuxinfo, " ", 6),
        attached = str_split_i(tmuxinfo, " ", 3),
        session_path = str_split_i(tmuxinfo, " ", 4)
      ) %>%
        mutate(
          runtime = difftime(now(), started, units = 'hours'),
          attached = if_else(as.numeric(attached) == 1, 'yes', 'no')
        ) %>%
        arrange(started)
    }
  })

  samplesheet <- reactive({
    file <- input$upload
  })
  
  selected <- reactive({
    getReactableState('table', 'selected')
  })
  
  session_selected <- reactive({
    tmux_sessions()[selected(), ]$session_id
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
      selectedFolder <- parseDirPath(volumes, input$fastq_folder)
      nbarcodes <- length(list.files(path = selectedFolder, pattern = "barcode*", recursive = F))
      htmlreport <- if_else(input$report, '-r', '')
      
      arguments <- c('-p', selectedFolder, '-c', samplesheet()$datapath, '-w', input$pipeline, htmlreport)  
      
      #:) remove empty strings
      #arguments <- arguments[arguments != ""] 
      cat(
        'Selected folder:\n', selectedFolder, '\n', '-------\n\n',
        'Number of barcodes:\n', nbarcodes, '\n',  '-------\n\n',
        'Command:\n',
        'process-ontseq.sh', arguments)
    }
  })
  
  
  
  # observers
  # samplesheet
  observeEvent(input$upload, {
    ext <- tools::file_ext(samplesheet()$datapath)
    notify_info(paste0(ext, ' file uploaded'), position = 'center-bottom')
  })
  
  #main call
  observeEvent(input$start, {
    if (is.integer(input$fastq_folder)) {
      notify_info("No fastq folder selected", position = 'center-bottom')
    } else if (is.null(samplesheet()$datapath)) {
      notify_info("No samplesheet uploaded", position = 'center-bottom')
    } else {
      new_session_name <- paste0(digest::digest(runif(1), algo = 'crc32'), '-', input$session_name, '-', input$pipeline)
      selectedFolder <- parseDirPath(volumes, input$fastq_folder)
      htmlreport <- if_else(input$report, '-r', '')
      
      # launch new session
      args1 <- c('new', '-d', '-s', new_session_name)
      system2('tmux', args = args1)
      
      # execute pipeline in the new session
      string <- paste(
        'ont-plasmid.sh', 'Space', '-p', 'Space', selectedFolder, 'Space',  
        '-c', 'Space', samplesheet()$datapath, 'Space', '-w', 'Space', input$pipeline, 'Space', htmlreport, sep = ' '
      )
      args2 <- c('send-keys', '-t', new_session_name, string, 'C-m')
      system2('tmux', args = args2)
      notify_success(text = paste0('Started session ', new_session_name), position = 'center-bottom')
     
    }
  })
  
  observeEvent(input$show_session, {
    #observe({
    #session_selected <- tmux_sessions()[selected(), ]$session_id
    withCallingHandlers({
      shinyjs::html(id = "stdout", "")
      #args <- paste0(' a', ' -t ', session_selected)
      args <- c('capture-pane', '-S', '-', '-E', '-', '-pt', session_selected())
      
      p <- processx::run(
        'tmux', args = args,
        stdout_callback = function(line, proc) {message(line)},
        #stdout_line_callback = function(line, proc) {message(line)},
        stderr_to_stdout = TRUE,
        error_on_status = FALSE
      )
    },
    message = function(m) {
      shinyjs::html(id = "stdout", html = m$message, add = T);
      #runjs("document.getElementById('stdout').parentElement.scrollTo(0,1e9);")
      runjs("document.getElementById('stdout').parentElement.scrollTo({ top: 1e9, behavior: 'smooth' });")
    }
    )
  })
  
  # close session
  observeEvent(input$kill, {
    #session_selected <- tmux_sessions()[selected(), ]$session_id
    args <- paste0('kill-session -t ', session_selected())
    if (!is.null(selected())) {
      system2('tmux', args = args)
      notify_success(text = paste0('Session ', session_selected(), ' killed!'), position = 'center-bottom')
    } else{
      notify_failure('Select session first!', timeout = 2000, position = 'center-bottom')
    }  
  })
  
  # send ctrl-c
  observeEvent(input$ctrlc, {
    #session_selected <- tmux_sessions()[selected(), ]$session_id
    args <- paste0('send-keys -t ', session_selected(), ' C-c')
    if (!is.null(selected())) {
      system2('tmux', args = args)
      notify_success(text = paste0('Ctrl-C sent to session ', session_selected()), timeout = 2000, position = 'center-bottom')
    } else {
      notify_failure('Select session first!', timeout = 2000, position = 'center-bottom')
    }
    
  })
  
  observe({
    if(!is.null(selected())) {
      shinyjs::enable('log')
    } else {
      shinyjs::disable('log')
    }
  })
  
  #save log
  logdata <- reactive({
    args <- c('capture-pane', '-S', '-', '-E', '-', '-pt', session_selected())
    if (!is.null(selected())) {
      system2('tmux', args = args, stdout = T)
    } else {
      ''
    }  
  })
  
  output$log <- downloadHandler(
    filename = function() {
      paste0(session_selected(), '.log')
      },
    content = function(con) {
      if(!is.null(selected())) {
        writeLines(logdata(), con)
      } else {
        notify_failure('Select session first!', timeout = 2000, position = 'center-bottom')
      }
    }
  )
  
  # outputs
  output$table <- renderReactable({
    reactable(
      empty_df,
      #tmux_sessions(), 
      pagination = FALSE, highlight = TRUE, height = 200, compact = T, 
      fullWidth = T, selection = 'single', onClick = 'select', defaultSelected = 1,
      theme = reactableTheme(
        rowSelectedStyle = list(backgroundColor = "#eee", boxShadow = "inset 2px 0 0 0 #7B241C")
      ),
      style = list(fontSize = '90%'),
      columns = list(
        started = colDef(format = colFormat(datetime = T, locales = 'en-GB')),
        runtime = colDef(format = colFormat(suffix = ' h', digits = 2))
      )
    )
  })
  
  observe({
    updateReactable('table', data = tmux_sessions(), selected = selected())
  })
  
}


if (!bin_on_path('tmux') || !(bin_on_path('ont-plasmid.sh'))) {
  shinyApp(ui, server_fallback)
} else {
  shinyApp(ui, server)
}


