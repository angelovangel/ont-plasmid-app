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
library(tools)


# Shiny app to run the ONT plasmid pipeline using pueue (BCL)

bin_on_path = function(bin) {
  exit_code = suppressWarnings(system2("command", args = c("-v", bin), stdout = FALSE))
  return(exit_code == 0)
}

sidebar <- sidebar(
  title = '',
  selectizeInput('pipeline', 'Assembly pipeline', choices = c('amplicon', 'bacterial genome' = 'genome', 'plasmid'), selected = 'plasmid'),
  fileInput('upload', 'Upload sample sheet', multiple = F, accept = c('.xlsx', '.csv'), placeholder = 'xlsx or csv file'),
  shinyDirButton("fastq_folder", "Select fastq_pass folder", title ='Please select a fastq_pass folder from a run', multiple = F),
  #tags$hr(),
  textInput('session_name', 'Name for new session', value = 'assembly', placeholder = 'to identify later'),
  checkboxInput('report', 'Faster html report', value = T),
  checkboxInput('map_reads', 'Map reads to assembly', value = T),
  checkboxInput('transfer', 'Transfer results to Wahadrive', value = T),
  checkboxInput('singularity', 'Run with Singularity', value = F),
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
      max_height = '300px',
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
    text = 'Please make sure tmux and process-ontseq.sh is in your path', 
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
    # make start active later in stdout
    shinyjs::disable('start')
  }
  
  
  empty_df <- data.frame(
    session_id = NA,
    started = NA,
    runtime = NA
    #command = NA,
    #active = NA,
    #attached = NA
    #session_path = NA
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
        runtime = NA
        #command = str_split_i(tmuxinfo, " ", 5),
        #active = str_split_i(tmuxinfo, " ", 6),
        #attached = str_split_i(tmuxinfo, " ", 3)
        #session_path = str_split_i(tmuxinfo, " ", 4)
      ) %>%
        mutate(
          runtime = difftime(now(), started, units = 'hours')
          #attached = if_else(as.numeric(attached) == 1, 'yes', 'no')
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
    req(samplesheet())
    ext <- tools::file_ext(samplesheet()$datapath)
    validate(need(ext == 'csv' | ext == 'xlsx', 'Please upload a csv or excel (xlsx) file'))
    if (ext == 'csv') {
      # deal with samplesheets lacking complete final line, e.g. CRLF
      # read 2 times, first time to capture warning
      x <- tryCatch(
        read.csv(samplesheet()$datapath, header = T), 
        warning = function(w) {w}
      )
      if (inherits(x, 'simpleWarning')) {
        notify_warning(x$message, position = 'center-center', timeout = 5000)
        #notify_warning('This samplesheet will work but the last sample may be omitted', position = 'center-center', timeout = 5000)
        x <- read.csv(samplesheet()$datapath, header = T)
      }
      # check if sample and barcode columns are present
      if (sum(c('user', 'sample', 'dna_size', 'barcode') %in% colnames(x)) != 4) {
        notify_failure('Samplesheet must have columns "user", "sample", "dna_size" and "barcode"', position = 'center-center', timeout = 5000)
        shinyjs::disable('start')
      } else {
        notify_success('Samplesheet OK', position = 'center-center', timeout = 1000)
        shinyjs::enable('start') 
        # check that selected folder ends in fastq_pass
      }
      x
    } else if (ext == 'xlsx') {
      y <- read_excel(samplesheet()$datapath)
      if (sum(c('sample', 'barcode') %in% colnames(y)) != 2) {
        notify_failure('Samplesheet must have columns "user", "sample", "dna_size" and "barcode"', position = 'center-center', timeout = 5000)
        shinyjs::disable('start')
      } else {
        notify_success('Samplesheet OK', position = 'center-center', timeout = 3000)
        shinyjs::enable('start')
      }
      y
    }
  })
  
  # observers
  observeEvent(input$fastq_folder, {
    # start checking if something is selected, initially it is integer
    if (!is.integer(input$fastq_folder)) {
      path <- parseDirPath(volumes, input$fastq_folder)
      if (str_ends(path, 'fastq_pass')) {
        notify_success(path, position = 'center-center', timeout = 3000)
        shinyjs::enable('start')
      } else {
        notify_failure('Select a fastq_pass folder!', position = 'center-center', timeout = 3000)
        #shinyjs::disable('start')
      }
    }
  })
  
  # disable mapping if amplicon selected
  observeEvent(input$pipeline, {
    if (input$pipeline == 'amplicon') {
      updateCheckboxInput('map_reads', value = F, session = session)
      shinyjs::disable('map_reads')
    } else {
      updateCheckboxInput('map_reads', value = T, session = session)
      shinyjs::enable('map_reads')
    }
  })
  #main call
  observeEvent(input$start, {
    if (is.integer(input$fastq_folder)) {
      notify_info("No fastq folder selected", position = 'center-bottom')
    } else if (is.null(samplesheet()$datapath)) {
      notify_info("No samplesheet uploaded", position = 'center-bottom')
    } else {
      # new_session_name <- paste0(
      #   digest::digest(runif(1), algo = 'crc32'), '-', 
      #   stringi::stri_replace_all_charclass(input$session_name, "\\p{WHITE_SPACE}", ""), 
      #   '-', input$pipeline
      # )
      new_session_name <- paste0(
        Sys.time() %>% format('%Y%m%d-%H%M%S'), "-", input$pipeline, "-",
        stringi::stri_replace_all_charclass(input$session_name, "\\p{WHITE_SPACE}", "")
      )
        
      selectedFolder <- parseDirPath(volumes, input$fastq_folder)
      htmlreport <- if_else(input$report, '-r', '')
      singularity <- if_else(input$singularity, '-s', '')
      mapping <- if_else(input$map_reads, '-m', '')
      transfer <- if_else(input$transfer, '-t', '')
      
      # launch new session
      args1 <- c('new', '-d', '-s', new_session_name)
      system2('tmux', args = args1)
      
      # execute pipeline in the new session
      string <- paste(
        'ont-plasmid.sh', 'Space', '-p', 'Space', selectedFolder, 'Space',  
        '-c', 'Space', samplesheet()$datapath, 'Space', '-w', 'Space', input$pipeline, 'Space', '-n', 'Space', new_session_name, 'Space', 
        htmlreport, 'Space', mapping, 'Space', singularity, 'Space', transfer, sep = ' '
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


