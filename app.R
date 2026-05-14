library(shiny)
library(bslib)

ui <- page_sidebar(
  title = "RNA-seq Workflow Runner",
  sidebar = sidebar(
    width = 350,
    textInput("input_dir", "Input FASTQ directory", placeholder = "/data/reads"),
    textInput("output_dir", "Output directory", placeholder = "/data/results"),
    textInput("genome_file", "Genome FASTA file", placeholder = "/refs/ecoli.fa"),
    textInput("annotation_file", "Annotation file (GFF/GTF)", placeholder = "/refs/ecoli.gff"),
    selectInput("alignment_type", "Alignment type", choices = c("bowtie2", "salmon", "kallisto")),
    actionButton("run_btn", "Run workflow", class = "btn-primary"),
    hr(),
    actionButton("refresh_btn", "Refresh output status")
  ),
  layout_column_wrap(
    width = 1,
    card(
      card_header("Command preview"),
      verbatimTextOutput("cmd_preview")
    ),
    card(
      card_header("Run log"),
      verbatimTextOutput("run_log")
    ),
    card(
      card_header("Detected output files"),
      tableOutput("output_files")
    )
  )
)

server <- function(input, output, session) {
  workflow_path <- normalizePath("workflow.sh", mustWork = FALSE)

  log_text <- reactiveVal("Workflow is ready.\n")

  workflow_cmd <- reactive({
    req(input$input_dir, input$output_dir, input$genome_file, input$annotation_file)

    sprintf(
      "bash %s --INPUT %s --OUTPUT %s --GENOME_FILE %s --ANNOTATION_FILE %s --ALIGNMENT_TYPE %s",
      shQuote(workflow_path),
      shQuote(input$input_dir),
      shQuote(input$output_dir),
      shQuote(input$genome_file),
      shQuote(input$annotation_file),
      shQuote(input$alignment_type)
    )
  })

  output$cmd_preview <- renderText({
    if (
      nzchar(input$input_dir) &&
      nzchar(input$output_dir) &&
      nzchar(input$genome_file) &&
      nzchar(input$annotation_file)
    ) {
      workflow_cmd()
    } else {
      "Fill all required fields to preview command."
    }
  })

  observeEvent(input$run_btn, {
    if (!file.exists(workflow_path)) {
      showNotification("workflow.sh was not found in the app directory.", type = "error")
      return()
    }

    required_values <- c(input$input_dir, input$output_dir, input$genome_file, input$annotation_file)
    if (any(!nzchar(required_values))) {
      showNotification("Please complete all required fields.", type = "error")
      return()
    }

    cmd <- workflow_cmd()
    log_text(paste0(log_text(), "\n---\nRunning:\n", cmd, "\n\n"))

    result <- system2("bash", c("-lc", cmd), stdout = TRUE, stderr = TRUE)
    status <- attr(result, "status")
    if (is.null(status)) status <- 0

    log_text(
      paste0(
        log_text(),
        paste(result, collapse = "\n"),
        "\n\nExit status: ",
        status,
        "\n"
      )
    )

    if (status == 0) {
      showNotification("Workflow finished successfully.", type = "message")
    } else {
      showNotification("Workflow failed. Check the run log.", type = "error")
    }
  })

  output$run_log <- renderText({
    log_text()
  })

  output_files <- reactiveVal(data.frame(File = character(), stringsAsFactors = FALSE))

  observeEvent(input$refresh_btn, {
    if (!nzchar(input$output_dir) || !dir.exists(input$output_dir)) {
      output_files(data.frame(File = "Output directory does not exist yet.", stringsAsFactors = FALSE))
      return()
    }

    files <- list.files(input$output_dir, recursive = TRUE, full.names = FALSE)
    if (length(files) == 0) {
      output_files(data.frame(File = "No files found in output directory.", stringsAsFactors = FALSE))
    } else {
      output_files(data.frame(File = files, stringsAsFactors = FALSE))
    }
  })

  output$output_files <- renderTable({
    output_files()
  })
}

shinyApp(ui, server)
