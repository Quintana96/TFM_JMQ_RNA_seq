#' server_tab_config.R
#' Logica server de la Tab 1 (Configuracion):
#'   - Seleccion de directorio FASTQ (shinyFiles)
#'   - Deteccion de muestras (cached) y preview
#'   - Validacion (errores, checklist)
#'   - Boton "Continuar al procesamiento"
#'   - Carga de matriz de conteos externa (modo "load")
#'
#' Devuelve invisible: lista de reactivos expuestos para los otros tabs
#' (notablemente `samples_eff`, `effective_*` y `workflow_cmd`).

server_tab_config <- function(input, output, session, state) {
  output$workflow_path_text <- renderText({ state$workflow_path })

  # ── UI dinamica: label genoma/transcriptoma ────────────────────────────────
  output$genome_label_ui <- renderUI({
    if (isTRUE(input$analysis_type == "alignment"))
      fileInput("genome_file_upload", "Genoma de referencia (FASTA)",
                accept = c(".fa", ".fasta", ".fna", ".gz"), multiple = FALSE)
    else
      fileInput("genome_file_upload", "Transcriptoma de referencia (FASTA)",
                accept = c(".fa", ".fasta", ".fna", ".gz"), multiple = FALSE)
  })

  # ── Seleccion del directorio de FASTQ ──────────────────────────────────────
  #
  # El dialogo del sistema es el camino principal: se abre donde el usuario
  # quiera, sin raices ni arbol web que aprender. Cuando no se puede abrir
  # —dentro de un contenedor o sin escritorio— se cae al selector de shinyFiles,
  # que sigue funcionando aunque el servidor este en otra maquina.
  #
  # Lo que se muestra y lo que consume el resto de la pestana es SIEMPRE el
  # campo de texto: el dialogo se limita a rellenarlo. Con dos fuentes de verdad
  # (el dialogo y el campo) habria que decidir cual gana cada vez que difieren,
  # y esa decision no tiene una respuesta buena.
  usar_dialogo_nativo <- isTRUE(dialogo_nativo_disponible())

  output$input_dir_btn_ui <- renderUI({
    if (usar_dialogo_nativo) {
      actionButton("input_dir_nativo_btn", "Examinar...", icon = icon("folder-open"),
                   class = "btn-picker")
    } else {
      shinyFiles::shinyDirButton("input_dir_btn", "Examinar...", "Seleccionar...",
                                 class = "btn-picker")
    }
  })

  if (!usar_dialogo_nativo) {
    shinyFiles::shinyDirChoose(input, "input_dir_btn", roots = state$roots,
                               session = session)
    observeEvent(input$input_dir_btn, {
      ruta <- as.character(shinyFiles::parseDirPath(state$roots, input$input_dir_btn))
      if (length(ruta) && nzchar(ruta)) {
        updateTextInput(session, "input_dir_texto", value = ruta)
      }
    })
  }

  observeEvent(input$input_dir_nativo_btn, {
    # Se parte de lo que ya haya escrito para no obligar a volver a navegar
    # desde el principio al corregir una seleccion.
    actual <- validar_directorio(input$input_dir_texto)$ruta
    ruta <- elegir_directorio_nativo(
      titulo = "Selecciona la carpeta con los FASTQ",
      inicio = if (nzchar(actual)) actual else path.expand("~"))
    # NULL significa que se cancelo: se deja lo que hubiera, sin avisos.
    if (!is.null(ruta)) updateTextInput(session, "input_dir_texto", value = ruta)
  })

  input_dir_val <- reactive({
    validar_directorio(input$input_dir_texto)$ruta
  })

  output$input_dir_aviso <- renderUI({
    v <- validar_directorio(input$input_dir_texto)
    if (is.null(v$error)) return(NULL)
    div(class = "text-danger small mb-2", icon("triangle-exclamation"), " ", v$error)
  })

  output_dir_val <- reactive(state$pending_output_dir() %||% "")

  output$input_dir_path <- renderText({ input_dir_val() })
  output$output_base_path <- renderText({
    paste0(state$outputs_dir, "/<fecha_hora>_<tipo>_<herramienta>")
  })

  # ── Deteccion de muestras: debounce + cache compartido ─────────────────────
  input_dir_debounced <- debounce(reactive(input_dir_val()), 600)
  samples_eff <- reactive({
    detect_samples(input_dir_debounced(), input$read_type %||% "pe")
  })

  output$sample_preview_ui <- renderUI({
    dir_val <- input_dir_val() %||% ""
    if (!nzchar(dir_val))
      return(div(class = "text-muted", icon("folder-open"),
                 " Introduce un directorio de entrada."))
    if (!dir.exists(dir_val))
      return(div(class = "text-danger", icon("circle-xmark"), " El directorio no existe."))
    samples <- samples_eff()
    if (length(samples) == 0)
      return(div(class = "alert alert-warning", icon("triangle-exclamation"),
                 if (isTRUE(input$read_type == "se"))
                   " No se encontraron FASTQ single-end (*.fastq.gz / *.fastq)."
                 else
                   " No se encontraron FASTQ paired-end (*_1.fastq.gz / *_R1.fastq.gz)."))

    bad <- bad_sample_chars(samples)
    r2_miss <- if (isTRUE(input$read_type == "se")) rep(FALSE, length(samples)) else missing_r2(dir_val, samples)

    tagList(
      div(class = "alert alert-info py-2",
          icon("circle-check"), sprintf(" %d muestra(s) %s detectadas.",
                                        length(samples), read_type_label(input$read_type %||% "pe"))),
      if (length(bad) > 0)
        div(class = "alert alert-warning py-2", icon("triangle-exclamation"),
            " Nombres con caracteres especiales (espacios, @, etc.) pueden causar",
            " problemas: ", tags$b(paste(bad, collapse = ", "))),
      if (any(r2_miss))
        div(class = "alert alert-danger py-2", icon("circle-xmark"),
            " Faltan R2 para: ", tags$b(paste(samples[r2_miss], collapse = ", "))),
      tags$table(
        class = "table table-sm table-striped",
        tags$thead(tags$tr(tags$th("Muestra"), tags$th("Lectura 1 / FASTQ"), tags$th("R2"))),
        tags$tbody(lapply(seq_along(samples), function(i) {
          tags$tr(
            tags$td(samples[i]),
            tags$td(tags$span(style = "color:#315342;font-weight:700;", "✓")),
            tags$td(if (isTRUE(input$read_type == "se")) tags$span(class = "text-muted", "No aplica")
                    else if (!r2_miss[[i]]) tags$span(style = "color:#315342;font-weight:700;", "✓")
                    else tags$span(style = "color:#8A2F2F;font-weight:700;", "✗ falta"))
          )
        }))
      )
    )
  })

  # ── Validaciones ──────────────────────────────────────────────────────────
  val_errors <- reactive({
    errs <- character(0)
    dir_in  <- input_dir_val() %||% ""
    gf      <- if (!is.null(input$genome_file_upload) && nrow(input$genome_file_upload) > 0)
                  input$genome_file_upload$datapath else ""
    aln     <- input$analysis_type %||% "alignment"

    if (!nzchar(dir_in))           errs <- c(errs, "Directorio de FASTQs: campo vacio.")
    else if (!dir.exists(dir_in))  errs <- c(errs, "Directorio de FASTQs: no existe.")
    else {
      samps <- samples_eff()
      if (length(samps) == 0)
        errs <- c(errs, if (isTRUE(input$read_type == "se"))
                          "No se encontraron archivos FASTQ single-end."
                        else "No se encontraron archivos FASTQ R1.")
      else if (!isTRUE(input$read_type == "se")) {
        bad <- samps[missing_r2(dir_in, samps)]
        if (length(bad)) errs <- c(errs, paste("Faltan R2 para:", paste(bad, collapse = ", ")))
      }
    }
    if (!nzchar(gf))              errs <- c(errs, "FASTA de referencia: campo vacio.")
    else if (!file.exists(gf))    errs <- c(errs, "FASTA de referencia: no existe.")
    if (aln == "alignment") {
      af <- if (!is.null(input$annotation_file_upload) && nrow(input$annotation_file_upload) > 0)
        input$annotation_file_upload$datapath
      else ""
      if (!nzchar(af))            errs <- c(errs, "Anotacion GFF/GTF: campo vacio.")
      else if (!file.exists(af))  errs <- c(errs, "Anotacion GFF/GTF: no existe.")
    }
    if (isTRUE(input$analysis_type == "pseudo") &&
        identical(input$pseudo_tool %||% "salmon", "kallisto") &&
        isTRUE(input$read_type == "se")) {
      fl <- input$fragment_length %||% NA_real_
      fsd <- input$fragment_sd %||% NA_real_
      if (is.na(fl) || fl <= 0)
        errs <- c(errs, "Kallisto single-end: longitud media de fragmento invalida.")
      if (is.na(fsd) || fsd <= 0)
        errs <- c(errs, "Kallisto single-end: desviacion estandar de fragmento invalida.")
    }
    if (!dir.exists(state$outputs_dir))
      errs <- c(errs, paste0("No se pudo acceder a la carpeta de salidas: ", state$outputs_dir))
    if (!file.exists(state$workflow_path))
      errs <- c(errs, paste0("workflow.sh no encontrado en: ", state$workflow_path))
    errs
  })

  # ── Checklist (circulos rojo/verde/amarillo) ──────────────────────────────
  checklist_status <- reactive({
    dir_in <- input_dir_val() %||% ""
    dir_ok <- nzchar(dir_in) && dir.exists(dir_in)
    samples_ok <- FALSE
    if (dir_ok) {
      samps <- samples_eff()
      samples_ok <- length(samps) > 0 &&
        (isTRUE(input$read_type == "se") || !any(missing_r2(dir_in, samps)))
    }
    genome_ok <- if (!is.null(input$genome_file_upload) && nrow(input$genome_file_upload) > 0) "ok" else "missing"
    if (isTRUE(input$analysis_type == "alignment")) {
      annot_status <- if (!is.null(input$annotation_file_upload) && nrow(input$annotation_file_upload) > 0) "ok" else "missing"
    } else {
      annot_status <- if (!is.null(input$annotation_file_pseudo_upload) && nrow(input$annotation_file_pseudo_upload) > 0) "ok" else "optional"
    }
    valid_ok <- if (length(val_errors()) == 0) "ok" else "missing"
    list(
      "Directorio FASTQ"      = if (dir_ok) "ok" else "missing",
      "Muestras detectadas"   = if (samples_ok) "ok" else "missing",
      "Genoma/transcriptoma"  = genome_ok,
      "Anotacion"             = annot_status,
      "Validacion"            = valid_ok
    )
  })

  output$checklist_ui <- renderUI({
    items <- checklist_status()
    # El estado se transmite por texto Y por color, no solo por color: un
    # circulo de color no lo distingue quien tiene daltonismo ni lo anuncia un
    # lector de pantalla. El aria-label lo deja explicito.
    tags$ul(class = "checklist",
      lapply(names(items), function(nm) {
        st  <- as.character(items[[nm]])
        cfg <- CHECKLIST_STATES[[st]] %||%
          list(nivel = "neutral", simbolo = "—", etiqueta = "desconocido")
        tags$li(
          tags$span(nm),
          status_pill(cfg$simbolo, cfg$nivel, paste0(nm, ": ", cfg$etiqueta))
        )
      })
    )
  })

  output$validation_ui <- renderUI({
    errs <- val_errors()
    if (length(errs) == 0)
      div(class = "alert alert-success mb-0", icon("circle-check"),
          " Todos los campos correctos. Puedes continuar.")
    else
      div(class = "alert alert-danger mb-0",
          tags$b(icon("circle-xmark"), " Corrige antes de continuar:"),
          tags$ul(class = "mb-0 mt-1", lapply(errs, tags$li)))
  })

  observe({
    if (length(val_errors()) > 0) shinyjs::disable("btn_to_processing")
    else                           shinyjs::enable("btn_to_processing")
  })

  # ── Herramienta y anotacion efectivas (compartidas) ────────────────────────
  effective_tool <- reactive({
    if (isTRUE(input$analysis_type == "alignment")) "bowtie2"
    else input$pseudo_tool %||% "salmon"
  })
  effective_read_type <- reactive({
    if (isTRUE(input$read_type == "se")) "se" else "pe"
  })
  effective_fragment_length <- reactive({ input$fragment_length %||% 200 })
  effective_fragment_sd     <- reactive({ input$fragment_sd %||% 20 })

  effective_genome_file <- reactive({
    if (!is.null(input$genome_file_upload) && nrow(input$genome_file_upload) > 0)
      input$genome_file_upload$datapath
    else ""
  })

  effective_annotation <- reactive({
    if (isTRUE(input$analysis_type == "alignment")) {
      if (!is.null(input$annotation_file_upload) && nrow(input$annotation_file_upload) > 0)
        input$annotation_file_upload$datapath
      else ""
    } else {
      if (!is.null(input$annotation_file_pseudo_upload) && nrow(input$annotation_file_pseudo_upload) > 0)
        input$annotation_file_pseudo_upload$datapath
      else "/dev/null"
    }
  })

  # ── Comando del workflow (preview) ─────────────────────────────────────────
  workflow_cmd <- reactive({
    req(input_dir_val(), output_dir_val())
    sprintf(
      "bash %s --INPUT %s --OUTPUT %s --GENOME_FILE %s --ANNOTATION_FILE %s --ALIGNMENT_TYPE %s --READ_TYPE %s --FRAGMENT_LENGTH %s --FRAGMENT_SD %s",
      shQuote(state$workflow_path),
      shQuote(input_dir_val()), shQuote(output_dir_val()),
      shQuote(effective_genome_file()), shQuote(effective_annotation()), shQuote(effective_tool()),
      shQuote(effective_read_type()), shQuote(effective_fragment_length()), shQuote(effective_fragment_sd())
    )
  })

  # ── Navegar Tab 1 -> Tab 2 ────────────────────────────────────────────────
  observeEvent(input$btn_to_processing, {
    req(length(val_errors()) == 0)
    state$process_unlocked(TRUE)
    samps <- samples_eff()
    run_output_dir <- create_run_output_dir(state$outputs_dir, input$analysis_type, effective_tool())
    state$pending_output_dir(run_output_dir)
    state$config_snap(list(
      analysis_type = input$analysis_type, tool = effective_tool(),
      input_dir  = input_dir_val(), output_dir = run_output_dir,
      genome_file = effective_genome_file(), annotation = effective_annotation(),
      n_samples = length(samps), read_type = read_type_label(effective_read_type()),
      fragment_length = effective_fragment_length(), fragment_sd = effective_fragment_sd()
    ))
    state$log_text(paste(
      ts_log("=== Configuracion validada ==="),
      ts_log(paste("Analisis:", input$analysis_type, "/", effective_tool())),
      ts_log(paste("Muestras:", length(samps))),
      ts_log(paste("Entrada:", input_dir_val())),
      ts_log(paste("Salida:",  run_output_dir)), "",
      sep = "\n"
    ))
    updateNavbarPage(session, "main_nav", selected = "tab_process")
  })

  observeEvent(input$btn_back, {
    updateNavbarPage(session, "main_nav", selected = "tab_config")
  })

  # Mejora UX: boton de acceso directo desde el aviso de Tab 2 bloqueada
  observeEvent(input$btn_goto_config, {
    updateNavbarPage(session, "main_nav", selected = "tab_config")
  })

  # ── Carga de matriz de conteos externa ────────────────────────────────────
  observeEvent(input$btn_load_existing, {
    counts <- NULL
    loaded_dir <- ""
    loaded_tool <- "matriz subida"
    if (!is.null(input$upload_counts)) {
      counts <- tryCatch(
        read.table(input$upload_counts$datapath, header = TRUE, row.names = 1,
                   sep = "\t", comment.char = "#"),
        error = function(e)
          tryCatch(read.csv(input$upload_counts$datapath, row.names = 1),
                   error = function(e2) NULL)
      )
      loaded_dir <- dirname(input$upload_counts$datapath)
    }
    if (is.null(counts)) {
      showNotification(
        "No se pudo cargar la matriz de conteos. Verifica el formato o usa la pestana Resultados para abrir ejecuciones previas.",
        type = "error")
      return()
    }
    counts <- round(counts)
    state$data_rv$count_matrix <- counts
    state$data_rv$counts_ready <- TRUE
    state$data_rv$source       <- "uploaded"

    files <- if (nzchar(loaded_dir) && dir.exists(loaded_dir)) {
      list.files(loaded_dir, recursive = TRUE, full.names = FALSE)
    } else character(0)
    if (length(files)) {
      state$output_files_rv(file_table_for_files(loaded_dir, files))
    } else {
      upload_size <- if (!is.null(input$upload_counts)) {
        file.info(input$upload_counts$datapath)$size
      } else NA_real_
      state$output_files_rv(data.frame(
        Archivo = "Matriz cargada sin directorio de resultados.",
        `Tamano` = fmt_bytes(upload_size),
        stringsAsFactors = FALSE, check.names = FALSE
      ))
    }
    state$run_params_rv(list(
      analysis_type = if (identical(loaded_tool, "bowtie2")) "alignment" else "pseudo",
      tool = loaded_tool,
      input_dir = "resultados previos",
      output_dir = loaded_dir,
      genome_file = "—",
      annotation_file = "—",
      n_samples = ncol(counts),
      read_type = "Paired-end",
      started_at = Sys.time(),
      r_version = paste(R.version$major, R.version$minor, sep = ".")
    ))
    state$analysis_done(TRUE)
    showNotification(
      sprintf("Datos cargados: %d genes x %d muestras.",
              nrow(counts), ncol(counts)), type = "message", duration = 8
    )
    state$results_refresh(Sys.time())
    updateNavbarPage(session, "main_nav", selected = "tab_results")
  })

  # Exposicion de reactivos para los otros tabs (via state$shared)
  state$shared <- list(
    input_dir_val             = input_dir_val,
    output_dir_val            = output_dir_val,
    samples_eff               = samples_eff,
    val_errors                = val_errors,
    effective_tool            = effective_tool,
    effective_read_type       = effective_read_type,
    effective_fragment_length = effective_fragment_length,
    effective_fragment_sd     = effective_fragment_sd,
    effective_genome_file     = effective_genome_file,
    effective_annotation      = effective_annotation,
    workflow_cmd              = workflow_cmd
  )

  # ── Aviso de anotacion con splicing frente a bowtie2 (B10) ─────────────────
  # `bowtie2` no es splice-aware: sobre un procariota es correcto porque no hay
  # intrones, pero sobre un eucariota pierde las lecturas que cruzan uniones
  # exon-exon y subestima los conteos de forma sistematica. En lugar de preguntar
  # al usuario por el organismo, se deduce de la anotacion que acaba de cargar.
  splice_check <- reactive({
    if (!isTRUE(input$analysis_type == "alignment")) return(NULL)
    af <- effective_annotation()
    if (is.null(af) || !nzchar(af %||% "") || !file.exists(af)) return(NULL)
    detect_spliced_annotation(af)
  })

  output$splice_warning_ui <- renderUI({
    s <- splice_check()
    if (is.null(s) || !isTRUE(s$spliced)) return(NULL)
    div(class = "alert alert-warning py-2 px-3 small mb-2",
        icon("triangle-exclamation"),
        tags$b(" La anotacion describe genes con splicing y bowtie2 no es splice-aware."),
        tags$div(
          class = "mt-1",
          paste0("El ", round(100 * s$frac, 1), " % de los genes con exones anotados ",
                 "tiene mas de uno (", fmt_int(s$n_genes_multiexon), " de ",
                 fmt_int(s$n_genes_with_exons), "). bowtie2 no alinea lecturas que ",
                 "cruzan uniones exon-exon, asi que los conteos quedaran ",
                 "subestimados de forma sistematica, y mas en los genes con mas ",
                 "intrones.")),
        tags$div(class = "mt-1",
                 tags$b("Que hacer: "),
                 paste("usa la ruta de pseudoalineamiento (salmon o kallisto) con un",
                       "transcriptoma, que si maneja isoformas. Si necesitas",
                       "alineamiento sobre el genoma, hace falta un alineador",
                       "splice-aware como STAR o HISAT2, que este pipeline no",
                       "incluye.")))
  })

  # ── Calculadora de potencia a priori (item 22) ─────────────────────────────
  pw_params <- reactive({
    list(n = max(2, input$pw_n %||% 3),
         effect = max(1.01, input$pw_effect %||% 2),
         cv = max(0.01, input$pw_cv %||% 0.4),
         depth = max(1, input$pw_depth %||% 20))
  })

  # Estimacion de los parametros a partir de la matriz cargada. Pedir el CV y la
  # profundidad "a ojo" es la parte mas fragil del calculo, y la diferencia no es
  # cosmetica: sobre datos reales, los valores por defecto (cv 0,4, profundidad
  # 20) y los medidos sobre la misma matriz dan potencias muy distintas.
  pw_estimated <- reactiveVal(NULL)
  observeEvent(input$pw_estimate_btn, {
    cm <- state$data_rv$count_matrix
    if (is.null(cm) || !length(cm)) {
      showNotification(paste0(
        "No hay ninguna matriz de conteos cargada en esta sesion. Carga una en ",
        "esta pestana o analiza una ejecucion guardada, y vuelve a pulsar."),
        type = "warning", duration = 12)
      return()
    }
    withProgress(message = "Estimando parametros de los datos...", value = 0.4, {
      est <- estimate_power_params(cm, NULL)
    })
    if (is.null(est) || !is.finite(est$cv %||% NA_real_)) {
      showNotification("No se han podido estimar los parametros de estos datos.",
                       type = "error", duration = 10)
      return()
    }
    updateNumericInput(session, "pw_cv", value = round(est$cv, 3))
    updateNumericInput(session, "pw_depth", value = round(est$depth))
    if (!is.na(est$n_por_grupo %||% NA)) {
      updateNumericInput(session, "pw_n", value = est$n_por_grupo)
    }
    pw_estimated(est)
  })

  output$pw_estimate_note <- renderUI({
    e <- pw_estimated()
    if (is.null(e)) {
      return(HTML(paste("Con una matriz cargada, el CV biologico y la profundidad",
                        "se miden en lugar de suponerse.")))
    }
    HTML(paste0("Medidos sobre la matriz cargada: CV = ", round(e$cv, 3),
                " (BCV de edgeR), profundidad = ", round(e$depth),
                " conteos por gen."))
  })

  output$pw_verdict <- renderUI({
    p <- pw_params()
    r <- power_for_n(p$n, cv = p$cv, effect = p$effect, depth = p$depth)
    if (!is.null(r$error)) {
      return(div(class = "alert alert-secondary py-2 px-2 small", r$error))
    }
    i <- interpret_power(p$n, r$power)
    need <- n_for_power(0.8, cv = p$cv, effect = p$effect, depth = p$depth)
    cls <- switch(i$level, "ok" = "alert-success", "aviso" = "alert-warning",
                  "bajo" = "alert-warning", "alert-secondary")
    div(class = paste("alert py-2 px-3 mb-2", cls),
        tags$b(i$label),
        tags$div(class = "small mt-1", i$detail),
        if (!is.na(need$n)) tags$div(
          class = "small mt-1",
          tags$b("Para potencia 0,8 harian falta "), ceiling(need$n),
          " replicas por grupo con estos parametros.") else NULL)
  })

  output$pw_curve <- plotly::renderPlotly({
    p <- pw_params()
    df <- power_curve(2:24, cv = p$cv, effect = p$effect, depth = p$depth)
    if (is.null(df)) return(plotly_message("No se ha podido calcular la curva."))
    plotly::plot_ly(df, x = ~n, y = ~power, type = "scatter", mode = "lines+markers",
                    line = list(color = "#244B34"),
                    marker = list(color = "#7BBF9A", size = 7),
                    text = ~paste0("n = ", n, "<br>potencia: ", round(100 * power, 1), " %"),
                    hoverinfo = "text") |>
      plotly::layout(
        xaxis = list(title = "replicas por grupo"),
        yaxis = list(title = "potencia", range = c(0, 1), tickformat = ".0%"),
        shapes = list(
          list(type = "line", xref = "paper", x0 = 0, x1 = 1, y0 = 0.8, y1 = 0.8,
               line = list(dash = "dash", color = "#7BBF9A")),
          # Referencia empirica de Schurch et al.: 6 replicas como minimo
          # razonable, con independencia de lo que diga el calculo.
          list(type = "line", x0 = 6, x1 = 6, yref = "paper", y0 = 0, y1 = 1,
               line = list(dash = "dot", color = "#F4A6A6"))),
        annotations = list(
          list(x = 6, y = 1, yref = "paper", text = "6 replicas (Schurch 2016)",
               showarrow = FALSE, xanchor = "left", font = list(size = 10)))
      )
  })

  invisible(state$shared)
}
