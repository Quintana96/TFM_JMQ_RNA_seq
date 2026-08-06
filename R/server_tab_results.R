#' server_tab_results.R
#' Logica server de la Tab 3 (Resultados):
#'   - renderUI("tab3_content") segun analysis_done o ejecuciones guardadas
#'   - Selector de ejecucion, parametros inferidos, summary cacheado
#'   - Tablas DT, plots plotly, descargas (CSV y plotly_download)

server_tab_results <- function(input, output, session, state) {
  outputs_dir <- state$outputs_dir
  workflow_path <- state$workflow_path

  # ── Carpetas disponibles ──────────────────────────────────────────────────
  available_result_choices <- reactive({
    state$results_refresh()
    result_choices(outputs_dir)
  })

  selected_result_dir <- reactive({
    choices <- available_result_choices()
    selected <- input$selected_result_dir %||% ""
    if (nzchar(selected) && selected %in% unname(choices)) return(selected)

    current <- state$run_params_rv()$output_dir %||% ""
    if (nzchar(current) && current %in% unname(choices)) return(current)

    if (length(choices)) return(unname(choices)[1])

    if (state$analysis_done() && nzchar(current)) return(current)
    ""
  })

  selected_result_params <- reactive({
    out_dir <- selected_result_dir()
    current <- state$run_params_rv()
    current_dir <- current$output_dir %||% ""
    if (length(current) && nzchar(current_dir) && nzchar(out_dir) &&
        identical(normalizePath(current_dir, mustWork = FALSE),
                  normalizePath(out_dir, mustWork = FALSE))) {
      return(current)
    }
    if (nzchar(out_dir) && dir.exists(out_dir)) return(infer_result_params(out_dir, workflow_path))
    list()
  })

  selected_result_files <- reactive({
    file_table_for_dir(selected_result_dir())
  })

  # Cache: summary, counts y general_table.
  # Forma idiomatica: bindCache(reactive_obj, key1, key2, ...). Las claves se
  # evaluan como expresiones reactivas, por lo que pasar selected_result_dir()
  # y state$results_refresh() captura ambas dependencias sin ambigüedad de NSE.
  selected_result_summary <- reactive({
    p <- selected_result_params()
    out_dir <- selected_result_dir()
    if (!length(p) || !nzchar(out_dir)) return(NULL)
    summarise_result(out_dir, p)
  })
  if (exists("bindCache", mode = "function")) {
    selected_result_summary <- bindCache(
      selected_result_summary,
      selected_result_dir(),
      state$results_refresh()
    )
  }

  # ── rRNA por muestra ───────────────────────────────────────────────────────
  selected_rrna <- reactive({
    p <- selected_result_params()
    out_dir <- selected_result_dir()
    if (!length(p) || !nzchar(out_dir)) return(NULL)
    counts <- tryCatch(
      load_counts_from_workflow(out_dir, p$tool %||% "",
                                annotation_file = annotation_file_for_run(out_dir)),
      error = function(e) NULL)
    if (is.null(counts)) return(NULL)
    af <- annotation_file_for_run(out_dir)
    ids <- if (!is.null(af)) rrna_ids_from_annotation(af) else character(0)
    rrna_fraction_per_sample(counts, ids)
  })

  output$rrna_note <- renderUI({
    r <- selected_rrna()
    if (is.null(r)) return(NULL)
    if (is.null(r$table)) {
      return(div(class = "small text-muted",
                 paste("No se han podido identificar genes de rRNA. Los",
                       "identificadores de la matriz no llevan informacion de tipo",
                       "(p. ej. locus tags), asi que hace falta el fichero de",
                       "anotacion de la ejecucion.")))
    }
    tagList(
      div(class = "small text-muted mb-1",
          tags$b("Origen: "), r$source, "  ·  ", fmt_int(r$n_rrna_genes),
          " genes de rRNA  ·  variacion entre muestras: ",
          if (is.finite(r$spread)) paste0(round(100 * r$spread, 1), " puntos") else "—"),
      if (!is.null(r$alert)) div(class = "alert alert-warning py-2 px-2 small mb-1",
                                 icon("triangle-exclamation"), " ", r$alert) else NULL
    )
  })

  make_rrna_plot <- function(r) {
    if (is.null(r) || is.null(r$table)) {
      return(plotly_message("Sin datos de rRNA para esta ejecucion."))
    }
    df <- r$table
    df$pct <- 100 * df$frac
    plotly::plot_ly(df, x = ~sample_id, y = ~pct, type = "bar",
                    marker = list(color = "#7BBF9A"),
                    text = ~paste0("Muestra: ", sample_id,
                                   "<br>rRNA: ", round(pct, 2), " %",
                                   "<br>lecturas en rRNA: ", fmt_int(rrna_reads),
                                   "<br>total asignadas: ", fmt_int(total_reads)),
                    hoverinfo = "text") |>
      plotly::layout(xaxis = list(title = ""),
                     yaxis = list(title = "% de lecturas asignadas a rRNA"))
  }

  output$rrna_plot <- plotly::renderPlotly(make_rrna_plot(selected_rrna()))
  output$download_rrna_plot <- plotly_download(
    "rrna_por_muestra", function() make_rrna_plot(selected_rrna()))

  selected_counts_tables <- reactive({
    p <- selected_result_params()
    out_dir <- selected_result_dir()
    if (!length(p) || !nzchar(out_dir)) {
      return(list(libs = data.frame(Mensaje = "Sin ejecucion seleccionada."),
                  top = data.frame(Mensaje = "Sin ejecucion seleccionada.")))
    }
    counts_tables(out_dir, p$tool %||% "")
  })
  if (exists("bindCache", mode = "function")) {
    selected_counts_tables <- bindCache(
      selected_counts_tables,
      selected_result_dir(),
      state$results_refresh()
    )
  }

  selected_general_table <- reactive({
    result_general_table(selected_result_dir())
  })
  if (exists("bindCache", mode = "function")) {
    selected_general_table <- bindCache(
      selected_general_table,
      selected_result_dir(),
      state$results_refresh()
    )
  }

  # ── MultiQC link ──────────────────────────────────────────────────────────
  selected_multiqc_href <- reactive({
    out_dir <- selected_result_dir()
    report <- file.path(out_dir, "multiqc_report.html")
    if (!nzchar(out_dir) || !file.exists(report)) return("")

    base <- normalizePath(outputs_dir, winslash = "/", mustWork = TRUE)
    target <- normalizePath(report, winslash = "/", mustWork = TRUE)
    if (!startsWith(target, paste0(base, "/"))) return("")

    rel <- substring(target, nchar(base) + 2L)
    rel_parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
    paste0("saved_outputs/", paste(vapply(rel_parts, URLencode, character(1), reserved = TRUE), collapse = "/"))
  })

  output$result_selector_ui <- renderUI({
    choices <- available_result_choices()
    if (!length(choices)) {
      return(div(class = "alert alert-info mb-0", icon("folder-open"),
                 " Aun no hay ejecuciones guardadas en outputs/."))
    }
    tagList(
      div(style = "display:flex;gap:10px;align-items:flex-end;width:100%;",
        div(style = "flex:1 1 auto;min-width:0;",
          selectizeInput(
            "selected_result_dir",
            "Ejecucion guardada",
            choices = choices,
            selected = selected_result_dir(),
            options = list(placeholder = 'Escribe para buscar...'),
            width = "100%"
          )
        ),
        div(style = "flex:0 0 auto;margin-bottom:16px;",
          shinyDirButton("select_result_dir_btn", "Seleccionar carpeta...", "Seleccionar...")
        )
      ),
      tags$div(style = "margin-top:6px;",
               tags$small(class = "text-muted",
                          "También puedes escribir para filtrar las ejecuciones guardadas."))
    )
  })

  output$multiqc_open_ui <- renderUI({
    href <- selected_multiqc_href()
    if (!nzchar(href)) {
      return(tags$small(class = "text-muted", "MultiQC no disponible para esta ejecucion."))
    }
    tags$a(
      href = href,
      target = "_blank",
      rel = "noopener noreferrer",
      class = "btn btn-sm btn-primary",
      style = "align-self:flex-start;",
      title = "Abrir el informe MultiQC de esta ejecución en una nueva pestaña",
      tagList(icon("up-right-from-square"), " Abrir MultiQC")
    )
  })

  # ── shinyFiles: seleccionar carpeta de outputs ────────────────────────────
  shinyFiles::shinyDirChoose(input, "select_result_dir_btn",
                             roots = state$roots_results, session = session)

  observeEvent(input$refresh_results_btn, {
    state$results_refresh(Sys.time())
    choices <- result_choices(outputs_dir)
    updateSelectInput(session, "selected_result_dir",
                      choices = choices,
                      selected = if (length(choices)) unname(choices)[1] else character(0))
  })

  observeEvent(input$select_result_dir_btn, {
    req(input$select_result_dir_btn)
    sel <- tryCatch(shinyFiles::parseDirPath(state$roots_results, input$select_result_dir_btn),
                    error = function(e) NULL)
    if (!is.null(sel)) {
      sel <- as.character(sel)
      if (nzchar(sel) && dir.exists(sel)) {
        state$results_refresh(Sys.time())
        updateSelectInput(session, "selected_result_dir",
                          choices = result_choices(outputs_dir), selected = sel)
      }
    }
  })

  # ── Contenido Tab 3 ───────────────────────────────────────────────────────
  output$tab3_content <- renderUI({
    if (!state$analysis_done() && !length(available_result_choices()))
      return(ui_tab_results_empty())
    p <- selected_result_params()
    if (!length(p)) return(NULL)
    ui_tab_results_content(selected_result_summary())
  })

  # ── Tablas DT ─────────────────────────────────────────────────────────────
  output$output_files_table <- renderDT({
    dt_table(selected_result_files(), page_length = 20)
  })

  output$result_interpretation_ui <- renderUI({
    s <- selected_result_summary()
    p <- selected_result_params()
    out_dir <- selected_result_dir()
    if (is.null(s) || !length(p)) return(NULL)

    items <- list(
      tags$li(tags$b("Estado: "), status_badge(s$status)),
      tags$li(tags$b("Carpeta: "), tags$code(out_dir)),
      tags$li(tags$b("Tamano total: "), s$total_size),
      tags$li(tags$b("MultiQC: "), if (isTRUE(s$has_multiqc)) "disponible" else "no encontrado")
    )

    if (!is.na(s$mean_mapped)) {
      items <- c(items, list(tags$li(tags$b("Mapeo medio: "), pct_label(s$mean_mapped),
                                if (s$mean_mapped < 70) tags$span(class = "text-danger", "  Revisar: bajo para RNA-seq bacteriano.") else NULL)))
    }
    if (!is.na(s$mean_trim_survival)) {
      items <- c(items, list(tags$li(tags$b("Lecturas retenidas tras fastp: "), pct_label(s$mean_trim_survival),
                                if (s$mean_trim_survival < 80) tags$span(class = "text-warning", "  Perdida alta durante trimming.") else NULL)))
    }
    if (!is.na(s$mean_q30)) {
      items <- c(items, list(tags$li(tags$b("Q30 post-trimming: "), pct_label(100 * s$mean_q30),
                                if (s$mean_q30 < 0.85) tags$span(class = "text-warning", "  Calidad post-trimming moderada/baja.") else NULL)))
    }
    if (!is.na(s$fastqc_fail)) {
      items <- c(items, list(tags$li(tags$b("Checks FastQC fallidos: "), s$fastqc_fail,
                                if (s$fastqc_fail > 0) tags$span(class = "text-warning", "  Revisar pestaña Calidad.") else NULL)))
    }

    tagList(
      tags$ul(class = "mb-0", items),
      if (file.exists(file.path(out_dir, "multiqc_report.html")))
        div(class = "alert alert-info mt-3 mb-0",
            icon("circle-info"),
            " Abre multiqc_report.html desde la carpeta de resultados para el informe interactivo completo.")
    )
  })

  output$run_stats_table <- renderDT({
    dt_table(selected_general_table())
  })

  output$fastqc_table <- renderDT({
    dt_table(fastqc_table(selected_result_dir()), page_length = 12)
  })

  output$alignment_table <- renderDT({
    p <- selected_result_params()
    dt_table(alignment_table(selected_result_dir(), p$tool %||% ""), page_length = 12)
  })

  output$count_lib_table <- renderDT({
    dt_table(selected_counts_tables()$libs, page_length = 12)
  })

  output$count_top_table <- renderDT({
    dt_table(selected_counts_tables()$top, page_length = 15)
  })

  output$artifact_table <- renderDT({
    dt_table(important_artifacts(selected_result_dir()), page_length = 12)
  })

  output$selected_log_tail <- renderText({
    log_tail_text(selected_result_dir())
  })

  # ── Descargas CSV ─────────────────────────────────────────────────────────
  output$download_run_stats <- csv_download(
    "multiqc_stats",
    function() selected_general_table()
  )
  output$download_fastqc <- csv_download(
    "fastqc",
    function() fastqc_table(selected_result_dir())
  )
  output$download_alignment <- csv_download(
    "alignment_metrics",
    function() {
      p <- selected_result_params()
      alignment_table(selected_result_dir(), p$tool %||% "")
    }
  )
  output$download_count_lib <- csv_download(
    "count_libraries",
    function() selected_counts_tables()$libs
  )
  output$download_count_top <- csv_download(
    "top_counts",
    function() selected_counts_tables()$top
  )
  output$download_artifacts <- csv_download(
    "artifacts",
    function() important_artifacts(selected_result_dir())
  )

  output$download_log <- downloadHandler(
    filename = function() paste0("rnaseq_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"),
    content = function(f) {
      log_file <- file.path(selected_result_dir(), "workflow_live.log")
      if (file.exists(log_file)) file.copy(log_file, f, overwrite = TRUE)
      else writeLines(state$log_text(), f)
    })
  output$download_filelist <- downloadHandler(
    filename = function() paste0("output_files_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(f) write.csv(selected_result_files(), f, row.names = FALSE))

  # ── QC adicional: tablas ──────────────────────────────────────────────────
  output$align_qc_summary_table  <- renderDT({ dt_table(align_qc_summary(selected_result_dir())) })
  output$align_qc_alerts_table   <- renderDT({ dt_table(align_qc_alerts(selected_result_dir())) })
  output$pseudo_qc_summary_table <- renderDT({ dt_table(pseudo_qc_summary(selected_result_dir())) })
  output$pseudo_qc_quant_table   <- renderDT({ dt_table(pseudo_qc_quant_table(selected_result_dir()), filter = "top") })
  output$pseudo_qc_alerts_table  <- renderDT({ dt_table(pseudo_qc_alerts(selected_result_dir())) })

  # ── QC adicional: plots plotly ────────────────────────────────────────────
  align_qc_mapping_plot_obj <- function() {
    df <- align_qc_summary(selected_result_dir())
    if (!has_real_rows(df)) return(plotly_message(df$Mensaje[1]))

    # Si tenemos unique + multimapped + unmapped, mostramos la descomposicion completa.
    # Si no (multimapping_rate desconocido), mostramos mapped vs unmapped sin partir.
    full_cols <- c("unique_reads", "multimapped_reads", "unmapped_reads")
    has_full <- all(full_cols %in% names(df)) && !all(is.na(as.matrix(df[, full_cols])))

    if (has_full) {
      df[, full_cols] <- lapply(df[, full_cols, drop = FALSE], function(x) ifelse(is.na(x), 0, x))
      return(
        plotly::plot_ly(df, x = ~sample_id) |>
          plotly::add_bars(y = ~unique_reads, name = "Unicas", marker = list(color = "#7BBF9A")) |>
          plotly::add_bars(y = ~multimapped_reads, name = "Multimapeadas", marker = list(color = "#F6D58A")) |>
          plotly::add_bars(y = ~unmapped_reads, name = "No alineadas", marker = list(color = "#F4A6A6")) |>
          plotly::layout(
            barmode = "stack",
            xaxis = list(title = "Muestra"),
            yaxis = list(title = "Lecturas"),
            legend = list(orientation = "h", x = 0, y = 1.12),
            margin = list(b = 90)
          )
      )
    }

    fallback_cols <- c("mapped_reads", "unmapped_reads")
    if (!all(fallback_cols %in% names(df)) || all(is.na(as.matrix(df[, fallback_cols]))))
      return(plotly_message("No se encontraron metricas suficientes para representar lecturas unicas, multimapeadas y no alineadas."))
    df[, fallback_cols] <- lapply(df[, fallback_cols, drop = FALSE], function(x) ifelse(is.na(x), 0, x))
    plotly::plot_ly(df, x = ~sample_id) |>
      plotly::add_bars(y = ~mapped_reads, name = "Mapeadas (sin desglose multi)", marker = list(color = "#8BC9A6")) |>
      plotly::add_bars(y = ~unmapped_reads, name = "No alineadas", marker = list(color = "#F4A6A6")) |>
      plotly::layout(
        barmode = "stack",
        xaxis = list(title = "Muestra"),
        yaxis = list(title = "Lecturas"),
        legend = list(orientation = "h", x = 0, y = 1.12),
        margin = list(b = 90),
        annotations = list(list(
          xref = "paper", yref = "paper", x = 0, y = 1.18, showarrow = FALSE,
          text = "Sin tasa de multimapping disponible: se muestra mapeo total.",
          font = list(size = 11, color = "#5C4A16")
        ))
      )
  }

  align_qc_assignment_plot_obj <- function() {
    df <- align_qc_summary(selected_result_dir())
    if (!has_real_rows(df)) return(plotly_message(df$Mensaje[1]))
    cols <- c("assigned_reads", "unassigned_reads")
    if (!all(cols %in% names(df)) || all(is.na(as.matrix(df[, cols]))))
      return(plotly_message("No se encontraron metricas de asignacion genica para este analisis."))
    df[, cols] <- lapply(df[, cols, drop = FALSE], function(x) ifelse(is.na(x), 0, x))
    plotly::plot_ly(df, x = ~sample_id) |>
      plotly::add_bars(y = ~assigned_reads, name = "Asignadas", marker = list(color = "#8BC9A6")) |>
      plotly::add_bars(y = ~unassigned_reads, name = "No asignadas", marker = list(color = "#D7EEF1")) |>
      plotly::layout(
        barmode = "stack",
        xaxis = list(title = "Muestra"),
        yaxis = list(title = "Lecturas"),
        legend = list(orientation = "h", x = 0, y = 1.12),
        margin = list(b = 90)
      )
  }

  align_qc_region_plot_obj <- function() {
    plotly_message("No se encontraron metricas exonicas, intronicas o intergenicas para este analisis.")
  }

  align_qc_gene_body_plot_obj <- function() {
    plotly_message("No se encontro cobertura 5'-3' del cuerpo genico para este analisis.")
  }

  pseudo_qc_rate_plot_obj <- function() {
    df <- pseudo_qc_summary(selected_result_dir())
    if (!has_real_rows(df)) return(plotly_message(df$Mensaje[1]))
    if (!"pseudoalignment_rate" %in% names(df) || all(is.na(df$pseudoalignment_rate)))
      return(plotly_message("No se encontro pseudoalignment_rate para este analisis."))
    df$pseudoalignment_rate_pct <- 100 * df$pseudoalignment_rate
    plotly::plot_ly(
      df, x = ~sample_id, y = ~pseudoalignment_rate_pct,
      type = "bar", marker = list(color = "#8BC9A6"),
      text = ~paste0(round(pseudoalignment_rate_pct, 2), "%"),
      hovertemplate = "Muestra: %{x}<br>Rate: %{y:.2f}%<extra></extra>"
    ) |>
      plotly::layout(
        xaxis = list(title = "Muestra"),
        yaxis = list(title = "Pseudoalignment rate (%)", range = c(0, 100)),
        shapes = list(list(
          type = "line", xref = "paper", x0 = 0, x1 = 1,
          y0 = 100 * qc_thresholds$pseudoalignment_rate_warning,
          y1 = 100 * qc_thresholds$pseudoalignment_rate_warning,
          line = list(color = "#F4A6A6", dash = "dash")
        )),
        margin = list(b = 90)
      )
  }

  pseudo_qc_tpm_plot_obj <- function() {
    q <- pseudo_qc_quant_table(selected_result_dir())
    if (!has_real_rows(q)) return(plotly_message(q$Mensaje[1]))
    if (!all(c("TPM", "sample_id") %in% names(q)))
      return(plotly_message("No se encontraron valores TPM para este analisis."))
    q$TPM <- num_or_na(q$TPM)
    q$log_tpm <- log10(q$TPM + 1)
    plotly::plot_ly(
      q, x = ~sample_id, y = ~log_tpm,
      type = "box", color = ~sample_id,
      boxpoints = "outliers",
      hovertemplate = "Muestra: %{x}<br>log10(TPM+1): %{y:.3f}<extra></extra>"
    ) |>
      plotly::layout(
        showlegend = FALSE,
        xaxis = list(title = "Muestra"),
        yaxis = list(title = "log10(TPM + 1)"),
        margin = list(b = 90)
      )
  }

  pseudo_qc_detected_plot_obj <- function() {
    df <- pseudo_qc_summary(selected_result_dir())
    if (!has_real_rows(df)) return(plotly_message(df$Mensaje[1]))
    if (!"transcripts_detected" %in% names(df) || all(is.na(df$transcripts_detected)))
      return(plotly_message("No se pudo calcular el numero de transcritos detectados."))
    plotly::plot_ly(
      df, x = ~sample_id, y = ~transcripts_detected,
      type = "bar", marker = list(color = "#7BBF9A"),
      hovertemplate = "Muestra: %{x}<br>Transcritos detectados: %{y}<extra></extra>"
    ) |>
      plotly::layout(
        xaxis = list(title = "Muestra"),
        yaxis = list(title = "Transcritos detectados"),
        margin = list(b = 90)
      )
  }

  pseudo_qc_scatter_plot_obj <- function() {
    q <- pseudo_qc_quant_table(selected_result_dir())
    if (!has_real_rows(q)) return(plotly_message(q$Mensaje[1]))
    if (!all(c("TPM", "NumReads", "sample_id") %in% names(q)))
      return(plotly_message("No se encontraron TPM y NumReads para generar el scatter plot."))
    q$TPM <- num_or_na(q$TPM)
    q$NumReads <- num_or_na(q$NumReads)
    keep <- is.finite(q$TPM) & is.finite(q$NumReads)
    if (!any(keep)) return(plotly_message("No hay valores numericos validos de TPM y NumReads."))
    q <- q[keep, , drop = FALSE]
    q$log_tpm <- log10(q$TPM + 1)
    q$log_reads <- log10(q$NumReads + 1)
    hover_text <- paste0(
      "Muestra: ", q$sample_id,
      "<br>Name: ", q$Name %||% "",
      if ("Type" %in% names(q)) paste0("<br>Type: ", q$Type) else "",
      "<br>TPM: ", round(q$TPM, 3),
      "<br>NumReads: ", round(q$NumReads, 3)
    )
    plotly::plot_ly(
      q, x = ~log_tpm, y = ~log_reads,
      type = "scatter", mode = "markers",
      color = ~sample_id, text = hover_text, hoverinfo = "text",
      marker = list(size = 6, opacity = 0.55)
    ) |>
      plotly::layout(
        xaxis = list(title = "log10(TPM + 1)"),
        yaxis = list(title = "log10(NumReads + 1)"),
        legend = list(orientation = "h", x = 0, y = 1.12)
      )
  }

  output$align_qc_mapping_plot    <- plotly::renderPlotly(align_qc_mapping_plot_obj())
  output$align_qc_assignment_plot <- plotly::renderPlotly(align_qc_assignment_plot_obj())
  output$align_qc_region_plot     <- plotly::renderPlotly(align_qc_region_plot_obj())
  output$align_qc_gene_body_plot  <- plotly::renderPlotly(align_qc_gene_body_plot_obj())
  output$pseudo_qc_rate_plot      <- plotly::renderPlotly(pseudo_qc_rate_plot_obj())
  output$pseudo_qc_tpm_plot       <- plotly::renderPlotly(pseudo_qc_tpm_plot_obj())
  output$pseudo_qc_detected_plot  <- plotly::renderPlotly(pseudo_qc_detected_plot_obj())
  output$pseudo_qc_scatter_plot   <- plotly::renderPlotly(pseudo_qc_scatter_plot_obj())

  # ── Descargas QC adicional ────────────────────────────────────────────────
  output$download_align_qc_summary <- downloadHandler(
    filename = function() paste0("align_qc_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(f) write.csv(align_qc_summary(selected_result_dir()), f, row.names = FALSE)
  )
  output$download_align_qc_alerts <- csv_download(
    "align_qc_alerts",
    function() align_qc_alerts(selected_result_dir())
  )
  output$download_pseudo_qc_summary <- downloadHandler(
    filename = function() paste0("pseudo_qc_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(f) write.csv(pseudo_qc_summary(selected_result_dir()), f, row.names = FALSE)
  )
  output$download_pseudo_quant <- downloadHandler(
    filename = function() paste0("pseudo_quantification_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(f) write.csv(pseudo_qc_quant_table(selected_result_dir()), f, row.names = FALSE)
  )
  output$download_pseudo_qc_alerts <- csv_download(
    "pseudo_qc_alerts",
    function() pseudo_qc_alerts(selected_result_dir())
  )

  output$download_align_qc_mapping_plot    <- plotly_download("align_qc_mapping",    align_qc_mapping_plot_obj)
  output$download_align_qc_assignment_plot <- plotly_download("align_qc_assignment", align_qc_assignment_plot_obj)
  output$download_align_qc_region_plot     <- plotly_download("align_qc_region",     align_qc_region_plot_obj)
  output$download_align_qc_gene_body_plot  <- plotly_download("align_qc_gene_body",  align_qc_gene_body_plot_obj)
  output$download_pseudo_qc_rate_plot      <- plotly_download("pseudo_qc_rate",      pseudo_qc_rate_plot_obj)
  output$download_pseudo_qc_tpm_plot       <- plotly_download("pseudo_qc_tpm",       pseudo_qc_tpm_plot_obj)
  output$download_pseudo_qc_detected_plot  <- plotly_download("pseudo_qc_detected",  pseudo_qc_detected_plot_obj)
  output$download_pseudo_qc_scatter_plot   <- plotly_download("pseudo_qc_scatter",   pseudo_qc_scatter_plot_obj)

  output$download_qc_alerts <- downloadHandler(
    filename = function() paste0("qc_alerts_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(f) {
      a <- align_qc_alerts(selected_result_dir())
      p <- pseudo_qc_alerts(selected_result_dir())
      if (has_real_rows(a)) a$bloque <- "alineamiento"
      if (has_real_rows(p)) p$bloque <- "pseudoalineamiento"
      out <- if (has_real_rows(a) && has_real_rows(p)) rbind(a, p)
             else if (has_real_rows(a)) a
             else if (has_real_rows(p)) p
             else message_df("No se detectaron alertas automaticas con los umbrales actuales.")
      write.csv(out, f, row.names = FALSE)
    }
  )

  invisible(NULL)
}
