#' server_tab_deg_reports.R
#' Replicabilidad, comparación de métodos y reproducibilidad (items 20, 21 y 23).
#'
#' Parte del modulo de la pestana 4, separado de server_tab_deg.R por tamaño: el
#' fichero original llego a 1.608 líneas en una única función. NO se usa
#' moduleServer(): igual que el resto de la aplicación, se conservan los IDs
#' Shiny originales y el estado se pasa explicitamente.
#'
#' `ctx` es el contexto compartido del modulo (un environment, como `state`, para
#' que las asignaciones sean visibles entre las partes). Lo crea server_tab_deg()
#' y contiene los reactivos que varias partes necesitan.

server_tab_deg_reports <- function(input, output, session, state, ctx) {
  # ── Sugerencia de método según el tamaño muestral (item 21) ────────────────

  # ── Replicabilidad por bootstrap (item 20) ─────────────────────────────────
  boot_rv <- reactiveVal(NULL)

  observeEvent(input$deg_run_boot_btn, {
    req(state$deg_rv$results, state$deg_rv$counts, state$deg_rv$meta)
    n_boot <- input$deg_boot_n %||% 20
    # Swish trabaja sobre replicas inferenciales y a nivel de transcrito, no
    # sobre la matriz de conteos: run_deg() ni siquiera lo acepta como motor, así
    # que el remuestreo moriria en match.arg con un error incomprensible.
    if (identical(state$deg_rv$method, "Swish")) {
      showNotification(paste0(
        "El bootstrap de replicabilidad no está disponible para Swish: remuestrea ",
        "la matriz de conteos, y Swish parte de las replicas inferenciales de la ",
        "cuantificacion."), type = "warning", duration = 14)
      boot_rv(NULL); return()
    }
    withProgress(message = "Estimando replicabilidad...", value = 0, {
      # TODOS los parámetros salen del estado del ajuste, no de la interfaz. Los
      # selectores siguen vivos después de ajustar, así que leerlos aquí podia
      # medir la replicabilidad de un contraste o un diseño distintos de los que
      # produjeron los resultados que se están evaluando.
      #
      # Límite conocido: si el diseño incluye variables sustitutas, el remuestreo
      # reutiliza las SV estimadas sobre los datos COMPLETOS en lugar de
      # reestimarlas en cada remuestreo. Reestimarlas sería lo correcto en
      # sentido estricto, pero sva falla a menudo sobre submuestras y el coste se
      # multiplica por n_boot. Reutilizarlas tiende a dar una replicabilidad algo
      # OPTIMISTA, que es la dirección menos peligrosa para un aviso.
      res <- bootstrap_replicability(
        state$deg_rv$counts, state$deg_rv$meta,
        method = state$deg_rv$method %||% "DESeq2",
        ref_level = state$deg_rv$ref_level, contrast_num = state$deg_rv$contrast_num,
        batch = state$deg_rv$batch,
        fdr = state$deg_rv$fdr %||% 0.05,
        lfc_threshold = state$deg_rv$lfc_threshold %||% 0,
        design_formula = state$deg_rv$design_formula,
        test_coef = state$deg_rv$test_coef,
        seed = ANALYSIS_SEED,
        n_boot = n_boot,
        progress = function(i, n) setProgress(value = i / n,
                                              detail = paste0(i, " / ", n)))
    })
    if (!is.null(res$error)) {
      showNotification(res$error, type = "error", duration = 16)
      boot_rv(NULL); return()
    }
    boot_rv(res)
    if (identical(res$interpretation$level, "baja")) {
      showNotification(res$interpretation$detail, type = "warning", duration = 20)
    }
  })

  output$deg_boot_verdict <- renderUI({
    r <- boot_rv()
    if (is.null(r)) return(div(class = "small text-muted",
                               "Pulsa 'Estimar' para medir la replicabilidad."))
    cls <- switch(r$interpretation$level, "alta" = "alert-success",
                  "baja" = "alert-danger", "intermedia" = "alert-warning",
                  "alert-secondary")
    div(class = paste("alert py-2 px-3 mb-2", cls),
        tags$b(r$interpretation$label),
        tags$div(class = "small mt-1", r$interpretation$detail),
        tags$div(class = "small mt-1 text-muted",
                 paste0(r$n_ok, " remuestreos válidos",
                        if (r$n_failed > 0) paste0(", ", r$n_failed, " descartados") else "",
                        "  ·  top-N = ", r$top_n)))
  })

  output$deg_boot_table <- renderDT({
    r <- boot_rv()
    if (is.null(r)) return(dt_table(message_df("Sin estimación de replicabilidad.")))
    dt_table(r$summary, page_length = 5)
  })

  output$deg_boot_plot <- plotly::renderPlotly({
    r <- boot_rv()
    if (is.null(r)) return(plotly_message("Sin estimación de replicabilidad."))
    pb <- r$per_boot
    plotly::plot_ly() |>
      plotly::add_markers(x = pb$boot, y = pb$spearman, name = "Spearman",
                          marker = list(color = "#244B34", size = 8)) |>
      plotly::add_markers(x = pb$boot, y = pb$jaccard_topn, name = "Jaccard top-N",
                          marker = list(color = "#7BBF9A", size = 8)) |>
      plotly::layout(
        xaxis = list(title = "remuestreo"),
        yaxis = list(title = "concordancia", range = c(0, 1)),
        shapes = list(
          list(type = "line", xref = "paper", x0 = 0, x1 = 1,
               y0 = REPLICABILITY_HIGH, y1 = REPLICABILITY_HIGH,
               line = list(dash = "dot", color = "#7BBF9A")),
          list(type = "line", xref = "paper", x0 = 0, x1 = 1,
               y0 = REPLICABILITY_LOW, y1 = REPLICABILITY_LOW,
               line = list(dash = "dot", color = "#D9534F")))
      )
  })

  # ── Comparación entre métodos (item 21) ────────────────────────────────────
  compare_rv <- reactiveVal(NULL)

  observe({
    cur <- input$deg_method %||% "DESeq2"
    all_m <- DEG_METHODS_PARAMETRIC
    updateSelectInput(session, "deg_compare_method",
                      choices = setdiff(all_m, cur),
                      selected = isolate(input$deg_compare_method))
  })

  observeEvent(input$deg_run_compare_btn, {
    req(state$deg_rv$results, state$deg_rv$counts, state$deg_rv$meta)
    other <- input$deg_compare_method %||% ""
    if (!nzchar(other)) {
      showNotification("Elige un método con el que comparar.", type = "warning"); return()
    }
    # Comparar contra Swish cruzaria identificadores de TRANSCRITO con
    # identificadores de GEN, dando un solapamiento de 0 que parece un resultado
    # y no lo es.
    if (identical(state$deg_rv$method, "Swish")) {
      showNotification(paste0(
        "La comparación entre métodos no está disponible con Swish: sus ",
        "resultados son por transcrito y los del resto por gen, así que el ",
        "solapamiento no sería interpretable."), type = "warning", duration = 14)
      compare_rv(NULL); return()
    }
    withProgress(message = paste0("Corriendo ", other, "..."), value = 0.4, {
      # Mismo criterio que el bootstrap: el contraste y el diseño son los del
      # ajuste, no los que tengan los selectores en este instante. De lo
      # contrario se compararian dos análisis distintos entre si.
      r2 <- tryCatch(run_deg(
        state$deg_rv$counts, state$deg_rv$meta, method = other,
        ref_level = state$deg_rv$ref_level, contrast_num = state$deg_rv$contrast_num,
        batch = state$deg_rv$batch,
        fdr = state$deg_rv$fdr %||% 0.05,
        lfc_threshold = state$deg_rv$lfc_threshold %||% 0, shrink = FALSE,
        design_formula = state$deg_rv$design_formula,
        test_coef = state$deg_rv$test_coef),
        error = function(e) list(table = NULL, error = conditionMessage(e)))
    })
    if (is.null(r2$table)) {
      showNotification(paste0("El motor ", other, " fallo: ", r2$error %||% "—"),
                       type = "error", duration = 14)
      compare_rv(NULL); return()
    }
    compare_rv(deg_method_overlap(state$deg_rv$results, r2$table,
                                 fdr = state$deg_rv$fdr %||% 0.05,
                                 name_a = state$deg_rv$method, name_b = other))
  })

  output$deg_compare_summary <- renderUI({
    o <- compare_rv()
    if (is.null(o)) return(div(class = "small text-muted",
                              "Elige un método y pulsa 'Comparar'."))
    div(class = "small",
        tags$b(paste0(o$name_a, " vs ", o$name_b, ": ")),
        paste0(fmt_int(o$n_common), " genes en comun; ", fmt_int(o$only_a),
               " solo en ", o$name_a, "; ", fmt_int(o$only_b), " solo en ",
               o$name_b, ". Jaccard = ",
               if (is.na(o$jaccard)) "—" else round(o$jaccard, 3), "."))
  })

  output$deg_compare_plot <- plotly::renderPlotly({
    o <- compare_rv()
    if (is.null(o)) return(plotly_message("Sin comparación."))
    df <- data.frame(
      cat = c(paste0("Solo ", o$name_a), "En ambos", paste0("Solo ", o$name_b)),
      n = c(o$only_a, o$n_common, o$only_b), stringsAsFactors = FALSE)
    plotly::plot_ly(df, x = ~cat, y = ~n, type = "bar",
                    marker = list(color = c("#A8DADC", "#7BBF9A", "#F4A6A6")),
                    text = ~fmt_int(n), textposition = "outside",
                    hoverinfo = "text") |>
      plotly::layout(xaxis = list(title = ""),
                     yaxis = list(title = "genes significativos"))
  })

  # ── Reproducibilidad: informe e script (item 23) ───────────────────────────
  report_diagnostics <- reactive({
    list(pi0 = ctx$diag_pi0()$pi0 %||% NA_real_,
         verdict = {
           p <- ctx$diag_pvalues()
           if (is.null(p) || !length(p)) NULL
           else diagnose_pvalue_shape(p, ctx$diag_pi0()$pi0 %||% NA_real_)
         },
         na_breakdown = ctx$deg_na_breakdown(),
         cooks_dominant = state$deg_rv$cooks$dominant %||% NA_character_,
         replicability = boot_rv())
  })

  output$deg_script_preview <- renderText({
    if (is.null(state$deg_rv$results)) return("Lanza primero un análisis DEG.")
    build_deg_r_script(state$deg_rv) %||% "—"
  })

  output$download_deg_report <- downloadHandler(
    filename = function() paste0("informe_deg_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html"),
    content = function(f) {
      h <- build_deg_report_html(state$deg_rv, report_diagnostics())
      writeLines(h %||% "<html><body><p>Sin análisis DEG.</p></body></html>", f)
    }
  )
  output$download_deg_script <- downloadHandler(
    filename = function() paste0("analisis_deg_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".R"),
    content = function(f) {
      writeLines(build_deg_r_script(state$deg_rv) %||% "# Sin análisis DEG.", f)
    }
  )

  # ── Tabla de correspondencia de la seudonimizacion ────────────────────────
  # Se descarga APARTE del informe y de los resultados a propósito: es lo único
  # que permite volver a los identificadores reales, así que exportarla tiene que
  # ser una decisión deliberada y no un efecto colateral de descargar el informe.
  output$deg_pseudonym_ui <- renderUI({
    m <- state$deg_rv$pseudonym_map
    if (is.null(m) || !nrow(m)) return(NULL)
    tagList(
      div(class = "alert alert-secondary py-2 px-2 small mb-2",
          icon("user-shield"),
          tags$b(" Identificadores seudonimizados: "), nrow(m), " muestras. ",
          "Los entregables llevan los alias. Esto es seudonimizacion, no ",
          "anonimizacion: existe una tabla de correspondencia, así que los datos ",
          "siguen siendo datos personales a efectos del RGPD. Ten en cuenta además ",
          "que los propios niveles de expresión permiten inferir genotipos ",
          "(Schadt et al., Nature Genetics 2012), de modo que renombrar las ",
          "columnas reduce la exposicion accidental pero no anonimiza la matriz."),
      downloadButton("download_pseudonym_map",
                     " Descargar tabla de correspondencia", class = "btn-sm")
    )
  })

  output$download_pseudonym_map <- downloadHandler(
    filename = function() paste0("correspondencia_muestras_",
                                 format(Sys.time(), "%Y%m%d_%H%M%S"), ".tsv"),
    content = function(f) {
      m <- state$deg_rv$pseudonym_map
      if (is.null(m)) m <- data.frame(original = character(0), alias = character(0))
      utils::write.table(m, f, sep = "\t", quote = FALSE, row.names = FALSE)
      append_audit_log("export_correspondencia",
                       list(n_muestras = nrow(m)), outputs_dir = state$outputs_dir)
    }
  )

  invisible(NULL)
}
