#' server_tab_deg_reports.R
#' Replicabilidad, comparacion de metodos y reproducibilidad (items 20, 21 y 23).
#'
#' Parte del modulo de la pestana 4, separado de server_tab_deg.R por tamano: el
#' fichero original llego a 1.608 lineas en una unica funcion. NO se usa
#' moduleServer(): igual que el resto de la aplicacion, se conservan los IDs
#' Shiny originales y el estado se pasa explicitamente.
#'
#' `ctx` es el contexto compartido del modulo (un environment, como `state`, para
#' que las asignaciones sean visibles entre las partes). Lo crea server_tab_deg()
#' y contiene los reactivos que varias partes necesitan.

server_tab_deg_reports <- function(input, output, session, state, ctx) {
  # ── Sugerencia de metodo segun el tamano muestral (item 21) ────────────────
  output$deg_method_hint <- renderUI({
    df <- ctx$meta_rv()
    if (is.null(df) || !nrow(df)) return(NULL)
    cc <- input$deg_condition_col %||% "condition"
    if (!cc %in% names(df)) return(NULL)
    d <- df; d$condition <- d[[cc]]
    s <- suggest_robust_comparison(d)
    m <- input$deg_method %||% "DESeq2"
    if (m %in% DEG_METHODS_ROBUST && is.null(s)) {
      return(div(class = "alert alert-warning py-2 px-2 small mb-2",
                 icon("triangle-exclamation"),
                 paste(" Con pocas replicas los metodos robustos pierden potencia:",
                       "por debajo de 8 muestras por grupo, los parametricos son la",
                       "eleccion correcta.")))
    }
    if (is.null(s)) return(NULL)
    div(class = "alert alert-info py-2 px-2 small mb-2", icon("lightbulb"), " ", s$message)
  })

  # ── Replicabilidad por bootstrap (item 20) ─────────────────────────────────
  boot_rv <- reactiveVal(NULL)

  observeEvent(input$deg_run_boot_btn, {
    req(state$deg_rv$results, state$deg_rv$counts, state$deg_rv$meta)
    n_boot <- input$deg_boot_n %||% 20
    withProgress(message = "Estimando replicabilidad...", value = 0, {
      res <- bootstrap_replicability(
        state$deg_rv$counts, state$deg_rv$meta,
        method = state$deg_rv$method %||% "DESeq2",
        ref_level = input$deg_contrast_den, contrast_num = input$deg_contrast_num,
        batch = if (isTRUE(input$deg_use_batch)) input$deg_batch_col else NULL,
        fdr = state$deg_rv$fdr %||% 0.05,
        lfc_threshold = state$deg_rv$lfc_threshold %||% 0,
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
                 paste0(r$n_ok, " remuestreos validos",
                        if (r$n_failed > 0) paste0(", ", r$n_failed, " descartados") else "",
                        "  ·  top-N = ", r$top_n)))
  })

  output$deg_boot_table <- renderDT({
    r <- boot_rv()
    if (is.null(r)) return(dt_table(message_df("Sin estimacion de replicabilidad.")))
    df <- r$summary
    for (nm in c("q1", "mediana", "q3")) df[[nm]] <- round(df[[nm]], 4)
    dt_table(df, page_length = 5)
  })

  output$deg_boot_plot <- plotly::renderPlotly({
    r <- boot_rv()
    if (is.null(r)) return(plotly_message("Sin estimacion de replicabilidad."))
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

  # ── Comparacion entre metodos (item 21) ────────────────────────────────────
  compare_rv <- reactiveVal(NULL)

  observe({
    cur <- input$deg_method %||% "DESeq2"
    all_m <- c(DEG_METHODS_PARAMETRIC, "Wilcoxon",
               if (isTRUE(HAS_DEARSEQ)) "dearseq")
    updateSelectInput(session, "deg_compare_method",
                      choices = setdiff(all_m, cur),
                      selected = isolate(input$deg_compare_method))
  })

  observeEvent(input$deg_run_compare_btn, {
    req(state$deg_rv$results, state$deg_rv$counts, state$deg_rv$meta)
    other <- input$deg_compare_method %||% ""
    if (!nzchar(other)) {
      showNotification("Elige un metodo con el que comparar.", type = "warning"); return()
    }
    withProgress(message = paste0("Corriendo ", other, "..."), value = 0.4, {
      r2 <- tryCatch(run_deg(
        state$deg_rv$counts, state$deg_rv$meta, method = other,
        ref_level = input$deg_contrast_den, contrast_num = input$deg_contrast_num,
        batch = if (isTRUE(input$deg_use_batch)) input$deg_batch_col else NULL,
        fdr = state$deg_rv$fdr %||% 0.05,
        lfc_threshold = state$deg_rv$lfc_threshold %||% 0, shrink = FALSE),
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
                              "Elige un metodo y pulsa 'Comparar'."))
    div(class = "small",
        tags$b(paste0(o$name_a, " vs ", o$name_b, ": ")),
        paste0(fmt_int(o$n_common), " genes en comun; ", fmt_int(o$only_a),
               " solo en ", o$name_a, "; ", fmt_int(o$only_b), " solo en ",
               o$name_b, ". Jaccard = ",
               if (is.na(o$jaccard)) "—" else round(o$jaccard, 3), "."))
  })

  output$deg_compare_plot <- plotly::renderPlotly({
    o <- compare_rv()
    if (is.null(o)) return(plotly_message("Sin comparacion."))
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
    if (is.null(state$deg_rv$results)) return("Lanza primero un analisis DEG.")
    build_deg_r_script(state$deg_rv) %||% "—"
  })

  output$download_deg_report <- downloadHandler(
    filename = function() paste0("informe_deg_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html"),
    content = function(f) {
      h <- build_deg_report_html(state$deg_rv, report_diagnostics())
      writeLines(h %||% "<html><body><p>Sin analisis DEG.</p></body></html>", f)
    }
  )
  output$download_deg_script <- downloadHandler(
    filename = function() paste0("analisis_deg_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".R"),
    content = function(f) {
      writeLines(build_deg_r_script(state$deg_rv) %||% "# Sin analisis DEG.", f)
    }
  )

  invisible(NULL)
}
