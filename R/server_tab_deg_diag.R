#' server_tab_deg_diag.R
#' Diagnosticos post-ajuste de la pestana 4 (item 10 y B3c).
#'
#' Parte del modulo de la pestana 4, separado de server_tab_deg.R por tamano: el
#' fichero original llego a 1.608 lineas en una unica funcion. NO se usa
#' moduleServer(): igual que el resto de la aplicacion, se conservan los IDs
#' Shiny originales y el estado se pasa explicitamente.
#'
#' `ctx` es el contexto compartido del modulo (un environment, como `state`, para
#' que las asignaciones sean visibles entre las partes). Lo crea server_tab_deg()
#' y contiene los reactivos que varias partes necesitan.

server_tab_deg_diag <- function(input, output, session, state, ctx) {
  # ── Diagnosticos post-ajuste ───────────────────────────────────────────────

  deg_na_breakdown <- reactive({
    padj_na_breakdown(state$deg_rv$results)
  })

  output$deg_na_breakdown <- renderUI({
    b <- deg_na_breakdown()
    txt <- padj_na_breakdown_text(b)
    if (is.null(txt)) return(NULL)
    div(class = "small text-muted mt-2 pt-2 border-top",
        tags$b("Genes sin p-valor ajustado: "), txt,
        if (isTRUE(b$n_outlier > 0)) tags$span(
          " Los marcados como outlier tienen un valor extremo en alguna muestra:",
          " mirar la pestana de diagnosticos antes de descartarlos.") else NULL)
  })

  # p-valores usados en los diagnosticos. Por defecto solo los genes que pasan el
  # filtrado independiente, que es el conjunto sobre el que se controla la FDR;
  # incluir los de conteo muy bajo distorsiona la forma del histograma.
  diag_pvalues <- reactive({
    tab <- state$deg_rv$results
    if (is.null(tab)) return(NULL)
    if (identical(input$deg_diag_pv_subset %||% "tested", "tested")) {
      tab$pvalue[!is.na(tab$padj)]
    } else {
      tab$pvalue[!is.na(tab$pvalue)]
    }
  })

  diag_pi0 <- reactive({
    p <- diag_pvalues()
    if (is.null(p)) return(NULL)
    estimate_pi0(p)
  })

  output$deg_diag_verdict <- renderUI({
    p <- diag_pvalues(); pi <- diag_pi0()
    if (is.null(p) || !length(p)) return(NULL)
    d <- diagnose_pvalue_shape(p, pi$pi0 %||% NA_real_)
    cls <- switch(d$verdict,
      "esperado" = "alert-success",
      "sospechoso" = "alert-danger",
      "bimodal" = "alert-warning",
      "conservador" = "alert-warning",
      "alert-secondary")
    pi_txt <- if (is.null(pi) || is.na(pi$pi0)) "no estimable" else paste0(
      round(pi$pi0, 3),
      if (!is.na(pi$pi0_qvalue)) paste0(" (qvalue: ", round(pi$pi0_qvalue, 3), ")") else "")
    div(class = paste("alert py-2 px-3 mb-2", cls),
        tags$b(d$label), tags$div(class = "small mt-1", d$detail),
        tags$div(class = "small mt-1",
                 tags$b("pi0 (proporcion de nulas ciertas): "), pi_txt,
                 tags$span(class = "text-muted",
                           paste0("  ·  ", fmt_int(length(p)), " p-valores"))))
  })

  make_pvalue_hist <- function(p, pi0 = NA_real_) {
    h <- pvalue_hist_data(p, bins = 40)
    if (is.null(h)) return(plotly_message("Sin p-valores para el histograma."))
    # Suelo esperado: si una fraccion pi0 de los genes es nula, sus p-valores se
    # reparten uniformemente y cada bin recibe n*pi0/bins.
    floor_h <- if (!is.na(pi0)) length(p) * pi0 / nrow(h) else NA_real_
    pl <- plotly::plot_ly(h, x = ~mid, y = ~count, type = "bar",
                          marker = list(color = "#7BBF9A"),
                          text = ~paste0("p en [", round(mid - 0.0125, 3), ", ",
                                         round(mid + 0.0125, 3), "]<br>genes: ", count),
                          hoverinfo = "text", name = "observado") |>
      plotly::layout(
        xaxis = list(title = "p-valor", range = c(0, 1)),
        yaxis = list(title = "numero de genes"),
        showlegend = !is.na(floor_h),
        bargap = 0.02
      )
    if (!is.na(floor_h)) {
      pl <- pl |> plotly::add_lines(
        x = c(0, 1), y = c(floor_h, floor_h), inherit = FALSE,
        line = list(dash = "dash", color = "#F4A6A6"),
        name = "suelo esperado bajo la nula"
      )
    }
    pl
  }

  output$deg_pvalue_hist <- plotly::renderPlotly({
    p <- diag_pvalues()
    if (is.null(p) || !length(p)) return(plotly_message("Lanza primero un analisis DEG."))
    make_pvalue_hist(p, diag_pi0()$pi0 %||% NA_real_)
  })

  output$deg_pvalue_reference <- renderPlot({
    shapes <- reference_pvalue_shapes()
    op <- graphics::par(mfrow = c(1, 3), mar = c(3.2, 3.2, 2.4, 0.8), mgp = c(2, 0.6, 0))
    on.exit(graphics::par(op), add = TRUE)
    for (nm in names(shapes)) {
      graphics::hist(shapes[[nm]], breaks = seq(0, 1, length.out = 41),
                     col = "#A8DADC", border = "white", main = nm,
                     xlab = "p-valor", ylab = "genes", cex.main = 0.95)
    }
  })

  make_disp_plot <- function(dd) {
    if (is.null(dd) || !nrow(dd)) {
      return(plotly_message(paste("Las dispersiones solo las produce DESeq2.",
                                  "Relanza el analisis con ese motor.")))
    }
    d <- dd[!is.na(dd$baseMean) & dd$baseMean > 0, , drop = FALSE]
    plotly::plot_ly() |>
      plotly::add_markers(data = d, x = ~baseMean, y = ~dispGeneEst,
                          name = "estimacion por gen",
                          marker = list(size = 4, opacity = 0.35, color = "#C0C0C0"),
                          hoverinfo = "skip") |>
      plotly::add_markers(data = d, x = ~baseMean, y = ~dispersion,
                          name = "dispersion final",
                          marker = list(size = 4, opacity = 0.5, color = "#7BBF9A"),
                          hoverinfo = "skip") |>
      plotly::add_markers(data = d[order(d$baseMean), ], x = ~baseMean, y = ~dispFit,
                          name = "curva ajustada",
                          marker = list(size = 3, color = "#F4A6A6"),
                          hoverinfo = "skip") |>
      plotly::layout(
        xaxis = list(title = "media de conteos normalizados", type = "log"),
        yaxis = list(title = "dispersion", type = "log")
      )
  }

  output$deg_disp_plot <- plotly::renderPlotly({
    make_disp_plot(state$deg_rv$disp_data)
  })

  make_rle_plot <- function(counts, meta = NULL) {
    rl <- rle_summary(counts)
    if (is.null(rl)) return(plotly_message("Sin datos para el RLE."))
    rl$color <- ifelse(rl$flag, "Desviada", "Normal")
    plotly::plot_ly(rl, x = ~sample_id) |>
      plotly::add_segments(y = ~p05, yend = ~p95, x = ~sample_id, xend = ~sample_id,
                           line = list(color = "#C0C0C0", width = 1),
                           showlegend = FALSE, hoverinfo = "skip") |>
      plotly::add_segments(y = ~q1, yend = ~q3, x = ~sample_id, xend = ~sample_id,
                           line = list(color = "#A8DADC", width = 9),
                           showlegend = FALSE, hoverinfo = "skip") |>
      plotly::add_markers(y = ~med, color = ~color,
                          colors = c("Normal" = "#244B34", "Desviada" = "#D9534F"),
                          marker = list(size = 9),
                          text = ~paste0("Muestra: ", sample_id,
                                         "<br>mediana RLE: ", round(med, 3),
                                         "<br>IQR: ", round(iqr, 3)),
                          hoverinfo = "text") |>
      plotly::layout(
        xaxis = list(title = ""),
        yaxis = list(title = "RLE (log2 respecto a la mediana del gen)"),
        shapes = list(list(type = "line", xref = "paper", x0 = 0, x1 = 1,
                           y0 = 0, y1 = 0,
                           line = list(dash = "dash", color = "#7BBF9A", width = 1)))
      )
  }

  output$deg_rle_plot <- plotly::renderPlotly({
    req(state$deg_rv$counts)
    make_rle_plot(state$deg_rv$counts, state$deg_rv$meta)
  })

  output$deg_cooks_warning <- renderUI({
    ck <- state$deg_rv$cooks
    if (is.null(ck) || is.na(ck$dominant)) return(NULL)
    div(class = "alert alert-warning py-2 px-3 mb-2",
        icon("triangle-exclamation"),
        tags$b(paste0(" La muestra '", ck$dominant, "' concentra los outliers. ")),
        paste0("Con ", nrow(ck$table), " muestras, lo esperado por azar seria ",
               round(100 * ck$expected_frac, 1), " % cada una."))
  })

  make_cooks_plot <- function(ck) {
    if (is.null(ck) || is.null(ck$table)) {
      return(plotly_message(paste("Las distancias de Cook solo las produce DESeq2.",
                                  "Relanza el analisis con ese motor.")))
    }
    df <- ck$table
    df$pct <- 100 * df$frac
    plotly::plot_ly(df, x = ~sample_id, y = ~pct, type = "bar",
                    marker = list(color = ifelse(
                      !is.na(ck$dominant) & df$sample_id == ck$dominant,
                      "#D9534F", "#7BBF9A")),
                    text = ~paste0("Muestra: ", sample_id,
                                   "<br>genes donde es el maximo: ", n_max,
                                   " (", round(pct, 1), " %)",
                                   "<br>Cook maximo: ", signif(max_cooks, 3)),
                    hoverinfo = "text") |>
      plotly::layout(
        xaxis = list(title = ""),
        yaxis = list(title = "% de genes donde la muestra es el maximo de Cook"),
        shapes = list(list(type = "line", xref = "paper", x0 = 0, x1 = 1,
                           y0 = 100 * ck$expected_frac, y1 = 100 * ck$expected_frac,
                           line = list(dash = "dot", color = "#60756A")))
      )
  }

  output$deg_cooks_plot <- plotly::renderPlotly({
    make_cooks_plot(state$deg_rv$cooks)
  })

  output$download_deg_pvalue_hist <- plotly_download(
    "deg_pvalue_hist",
    function() {
      p <- diag_pvalues()
      if (is.null(p)) return(plotly_message("Sin p-valores."))
      make_pvalue_hist(p, diag_pi0()$pi0 %||% NA_real_)
    }
  )
  output$download_deg_disp_plot <- plotly_download(
    "deg_dispersiones", function() make_disp_plot(state$deg_rv$disp_data))
  output$download_deg_rle_plot <- plotly_download(
    "deg_rle", function() {
      if (is.null(state$deg_rv$counts)) return(plotly_message("Sin conteos."))
      make_rle_plot(state$deg_rv$counts, state$deg_rv$meta)
    })
  output$download_deg_cooks_plot <- plotly_download(
    "deg_cooks", function() make_cooks_plot(state$deg_rv$cooks))

  # ── Diagnostico de sesgo de longitud (B3c) ─────────────────────────────────
  # Las longitudes salen de la anotacion de la ejecucion. Con una matriz subida no
  # hay anotacion asociada, asi que el diagnostico no es calculable y se dice.
  deg_gene_lengths <- reactive({
    src <- input$deg_source %||% "current"
    af <- switch(
      src,
      "saved"   = annotation_file_for_run(input$selected_deg_run_dir %||% ""),
      "current" = {
        p <- state$run_params_rv()
        p$annotation_file %||% annotation_file_for_run(p$output_dir %||% "")
      },
      # Una matriz subida no tiene anotacion asociada. Antes caia en la rama de
      # la ejecucion actual y usaba SU anotacion, que no tiene por que
      # corresponder a los genes de la matriz: si los identificadores
      # solapaban parcialmente, el diagnostico de sesgo de longitud salia
      # calculado sobre longitudes de otro organismo.
      NULL
    )
    if (is.null(af) || !nzchar(af %||% "")) return(NULL)
    gene_lengths_from_annotation(af)
  })

  deg_length_bias <- reactive({
    req(state$deg_rv$results)
    len <- deg_gene_lengths()
    if (is.null(len)) return(NULL)
    length_bias_diagnostic(state$deg_rv$results, len,
                           fdr = state$deg_rv$fdr %||% 0.05)
  })

  output$deg_lenbias_verdict <- renderUI({
    if (is.null(state$deg_rv$results)) {
      return(div(class = "small text-muted", "Lanza primero un analisis DEG."))
    }
    if (is.null(deg_gene_lengths())) {
      return(div(class = "alert alert-secondary py-2 px-2 small",
                 icon("circle-info"),
                 paste(" No se puede evaluar: hacen falta las longitudes de gen, que",
                       "vienen del fichero de anotacion de la ejecucion. Con una matriz",
                       "de conteos subida no hay anotacion asociada.")))
    }
    lb <- deg_length_bias()
    if (is.null(lb) || is.null(lb$verdict)) {
      return(div(class = "alert alert-secondary py-2 px-2 small",
                 "No hay suficientes genes con longitud conocida para evaluarlo."))
    }
    cls <- switch(lb$verdict$level, "aviso" = "alert-warning", "leve" = "alert-info",
                  "ok" = "alert-success", "alert-secondary")
    div(class = paste("alert py-2 px-3 mb-2", cls),
        tags$b(lb$verdict$label),
        tags$div(class = "small mt-1", lb$verdict$detail),
        tags$div(class = "small mt-1 text-muted",
                 paste0(fmt_int(lb$n_used), " genes con longitud conocida")))
  })

  make_deg_lenbias_plot <- function() {
    lb <- tryCatch(deg_length_bias(), error = function(e) NULL)
    if (is.null(lb) || is.null(lb$table)) {
      return(plotly_message("Sin diagnostico de sesgo de longitud."))
    }
    df <- lb$table
    df$pct <- 100 * df$prop_de
    overall <- 100 * sum(df$n_de) / sum(df$n_genes)
    plotly::plot_ly(
      df, x = ~len_median, y = ~pct, type = "scatter", mode = "lines+markers",
      line = list(color = "#244B34", shape = "spline"),
      marker = list(color = "#7BBF9A", size = 9),
      text = ~paste0("longitud mediana: ", fmt_int(len_median), " pb",
                     "<br>genes: ", fmt_int(n_genes),
                     "<br>diferenciales: ", fmt_int(n_de), " (", round(pct, 1), " %)"),
      hoverinfo = "text"
    ) |>
      plotly::layout(
        xaxis = list(title = "longitud del gen (pb, mediana del bin)", type = "log"),
        yaxis = list(title = "% de genes diferenciales"),
        shapes = list(
          # Referencia: la proporcion global. Una curva plana sobre esta linea
          # significa que la longitud no influye.
          list(type = "line", xref = "paper", x0 = 0, x1 = 1,
               y0 = overall, y1 = overall,
               line = list(dash = "dash", color = "#F4A6A6"))),
        annotations = list(
          list(x = 1, y = overall, xref = "paper", yanchor = "bottom",
               text = paste0("proporcion global: ", round(overall, 1), " %"),
               showarrow = FALSE, xanchor = "right", font = list(size = 10)))
      )
  }

  output$deg_lenbias_plot <- plotly::renderPlotly(make_deg_lenbias_plot())
  output$download_deg_lenbias_plot <- plotly_download(
    "deg_sesgo_longitud", function() make_deg_lenbias_plot())

  # El informe reproducible necesita estos tres diagnosticos. Se publican en el
  # contexto en lugar de duplicar su calculo.
  ctx$diag_pvalues      <- diag_pvalues
  ctx$diag_pi0          <- diag_pi0
  ctx$deg_na_breakdown  <- deg_na_breakdown

  invisible(NULL)
}
