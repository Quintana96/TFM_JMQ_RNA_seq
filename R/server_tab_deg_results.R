#' server_tab_deg_results.R
#' Tabla de resultados, volcano, MA, PCA y heatmaps de la pestana 4.
#'
#' Parte del modulo de la pestana 4, separado de server_tab_deg.R por tamaño: el
#' fichero original llego a 1.608 líneas en una única función. NO se usa
#' moduleServer(): igual que el resto de la aplicación, se conservan los IDs
#' Shiny originales y el estado se pasa explicitamente.
#'
#' `ctx` es el contexto compartido del modulo (un environment, como `state`, para
#' que las asignaciones sean visibles entre las partes). Lo crea server_tab_deg()
#' y contiene los reactivos que varias partes necesitan.

server_tab_deg_results <- function(input, output, session, state, ctx) {
  # ── Tabla DT ───────────────────────────────────────────────────────────────
  output$deg_table <- renderDT({
    df <- ctx$deg_filtered()
    if (is.null(df)) return(dt_table(message_df("Sin resultados DEG.")))
    df_r <- df
    # Las columnas que un motor no rellena (log2FC_shrunk solo lo produce
    # DESeq2, lfcSE no lo produce edgeR) se ocultan en vez de mostrar una
    # columna entera de NA.
    for (nm in intersect(c("log2FC_shrunk", "lfcSE", "stat",
                           "log2FC_lower", "log2FC_upper"), names(df_r))) {
      if (all(is.na(df_r[[nm]]))) df_r[[nm]] <- NULL
    }
    if ("log2FC" %in% names(df_r)) {
      df_r$direction <- ifelse(df_r$log2FC >= 0, "Up", "Down")
      df_r <- df_r[, c("gene", "direction", setdiff(names(df_r), c("gene", "direction"))), drop = FALSE]
    }
    dt_table_num(df_r, page_length = 15, filter = "top")
  })

  # ── Helpers de ploteo (reutilizados en render y en descarga) ──────────────

  #' Elige el eje de fold-change: el encogido si existe. Los estimadores de
  #' máxima verosimilitud están sesgados hacía valores exagerados en genes de
  #' baja expresión, así que un volcano construido sobre ellos destaca
  #' visualmente los genes peor estimados. Se ordena y se dibuja con el
  #' encogido; se testea con el MLE.
  lfc_axis <- function(df) {
    if ("log2FC_shrunk" %in% names(df) && !all(is.na(df$log2FC_shrunk))) {
      list(values = df$log2FC_shrunk, label = "log2 Fold Change (encogido)")
    } else {
      list(values = df$log2FC, label = "log2 Fold Change (MLE)")
    }
  }

  plot_title <- function(contrast, extra = NULL) {
    if (is.null(contrast) || is.na(contrast)) return(extra)
    paste0(c(paste0("Contraste: ", contrast), extra), collapse = "  ·  ")
  }

  make_deg_volcano_plot <- function(df, fdr_thr = 0.05, lfc_thr = 0,
                                    contrast = NULL) {
    ax <- lfc_axis(df)
    df$x <- ax$values
    df$minus_log10_p <- -log10(pmax(df$pvalue, .Machine$double.xmin))
    # La significacion es padj <= FDR y nada más: cuando hay umbral de
    # fold-change ya está dentro del test, así que volver a cortar por |log2FC|
    # aquí sería el filtro post-hoc que estamos eliminando.
    sig <- !is.na(df$padj) & df$padj <= fdr_thr
    df$significant <- ifelse(is.na(sig), FALSE, sig)
    df$color <- ifelse(df$significant, "Significativo", "No significativo")
    fdr_y <- -log10(pmax(fdr_thr, .Machine$double.xmin))
    shapes <- list(
      list(type = "line", xref = "paper", x0 = 0, x1 = 1,
           y0 = fdr_y, y1 = fdr_y,
           line = list(dash = "dot", color = "#A8DADC"))
    )
    if (has_lfc_threshold(lfc_thr)) {
      shapes <- c(shapes, list(
        list(type = "line", x0 = -lfc_thr, x1 = -lfc_thr, yref = "paper",
             y0 = 0, y1 = 1, line = list(dash = "dot", color = "#F4A6A6")),
        list(type = "line", x0 = lfc_thr, x1 = lfc_thr, yref = "paper",
             y0 = 0, y1 = 1, line = list(dash = "dot", color = "#F4A6A6"))
      ))
    }
    p <- plotly::plot_ly(
      df, x = ~x, y = ~minus_log10_p, color = ~color,
      colors = DEG_SIG_COLORS,
      type = "scatter", mode = "markers",
      text = ~paste0("Gen: ", gene, "<br>", ax$label, ": ", round(x, 3),
                     "<br>padj: ", signif(padj, 3)),
      hoverinfo = "text",
      marker = list(size = 6, opacity = 0.7)
    ) |>
      plotly::layout(
        title = list(text = plot_title(contrast,
                       if (has_lfc_threshold(lfc_thr))
                         paste0("umbral del test |log2FC| > ", lfc_thr)),
                     font = list(size = 12)),
        xaxis = list(title = ax$label),
        yaxis = list(title = "-log10(pvalue)"),
        shapes = shapes
      )
    top_sig <- df[df$significant & !is.na(df$gene), ]
    if (nrow(top_sig) > 0) {
      top_sig <- top_sig[order(top_sig$padj, na.last = TRUE), ]
      top_sig <- head(top_sig, 10)
      p <- p |> plotly::add_annotations(
        data = top_sig,
        x = ~x, y = ~minus_log10_p,
        text = ~gene, showarrow = TRUE,
        arrowhead = 2, arrowsize = 0.5, arrowwidth = 1,
        arrowcolor = "#60756A",
        font = list(size = 9, color = "#20332A"),
        ax = 20, ay = -20
      )
    }
    p
  }

  make_deg_ma_plot <- function(df, fdr_thr = 0.05, contrast = NULL) {
    ax <- lfc_axis(df)
    df$y <- ax$values
    df$significant <- !is.na(df$padj) & df$padj <= fdr_thr
    df$color <- ifelse(df$significant, "Significativo", "No significativo")
    df$log_base <- log10(pmax(df$baseMean, 1))
    plotly::plot_ly(
      df, x = ~log_base, y = ~y, color = ~color,
      colors = DEG_SIG_COLORS,
      type = "scatter", mode = "markers",
      text = ~paste0("Gen: ", gene, "<br>baseMean: ", signif(baseMean, 3),
                     "<br>", ax$label, ": ", round(y, 3),
                     "<br>padj: ", signif(padj, 3)),
      hoverinfo = "text",
      marker = list(size = 6, opacity = 0.7)
    ) |>
      plotly::layout(
        title = list(text = plot_title(contrast), font = list(size = 12)),
        xaxis = list(title = "log10(baseMean)"),
        yaxis = list(title = ax$label),
        shapes = list(
          list(type = "line", xref = "paper", x0 = 0, x1 = 1,
               y0 = 0, y1 = 0,
               line = list(dash = "dash", color = "#7BBF9A", width = 1))
        )
      )
  }

  # ── Volcano ────────────────────────────────────────────────────────────────
  output$deg_volcano_plot <- plotly::renderPlotly({
    req(state$deg_rv$results)
    make_deg_volcano_plot(
      state$deg_rv$results,
      fdr_thr  = state$deg_rv$fdr %||% 0.05,
      lfc_thr  = state$deg_rv$lfc_threshold %||% 0,
      contrast = state$deg_rv$contrast
    )
  })

  # ── MA plot ────────────────────────────────────────────────────────────────
  output$deg_ma_plot <- plotly::renderPlotly({
    req(state$deg_rv$results)
    make_deg_ma_plot(state$deg_rv$results,
                     fdr_thr  = state$deg_rv$fdr %||% 0.05,
                     contrast = state$deg_rv$contrast)
  })

  # ── PCA ────────────────────────────────────────────────────────────────────
  # El selector se rellena con las columnas del samplesheet del ajuste, para que
  # se pueda colorear por cualquier covariable y no solo por la condición.
  observe({
    m <- state$deg_rv$meta
    if (is.null(m) || !nrow(m)) return()
    # Se descartan los identificadores únicos por muestra (sample_id y similares):
    # colorear por ellos da un color por punto, no informa de ningun agrupamiento
    # y además agota la paleta.
    is_id <- vapply(names(m), function(v) {
      length(unique(as.character(m[[v]]))) >= nrow(m) && !is.numeric(m[[v]])
    }, logical(1))
    ch <- names(m)[!is_id]
    if (!length(ch)) ch <- names(m)
    sel <- isolate(input$deg_pca_color)
    updateSelectInput(session, "deg_pca_color", choices = ch,
                      selected = if (!is.null(sel) && sel %in% ch) sel
                                 else if ("condition" %in% ch) "condition" else ch[1])
  })

  make_deg_pca_plot <- function(color_col = NULL) {
    if (is.null(state$deg_rv$vst_mat) || is.null(state$deg_rv$meta)) {
      return(plotly_message("Sin datos para el PCA."))
    }
    pcd <- pca_data(state$deg_rv$vst_mat, state$deg_rv$meta, ntop = 500)
    if (is.null(pcd) || !nrow(pcd)) return(plotly_message("No se pudo calcular el PCA."))
    ve <- attr(pcd, "var_explained")
    cc <- if (!is.null(color_col) && nzchar(color_col %||% "") &&
              color_col %in% names(pcd)) color_col
          else if ("condition" %in% names(pcd)) "condition" else "sample_id"
    # Una covariable continua (edad, dosis, tiempo, una variable sustituta) se deja
    # numérica para que plotly use una escala de color continua; forzarla a
    # discreta daría un color por valor y no se leeria nada.
    raw <- pcd[[cc]]
    pcd$color <- if (is.numeric(raw) && length(unique(raw)) > 5) raw else as.character(raw)
    p <- plotly::plot_ly(
      pcd, x = ~PC1, y = ~PC2, color = ~color,
      type = "scatter", mode = "markers",
      text = ~paste0("Muestra: ", sample_id, "<br>", cc, ": ", color),
      hoverinfo = "text",
      marker = list(size = 11, opacity = 0.85)
    ) |>
      plotly::layout(
        title = list(text = paste0("Color: ", cc), font = list(size = 12)),
        xaxis = list(title = paste0("PC1 (", round(100 * (ve[1] %||% NA), 1), "%)")),
        yaxis = list(title = paste0("PC2 (", round(100 * (ve[2] %||% NA), 1), "%)"))
      )
    p
  }

  output$deg_pca_plot <- plotly::renderPlotly({
    req(state$deg_rv$vst_mat, state$deg_rv$meta)
    make_deg_pca_plot(input$deg_pca_color)
  })

  # ── Heatmaps (gráficos base) ───────────────────────────────────────────────
  # Estos dos se dibujan con gráficos base, no con plotly, así que no pueden usar
  # `plotly_download()`. Se factorizan en helpers para que el render y la descarga
  # produzcan exactamente el mismo gráfico: antes el bloque de pheatmap estaba
  # duplicado literalmente en los cuatro sitios, y cualquier cambio de paleta o de
  # anotación había que hacerlo cuatro veces.
  draw_deg_heatmap <- function(n = 30) {
    if (is.null(state$deg_rv$vst_mat)) {
      plot.new(); title("Sin datos para el heatmap."); return(invisible())
    }
    genes <- top_var_genes(state$deg_rv$vst_mat, n = n)
    if (!length(genes)) {
      plot.new(); title("No hay genes para el heatmap."); return(invisible())
    }
    m <- state$deg_rv$vst_mat[genes, , drop = FALSE]
    pal <- grDevices::colorRampPalette(c("#A8DADC", "#FFFFFF", "#F4A6A6"))(50)
    if (isTRUE(HAS_PHEATMAP)) {
      meta <- state$deg_rv$meta
      ann <- if (!is.null(meta) && "condition" %in% names(meta)) {
        data.frame(condition = meta$condition, row.names = meta$sample_id,
                   stringsAsFactors = FALSE)
      } else NULL
      pheatmap::pheatmap(m, scale = "row", annotation_col = ann,
                         show_rownames = TRUE, show_colnames = TRUE,
                         color = pal, silent = FALSE)
    } else {
      stats::heatmap(m, scale = "row", col = pal)
    }
  }

  draw_deg_dist_heatmap <- function() {
    dm <- if (is.null(state$deg_rv$vst_mat)) NULL
          else sample_distance_matrix(state$deg_rv$vst_mat)
    if (is.null(dm)) {
      plot.new(); title("Sin matriz de distancias."); return(invisible())
    }
    pal <- grDevices::colorRampPalette(c("#244B34", "#A8DDB8", "#FFFFFF"))(50)
    if (isTRUE(HAS_PHEATMAP)) {
      pheatmap::pheatmap(dm, color = pal,
                         clustering_distance_rows = stats::as.dist(dm),
                         clustering_distance_cols = stats::as.dist(dm),
                         silent = FALSE)
    } else {
      stats::heatmap(dm, symm = TRUE, col = pal)
    }
  }

  #' Descarga PNG de un gráfico de gráficos base.
  png_download <- function(prefix, draw_fun, width = 1200, height = 900) {
    downloadHandler(
      filename = function() paste0(prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content = function(f) {
        grDevices::png(f, width = width, height = height, res = 120)
        on.exit(grDevices::dev.off(), add = TRUE)
        draw_fun()
      }
    )
  }

  output$deg_heatmap <- renderPlot({
    req(state$deg_rv$vst_mat)
    draw_deg_heatmap(input$deg_heatmap_topn %||% 30)
  })

  output$deg_dist_heatmap <- renderPlot({
    req(state$deg_rv$vst_mat)
    draw_deg_dist_heatmap()
  })

  # ── Descargas de esta sección ──────────────────────────────────────────────
  # Cada descarga vive junto al render que reutiliza, en lugar de en un bloque
  # "Descargas" al final del fichero: así es imposible que se desincronicen, que
  # es lo que le había pasado al PCA.
  output$download_deg_table <- csv_download(
    "deg_table_filtered",
    function() {
      df <- ctx$deg_filtered()
      if (is.null(df)) return(message_df("Sin resultados DEG."))
      df
    }
  )

  output$download_deg_volcano_plot <- plotly_download(
    "deg_volcano",
    function() {
      df <- state$deg_rv$results
      if (is.null(df)) return(plotly_message("Sin resultados DEG."))
      make_deg_volcano_plot(df,
                            fdr_thr  = state$deg_rv$fdr %||% 0.05,
                            lfc_thr  = state$deg_rv$lfc_threshold %||% 0,
                            contrast = state$deg_rv$contrast)
    }
  )

  output$download_deg_ma_plot <- plotly_download(
    "deg_ma",
    function() {
      df <- state$deg_rv$results
      if (is.null(df)) return(plotly_message("Sin resultados DEG."))
      make_deg_ma_plot(df, fdr_thr = state$deg_rv$fdr %||% 0.05,
                       contrast = state$deg_rv$contrast)
    }
  )

  output$download_deg_pca_plot <- plotly_download(
    "deg_pca", function() make_deg_pca_plot(input$deg_pca_color))

  output$download_deg_heatmap <- png_download(
    "deg_heatmap", function() draw_deg_heatmap(input$deg_heatmap_topn %||% 30))

  output$download_deg_dist_heatmap <- png_download(
    "deg_dist_heatmap", function() draw_deg_dist_heatmap(), width = 900, height = 800)

  invisible(NULL)
}
