#' server_tab_deg.R
#' Logica server de la Tab 4 (Expresion diferencial).
#'   - Carga counts desde 3 fuentes (run actual, run guardada, upload).
#'   - Editor inline de metadatos con autocompletado.
#'   - run_deg() cacheado en state$deg_rv$results.
#'   - Filtros aplicados al vuelo en reactivos derivados.
#'   - Enriquecimiento GO/KEGG bajo demanda.

server_tab_deg <- function(input, output, session, state) {

  outputs_dir <- state$outputs_dir

  # ── Contenido del nav_panel "4. Expresion diferencial" ────────────────────
  output$tab_deg_content <- renderUI({
    ui_tab_deg()
  })

  # ── Helpers de lectura de matriz subida ────────────────────────────────────
  read_uploaded_counts <- function(path) {
    if (!nzchar(path) || !file.exists(path)) return(NULL)
    sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
    df <- tryCatch(
      read.delim(path, sep = sep, check.names = FALSE, stringsAsFactors = FALSE,
                 row.names = 1, comment.char = ""),
      error = function(e) NULL
    )
    if (is.null(df) || !nrow(df) || !ncol(df)) return(NULL)
    cm <- as.matrix(df)
    storage.mode(cm) <- "numeric"
    cm
  }

  read_uploaded_meta <- function(path) {
    if (!nzchar(path) || !file.exists(path)) return(NULL)
    sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
    tryCatch(
      read.delim(path, sep = sep, check.names = FALSE, stringsAsFactors = FALSE,
                 comment.char = ""),
      error = function(e) NULL
    )
  }

  # ── Run guardada selector ─────────────────────────────────────────────────
  output$deg_saved_run_ui <- renderUI({
    state$results_refresh()
    choices <- result_choices(outputs_dir)
    if (!length(choices)) {
      return(div(class = "alert alert-info mb-0",
                 icon("folder-open"), " No hay ejecuciones guardadas en outputs/."))
    }
    selectizeInput("selected_deg_run_dir", "Ejecucion guardada",
                   choices = choices, selected = unname(choices)[1],
                   options = list(placeholder = "Escribe para buscar..."),
                   width = "100%")
  })

  # ── reactive: matriz de conteos segun fuente ──────────────────────────────
  deg_counts_source <- reactive({
    src <- input$deg_source %||% "current"
    if (identical(src, "current")) {
      cm <- state$data_rv$count_matrix
      if (!is.null(cm) && length(cm)) return(as.matrix(cm))
      # Fallback: intentamos cargar desde la run actual
      p <- state$run_params_rv()
      if (length(p) && nzchar(p$output_dir %||% "") && nzchar(p$tool %||% "")) {
        return(tryCatch(load_counts_from_workflow(p$output_dir, p$tool),
                        error = function(e) NULL))
      }
      return(NULL)
    }
    if (identical(src, "saved")) {
      sel <- input$selected_deg_run_dir %||% ""
      if (!nzchar(sel) || !dir.exists(sel)) return(NULL)
      p <- infer_result_params(sel, state$workflow_path)
      return(tryCatch(load_counts_from_workflow(sel, p$tool %||% ""),
                      error = function(e) NULL))
    }
    if (identical(src, "upload")) {
      up <- input$deg_counts_upload
      if (is.null(up) || !nrow(up)) return(NULL)
      return(read_uploaded_counts(up$datapath))
    }
    NULL
  })

  # ── reactiveValues: edicion inline del samplesheet ────────────────────────
  meta_rv <- reactiveVal(NULL)

  # Cuando llega un upload de samplesheet, lo cargamos
  observeEvent(input$deg_meta_upload, {
    df <- read_uploaded_meta(input$deg_meta_upload$datapath)
    if (is.null(df)) {
      showNotification("No se pudo leer el samplesheet.", type = "error")
      return()
    }
    if (!"sample_id" %in% names(df)) {
      # tolerar variantes habituales
      cand <- intersect(c("Sample", "sample", "Muestra", "sample_name"), names(df))
      if (length(cand)) names(df)[names(df) == cand[1]] <- "sample_id"
    }
    meta_rv(df)
  })

  # Sincronizar con muestras detectadas en counts (rellena sample_id, deja condition vacia)
  observeEvent(input$deg_meta_sync_btn, {
    cm <- deg_counts_source()
    if (is.null(cm) || !ncol(cm)) {
      showNotification("Carga primero una matriz de conteos.", type = "warning")
      return()
    }
    cur <- meta_rv()
    new_df <- data.frame(
      sample_id = colnames(cm),
      condition = "",
      batch     = "",
      stringsAsFactors = FALSE
    )
    if (!is.null(cur) && nrow(cur)) {
      # preservamos condition/batch ya editados para sample_ids comunes
      common <- intersect(new_df$sample_id, cur$sample_id)
      if (length(common)) {
        m_idx <- match(common, new_df$sample_id)
        c_idx <- match(common, cur$sample_id)
        for (col in c("condition", "batch")) {
          if (col %in% names(cur)) new_df[[col]][m_idx] <- as.character(cur[[col]][c_idx])
        }
      }
    }
    meta_rv(new_df)
  })

  observeEvent(input$deg_meta_add_row_btn, {
    cur <- meta_rv()
    if (is.null(cur)) cur <- data.frame(sample_id = character(0),
                                        condition = character(0),
                                        batch = character(0),
                                        stringsAsFactors = FALSE)
    cur <- rbind(cur, setNames(as.data.frame(as.list(rep("", ncol(cur))),
                                             stringsAsFactors = FALSE), names(cur)))
    meta_rv(cur)
  })

  # Renderizado editable
  output$deg_meta_table <- renderDT({
    df <- meta_rv()
    if (is.null(df) || !nrow(df)) {
      df <- data.frame(sample_id = character(0),
                       condition = character(0),
                       batch = character(0),
                       stringsAsFactors = FALSE)
    }
    datatable(df, editable = TRUE, rownames = FALSE,
              options = list(pageLength = 8, scrollX = TRUE, dom = "ftip"))
  }, server = TRUE)

  # Aplicar ediciones inline
  observeEvent(input$deg_meta_table_cell_edit, {
    info <- input$deg_meta_table_cell_edit
    df <- meta_rv()
    if (is.null(df)) return()
    r <- info$row
    c <- info$col + 1L  # DT col is 0-indexed when rownames=FALSE
    if (r > 0 && r <= nrow(df) && c > 0 && c <= ncol(df)) {
      df[r, c] <- DT::coerceValue(info$value, df[r, c])
      meta_rv(df)
    }
  })

  # Actualiza selectInputs de condition/batch/ref_level segun columnas del meta
  observe({
    df <- meta_rv()
    if (is.null(df) || !nrow(df)) return()
    cond_cols <- intersect(names(df), c("condition", "Condition", "group", "Group"))
    if (!length(cond_cols)) cond_cols <- setdiff(names(df), "sample_id")
    if (!length(cond_cols)) cond_cols <- "condition"
    updateSelectInput(session, "deg_condition_col",
                      choices = cond_cols,
                      selected = isolate(input$deg_condition_col) %||% cond_cols[1])

    batch_candidates <- setdiff(names(df), c("sample_id", isolate(input$deg_condition_col) %||% "condition"))
    if (!length(batch_candidates)) batch_candidates <- "batch"
    updateSelectInput(session, "deg_batch_col",
                      choices = batch_candidates,
                      selected = isolate(input$deg_batch_col) %||% batch_candidates[1])

    cc <- isolate(input$deg_condition_col) %||% cond_cols[1]
    if (cc %in% names(df)) {
      lvls <- unique(df[[cc]][!is.na(df[[cc]]) & nzchar(as.character(df[[cc]]))])
      if (length(lvls)) {
        updateSelectInput(session, "deg_ref_level",
                          choices = lvls,
                          selected = isolate(input$deg_ref_level) %||% lvls[1])
      }
    }
  })

  # ── Lanzar DEG ─────────────────────────────────────────────────────────────
  observeEvent(input$run_deg_btn, {
    cm <- deg_counts_source()
    df <- meta_rv()
    method <- input$deg_method %||% "DESeq2"

    if (is.null(cm) || !ncol(cm)) {
      showNotification("No hay matriz de conteos cargada.", type = "error"); return()
    }
    if (is.null(df) || !nrow(df)) {
      showNotification("Carga o crea un samplesheet.", type = "error"); return()
    }

    # Renombrar columna de condicion a 'condition' para el motor
    cond_col <- input$deg_condition_col %||% "condition"
    if (!cond_col %in% names(df)) {
      showNotification(paste0("La columna '", cond_col, "' no esta en el samplesheet."),
                       type = "error"); return()
    }
    df_work <- df
    if (cond_col != "condition") {
      df_work$condition <- df_work[[cond_col]]
    }

    val <- validate_samplesheet(df_work, colnames(cm))
    if (!val$ok) {
      showNotification(paste(val$errors, collapse = " | "), type = "error", duration = 12)
      return()
    }

    aligned <- align_counts_to_metadata(cm, df_work)
    if (length(aligned$warnings)) {
      for (w in aligned$warnings) showNotification(w, type = "warning", duration = 8)
    }
    cm_aln <- aligned$counts
    meta_aln <- aligned$meta

    if (is.null(cm_aln) || ncol(cm_aln) < 2) {
      showNotification("Tras alinear muestras quedan menos de 2 columnas.", type = "error"); return()
    }

    # Prefiltrado
    mc <- input$deg_min_count %||% 10
    ms <- input$deg_min_samples
    if (is.null(ms) || !is.finite(ms) || ms < 1) ms <- NULL
    cm_f <- prefilter_counts(cm_aln, min_count = mc, min_samples = ms)
    if (is.null(cm_f) || !nrow(cm_f)) {
      showNotification("Tras prefiltrar no quedan filas. Reduce los umbrales.", type = "error"); return()
    }

    ref <- input$deg_ref_level
    batch <- if (isTRUE(input$deg_use_batch)) input$deg_batch_col else NULL

    showNotification(paste0("Corriendo DEG (", method, "). Esto puede tardar..."),
                     type = "message", duration = 4)

    res <- tryCatch(
      run_deg(cm_f, meta_aln, method = method, ref_level = ref, batch = batch),
      error = function(e) list(table = NULL, error = conditionMessage(e), method = method)
    )

    if (is.null(res$table)) {
      showNotification(paste0("Error en ", method, ": ", res$error %||% "fallo desconocido"),
                       type = "error", duration = 12)
      return()
    }

    # Cache de transformacion para visualizacion
    vst_mat <- tryCatch(vst_or_rlog(cm_f, meta_aln), error = function(e) NULL)

    state$deg_rv$counts  <- cm_f
    state$deg_rv$meta    <- meta_aln
    state$deg_rv$method  <- method
    state$deg_rv$results <- res$table
    state$deg_rv$vst_mat <- vst_mat
    state$deg_rv$run_at  <- Sys.time()

    showNotification(paste0("DEG completado (", method, "): ",
                            nrow(res$table), " filas."),
                     type = "default", duration = 6)
  })

  output$deg_status_text <- renderText({
    if (is.null(state$deg_rv$results)) return("Sin ejecucion DEG. Pulsa 'Lanzar DEG'.")
    tab <- state$deg_rv$results
    fdr_thr <- input$deg_fdr_cutoff %||% 0.05
    lfc_thr <- input$deg_log2fc_cutoff %||% 1
    sig <- !is.na(tab$padj) & tab$padj <= fdr_thr & abs(tab$log2FC) >= lfc_thr
    n_up   <- sum(sig & !is.na(tab$log2FC) & tab$log2FC > 0, na.rm = TRUE)
    n_down <- sum(sig & !is.na(tab$log2FC) & tab$log2FC < 0, na.rm = TRUE)
    paste0(
      "Motor: ", state$deg_rv$method, "\n",
      "Total: ", nrow(tab), " genes\n",
      "FDR<", fdr_thr, ", |LFC|>", lfc_thr, ": ",
      sum(sig, na.rm = TRUE), " (", n_up, " up / ", n_down, " down)\n",
      "Ultima ejecucion: ", format(state$deg_rv$run_at, "%Y-%m-%d %H:%M:%S")
    )
  })

  # ── Tabla filtrada (reactivo derivado, rapido) ─────────────────────────────
  deg_filtered <- reactive({
    tab <- state$deg_rv$results
    if (is.null(tab)) return(NULL)
    apply_deg_filters(
      tab,
      fdr        = input$deg_fdr_cutoff %||% 0.05,
      abs_log2fc = input$deg_log2fc_cutoff %||% 1,
      base_mean  = input$deg_basemean_cutoff %||% 0
    )
  })

  # ── UI de resultados ───────────────────────────────────────────────────────
  output$deg_results_ui <- renderUI({
    if (is.null(state$deg_rv$results)) {
      return(div(class = "alert alert-info mt-3",
                 icon("circle-info"),
                 " Configura datos y metadatos, y pulsa 'Lanzar DEG' para ver los resultados."))
    }
    ui_tab_deg_results()
  })

  # ── Tabla DT ───────────────────────────────────────────────────────────────
  output$deg_table <- renderDT({
    df <- deg_filtered()
    if (is.null(df)) return(dt_table(message_df("Sin resultados DEG.")))
    df_r <- df
    if ("log2FC" %in% names(df_r)) {
      df_r$direction <- ifelse(df_r$log2FC >= 0, "Up", "Down")
      df_r <- df_r[, c("gene", "direction", setdiff(names(df_r), c("gene", "direction"))), drop = FALSE]
    }
    num_cols <- intersect(c("baseMean", "log2FC", "lfcSE", "stat", "pvalue", "padj"), names(df_r))
    for (nm in num_cols) df_r[[nm]] <- signif(df_r[[nm]], 4)
    dt_table(df_r, page_length = 15, filter = "top")
  })

  # ── Helpers de ploteo (reutilizados en render y en descarga) ──────────────
  make_deg_volcano_plot <- function(df, fdr_thr = 0.05, lfc_thr = 1) {
    df$minus_log10_p <- -log10(pmax(df$pvalue, .Machine$double.xmin))
    sig <- !is.na(df$padj) & df$padj <= fdr_thr & abs(df$log2FC) >= lfc_thr
    df$significant <- ifelse(is.na(sig), FALSE, sig)
    df$color <- ifelse(df$significant, "Significativo", "No significativo")
    fdr_y <- -log10(pmax(fdr_thr, .Machine$double.xmin))
    p <- plotly::plot_ly(
      df, x = ~log2FC, y = ~minus_log10_p, color = ~color,
      colors = c("Significativo" = "#7BBF9A", "No significativo" = "#C0C0C0"),
      type = "scatter", mode = "markers",
      text = ~paste0("Gen: ", gene, "<br>log2FC: ", round(log2FC, 3),
                     "<br>padj: ", signif(padj, 3)),
      hoverinfo = "text",
      marker = list(size = 6, opacity = 0.7)
    ) |>
      plotly::layout(
        xaxis = list(title = "log2 Fold Change"),
        yaxis = list(title = "-log10(pvalue)"),
        shapes = list(
          list(type = "line", x0 = -lfc_thr, x1 = -lfc_thr, yref = "paper",
               y0 = 0, y1 = 1, line = list(dash = "dot", color = "#F4A6A6")),
          list(type = "line", x0 = lfc_thr, x1 = lfc_thr, yref = "paper",
               y0 = 0, y1 = 1, line = list(dash = "dot", color = "#F4A6A6")),
          list(type = "line", xref = "paper", x0 = 0, x1 = 1,
               y0 = fdr_y, y1 = fdr_y,
               line = list(dash = "dot", color = "#A8DADC"))
        )
      )
    top_sig <- df[df$significant & !is.na(df$gene), ]
    if (nrow(top_sig) > 0) {
      top_sig <- top_sig[order(top_sig$padj, na.last = TRUE), ]
      top_sig <- head(top_sig, 10)
      p <- p |> plotly::add_annotations(
        data = top_sig,
        x = ~log2FC, y = ~minus_log10_p,
        text = ~gene, showarrow = TRUE,
        arrowhead = 2, arrowsize = 0.5, arrowwidth = 1,
        arrowcolor = "#60756A",
        font = list(size = 9, color = "#20332A"),
        ax = 20, ay = -20
      )
    }
    p
  }

  make_deg_ma_plot <- function(df, fdr_thr = 0.05) {
    df$significant <- !is.na(df$padj) & df$padj <= fdr_thr
    df$color <- ifelse(df$significant, "Significativo", "No significativo")
    df$log_base <- log10(pmax(df$baseMean, 1))
    plotly::plot_ly(
      df, x = ~log_base, y = ~log2FC, color = ~color,
      colors = c("Significativo" = "#7BBF9A", "No significativo" = "#C0C0C0"),
      type = "scatter", mode = "markers",
      text = ~paste0("Gen: ", gene, "<br>baseMean: ", signif(baseMean, 3),
                     "<br>log2FC: ", round(log2FC, 3),
                     "<br>padj: ", signif(padj, 3)),
      hoverinfo = "text",
      marker = list(size = 6, opacity = 0.7)
    ) |>
      plotly::layout(
        xaxis = list(title = "log10(baseMean)"),
        yaxis = list(title = "log2 Fold Change"),
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
      fdr_thr = input$deg_fdr_cutoff %||% 0.05,
      lfc_thr = input$deg_log2fc_cutoff %||% 1
    )
  })

  # ── MA plot ────────────────────────────────────────────────────────────────
  output$deg_ma_plot <- plotly::renderPlotly({
    req(state$deg_rv$results)
    make_deg_ma_plot(state$deg_rv$results, fdr_thr = input$deg_fdr_cutoff %||% 0.05)
  })

  # ── PCA ────────────────────────────────────────────────────────────────────
  output$deg_pca_plot <- plotly::renderPlotly({
    req(state$deg_rv$vst_mat, state$deg_rv$meta)
    pcd <- pca_data(state$deg_rv$vst_mat, state$deg_rv$meta, ntop = 500)
    if (is.null(pcd) || !nrow(pcd)) return(plotly_message("No se pudo calcular el PCA."))
    ve <- attr(pcd, "var_explained")
    color_col <- if ("condition" %in% names(pcd)) "condition" else "sample_id"
    pcd$color <- pcd[[color_col]]
    plotly::plot_ly(
      pcd, x = ~PC1, y = ~PC2, color = ~color,
      type = "scatter", mode = "markers",
      text = ~paste0("Muestra: ", sample_id,
                     "<br>", color_col, ": ", color),
      hoverinfo = "text",
      marker = list(size = 11, opacity = 0.85)
    ) |>
      plotly::layout(
        xaxis = list(title = paste0("PC1 (", round(100 * (ve[1] %||% NA), 1), "%)")),
        yaxis = list(title = paste0("PC2 (", round(100 * (ve[2] %||% NA), 1), "%)"))
      )
  })

  # ── Heatmap top-N ──────────────────────────────────────────────────────────
  output$deg_heatmap <- renderPlot({
    req(state$deg_rv$vst_mat)
    n <- input$deg_heatmap_topn %||% 30
    genes <- top_var_genes(state$deg_rv$vst_mat, n = n)
    if (!length(genes)) {
      plot.new(); title("No hay genes para el heatmap."); return(invisible())
    }
    m <- state$deg_rv$vst_mat[genes, , drop = FALSE]
    if (requireNamespace("pheatmap", quietly = TRUE)) {
      ann <- NULL
      meta <- state$deg_rv$meta
      if (!is.null(meta) && "condition" %in% names(meta)) {
        ann <- data.frame(condition = meta$condition,
                          row.names = meta$sample_id,
                          stringsAsFactors = FALSE)
      }
      pheatmap::pheatmap(m, scale = "row", annotation_col = ann,
                         show_rownames = TRUE, show_colnames = TRUE,
                         color = grDevices::colorRampPalette(c("#A8DADC", "#FFFFFF", "#F4A6A6"))(50),
                         silent = FALSE)
    } else {
      stats::heatmap(m, scale = "row",
                     col = grDevices::colorRampPalette(c("#A8DADC", "#FFFFFF", "#F4A6A6"))(50))
    }
  })

  # ── Distancia entre muestras ───────────────────────────────────────────────
  output$deg_dist_heatmap <- renderPlot({
    req(state$deg_rv$vst_mat)
    dm <- sample_distance_matrix(state$deg_rv$vst_mat)
    if (is.null(dm)) { plot.new(); title("Sin matriz de distancias."); return(invisible()) }
    if (requireNamespace("pheatmap", quietly = TRUE)) {
      pheatmap::pheatmap(dm,
                         color = grDevices::colorRampPalette(c("#244B34", "#A8DDB8", "#FFFFFF"))(50),
                         clustering_distance_rows = stats::as.dist(dm),
                         clustering_distance_cols = stats::as.dist(dm),
                         silent = FALSE)
    } else {
      stats::heatmap(dm, symm = TRUE,
                     col = grDevices::colorRampPalette(c("#244B34", "#A8DDB8", "#FFFFFF"))(50))
    }
  })

  # ── Enriquecimiento ────────────────────────────────────────────────────────
  enrich_rv <- reactiveVal(NULL)

  observeEvent(input$deg_run_enrich_btn, {
    req(state$deg_rv$results)
    df <- deg_filtered()
    if (is.null(df) || !nrow(df)) {
      showNotification("No hay genes significativos con los filtros actuales.",
                       type = "warning"); return()
    }
    ont <- input$deg_ontology %||% "BP"
    universe <- state$deg_rv$results$gene
    genes <- df$gene
    if (identical(ont, "KEGG")) {
      org_code <- trimws(input$deg_kegg_organism %||% "eco")
      if (!nzchar(org_code)) org_code <- "eco"
      res <- run_enrichment_kegg(genes, organism = org_code)
    } else {
      org <- if (isTRUE(HAS_ORGECDB)) "org.EcK12.eg.db" else NULL
      res <- run_enrichment_go(genes, universe = universe,
                               OrgDb = org, ont = ont, keyType = "SYMBOL")
    }
    if (is.null(res$table)) {
      showNotification(paste0("Enriquecimiento sin resultados: ", res$error %||% "—"),
                       type = "warning", duration = 8)
      enrich_rv(NULL); return()
    }
    enrich_rv(res$table)
  })

  output$deg_enrich_dotplot <- plotly::renderPlotly({
    df <- enrich_rv()
    if (is.null(df) || !nrow(df))
      return(plotly_message("Pulsa 'Calcular' para correr el enriquecimiento."))
    top <- enrichment_dotplot_data(df, top_n = 15)
    if (is.null(top) || !nrow(top))
      return(plotly_message("Sin terminos enriquecidos."))
    plotly::plot_ly(
      top,
      x = ~ -log10(p.adjust),
      y = ~stats::reorder(Description, -log10(p.adjust)),
      type = "scatter", mode = "markers",
      size = ~Count,
      marker = list(color = ~ -log10(p.adjust),
                    colorscale = list(c(0, "#A8DADC"), c(1, "#7BBF9A")),
                    sizemode = "area"),
      text = ~paste0("Termino: ", Description,
                     "<br>Count: ", Count,
                     "<br>p.adjust: ", signif(p.adjust, 3)),
      hoverinfo = "text"
    ) |>
      plotly::layout(
        xaxis = list(title = "-log10(p.adjust)"),
        yaxis = list(title = "", automargin = TRUE),
        margin = list(l = 220)
      )
  })

  output$deg_enrich_table <- renderDT({
    df <- enrich_rv()
    if (is.null(df) || !nrow(df)) return(dt_table(message_df("Sin enriquecimiento calculado.")))
    df_r <- df
    for (nm in intersect(c("pvalue", "p.adjust", "qvalue"), names(df_r)))
      df_r[[nm]] <- signif(df_r[[nm]], 3)
    dt_table(df_r, page_length = 15, filter = "top")
  })

  # ── Descargas ──────────────────────────────────────────────────────────────
  output$download_deg_table <- csv_download(
    "deg_table_filtered",
    function() {
      df <- deg_filtered()
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
                            fdr_thr = input$deg_fdr_cutoff %||% 0.05,
                            lfc_thr = input$deg_log2fc_cutoff %||% 1)
    }
  )
  output$download_deg_ma_plot <- plotly_download(
    "deg_ma",
    function() {
      df <- state$deg_rv$results
      if (is.null(df)) return(plotly_message("Sin resultados DEG."))
      make_deg_ma_plot(df, fdr_thr = input$deg_fdr_cutoff %||% 0.05)
    }
  )
  output$download_deg_pca_plot <- plotly_download(
    "deg_pca",
    function() {
      pcd <- pca_data(state$deg_rv$vst_mat, state$deg_rv$meta, ntop = 500)
      if (is.null(pcd)) return(plotly_message("Sin PCA."))
      plotly::plot_ly(pcd, x = ~PC1, y = ~PC2, type = "scatter", mode = "markers")
    }
  )

  output$download_deg_heatmap <- downloadHandler(
    filename = function() paste0("deg_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(f) {
      req(state$deg_rv$vst_mat)
      n <- input$deg_heatmap_topn %||% 30
      genes <- top_var_genes(state$deg_rv$vst_mat, n = n)
      if (!length(genes)) return(invisible())
      m <- state$deg_rv$vst_mat[genes, , drop = FALSE]
      grDevices::png(f, width = 1200, height = 900, res = 120)
      on.exit(grDevices::dev.off(), add = TRUE)
      if (requireNamespace("pheatmap", quietly = TRUE)) {
        ann <- NULL
        meta <- state$deg_rv$meta
        if (!is.null(meta) && "condition" %in% names(meta)) {
          ann <- data.frame(condition = meta$condition, row.names = meta$sample_id, stringsAsFactors = FALSE)
        }
        pheatmap::pheatmap(m, scale = "row", annotation_col = ann,
                           show_rownames = TRUE, show_colnames = TRUE,
                           color = grDevices::colorRampPalette(c("#A8DADC", "#FFFFFF", "#F4A6A6"))(50),
                           silent = FALSE)
      } else {
        stats::heatmap(m, scale = "row",
                       col = grDevices::colorRampPalette(c("#A8DADC", "#FFFFFF", "#F4A6A6"))(50))
      }
    }
  )

  output$download_deg_dist_heatmap <- downloadHandler(
    filename = function() paste0("deg_dist_heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(f) {
      req(state$deg_rv$vst_mat)
      dm <- sample_distance_matrix(state$deg_rv$vst_mat)
      req(!is.null(dm))
      grDevices::png(f, width = 900, height = 800, res = 120)
      on.exit(grDevices::dev.off(), add = TRUE)
      if (requireNamespace("pheatmap", quietly = TRUE)) {
        pheatmap::pheatmap(dm,
                           color = grDevices::colorRampPalette(c("#244B34", "#A8DDB8", "#FFFFFF"))(50),
                           clustering_distance_rows = stats::as.dist(dm),
                           clustering_distance_cols = stats::as.dist(dm),
                           silent = FALSE)
      } else {
        stats::heatmap(dm, symm = TRUE,
                       col = grDevices::colorRampPalette(c("#244B34", "#A8DDB8", "#FFFFFF"))(50))
      }
    }
  )

  output$download_enrich_table <- csv_download(
    "enriquecimiento",
    function() {
      df <- enrich_rv()
      if (is.null(df)) return(message_df("Sin enriquecimiento calculado."))
      df
    }
  )

  invisible(NULL)
}
