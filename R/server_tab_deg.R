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
        # Por defecto, el ultimo nivel contra el primero: reproduce lo que hacia
        # la app antes, pero ahora dicho explicitamente.
        prev_num <- isolate(input$deg_contrast_num)
        prev_den <- isolate(input$deg_contrast_den)
        updateSelectInput(session, "deg_contrast_num", choices = lvls,
                          selected = if (!is.null(prev_num) && prev_num %in% lvls) prev_num
                                     else lvls[length(lvls)])
        updateSelectInput(session, "deg_contrast_den", choices = lvls,
                          selected = if (!is.null(prev_den) && prev_den %in% lvls) prev_den
                                     else lvls[1])
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

    # El denominador del contraste es el nivel de referencia del factor.
    num   <- input$deg_contrast_num
    ref   <- input$deg_contrast_den
    batch <- if (isTRUE(input$deg_use_batch)) input$deg_batch_col else NULL

    if (!is.null(num) && !is.null(ref) && identical(num, ref)) {
      showNotification("El numerador y el denominador del contraste son el mismo nivel.",
                       type = "error"); return()
    }

    # Prefiltrado. En modo automatico pasamos el diseno para que filterByExpr
    # use el tamano del grupo mas pequeno; en manual, el grupo permite derivar
    # min_samples cuando el usuario lo deja vacio.
    pf_mode <- input$deg_prefilter_mode %||% "auto"
    mc <- input$deg_min_count %||% 10
    ms <- input$deg_min_samples
    if (is.null(ms) || !is.finite(ms) || ms < 1) ms <- NULL
    design_for_filter <- tryCatch({
      d <- build_design(meta_aln, ref, batch)
      stats::model.matrix(d$formula, data = d$meta)
    }, error = function(e) NULL)

    cm_f <- prefilter_counts(cm_aln, min_count = mc, min_samples = ms,
                             mode = pf_mode, design = design_for_filter,
                             group = meta_aln$condition)
    if (is.null(cm_f) || !nrow(cm_f)) {
      showNotification("Tras prefiltrar no quedan filas. Reduce los umbrales.", type = "error"); return()
    }
    pf_info <- attr(cm_f, "prefilter")

    fdr_target <- input$deg_fdr_target %||% 0.05
    lfc_thr    <- input$deg_lfc_threshold
    if (is.null(lfc_thr) || !is.finite(lfc_thr) || lfc_thr < 0) lfc_thr <- 0
    do_shrink  <- isTRUE(input$deg_shrink)

    showNotification(paste0("Corriendo DEG (", method, "). Esto puede tardar..."),
                     type = "message", duration = 4)

    res <- tryCatch(
      run_deg(cm_f, meta_aln, method = method, ref_level = ref, batch = batch,
              fdr = fdr_target, lfc_threshold = lfc_thr, shrink = do_shrink,
              contrast_num = num, use_ihw = isTRUE(input$deg_use_ihw)),
      error = function(e) list(table = NULL, error = conditionMessage(e), method = method)
    )

    if (is.null(res$table)) {
      showNotification(paste0("Error en ", method, ": ", res$error %||% "fallo desconocido"),
                       type = "error", duration = 12)
      return()
    }

    # Cache de transformacion para visualizacion
    vst_mat <- tryCatch(vst_or_rlog(cm_f, meta_aln), error = function(e) NULL)

    state$deg_rv$counts        <- cm_f
    state$deg_rv$meta          <- meta_aln
    state$deg_rv$method        <- method
    state$deg_rv$results       <- res$table
    state$deg_rv$vst_mat       <- vst_mat
    state$deg_rv$run_at        <- Sys.time()
    state$deg_rv$fdr           <- fdr_target
    state$deg_rv$lfc_threshold <- lfc_thr
    state$deg_rv$contrast      <- res$contrast
    state$deg_rv$n_levels      <- res$n_levels %||% NA_integer_
    state$deg_rv$shrink        <- res$shrink %||% "ninguno"
    state$deg_rv$prefilter     <- pf_info
    state$deg_rv$padj_method   <- res$padj_method %||% "BH"
    state$deg_rv$disp_data     <- res$disp_data
    state$deg_rv$cooks         <- res$cooks

    if (do_shrink && identical(method, "DESeq2") &&
        identical(state$deg_rv$shrink, "ninguno")) {
      showNotification(
        "No se pudo encoger el log2FC (apeglm/ashr no disponibles).",
        type = "warning", duration = 8
      )
    }
    if (isTRUE(input$deg_use_ihw) && identical(method, "DESeq2") &&
        !identical(state$deg_rv$padj_method, "IHW")) {
      showNotification(
        "No se pudo usar IHW; se ha corregido con Benjamini-Hochberg.",
        type = "warning", duration = 8
      )
    }
    # Si una muestra concentra los outliers de Cook, el problema es de la
    # muestra y no de los genes.
    dom <- res$cooks$dominant %||% NA_character_
    if (!is.na(dom)) {
      showNotification(
        paste0("La muestra '", dom, "' concentra la mayoria de los outliers de ",
               "Cook. Revisala en la pestana de diagnosticos antes de interpretar ",
               "los resultados."),
        type = "warning", duration = 16
      )
    }

    showNotification(paste0("DEG completado (", method, "): ",
                            nrow(res$table), " filas."),
                     type = "default", duration = 6)
  })

  # ── Banner del contraste testeado ──────────────────────────────────────────
  output$deg_contrast_banner <- renderUI({
    if (is.null(state$deg_rv$results)) return(NULL)
    ct <- state$deg_rv$contrast
    nl <- state$deg_rv$n_levels
    lfc_thr <- state$deg_rv$lfc_threshold %||% 0
    test_txt <- if (is.finite(lfc_thr) && lfc_thr > 0) {
      paste0("H0: |log2FC| <= ", lfc_thr, " (umbral dentro del test)")
    } else {
      "H0: log2FC = 0"
    }
    bits <- tagList(
      tags$b("Contraste: "), tags$span(ct %||% "no determinado"),
      tags$span(class = "text-muted",
                paste0("  ·  motor ", state$deg_rv$method,
                       "  ·  FDR objetivo ", state$deg_rv$fdr,
                       " (", state$deg_rv$padj_method %||% "BH", ")",
                       "  ·  ", test_txt))
    )
    # Con el contraste elegido explicitamente ya no hay comparacion oculta, asi
    # que con >2 niveles solo se recuerda que quedan otras por explorar.
    extra <- if (isTRUE(nl > 2)) {
      n_pairs <- nl * (nl - 1) / 2
      tags$div(class = "small mt-1 text-muted",
               paste0("condition tiene ", nl, " niveles: hay ", n_pairs,
                      " contrastes por pares posibles. Cambia numerador o ",
                      "denominador en la tarjeta 3 para testear otro."))
    } else NULL
    div(class = "alert alert-light border py-2 px-3 mb-2",
        icon("circle-info"), " ", bits, extra)
  })

  output$deg_status_text <- renderText({
    if (is.null(state$deg_rv$results)) return("Sin ejecucion DEG. Pulsa 'Lanzar DEG'.")
    tab <- state$deg_rv$results
    # La significacion se lee con el FDR del AJUSTE, no con un deslizador: es el
    # unico nivel para el que el control de FDR calculado es valido.
    fdr_thr <- state$deg_rv$fdr %||% 0.05
    lfc_thr <- state$deg_rv$lfc_threshold %||% 0
    sig <- !is.na(tab$padj) & tab$padj <= fdr_thr
    n_up   <- sum(sig & !is.na(tab$log2FC) & tab$log2FC > 0, na.rm = TRUE)
    n_down <- sum(sig & !is.na(tab$log2FC) & tab$log2FC < 0, na.rm = TRUE)
    pf <- state$deg_rv$prefilter
    pf_txt <- if (is.null(pf)) "" else paste0(
      "Prefiltrado (", pf$mode, "): ", fmt_int(pf$n_before), " -> ",
      fmt_int(pf$n_after), " genes\n"
    )
    paste0(
      "Motor: ", state$deg_rv$method, "\n",
      "Contraste: ", state$deg_rv$contrast %||% "no determinado", "\n",
      pf_txt,
      "Test: ", if (lfc_thr > 0) paste0("|log2FC| > ", lfc_thr, " dentro del modelo")
                else "log2FC != 0", "\n",
      "Encogido log2FC: ", state$deg_rv$shrink %||% "ninguno", "\n",
      "Significativos a FDR <= ", fdr_thr, ": ",
      sum(sig, na.rm = TRUE), " (", n_up, " up / ", n_down, " down)\n",
      "Ultima ejecucion: ", format(state$deg_rv$run_at, "%Y-%m-%d %H:%M:%S")
    )
  })

  # ── Tabla filtrada (reactivo derivado, rapido) ─────────────────────────────
  # El FDR viene del ajuste; |log2FC| y baseMean son filtros de visualizacion.
  deg_filtered <- reactive({
    tab <- state$deg_rv$results
    if (is.null(tab)) return(NULL)
    apply_deg_filters(
      tab,
      fdr        = state$deg_rv$fdr %||% 0.05,
      abs_log2fc = input$deg_log2fc_cutoff %||% 0,
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
    # Las columnas que un motor no rellena (log2FC_shrunk solo lo produce
    # DESeq2, lfcSE no lo produce edgeR) se ocultan en vez de mostrar una
    # columna entera de NA.
    for (nm in intersect(c("log2FC_shrunk", "lfcSE", "stat"), names(df_r))) {
      if (all(is.na(df_r[[nm]]))) df_r[[nm]] <- NULL
    }
    if ("log2FC" %in% names(df_r)) {
      df_r$direction <- ifelse(df_r$log2FC >= 0, "Up", "Down")
      df_r <- df_r[, c("gene", "direction", setdiff(names(df_r), c("gene", "direction"))), drop = FALSE]
    }
    num_cols <- intersect(c("baseMean", "log2FC", "log2FC_shrunk", "lfcSE",
                            "stat", "pvalue", "padj"), names(df_r))
    for (nm in num_cols) df_r[[nm]] <- signif(df_r[[nm]], 4)
    dt_table(df_r, page_length = 15, filter = "top")
  })

  # ── Helpers de ploteo (reutilizados en render y en descarga) ──────────────

  #' Elige el eje de fold-change: el encogido si existe. Los estimadores de
  #' maxima verosimilitud estan sesgados hacia valores exagerados en genes de
  #' baja expresion, asi que un volcano construido sobre ellos destaca
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
    # La significacion es padj <= FDR y nada mas: cuando hay umbral de
    # fold-change ya esta dentro del test, asi que volver a cortar por |log2FC|
    # aqui seria el filtro post-hoc que estamos eliminando.
    sig <- !is.na(df$padj) & df$padj <= fdr_thr
    df$significant <- ifelse(is.na(sig), FALSE, sig)
    df$color <- ifelse(df$significant, "Significativo", "No significativo")
    fdr_y <- -log10(pmax(fdr_thr, .Machine$double.xmin))
    shapes <- list(
      list(type = "line", xref = "paper", x0 = 0, x1 = 1,
           y0 = fdr_y, y1 = fdr_y,
           line = list(dash = "dot", color = "#A8DADC"))
    )
    if (is.finite(lfc_thr) && lfc_thr > 0) {
      shapes <- c(shapes, list(
        list(type = "line", x0 = -lfc_thr, x1 = -lfc_thr, yref = "paper",
             y0 = 0, y1 = 1, line = list(dash = "dot", color = "#F4A6A6")),
        list(type = "line", x0 = lfc_thr, x1 = lfc_thr, yref = "paper",
             y0 = 0, y1 = 1, line = list(dash = "dot", color = "#F4A6A6"))
      ))
    }
    p <- plotly::plot_ly(
      df, x = ~x, y = ~minus_log10_p, color = ~color,
      colors = c("Significativo" = "#7BBF9A", "No significativo" = "#C0C0C0"),
      type = "scatter", mode = "markers",
      text = ~paste0("Gen: ", gene, "<br>", ax$label, ": ", round(x, 3),
                     "<br>padj: ", signif(padj, 3)),
      hoverinfo = "text",
      marker = list(size = 6, opacity = 0.7)
    ) |>
      plotly::layout(
        title = list(text = plot_title(contrast,
                       if (lfc_thr > 0) paste0("umbral del test |log2FC| > ", lfc_thr)),
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
      colors = c("Significativo" = "#7BBF9A", "No significativo" = "#C0C0C0"),
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

  # ── Enriquecimiento ────────────────────────────────────────────────────────
  enrich_rv <- reactiveVal(NULL)
  enrich_mapping_rv <- reactiveVal(NULL)

  # El OrgDb disponible hoy es solo el de E. coli K12; se centraliza aqui para
  # que el selector de keyType y el enriquecimiento no se desincronicen.
  deg_orgdb <- reactive({
    if (isTRUE(HAS_ORGECDB)) "org.EcK12.eg.db" else NULL
  })

  # keyType seleccionable, poblado con los keytypes reales del OrgDb. Antes
  # estaba fijado a "SYMBOL" en el codigo, lo que con locus tags de E. coli
  # (b0001) daba un mapeo casi nulo sin ningun aviso.
  observe({
    kt <- orgdb_keytypes(deg_orgdb())
    if (!length(kt)) kt <- c("SYMBOL", "ENTREZID", "ENSEMBL", "ALIAS", "REFSEQ")
    preferred <- if ("SYMBOL" %in% kt) "SYMBOL" else kt[1]
    updateSelectInput(session, "deg_go_keytype", choices = kt,
                      selected = isolate(input$deg_go_keytype) %||% preferred)
  })

  observeEvent(input$deg_run_enrich_btn, {
    req(state$deg_rv$results)
    ont <- input$deg_ontology %||% "BP"
    approach <- input$deg_enrich_approach %||% "ora"
    org_code <- trimws(input$deg_kegg_organism %||% "eco")
    if (!nzchar(org_code)) org_code <- "eco"
    # Fondo = los genes efectivamente testeados, no el genoma completo. Aplica
    # tanto a GO como a KEGG.
    universe <- state$deg_rv$results$gene

    if (identical(approach, "gsea")) {
      # GSEA parte del ranking completo, sin umbralizar: por eso usa la tabla
      # entera y no la lista filtrada.
      metric <- input$deg_gsea_metric %||% "stat"
      rk <- deg_ranking_metric(state$deg_rv$results, metric)
      if (is.null(rk$ranked)) {
        showNotification(paste0("No se pudo construir el ranking: ",
                                rk$error %||% "—"), type = "error", duration = 10)
        enrich_rv(NULL); return()
      }
      if (!is.na(rk$ties_frac) && rk$ties_frac > 0.01) {
        showNotification(
          paste0("El ", round(100 * rk$ties_frac, 1), " % de los genes comparte ",
                 "valor en el ranking. GSEA no resuelve los empates, asi que su ",
                 "orden interno es arbitrario; considera usar la metrica 'stat'."),
          type = "warning", duration = 16
        )
      }
      res <- run_gsea(rk$ranked, ont = ont, OrgDb = deg_orgdb(),
                      organism = org_code,
                      keyType = if (identical(ont, "KEGG"))
                        input$deg_kegg_keytype %||% "kegg"
                      else input$deg_go_keytype %||% "SYMBOL",
                      exponent = 0)
    } else {
      df <- deg_filtered()
      if (is.null(df) || !nrow(df)) {
        showNotification("No hay genes significativos con los filtros actuales.",
                         type = "warning"); return()
      }
      genes <- df$gene
      if (identical(ont, "KEGG")) {
        res <- run_enrichment_kegg(genes, universe = universe, organism = org_code,
                                   keyType = input$deg_kegg_keytype %||% "kegg")
      } else {
        res <- run_enrichment_go(genes, universe = universe,
                                 OrgDb = deg_orgdb(), ont = ont,
                                 keyType = input$deg_go_keytype %||% "SYMBOL",
                                 simplify_terms = isTRUE(input$deg_go_simplify))
      }
    }

    enrich_mapping_rv(res$mapping)
    # Un mapeo bajo hace el resultado no interpretable, asi que se avisa aunque
    # el enriquecimiento haya devuelto terminos.
    mp <- res$mapping
    if (!is.null(mp) && !is.na(mp$rate %||% NA) && mp$rate < 0.5) {
      showNotification(
        paste0("Atencion: solo se ha mapeado el ", round(100 * mp$rate, 1),
               " % de los genes (", mp$n_mapped, "/", mp$n_input,
               "). Revisa el keyType: los IDs de featureCounts suelen ser locus ",
               "tags, no simbolos."),
        type = "warning", duration = 16
      )
    }

    if (is.null(res$table)) {
      showNotification(paste0("Enriquecimiento sin resultados: ", res$error %||% "—"),
                       type = "warning", duration = 8)
      enrich_rv(NULL); return()
    }
    enrich_rv(res$table)
  })

  output$deg_enrich_mapping <- renderUI({
    mp <- enrich_mapping_rv()
    txt <- mapping_rate_text(mp)
    if (is.null(txt)) return(NULL)
    low <- !is.null(mp) && !is.na(mp$rate %||% NA) && mp$rate < 0.5
    div(class = paste("small py-1 px-2 mb-1 rounded",
                      if (low) "alert alert-warning mb-2" else "text-muted"),
        if (low) tagList(icon("triangle-exclamation"), " ") else NULL,
        tags$b("Tasa de mapeo: "), txt)
  })

  #' El dotplot sirve a dos tipos de tabla: el ORA trae GeneRatio/Count y se
  #' ordena por significacion; GSEA trae NES/setSize, y ahi lo informativo es el
  #' signo del NES (si el conjunto sube o baja), no solo el p-valor.
  make_enrich_dotplot <- function(df) {
    if (is.null(df) || !nrow(df))
      return(plotly_message("Pulsa 'Calcular' para correr el enriquecimiento."))
    top <- enrichment_dotplot_data(df, top_n = 15)
    if (is.null(top) || !nrow(top))
      return(plotly_message("Sin terminos enriquecidos."))

    if ("NES" %in% names(top)) {
      top$size_val <- top$setSize %||% rep(10, nrow(top))
      return(plotly::plot_ly(
        top,
        x = ~NES,
        y = ~stats::reorder(Description, NES),
        type = "scatter", mode = "markers",
        size = ~size_val,
        marker = list(color = ~ifelse(NES > 0, "#7BBF9A", "#F4A6A6"),
                      sizemode = "area"),
        text = ~paste0("Conjunto: ", Description,
                       "<br>NES: ", round(NES, 3),
                       "<br>genes en el conjunto: ", size_val,
                       "<br>p.adjust: ", signif(p.adjust, 3)),
        hoverinfo = "text"
      ) |>
        plotly::layout(
          xaxis = list(title = "NES (score de enriquecimiento normalizado)"),
          yaxis = list(title = "", automargin = TRUE),
          margin = list(l = 220),
          shapes = list(list(type = "line", x0 = 0, x1 = 0, yref = "paper",
                             y0 = 0, y1 = 1,
                             line = list(dash = "dot", color = "#60756A")))
        ))
    }

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
  }

  output$deg_enrich_dotplot <- plotly::renderPlotly({
    make_enrich_dotplot(enrich_rv())
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
      make_deg_ma_plot(df, fdr_thr  = state$deg_rv$fdr %||% 0.05,
                       contrast = state$deg_rv$contrast)
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
