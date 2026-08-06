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
        return(tryCatch(load_counts_from_workflow(
                          p$output_dir, p$tool,
                          annotation_file = p$annotation_file %||%
                            annotation_file_for_run(p$output_dir)),
                        error = function(e) NULL))
      }
      return(NULL)
    }
    if (identical(src, "saved")) {
      sel <- input$selected_deg_run_dir %||% ""
      if (!nzchar(sel) || !dir.exists(sel)) return(NULL)
      p <- infer_result_params(sel, state$workflow_path)
      return(tryCatch(load_counts_from_workflow(
                        sel, p$tool %||% "",
                        annotation_file = annotation_file_for_run(sel)),
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

  # ── Diseno avanzado: validacion en vivo ────────────────────────────────────
  # Se valida mientras se escribe, para que el error no llegue como un mensaje
  # criptico de DESeq2 despues de esperar el ajuste.
  design_validation <- reactive({
    if (!isTRUE(input$deg_advanced_design)) return(NULL)
    df <- meta_rv()
    if (is.null(df) || !nrow(df)) return(NULL)
    validate_design_formula(input$deg_design_formula %||% "~ condition", df)
  })

  output$deg_design_feedback <- renderUI({
    v <- design_validation()
    if (is.null(v)) return(NULL)
    if (!isTRUE(v$ok)) {
      return(div(class = "alert alert-danger py-2 px-2 small mb-2",
                 icon("circle-xmark"), " ", paste(v$errors, collapse = " ")))
    }
    tagList(
      div(class = "alert alert-success py-2 px-2 small mb-1",
          icon("circle-check"), " ", design_summary_text(v)),
      if (length(v$warnings)) div(class = "alert alert-warning py-2 px-2 small mb-1",
                                  icon("triangle-exclamation"), " ",
                                  paste(v$warnings, collapse = " ")) else NULL
    )
  })

  # El selector de coeficiente se rellena con los nombres REALES del ultimo
  # ajuste: DESeq2 y model.matrix nombran las interacciones distinto
  # ("a.b" vs "a:b"), asi que adivinarlos daria opciones invalidas.
  observe({
    avail <- state$deg_rv$coef_available
    if (is.null(avail) || !length(avail)) {
      updateSelectInput(session, "deg_test_coef",
                        choices = c("Automatico (contraste de la condicion)" = ""))
      return()
    }
    ch <- setdiff(avail, c("Intercept", "(Intercept)"))
    updateSelectInput(
      session, "deg_test_coef",
      choices = c(c("Automatico (contraste de la condicion)" = ""), ch),
      selected = isolate(input$deg_test_coef) %||% "")
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

    # Formula del diseno: la libre si esta activada, y si no la construye el motor
    # a partir de ref_level/batch.
    dsg_formula <- if (isTRUE(input$deg_advanced_design)) {
      input$deg_design_formula %||% "~ condition"
    } else NULL
    test_coef <- if (isTRUE(input$deg_advanced_design)) {
      tc <- input$deg_test_coef %||% ""
      if (nzchar(tc)) tc else NULL
    } else NULL

    # Variables sustitutas: se anaden al DISENO (no se corrigen los conteos),
    # que es la forma correcta de tratarlas para testear.
    if (isTRUE(input$deg_use_sva)) {
      base_f <- stats::as.formula(dsg_formula %||%
        if (!is.null(batch) && nzchar(batch)) paste0("~ ", batch, " + condition")
        else "~ condition")
      n_req <- input$deg_n_sv %||% 0
      svres <- estimate_surrogate_vars(
        cm_f, meta_aln, base_f,
        n_sv = if (is.finite(n_req) && n_req >= 1) n_req else NULL)
      if (is.null(svres$sv)) {
        showNotification(paste0("No se han podido estimar variables sustitutas: ",
                                svres$error %||% "—"),
                         type = "error", duration = 14)
        return()
      }
      meta_aln <- cbind(meta_aln, as.data.frame(svres$sv))
      dsg_formula <- paste(deparse1(base_f), "+",
                           paste(colnames(svres$sv), collapse = " + "))
      msg <- paste0(svres$n_sv, " variable(s) sustituta(s) anadidas al diseno")
      if (!is.na(svres$n_sv_estimated %||% NA) && svres$n_sv_estimated > svres$n_sv) {
        msg <- paste0(msg, " (sva propuso ", svres$n_sv_estimated,
                      ", recortadas para conservar grados de libertad)")
      }
      showNotification(paste0(msg, "."), type = "message", duration = 10)
    }

    showNotification(paste0("Corriendo DEG (", method, "). Esto puede tardar..."),
                     type = "message", duration = 4)

    res <- if (identical(method, "Swish")) {
      # Swish es el unico motor que no parte de la matriz de conteos: la
      # incertidumbre de la cuantificacion vive en las replicas inferenciales de
      # los ficheros de salmon/kallisto, no en la matriz ya resumida. Por eso
      # necesita el directorio de la ejecucion y no funciona con una matriz subida.
      run_info <- if (identical(input$deg_source %||% "current", "saved")) {
        list(dir = input$selected_deg_run_dir %||% "",
             tool = (infer_result_params(input$selected_deg_run_dir %||% "",
                                         state$workflow_path))$tool %||% "")
      } else {
        p <- state$run_params_rv()
        list(dir = p$output_dir %||% "", tool = p$tool %||% "")
      }
      if (!nzchar(run_info$dir) || !dir.exists(run_info$dir) ||
          !run_info$tool %in% c("salmon", "kallisto")) {
        list(table = NULL, method = method, error = paste0(
          "Swish necesita una ejecucion de salmon o kallisto con replicas ",
          "inferenciales. Selecciona la ejecucion actual o una guardada como ",
          "fuente de datos (no una matriz subida)."))
      } else {
        tryCatch(
          run_deg_swish(meta_aln, run_info$dir, run_info$tool,
                        annotation_file = annotation_file_for_run(run_info$dir),
                        ref_level = ref, contrast_num = num, batch = batch,
                        fdr = fdr_target),
          error = function(e) list(table = NULL, error = conditionMessage(e),
                                   method = method))
      }
    } else {
      tryCatch(
        run_deg(cm_f, meta_aln, method = method, ref_level = ref, batch = batch,
                fdr = fdr_target, lfc_threshold = lfc_thr, shrink = do_shrink,
                contrast_num = num, use_ihw = isTRUE(input$deg_use_ihw),
                design_formula = dsg_formula, test_coef = test_coef),
        error = function(e) list(table = NULL, error = conditionMessage(e), method = method)
      )
    }
    if (identical(method, "Swish")) res$method <- method

    if (is.null(res$table)) {
      showNotification(paste0("Error en ", method, ": ", res$error %||% "fallo desconocido"),
                       type = "error", duration = 12)
      return()
    }

    # Cache de transformacion para visualizacion. Aqui SI se corrige la matriz si
    # el usuario lo ha pedido, porque afecta solo a los graficos: el test ya se
    # ha hecho sobre los conteos sin tocar.
    viz_corr <- input$deg_viz_correction %||% "none"
    counts_for_viz <- cm_f
    viz_note <- NULL
    if (identical(viz_corr, "combat") && !is.null(batch) && nzchar(batch %||% "") &&
        batch %in% names(meta_aln)) {
      cb <- combat_seq_counts(cm_f, meta_aln[[batch]], meta_aln$condition)
      if (!is.null(cb$counts)) {
        counts_for_viz <- cb$counts
        viz_note <- paste0("Graficos sobre conteos ajustados con ComBat-seq por '", batch, "'.")
      } else {
        showNotification(paste0("ComBat-seq fallo: ", cb$error %||% "—",
                                ". Los graficos usan los conteos sin corregir."),
                         type = "warning", duration = 12)
      }
    } else if (identical(viz_corr, "combat")) {
      showNotification(paste0("ComBat-seq necesita una columna de batch: activa ",
                              "'Incluir variable batch' en la tarjeta 3."),
                       type = "warning", duration = 12)
    }

    vst_mat <- tryCatch(vst_or_rlog(counts_for_viz, meta_aln), error = function(e) NULL)

    if (identical(viz_corr, "rbe") && !is.null(vst_mat)) {
      if (!is.null(batch) && nzchar(batch %||% "") && batch %in% names(meta_aln)) {
        dm <- tryCatch(stats::model.matrix(~ condition, data = meta_aln),
                       error = function(e) NULL)
        rb <- remove_batch_for_plots(vst_mat, meta_aln[[batch]], design = dm)
        if (is.null(rb$error)) {
          vst_mat <- rb$mat
          viz_note <- paste0("Graficos con el efecto de '", batch,
                             "' eliminado (removeBatchEffect).")
        } else {
          showNotification(paste0("removeBatchEffect fallo: ", rb$error),
                           type = "warning", duration = 12)
        }
      } else {
        showNotification(paste0("removeBatchEffect necesita una columna de batch: ",
                                "activa 'Incluir variable batch' en la tarjeta 3."),
                         type = "warning", duration = 12)
      }
    }

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
    state$deg_rv$design        <- res$design %||% "~ condition"
    state$deg_rv$coef          <- res$coef
    state$deg_rv$coef_available <- res$coef_available %||% character(0)
    state$deg_rv$viz_note      <- viz_note

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
                       "  ·  diseno ", state$deg_rv$design %||% "~ condition",
                       "  ·  coef ", state$deg_rv$coef %||% "—",
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
    # Si los graficos llevan una correccion que el test no lleva, hay que decirlo
    # justo al lado del contraste, no en una ayuda escondida.
    vn <- state$deg_rv$viz_note
    viz <- if (!is.null(vn)) tags$div(
      class = "small mt-1", icon("eye"), " ", tags$b("Solo en los graficos: "), vn,
      tags$span(class = "text-muted",
                " El test se ha hecho sobre los conteos sin corregir.")) else NULL
    div(class = "alert alert-light border py-2 px-3 mb-2",
        icon("circle-info"), " ", bits, extra, viz)
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

  # ── Sugerencia de metodo segun el tamano muestral (item 21) ────────────────
  output$deg_method_hint <- renderUI({
    df <- meta_rv()
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
    list(pi0 = diag_pi0()$pi0 %||% NA_real_,
         verdict = {
           p <- diag_pvalues()
           if (is.null(p) || !length(p)) NULL
           else diagnose_pvalue_shape(p, diag_pi0()$pi0 %||% NA_real_)
         },
         na_breakdown = deg_na_breakdown(),
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
