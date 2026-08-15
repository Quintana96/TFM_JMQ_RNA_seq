#' server_tab_deg.R
#' Nucleo del modulo de la pestana 4 (Expresion diferencial).
#'
#' Responsabilidades de ESTE fichero:
#'   - carga de la matriz de conteos desde las tres fuentes;
#'   - editor inline del samplesheet y los selectores que dependen de el;
#'   - validacion en vivo de la formula de diseno;
#'   - el observer que ajusta el modelo, y el banner y el estado que lo resumen;
#'   - `deg_filtered`, la seleccion de significativos con los filtros visuales.
#'
#' El resto vive en cuatro ficheros hermanos, porque este llego a 1.608 lineas en
#' una unica funcion:
#'   - server_tab_deg_results.R  tabla, volcano, MA, PCA, heatmaps
#'   - server_tab_deg_diag.R     diagnosticos post-ajuste y sesgo de longitud
#'   - server_tab_deg_enrich.R   enriquecimiento funcional y GSEA
#'   - server_tab_deg_reports.R  replicabilidad, comparacion, reproducibilidad
#'
#' No se usa moduleServer() a proposito, igual que en el resto de la aplicacion:
#' obligaria a renombrar todos los IDs Shiny. La comunicacion entre partes es
#' explicita a traves de `ctx` (ver el final del fichero).

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
      if (!is.null(cm) && length(cm)) {
        m <- as.matrix(cm)
        attr(m, "counts_origin") <- list(
          tipo = "Ejecucion actual de la sesion",
          ruta = state$run_params_rv()$output_dir %||% "—",
          detalle = paste0("Fuente: ", state$data_rv$source %||% "—"))
        return(m)
      }
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
      cm <- tryCatch(load_counts_from_workflow(
                       sel, p$tool %||% "",
                       annotation_file = annotation_file_for_run(sel)),
                     error = function(e) NULL)
      if (!is.null(cm)) attr(cm, "counts_origin") <- list(
        tipo = "Ejecucion guardada", ruta = sel,
        detalle = paste0("Cuantificador: ", p$tool %||% "—"))
      return(cm)
    }
    if (identical(src, "upload")) {
      up <- input$deg_counts_upload
      if (is.null(up) || !nrow(up)) return(NULL)
      cm <- read_uploaded_counts(up$datapath)
      # El md5 del fichero subido es lo unico que permite despues demostrar que
      # dos analisis partieron de la misma matriz.
      if (!is.null(cm)) attr(cm, "counts_origin") <- list(
        tipo = "Matriz subida", ruta = up$name,
        md5 = unname(tryCatch(tools::md5sum(up$datapath), error = function(e) NA_character_)),
        detalle = paste0(fmt_bytes(up$size %||% NA_real_)))
      return(cm)
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
    # La procedencia y el metodo de resumen a gen se leen AQUI: los atributos se
    # pierden al alinear y prefiltrar, y el informe los necesita para declarar
    # de donde salio la matriz y si se uso tximport o el respaldo.
    counts_origin_info <- attr(cm, "counts_origin")
    counts_source_info <- attr(cm, "counts_source")

    if (is.null(cm) || !ncol(cm)) {
      # Un "no hay matriz" a secas no dice donde mirar. El fallo depende de la
      # fuente elegida en la tarjeta 1, asi que el mensaje la nombra y explica
      # que falta exactamente en cada caso.
      src <- input$deg_source %||% "current"
      detalle <- switch(src,
        "current" = paste0(
          "Fuente: 'Ejecucion actual', pero no hay ninguna cargada en esta ",
          "sesion. Si has subido la matriz en la pestana 1, pulsa alli ",
          "'Analisis a partir de matriz de conteos' para cargarla; si vienes de ",
          "una ejecucion guardada, elige 'Ejecucion guardada' aqui."),
        "saved" = {
          sel <- input$selected_deg_run_dir %||% ""
          if (!nzchar(sel)) {
            "Fuente: 'Ejecucion guardada', pero no hay ninguna seleccionada."
          } else if (!dir.exists(sel)) {
            paste0("Fuente: 'Ejecucion guardada', pero el directorio ya no ",
                   "existe: ", sel)
          } else {
            paste0("Fuente: 'Ejecucion guardada' (", basename(sel),
                   "), pero no se ha podido leer su matriz de conteos. ",
                   "Comprueba que existe 04_counts/count_matrix.tsv, o las ",
                   "cuantificaciones en 03_alignments/.")
          }
        },
        "upload" = paste0(
          "Fuente: 'Matriz subida', pero no se ha podido leer el fichero. ",
          "Comprueba que la subida ha terminado y que es un TSV o CSV con los ",
          "genes en la primera columna y una muestra por columna."),
        "No se ha podido obtener la matriz de conteos."
      )
      showNotification(detalle, type = "error", duration = 20)
      return()
    }
    if (is.null(df) || !nrow(df)) {
      showNotification("Carga o crea un samplesheet.", type = "error"); return()
    }

    # Seudonimizacion opcional de los identificadores de muestra. Se aplica ANTES
    # de alinear y ajustar, para que los alias viajen a todo lo que se genere
    # despues: graficos, informe, script y ficheros persistidos.
    pseudo_map <- NULL
    if (isTRUE(input$deg_pseudonymize)) {
      ps <- pseudonymize_dataset(cm, df)
      cm <- ps$counts; df <- ps$meta; pseudo_map <- ps$map
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

    # Semillas usadas en este ajuste. Se registran para poder declararlas en el
    # informe y reproducir el resultado: sin esto, un analisis con sva o con
    # Swish no es reproducible aunque el resto de parametros coincida.
    seeds_used <- list()

    # Variables sustitutas: se anaden al DISENO (no se corrigen los conteos),
    # que es la forma correcta de tratarlas para testear.
    if (isTRUE(input$deg_use_sva)) {
      base_f <- stats::as.formula(dsg_formula %||%
        if (!is.null(batch) && nzchar(batch)) paste0("~ ", batch, " + condition")
        else "~ condition")
      n_req <- input$deg_n_sv %||% 0
      svres <- estimate_surrogate_vars(
        cm_f, meta_aln, base_f,
        n_sv = if (is.finite(n_req) && n_req >= 1) n_req else NULL,
        seed = ANALYSIS_SEED)
      seeds_used$sva <- ANALYSIS_SEED
      if (!is.null(svres$error)) {
        showNotification(paste0("No se han podido estimar variables sustitutas: ",
                                svres$error),
                         type = "error", duration = 14)
        return()
      }
      # `num.sv` puede estimar 0: no hay estructura latente que modelar. Antes se
      # forzaba a 1 y se metia una covariable espuria que gastaba un grado de
      # libertad y podia absorber senal de la condicion. Ahora se sigue sin SV.
      if (is.null(svres$sv) || !svres$n_sv) {
        showNotification(paste0("sva no ha encontrado variacion latente que ",
                                "modelar (0 variables sustitutas): el analisis ",
                                "continua con el diseno sin SV."),
                         type = "message", duration = 12)
      } else {
        meta_aln <- cbind(meta_aln, as.data.frame(svres$sv))
        dsg_formula <- paste(deparse1(base_f), "+",
                             paste(colnames(svres$sv), collapse = " + "))
        seeds_used$n_sv <- svres$n_sv
        msg <- paste0(svres$n_sv, " variable(s) sustituta(s) anadidas al diseno")
        if (!is.na(svres$n_sv_estimated %||% NA) && svres$n_sv_estimated > svres$n_sv) {
          msg <- paste0(msg, " (sva propuso ", svres$n_sv_estimated,
                        ", recortadas para conservar grados de libertad)")
        }
        showNotification(paste0(msg, "."), type = "message", duration = 10)
      }
    }

    # Barra de progreso en lugar de un aviso que desaparece a los 4 segundos: el
    # ajuste tarda decenas de segundos (mas si hay encogido o IHW) y sin feedback
    # persistente la aplicacion parece colgada.
    res <- withProgress(message = paste0("Ajustando el modelo (", method, ")"),
                        value = 0.15, {
      setProgress(value = 0.25, detail = "estimando dispersiones y ajustando")
      out <- if (identical(method, "Swish")) {
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
                        fdr = fdr_target, seed = ANALYSIS_SEED),
          error = function(e) list(table = NULL, error = conditionMessage(e),
                                   method = method))
      }
    } else {
      tryCatch(
        run_deg(cm_f, meta_aln, method = method, ref_level = ref, batch = batch,
                fdr = fdr_target, lfc_threshold = lfc_thr, shrink = do_shrink,
                contrast_num = num, use_ihw = isTRUE(input$deg_use_ihw),
                design_formula = dsg_formula, test_coef = test_coef,
                outliers = input$deg_outliers %||% "na"),
        error = function(e) list(table = NULL, error = conditionMessage(e), method = method)
      )
      }
      setProgress(value = 0.85, detail = "preparando visualizaciones")
      out
    })
    quant_tool_used <- NULL
    if (identical(method, "Swish")) {
      res$method <- method
      # Swish permuta etiquetas: sin semilla, el mismo input da q-valores
      # distintos en cada ejecucion.
      seeds_used$swish <- ANALYSIS_SEED
      seeds_used$swish_nperms <- res$n_perms %||% SWISH_NPERMS
      quant_tool_used <- run_info$tool
    }

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
    # Contraste y diseno con los que se ha ajustado REALMENTE. Los selectores de
    # la interfaz pueden cambiar despues sin relanzar el analisis, asi que el
    # informe, el script exportado, el bootstrap y la comparacion de metodos
    # leen de aqui y no de `input$...`.
    state$deg_rv$ref_level      <- ref
    state$deg_rv$contrast_num   <- num
    state$deg_rv$batch          <- batch
    state$deg_rv$design_formula <- dsg_formula
    state$deg_rv$test_coef      <- test_coef
    state$deg_rv$quant_tool     <- quant_tool_used
    state$deg_rv$seeds          <- seeds_used
    state$deg_rv$pseudonym_map  <- pseudo_map
    state$deg_rv$counts_origin  <- counts_origin_info
    state$deg_rv$counts_source  <- counts_source_info
    # Directorio de la ejecucion de origen, si lo hay: da acceso a versions.tsv
    # y checksums.tsv del pipeline para cerrar el ciclo de trazabilidad.
    state$deg_rv$run_dir <- if (identical(input$deg_source %||% "current", "saved")) {
      input$selected_deg_run_dir %||% NULL
    } else state$run_params_rv()$output_dir %||% NULL

    if (do_shrink && identical(method, "DESeq2") &&
        identical(state$deg_rv$shrink, "ninguno")) {
      showNotification(
        "No se pudo encoger el log2FC (apeglm/ashr no disponibles).",
        type = "warning", duration = 8
      )
    }
    if (!is.null(res$design_warning)) {
      showNotification(res$design_warning, type = "warning", duration = 16)
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

    # Persistencia del analisis. Hasta ahora el informe y el script solo
    # existian si el usuario los descargaba, de modo que una ejecucion del
    # pipeline dejaba rastro en disco pero los analisis hechos sobre ella no:
    # no habia forma de saber cuantos se lanzaron ni con que parametros.
    dest <- persist_deg_analysis(state$deg_rv, base_dir = state$deg_rv$run_dir,
                                 outputs_dir = state$outputs_dir)
    append_audit_log("deg_run", list(
      motor = method, contraste = res$contrast, fdr = fdr_target,
      lfc_umbral = lfc_thr, diseno = state$deg_rv$design,
      genes = nrow(res$table),
      significativos = sum(!is.na(res$table$padj) & res$table$padj <= fdr_target),
      origen = counts_origin_info$tipo %||% "—",
      guardado_en = dest %||% "no guardado"
    ), outputs_dir = state$outputs_dir)

    showNotification(
      paste0("DEG completado (", method, "): ", nrow(res$table), " filas.",
             if (!is.null(dest)) paste0(" Guardado en ", basename(dirname(dest)),
                                        "/", basename(dest), ".") else ""),
      type = "default", duration = 8)
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

  # Lista de significativos SIN los filtros de visualizacion: unicamente el FDR
  # con el que se ajusto el modelo.
  #
  # Es la que debe alimentar el enriquecimiento. Usar la tabla filtrada hacia que
  # mover un deslizador declarado "solo visual" cambiase el resultado del ORA, y
  # ademas recortaba la lista sin recortar el universo, que es exactamente el
  # sesgo de fondo mal definido que la Fase 1 elimino de las listas DEG.
  deg_significant <- reactive({
    tab <- state$deg_rv$results
    if (is.null(tab)) return(NULL)
    deg_significant_genes(tab, fdr = state$deg_rv$fdr %||% 0.05)
  })

  # Universo del enriquecimiento: los genes efectivamente EVALUABLES, es decir
  # los que tienen padj. Los descartados por el filtrado independiente nunca
  # podrian haber entrado en la lista de significativos, asi que incluirlos en el
  # fondo infla artificialmente el enriquecimiento.
  deg_universe <- reactive({
    tab <- state$deg_rv$results
    if (is.null(tab)) return(NULL)
    deg_testable_universe(tab)
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

  # ── Contexto compartido y delegacion en las partes del modulo ──────────────
  #
  # La pestana 4 llego a 1.608 lineas en una sola funcion, asi que se ha partido
  # en cuatro ficheros por area de responsabilidad. `ctx` es un environment (no
  # una lista) por el mismo motivo que `state`: en R los environments son
  # pass-by-reference, de modo que lo que una parte publica en `ctx` lo ven las
  # demas. Las listas serian copy-on-modify y la publicacion se perderia.
  #
  # Lo que va en `ctx` es exactamente lo que MAS de una parte necesita; todo lo
  # demas queda local a su fichero.
  ctx <- new.env(parent = emptyenv())
  ctx$meta_rv           <- meta_rv            # lo lee la sugerencia de metodo
  ctx$deg_counts_source <- deg_counts_source
  ctx$deg_filtered      <- deg_filtered       # tabla y graficos (filtros visuales)
  ctx$deg_significant   <- deg_significant    # enriquecimiento: solo el FDR del ajuste
  ctx$deg_universe      <- deg_universe       # fondo del enriquecimiento

  # El orden importa: los diagnosticos publican en `ctx` los reactivos que el
  # informe reproducible consume, asi que tienen que registrarse antes.
  server_tab_deg_results(input, output, session, state, ctx)
  server_tab_deg_diag(input, output, session, state, ctx)
  server_tab_deg_enrich(input, output, session, state, ctx)
  server_tab_deg_reports(input, output, session, state, ctx)

  invisible(NULL)
}
