#' server_tab_processing.R
#' Logica server de la Tab 2 (Procesamiento):
#'   - renderUI("tab2_content") segun process_unlocked
#'   - Lanzamiento del workflow (processx no bloqueante o system2 bloqueante)
#'   - Polling periodico, parsing de progreso, heartbeat
#'   - Boton stop, refresh
#'
#' Comparte estado con tab_config y tab_results via `state`.

server_tab_processing <- function(input, output, session, state) {
  proc_rv <- state$proc_rv
  shared <- state$shared
  outputs_dir <- state$outputs_dir
  workflow_path <- state$workflow_path

  # Timer de polling
  poll_timer <- reactiveTimer(500)

  progress_tick <- function() {
    if (isTRUE(proc_rv$running)) poll_timer()
  }

  #' Anade un fragmento al log en disco y al reactivo
  append_run_log <- function(chunk) {
    chunk <- terminal_text(chunk)
    if (!nzchar(chunk)) return(invisible(FALSE))
    if (!is.null(proc_rv$log_file) && nzchar(proc_rv$log_file)) {
      cat(chunk, file = proc_rv$log_file, append = TRUE)
      proc_rv$log_seen_size <- file.info(proc_rv$log_file)$size %||% 0
    }
    state$log_text(trim_log_text(paste0(state$log_text(), chunk)))
    proc_rv$last_output_time <- Sys.time()
    invisible(TRUE)
  }

  #' Linea de log con timestamp y newline
  log_line <- function(msg) append_run_log(paste0(ts_log(msg), "\n"))

  #' Sincroniza el log de disco al reactivo (tail-only)
  sync_run_log_file <- function() {
    log_file <- proc_rv$log_file
    if (is.null(log_file) || !nzchar(log_file) || !file.exists(log_file)) return("")
    size <- file.info(log_file)$size %||% 0
    old_size <- proc_rv$log_seen_size %||% 0
    if (is.na(size) || size <= old_size) return("")
    if (old_size < 0 || old_size > size) old_size <- 0

    con <- file(log_file, open = "rb")
    on.exit(close(con), add = TRUE)
    seek(con, where = old_size, origin = "start")
    raw <- readBin(con, what = "raw", n = size - old_size)
    proc_rv$log_seen_size <- size

    chunk <- terminal_text(rawToChar(raw))
    if (nzchar(chunk)) {
      state$log_text(trim_log_text(paste0(state$log_text(), chunk)))
      proc_rv$last_output_time <- Sys.time()
    }
    chunk
  }

  #' Actualiza checkpoints / muestras a partir de las nuevas lineas
  update_progress_from_log <- function(chunk) {
    lines <- strsplit(terminal_text(chunk), "\n", fixed = TRUE)[[1]]
    lines <- lines[nzchar(lines)]
    for (line in lines) {
      if (grepl("fastqc", line, ignore.case = TRUE) && proc_rv$cp_idx < 1L)
        proc_rv$cp_idx <- 1L
      if (grepl("Processing sample:", line, ignore.case = TRUE) ||
          grepl("bowtie2|salmon quant|kallisto quant|fastp", line, ignore.case = TRUE)) {
        proc_rv$cp_idx <- max(proc_rv$cp_idx, 2L)
      }
      if (grepl("samtools sort|samtools index|featureCounts|abundance\\.tsv|abundance\\.h5|quant(?:\\.sf)?", line, ignore.case = TRUE))
        proc_rv$cp_idx <- max(proc_rv$cp_idx, 3L)
      if (grepl("Generating count matrix", line, ignore.case = TRUE)) {
        if (!is.null(proc_rv$cur_sample)) {
          ss <- proc_rv$samp_stat
          ss[[proc_rv$cur_sample]] <- "done"
          proc_rv$samp_stat <- ss
          proc_rv$cur_sample <- NULL
        }
        proc_rv$cp_idx <- max(proc_rv$cp_idx, 4L)
      }
      if (grepl("multiqc", line, ignore.case = TRUE))
        proc_rv$cp_idx <- max(proc_rv$cp_idx, 6L)
      if (grepl("Analysis completed|Analysis completed successfully", line, ignore.case = TRUE))
        proc_rv$cp_idx <- length(proc_rv$checkpoints)
      if (grepl("Processing sample:", line, ignore.case = TRUE)) {
        sname <- trimws(sub(".*Processing sample:\\s*", "", line, ignore.case = TRUE))
        cur   <- proc_rv$cur_sample
        if (!is.null(cur)) {
          ss <- proc_rv$samp_stat; ss[[cur]] <- "done"; proc_rv$samp_stat <- ss
          bytes <- proc_rv$sample_sizes[[cur]] %||% 0
          proc_rv$bytes_done <- (proc_rv$bytes_done %||% 0) + bytes
        }
        proc_rv$cur_sample <- sname
        ss <- proc_rv$samp_stat; ss[[sname]] <- "running"; proc_rv$samp_stat <- ss
      }
    }
  }

  #' Post-procesado comun tras una run (exit_code, output_dir).
  #' Carga matriz, actualiza tabla de archivos y navega a Tab 3 si OK.
  finalize_run <- function(exit_code, output_dir) {
    if (exit_code == 0) {
      p <- state$run_params_rv()
      counts <- tryCatch(
        load_counts_from_workflow(
          output_dir, p$tool,
          annotation_file = p$annotation_file %||% annotation_file_for_run(output_dir)),
        error = function(e) NULL)
      if (!is.null(counts)) {
        state$data_rv$count_matrix <- counts
        state$data_rv$counts_ready <- TRUE
        log_line(sprintf("Matriz cargada: %d genes x %d muestras",
                         nrow(counts), ncol(counts)))
        # Si se ha tenido que degradar a est_counts crudos, decirlo: antes
        # fallaba en silencio y parecia que se habia usado tximport.
        src <- attr(counts, "counts_source")
        if (!is.null(src)) {
          log_line(paste0("Origen de los conteos: ", src$method,
                          if (!isTRUE(src$ok)) paste0(" — ", src$detail) else ""))
        }
      }
      files <- list.files(output_dir, recursive = TRUE, full.names = FALSE)
      if (length(files) > 0) {
        state$output_files_rv(file_table_for_files(output_dir, files))
      }
      log_line("=== Analisis finalizado OK ===")
      state$analysis_done(TRUE)
      state$results_refresh(Sys.time())
      updateSelectInput(session, "selected_result_dir",
                        choices = result_choices(outputs_dir),
                        selected = output_dir)
      proc_rv$cp_idx <- length(proc_rv$checkpoints)
      showNotification("Workflow finalizado correctamente.", type = "message")
      updateNavbarPage(session, "main_nav", selected = "tab_results")
    } else {
      log_line("=== ERROR en el workflow ===")
      showNotification(
        paste0("Error (codigo ", exit_code, "). Revisa el log de ejecucion."),
        type = "error", duration = 12
      )
    }
  }

  # ── Observer de polling (processx) ─────────────────────────────────────────
  observe({
    req(HAS_PROCESSX, proc_rv$running, !is.null(proc_rv$proc))
    poll_timer()

    proc <- proc_rv$proc
    if (!proc$is_alive()) {
      remaining <- sync_run_log_file()
      if (nzchar(remaining)) update_progress_from_log(remaining)

      exit_code <- proc$get_exit_status() %||% 0L
      proc_rv$running <- FALSE
      proc_rv$end_time <- Sys.time()
      shinyjs::disable("stop_btn")
      shinyjs::enable("run_btn")
      log_line(paste0("Codigo de salida: ", exit_code))

      if (!is.null(proc_rv$cur_sample)) {
        ss <- proc_rv$samp_stat
        ss[[proc_rv$cur_sample]] <- "done"
        proc_rv$samp_stat <- ss
        bytes <- proc_rv$sample_sizes[[proc_rv$cur_sample]] %||% 0
        proc_rv$bytes_done <- (proc_rv$bytes_done %||% 0) + bytes
        proc_rv$cur_sample <- NULL
      }
      # Solo se dan por completados TODOS los pasos si el proceso termino bien.
      # Con un codigo de salida distinto de 0, el paso que estaba en curso queda
      # marcado como fallido y los posteriores siguen pendientes: marcar la lista
      # entera en verde contradecia la notificacion de error que se muestra a
      # continuacion.
      if (identical(exit_code, 0L) || identical(exit_code, 0)) {
        proc_rv$cp_idx <- length(proc_rv$checkpoints)
        proc_rv$cp_failed <- NA_integer_
        proc_rv$cp_failed_kind <- NA_character_
      } else {
        proc_rv$cp_failed <- min(proc_rv$cp_idx + 1L, length(proc_rv$checkpoints))
        proc_rv$cp_failed_kind <- "error"
      }

      finalize_run(exit_code, state$run_params_rv()$output_dir %||% "")
    } else {
      new_chunk <- sync_run_log_file()
      if (nzchar(new_chunk)) {
        update_progress_from_log(new_chunk)
      } else if (!is.null(proc_rv$last_output_time)) {
        now <- Sys.time()
        quiet_for <- as.numeric(difftime(now, proc_rv$last_output_time, units = "secs"))
        heartbeat_for <- if (is.null(proc_rv$last_heartbeat)) Inf else
          as.numeric(difftime(now, proc_rv$last_heartbeat, units = "secs"))
        if (quiet_for >= 30 && heartbeat_for >= 30) {
          log_line(sprintf("Proceso activo sin nueva salida desde hace %s.", fmt_elapsed(quiet_for)))
          proc_rv$last_heartbeat <- now
        }
      }
    }
  })

  # ── Contenido renderUI("tab2_content") ────────────────────────────────────
  output$tab2_content <- renderUI({
    if (!state$process_unlocked())
      return(ui_tab_processing_locked())
    cfg <- state$config_snap()
    if (length(cfg) == 0) return(NULL)

    total_sz <- {
      cfg_read_type <- if (identical(cfg$read_type, "Single-end")) "se" else "pe"
      cfg_samples <- detect_samples(cfg$input_dir, cfg_read_type)
      sizes <- sample_fastq_sizes(cfg$input_dir, cfg_samples, cfg_read_type)
      if (length(sizes)) fmt_bytes(sum(sizes, na.rm = TRUE)) else "—"
    }

    ui_tab_processing_content(cfg, total_sz)
  })

  # ── Outputs de progreso ───────────────────────────────────────────────────
  output$elapsed_text <- renderText({
    progress_tick()
    if (is.null(proc_rv$start_time)) return("—")
    end_time <- proc_rv$end_time %||% Sys.time()
    fmt_elapsed(as.numeric(difftime(end_time, proc_rv$start_time, units = "secs")))
  })

  output$eta_text <- renderText({
    progress_tick()
    if (is.null(proc_rv$start_time) || !proc_rv$running) return("—")
    elapsed <- as.numeric(difftime(Sys.time(), proc_rv$start_time, units = "secs"))
    total_bytes <- proc_rv$total_bytes %||% 0
    bytes_done <- proc_rv$bytes_done %||% 0
    if (total_bytes > 0 && bytes_done > 0) {
      speed <- bytes_done / elapsed
      if (speed > 0) {
        remaining <- (total_bytes - bytes_done) / speed
        if (remaining < 0) return("finalizando...")
        return(paste0("~", fmt_elapsed(remaining)))
      }
    }
    n_done  <- sum(vapply(proc_rv$samp_stat, identical, logical(1), "done"))
    n_total <- proc_rv$n_total
    if (n_total == 0) return("—")
    if (n_done == 0) {
      cur <- proc_rv$cur_sample %||% ""
      if (nzchar(cur)) {
        remaining <- elapsed * (n_total - 1)
        if (remaining < 0) return("finalizando...")
        return(paste0("~", fmt_elapsed(remaining)))
      }
      return("calculando...")
    }
    avg <- elapsed / n_done
    remaining <- avg * (n_total - n_done)
    if (remaining < 0) return("finalizando...")
    paste0("~", fmt_elapsed(remaining))
  })

  output$remaining_pct <- renderText({
    progress_tick()
    total_bytes <- proc_rv$total_bytes %||% 0
    bytes_done <- proc_rv$bytes_done %||% 0
    if (total_bytes > 0 && bytes_done > 0) {
      pct_done <- min(100, bytes_done / total_bytes * 100)
      return(paste0(round(pct_done, 1), "%"))
    }
    n_done <- sum(vapply(proc_rv$samp_stat, identical, logical(1), "done"))
    n_total <- proc_rv$n_total
    if (n_total > 0) {
      cur <- proc_rv$cur_sample %||% ""
      progress <- if (nzchar(cur) && n_done < n_total) {
        (n_done + 0.5) / n_total
      } else {
        n_done / n_total
      }
      pct_done <- min(100, max(0, progress * 100))
      return(paste0(round(pct_done, 1), "%"))
    }
    "—"
  })

  output$samp_prog_text <- renderText({
    progress_tick()
    n_done  <- sum(vapply(proc_rv$samp_stat, identical, logical(1), "done"))
    n_total <- proc_rv$n_total
    if (n_total == 0) return("—")
    cur <- proc_rv$cur_sample %||% ""
    if (nzchar(cur)) {
      sprintf("%d / %d completadas (%s en curso)", n_done, n_total, cur)
    } else {
      sprintf("%d / %d completadas", n_done, n_total)
    }
  })

  output$checkpoints_ui <- renderUI({
    progress_tick()
    cps    <- proc_rv$checkpoints
    cp_idx <- proc_rv$cp_idx
    if (length(cps) == 0) return(NULL)
    failed <- proc_rv$cp_failed
    kind   <- proc_rv$cp_failed_kind
    # El simbolo va acompanado de texto porque el color por si solo no es un
    # canal accesible: sin el, un usuario con daltonismo o un lector de pantalla
    # no distinguen "completado" de "fallido".
    tags$ul(class = "list-unstyled mb-0",
      lapply(seq_along(cps), function(i) {
        if (!is.na(failed) && i == failed) {
          etiqueta <- if (identical(kind, "cancelled")) "cancelado" else "fallido"
          tags$li(style = "color:#8A1F1F;font-weight:650;",
                  tags$span(if (identical(kind, "cancelled")) "■ " else "✗ "),
                  cps[i],
                  tags$span(class = "small ms-1", paste0("(", etiqueta, ")")))
        } else if (i <= cp_idx) {
          tags$li(style = "color:#315342;font-weight:600;", tags$span("✓ "), cps[i])
        } else if (proc_rv$running && i == cp_idx + 1L) {
          tags$li(style = "color:#8A6D1C; font-weight:650;",
                  tags$span("⟳ "), cps[i])
        } else {
          tags$li(style = "color:#7C9185;", tags$span("◷ "), cps[i])
        }
      })
    )
  })

  output$sample_status_ui <- renderUI({
    progress_tick()
    ss <- proc_rv$samp_stat
    if (length(ss) == 0) return(NULL)
    tags$table(class = "table table-sm",
      tags$thead(tags$tr(tags$th("Muestra"), tags$th("Estado"))),
      tags$tbody(lapply(names(ss), function(s) {
        st <- ss[[s]]
        icon_el <- switch(st,
          done      = tags$span(style = "color:#315342;font-weight:600;",  "✓ Completada"),
          running   = tags$span(style = "color:#8A6D1C;font-weight:600;", "⟳ Procesando..."),
          cancelled = tags$span(style = "color:#60756A;", "■ Cancelada"),
          tags$span(style = "color:#7C9185;", "◷ Pendiente")
        )
        tags$tr(tags$td(s), tags$td(icon_el))
      }))
    )
  })

  output$cmd_preview_text <- renderText({
    tryCatch(shared$workflow_cmd(),
             error = function(e) "Faltan campos para generar el comando.")
  })

  # ── Ejecutar workflow ──────────────────────────────────────────────────────
  observeEvent(input$run_btn, {
    errs <- shared$val_errors()
    if (length(errs) > 0) {
      showNotification(paste(errs, collapse = "\n"), type = "error", duration = 8); return()
    }
    if (!file.exists(workflow_path)) {
      showNotification("workflow.sh no encontrado.", type = "error"); return()
    }
    if (!nzchar(shared$output_dir_val())) {
      state$pending_output_dir(create_run_output_dir(outputs_dir, input$analysis_type, shared$effective_tool()))
    }

    genome_path <- tryCatch(
      prepare_uploaded_input_file(input$genome_file_upload, shared$output_dir_val(), prefix = "genome"),
      error = function(e) {
        showNotification(paste0("Error al preparar el archivo FASTA: ", conditionMessage(e)), type = "error")
        NULL
      }
    )
    if (is.null(genome_path)) return()

    annotation_path <- tryCatch({
      if (isTRUE(input$analysis_type == "alignment")) {
        prepare_uploaded_input_file(input$annotation_file_upload, shared$output_dir_val(), prefix = "annotation")
      } else {
        prepare_uploaded_input_file(input$annotation_file_pseudo_upload, shared$output_dir_val(), prefix = "annotation", optional = TRUE)
      }
    }, error = function(e) {
      showNotification(paste0("Error al preparar el archivo de anotación: ", conditionMessage(e)), type = "error")
      NULL
    })
    if (is.null(annotation_path)) return()

    cmd <- sprintf(
      "bash %s --INPUT %s --OUTPUT %s --GENOME_FILE %s --ANNOTATION_FILE %s --ALIGNMENT_TYPE %s --READ_TYPE %s --FRAGMENT_LENGTH %s --FRAGMENT_SD %s",
      shQuote(workflow_path), shQuote(shared$input_dir_val()), shQuote(shared$output_dir_val()),
      shQuote(genome_path), shQuote(annotation_path), shQuote(shared$effective_tool()), shQuote(shared$effective_read_type()),
      shQuote(shared$effective_fragment_length()), shQuote(shared$effective_fragment_sd())
    )

    samps <- shared$samples_eff()
    state$run_params_rv(list(
      analysis_type = input$analysis_type, tool = shared$effective_tool(),
      input_dir = shared$input_dir_val(), output_dir = shared$output_dir_val(),
      genome_file = genome_path, annotation_file = annotation_path,
      n_samples = length(samps), read_type = read_type_label(shared$effective_read_type()),
      fragment_length = shared$effective_fragment_length(), fragment_sd = shared$effective_fragment_sd(),
      started_at = Sys.time(),
      r_version = paste(R.version$major, R.version$minor, sep = ".")
    ))

    state$analysis_done(FALSE)
    state$data_rv$counts_ready <- FALSE
    shinyjs::disable("run_btn")
    shinyjs::disable("stop_btn")
    on.exit({ if (!proc_rv$running) shinyjs::enable("run_btn") }, add = TRUE)

    # Checkpoints segun tipo de analisis
    cps <- if (input$analysis_type == "alignment")
      c("Construyendo indice Bowtie2",
        "Control de calidad inicial (FastQC)",
        "Alineamiento de muestras (Bowtie2 + fastp)",
        "Procesando BAM (samtools sort + index)",
        "Conteos por gen (featureCounts)",
        "Control de calidad post-trimming (FastQC)",
        "Informe global (MultiQC)")
    else
      c(paste0("Construyendo indice (", shared$effective_tool(), ")"),
        "Control de calidad inicial (FastQC)",
        paste0("Cuantificacion de muestras (", shared$effective_tool(), " + fastp)"),
        "Importando cuantificaciones",
        "Matriz de conteos",
        "Control de calidad post-trimming (FastQC)",
        "Informe global (MultiQC)")

    proc_rv$checkpoints <- cps
    proc_rv$cp_idx      <- 0L
    # Una ejecucion nueva parte sin el fallo de la anterior.
    proc_rv$cp_failed      <- NA_integer_
    proc_rv$cp_failed_kind <- NA_character_
    proc_rv$n_total     <- length(samps)
    proc_rv$samp_stat   <- setNames(as.list(rep("pending", length(samps))), samps)
    proc_rv$cur_sample  <- NULL
    proc_rv$sample_sizes <- sample_fastq_sizes(shared$input_dir_val(), samps, shared$effective_read_type())
    proc_rv$total_bytes <- sum(proc_rv$sample_sizes)
    proc_rv$bytes_done <- 0
    proc_rv$start_time  <- Sys.time()
    proc_rv$end_time    <- NULL
    proc_rv$last_output_time <- Sys.time()
    proc_rv$last_heartbeat   <- NULL
    proc_rv$log_file <- file.path(shared$output_dir_val(), "workflow_live.log")
    proc_rv$log_seen_size <- 0
    cat("", file = proc_rv$log_file, append = FALSE)

    append_run_log(paste0(
      ts_log("=== Iniciando analisis ==="), "\n",
      ts_log("Log completo en: "), proc_rv$log_file, "\n",
      ts_log(paste0("Lanzando: ", cmd)), "\n"))

    if (HAS_PROCESSX) {
      # Ejecucion no bloqueante con processx
      proc <- tryCatch(
        processx::process$new("bash", c("-lc", cmd),
                              stdout = proc_rv$log_file, stderr = "2>&1",
                              env = c("PATH" = Sys.getenv("PATH"))),
        error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL }
      )
      if (!is.null(proc)) {
        proc_rv$proc    <- proc
        proc_rv$running <- TRUE
        shinyjs::enable("stop_btn")
      } else {
        shinyjs::enable("run_btn")
      }
    } else {
      # Fallback bloqueante con system2 — reusa finalize_run para la post-ejecucion
      withProgress(message = "Ejecutando analisis RNA-seq...", value = 0, {
        setProgress(0.05, detail = "Preparando comando...")
        setProgress(0.15, detail = cps[1])
        setProgress(0.25, detail = cps[2])
        result <- tryCatch(
          system2("bash", c("-lc", cmd), stdout = TRUE, stderr = TRUE),
          error = function(e) structure(conditionMessage(e), status = 1L)
        )
        exit_status <- attr(result, "status") %||% 0L
        append_run_log(paste0(terminal_text(result), "\n",
                              ts_log(paste0("Codigo de salida: ", exit_status)), "\n"))
        setProgress(0.85, detail = "Cargando resultados...")
        finalize_run(exit_status, state$run_params_rv()$output_dir %||% "")
        setProgress(1, detail = if (exit_status == 0) "Completado." else "Finalizado con error.")
      })
    }
  })

  # ── Boton de detener proceso ───────────────────────────────────────────────
  observeEvent(input$stop_btn, {
    req(HAS_PROCESSX, proc_rv$running, !is.null(proc_rv$proc))
    proc <- proc_rv$proc
    if (proc$is_alive()) {
      proc$kill()
      proc$wait(2000)
      proc_rv$running <- FALSE
      proc_rv$end_time <- Sys.time()
      proc_rv$proc <- NULL
      # Cancelar no completa el pipeline: el paso en curso queda marcado como
      # cancelado y los siguientes, pendientes.
      proc_rv$cp_failed <- min(proc_rv$cp_idx + 1L, length(proc_rv$checkpoints))
      proc_rv$cp_failed_kind <- "cancelled"
      proc_rv$cur_sample <- NULL
      ss <- proc_rv$samp_stat
      if (length(ss) > 0L) {
        proc_rv$samp_stat <- lapply(ss, function(x) if (x %in% c("running", "pending")) "cancelled" else x)
      }
      log_line("=== Proceso detenido por el usuario ===")
      shinyjs::disable("stop_btn")
      shinyjs::enable("run_btn")
      showNotification("Proceso detenido.", type = "warning")
    }
  })

  # ── Log y refresh de archivos ──────────────────────────────────────────────
  output$run_log <- renderText({ state$log_text() })

  observeEvent(input$refresh_btn, {
    p <- state$run_params_rv()
    out_dir <- p$output_dir %||% shared$output_dir_val() %||% ""
    state$output_files_rv(file_table_for_dir(out_dir))
  })

  invisible(NULL)
}
