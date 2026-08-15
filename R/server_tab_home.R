#' server_tab_home.R
#' Logica server de la portada (Tab 0):
#'   - estado de cada uno de los cuatro pasos;
#'   - resumen de la sesion;
#'   - navegacion desde las tarjetas.
#'
#' Se registra DESPUES de server_tab_config, porque lee `state$shared`, que esa
#' funcion rellena al final. Aun asi todos los accesos van protegidos: el
#' contenido de la portada se renderiza antes de la primera invalidacion de
#' cualquier reactivo de configuracion.

server_tab_home <- function(input, output, session, state) {

  #' Ejecuciones guardadas en outputs/, para el paso 3.
  saved_runs <- reactive({
    state$results_refresh()
    result_choices(state$outputs_dir)
  })

  #' Numero de muestras detectadas en la configuracion actual, o NA.
  detected_samples <- reactive({
    f <- state$shared$samples_eff
    if (is.null(f)) return(NA_integer_)
    length(tryCatch(f(), error = function(e) character(0)))
  })

  #' Estado del paso 1. `val_errors()` es la misma fuente que usan el checklist
  #' y el boton de continuar, para que los tres no puedan discrepar.
  config_state <- reactive({
    f <- state$shared$val_errors
    if (is.null(f)) return(list(level = "neutral", label = "SIN EMPEZAR", cta = "Configurar"))
    errs <- tryCatch(f(), error = function(e) "pendiente")
    if (!length(errs)) list(level = "ok", label = "LISTO", cta = "Revisar")
    else list(level = "warn", label = paste0(length(errs), " PENDIENTE",
                                             if (length(errs) == 1) "" else "S"),
              cta = "Configurar")
  })

  process_state <- reactive({
    if (isTRUE(state$proc_rv$running))
      return(list(level = "info", label = "EJECUTANDO", cta = "Ver progreso"))
    if (isTRUE(state$analysis_done()) &&
        !identical(state$data_rv$source, "uploaded"))
      return(list(level = "ok", label = "COMPLETADO", cta = "Ver log"))
    if (isTRUE(state$process_unlocked()))
      return(list(level = "warn", label = "LISTO", cta = "Ejecutar"))
    list(level = "neutral", label = "BLOQUEADO", cta = "Requiere el paso 1")
  })

  results_state <- reactive({
    n <- length(saved_runs())
    if (isTRUE(state$data_rv$counts_ready) && identical(state$data_rv$source, "uploaded"))
      return(list(level = "ok", label = "MATRIZ CARGADA", cta = "Explorar"))
    if (n > 0) list(level = "ok",
                    label = paste0(n, " EJECUCION", if (n == 1) "" else "ES"),
                    cta = "Explorar")
    else list(level = "neutral", label = "SIN RESULTADOS", cta = "Explorar")
  })

  deg_state <- reactive({
    res <- state$deg_rv$results
    if (!is.null(res)) {
      return(list(level = "ok",
                  label = toupper(state$deg_rv$method %||% "AJUSTADO"),
                  cta = "Ver resultados"))
    }
    if (isTRUE(state$data_rv$counts_ready) || length(saved_runs()))
      list(level = "warn", label = "CON DATOS", cta = "Analizar")
    else list(level = "neutral", label = "SIN DATOS", cta = "Analizar")
  })

  #' Aviso contextual: una sola frase que dice cual es el siguiente paso util.
  #' Sustituye a que el usuario deduzca por si mismo por que un boton esta gris.
  home_hint <- reactive({
    if (isTRUE(state$proc_rv$running))
      return(div(class = "alert alert-info mb-0 py-2", icon("spinner"),
                 " Hay un workflow en ejecucion. Puedes seguirlo en el paso 2."))
    if (!is.null(state$deg_rv$results)) return(NULL)
    if (isTRUE(state$data_rv$counts_ready) || length(saved_runs()))
      return(div(class = "alert alert-secondary mb-0 py-2", icon("circle-info"),
                 " Ya hay conteos disponibles: puedes ir directamente al paso 4",
                 " sin volver a ejecutar el pipeline."))
    NULL
  })

  output$home_content <- renderUI({
    cfg <- config_state(); prc <- process_state()
    res <- results_state(); deg <- deg_state()

    n_samp <- detected_samples()
    n_runs <- length(saved_runs())

    # La matriz puede haber entrado por dos caminos distintos: la pestana 1
    # (state$data_rv) o directamente la pestana 4 (state$deg_rv$counts). Mirar
    # solo el primero hacia que el resumen dijera "sin matriz" con un analisis
    # DEG ya ajustado en pantalla.
    cm <- state$data_rv$count_matrix %||% state$deg_rv$counts

    fuente <- if (isTRUE(state$data_rv$counts_ready)) {
      switch(state$data_rv$source %||% "workflow",
             uploaded = "Matriz subida", workflow = "Workflow", "Workflow")
    } else if (!is.null(state$deg_rv$counts)) {
      switch(input$deg_source %||% "current",
             saved = "Ejecucion guardada", upload = "Matriz subida", "Ejecucion actual")
    } else "—"

    tiles <- list(
      stat_tile(fuente, "Fuente de los conteos",
                if (identical(fuente, "—")) "neutral" else "ok"),
      stat_tile(if (is.na(n_samp) || !n_samp) "—" else fmt_int(n_samp),
                "Muestras detectadas",
                if (!is.na(n_samp) && n_samp > 0) "ok" else "neutral"),
      stat_tile(if (is.null(cm)) "—" else paste0(fmt_int(nrow(cm)), " x ", fmt_int(ncol(cm))),
                "Matriz en memoria",
                if (is.null(cm)) "neutral" else "ok"),
      stat_tile(fmt_int(n_runs), "Ejecuciones en outputs/",
                if (n_runs > 0) "info" else "neutral")
    )

    ui_tab_home_content(
      steps = list(
        config  = list(pill = status_pill(cfg$label, cfg$level, paste("Configuracion:", cfg$label)),
                       cta = cfg$cta),
        process = list(pill = status_pill(prc$label, prc$level, paste("Procesamiento:", prc$label)),
                       cta = prc$cta),
        results = list(pill = status_pill(res$label, res$level, paste("Resultados:", res$label)),
                       cta = res$cta),
        deg     = list(pill = status_pill(deg$label, deg$level, paste("Expresion diferencial:", deg$label)),
                       cta = deg$cta)
      ),
      tiles = tiles,
      hint  = home_hint()
    )
  })

  # ── Navegacion desde las tarjetas ───────────────────────────────────────────
  nav_to <- function(btn, target) {
    observeEvent(input[[btn]], {
      updateNavbarPage(session, "main_nav", selected = target)
    }, ignoreInit = TRUE)
  }
  nav_to("home_go_config",  "tab_config")
  nav_to("home_go_process", "tab_process")
  nav_to("home_go_results", "tab_results")
  nav_to("home_go_deg",     "tab_deg")

  # Enlace de vuelta a la portada desde cualquier pestana.
  observeEvent(input$btn_home, {
    updateNavbarPage(session, "main_nav", selected = "tab_home")
  }, ignoreInit = TRUE)

  invisible(NULL)
}
