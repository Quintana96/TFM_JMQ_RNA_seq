#' state.R
#' Fabrica del estado compartido entre las funciones server de cada pestana.
#'
#' La aplicacion comparte mucho estado entre las tres pestanas (logs, snapshots
#' de configuracion, resultados pendientes, etc.). En lugar de usar moduleServer
#' con namespaces (que obligaria a renombrar todos los IDs Shiny), exponemos un
#' objeto `state` con todos los reactiveVal/reactiveValues. Cada funcion server
#' lo recibe como argumento y opera sobre el.

#' Crea el estado compartido. Debe llamarse dentro de la funcion server,
#' porque depende del scope de session (para addResourcePath).
#'
#' El estado se devuelve como `environment` (no `list`) para que las
#' asignaciones `state$shared <- list(...)` realizadas en server_tab_config
#' sean visibles desde server_tab_processing y server_tab_results.
#' En R, los environments son pass-by-reference; las listas, copy-on-modify.
create_app_state <- function(session) {
  workflow_path <- normalizePath("workflow.sh", mustWork = FALSE)
  outputs_dir <- outputs_base_dir()
  if (!dir.exists(outputs_dir)) dir.create(outputs_dir, recursive = TRUE, showWarnings = FALSE)
  outputs_dir_norm <- normalizePath(outputs_dir, mustWork = TRUE)
  addResourcePath("saved_outputs", outputs_dir_norm)

  # Roots para shinyFiles (compartidos entre las pestanas 1 y 3)
  roots <- c(wd = normalizePath(getwd()), home = normalizePath("~"))
  roots_results <- c(outputs = normalizePath(outputs_dir, mustWork = FALSE))

  state <- new.env(parent = emptyenv())

  # Rutas y constantes
  state$workflow_path <- workflow_path
  state$outputs_dir   <- outputs_dir_norm
  state$roots         <- roots
  state$roots_results <- roots_results

  # Log de ejecucion
  state$log_text <- reactiveVal(paste0(ts_log("Workflow listo.\n")))

  # Tabla de archivos del output_dir actual
  state$output_files_rv <- reactiveVal(
    data.frame(Archivo = character(), `Tamano` = character(),
               stringsAsFactors = FALSE, check.names = FALSE)
  )

  # Navegacion entre pestanas
  state$process_unlocked   <- reactiveVal(FALSE)
  state$analysis_done      <- reactiveVal(FALSE)
  state$config_snap        <- reactiveVal(list())
  state$run_params_rv      <- reactiveVal(list())
  state$pending_output_dir <- reactiveVal("")
  state$results_refresh    <- reactiveVal(Sys.time())

  # Estado del proceso lanzado con processx
  state$proc_rv <- reactiveValues(
    proc            = NULL,
    running         = FALSE,
    start_time      = NULL,
    end_time        = NULL,
    checkpoints     = character(0),
    cp_idx          = 0,
    samp_stat       = list(),
    cur_sample      = NULL,
    n_total         = 0,
    log_file        = NULL,
    log_seen_size   = 0,
    last_output_time = NULL,
    last_heartbeat   = NULL,
    sample_sizes    = list(),
    total_bytes     = 0,
    bytes_done      = 0
  )

  # Matriz de conteos cargada (workflow o uploaded)
  state$data_rv <- reactiveValues(
    count_matrix = NULL,
    counts_ready = FALSE,
    source       = "workflow"
  )

  # Resultado de la pestana DEG (Tab 4). Reactivos para que las renders
  # cuelguen automaticamente al lanzar un nuevo run_deg().
  state$deg_rv <- reactiveValues(
    counts  = NULL,
    meta    = NULL,
    method  = NULL,
    results = NULL,
    vst_mat = NULL,
    run_at  = NULL
  )

  # Slot para reactivos expuestos por server_tab_config (rellenado al final
  # de server_tab_config). Se inicializa a NULL para que `state$shared` sea
  # consultable sin error antes de la primera asignacion.
  state$shared <- NULL

  state
}
