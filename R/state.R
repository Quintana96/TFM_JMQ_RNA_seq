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

  # Roots para shinyFiles (compartidos entre las pestanas 1 y 3).
  #
  # Se pasan como FUNCION, no como vector: shinyFiles la reevalua en cada
  # peticion, de modo que un disco externo montado con la app ya abierta aparece
  # sin reiniciarla. Antes eran solo el directorio del proyecto y el home, asi
  # que no habia forma de llegar a FASTQ, genomas o anotaciones guardados en
  # otro volumen o fuera de la carpeta del usuario.
  #
  # getVolumes() enumera los volumenes montados (en macOS incluye
  # /Volumes/<disco>, y los externos segun se conecten). Se añade la raiz del
  # sistema para poder navegar a cualquier ruta absoluta.
  roots <- function() {
    vols <- tryCatch(shinyFiles::getVolumes()(), error = function(e) character(0))
    base <- c("Proyecto" = normalizePath(getwd(), mustWork = FALSE),
              "Inicio"   = normalizePath("~", mustWork = FALSE),
              "Sistema (/)" = "/")
    # Sin duplicar rutas: si un volumen coincide con una de las anteriores, se
    # queda la primera (la de nombre mas descriptivo).
    todo <- c(base, vols)
    todo[!duplicated(normalizePath(todo, mustWork = FALSE))]
  }
  # El selector de resultados arranca en outputs/ por comodidad, pero tambien
  # permite salir de ahi: una ejecucion guardada puede vivir en otro disco.
  roots_results <- function() {
    c("Resultados" = normalizePath(outputs_dir, mustWork = FALSE), roots())
  }

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
    data.frame(Archivo = character(), `Tamaño` = character(),
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
    # Paso en el que se interrumpio la ejecucion, y por que. Sin esto, al fallar
    # o cancelar se marcaba la lista ENTERA como completada y el usuario veia
    # todos los pasos en verde justo despues del aviso de error.
    cp_failed       = NA_integer_,
    cp_failed_kind  = NA_character_,   # "error" | "cancelled"
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
  # `fdr` y `lfc_threshold` guardan los parametros con los que se calculo la
  # tabla que hay ahora en `results`, no los de la interfaz en este instante:
  # son los unicos validos para declarar significacion.
  #
  # Ya no exigen relanzar. Ajustar el modelo y extraer la tabla son dos
  # operaciones separadas (ver `deg_reextract()`), y la segunda cuesta un 4 % de
  # la primera, asi que el FDR y el umbral del test se recalculan en vivo sobre
  # el mismo ajuste. Lo que sigue exigiendo relanzar es lo que cambia el AJUSTE:
  # motor, diseño, batch, variables sustitutas, prefiltrado y encogido.
  state$deg_rv <- reactiveValues(
    counts        = NULL,
    meta          = NULL,
    method        = NULL,
    results       = NULL,
    vst_mat       = NULL,
    run_at        = NULL,
    fdr           = 0.05,
    lfc_threshold = 0,
    # Objeto ajustado reutilizable (`dds`, `glmQLFit` o `lmFit` segun el motor).
    # Es el estado pesado de la sesion —decenas de MB con un dataset humano— y
    # es lo que se cambia por no repetir cinco segundos de ajuste cada vez que
    # se mueve un deslizador.
    fit           = NULL,
    # Parametros con los que se extrajo `results`, para no reextraer en balde.
    extract_params = NULL,
    # Huella de los parametros que definen el AJUSTE. Si la interfaz deja de
    # coincidir con ella, lo que se esta viendo corresponde a otro modelo y hay
    # que decirlo: sin esta comparacion, mover un selector que exige reajuste no
    # produce ningun efecto visible y se lee como que la aplicacion esta rota.
    fit_signature  = NULL,
    reextracted_at = NULL,
    contrast      = NULL,
    n_levels      = NA_integer_,
    shrink        = "ninguno",
    prefilter     = NULL,
    padj_method   = "BH",
    # Diagnosticos post-ajuste (solo los produce DESeq2): tabla de dispersiones
    # para el equivalente de plotDispEsts y reparto de los maximos de Cook.
    disp_data     = NULL,
    cooks         = NULL,
    # Diseño efectivamente ajustado, coeficiente testeado y coeficientes
    # disponibles (para el selector del modo avanzado).
    #
    # Hay TRES formas del diseño y no son intercambiables:
    #   `design`       etiqueta legible (banner, informe). Para Wilcoxon es prosa.
    #   `design_code`  la misma formula como CODIGO, o NULL si el motor no ajusta
    #                  modelo. La consume el script exportado.
    #   `design_base`  el diseño ANTES de añadir las variables sustitutas, o NULL
    #                  si no se uso sva. Es el modelo con el que hay que
    #                  reestimarlas para reproducir el ajuste; usar `design_code`
    #                  (que ya las incluye) seria circular.
    design        = "~ condition",
    design_code   = NULL,
    design_base   = NULL,
    # Tratamiento de los outliers de Cook con el que se ajusto ("na" | "refit" |
    # "keep"). Cambia que genes tienen padj, asi que es un parametro del test y
    # tiene que viajar al informe y al script.
    outliers      = "na",
    coef          = NULL,
    coef_available = character(0),
    # Contraste y diseño CONGELADOS en el momento del ajuste. Los selectores de
    # la interfaz siguen vivos despues de ajustar, asi que leerlos mas tarde
    # (informe, bootstrap, comparacion de metodos) puede describir un analisis
    # distinto del que se ejecuto. Todo lo que reproduzca o evalue el ajuste
    # debe consumir estos campos y NO `input$...`.
    ref_level      = NULL,   # denominador del contraste
    contrast_num   = NULL,   # numerador del contraste
    batch          = NULL,   # columna de batch incluida en el diseño, si la hay
    design_formula = NULL,   # formula libre usada (incluye las SV), o NULL
    test_coef      = NULL,   # coeficiente elegido en el modo avanzado, o NULL
    quant_tool     = NULL,   # salmon/kallisto: lo necesita el motor Swish
    # Semillas y parametros estocasticos, para poder declararlos en el informe
    # y reproducir el resultado exactamente.
    seeds         = NULL,
    # Procedencia de la matriz: de donde salio (con md5 si es un fichero) y como
    # se resumio a gen (tximport con offsets de longitud, o el respaldo de
    # est_counts redondeados). Lo segundo es una degradacion que el informe debe
    # declarar, porque cambia la validez de los resultados.
    counts_origin = NULL,
    counts_source = NULL,
    run_dir       = NULL,
    # Tabla de correspondencia de la seudonimizacion, si se aplico. NO viaja al
    # informe: se descarga aparte, para que exportar los identificadores reales
    # sea una decision explicita.
    pseudonym_map = NULL,
    # Ultimo enriquecimiento ejecutado: parametros y resumen del resultado.
    enrich        = NULL,
    # Nota sobre la correccion aplicada SOLO a los graficos, si la hay.
    viz_note      = NULL
  )

  # Slot para reactivos expuestos por server_tab_config (rellenado al final
  # de server_tab_config). Se inicializa a NULL para que `state$shared` sea
  # consultable sin error antes de la primera asignacion.
  state$shared <- NULL

  state
}
