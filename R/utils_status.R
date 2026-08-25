#' utils_status.R
#' Estado de una run a partir de su log y badges visuales.

#' Lee el estado de salida que el workflow deja como fichero.
#'
#' @return list(exit_code, status, finished_at) o NULL si no existe.
read_exit_status <- function(out_dir) {
  f <- file.path(out_dir, "exit_status.tsv")
  if (!file.exists(f)) return(NULL)
  df <- tryCatch(
    utils::read.delim(f, header = FALSE, col.names = c("key", "value"),
                      stringsAsFactors = FALSE, quote = ""),
    error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(NULL)
  as.list(stats::setNames(trimws(df$value), trimws(df$key)))
}

#' Infiere el estado de una ejecución: "completado" / "error" / "incompleto" / "sin log".
#'
#' Prioriza el fichero `exit_status.tsv` que escribe el workflow con un `trap
#' EXIT`, porque es un dato explícito. El respaldo por texto del log queda para
#' ejecuciones anteriores a que ese fichero existiera, y es fragil por dos
#' motivos: depende de una frase concreta en ingles, y clasifica como fallida
#' cualquier ejecución cuyo log contenga la palabra "Error", aunque venga de un
#' aviso inocuo de una herramienta.
status_from_log <- function(out_dir) {
  st <- read_exit_status(out_dir)
  if (!is.null(st) && !is.null(st$status)) {
    return(switch(st$status, success = "completado", error = "error", "incompleto"))
  }
  log_file <- file.path(out_dir, "workflow_live.log")
  if (!file.exists(log_file)) return("sin log")
  txt <- read_tail_text(log_file, max_bytes = 512000L)
  if (grepl("Analysis completed successfully|Análisis finalizado OK", txt, ignore.case = TRUE))
    return("completado")
  if (grepl("ERROR|Error \\(código|fallo en la línea", txt, ignore.case = TRUE))
    return("error")
  "incompleto"
}

#' Versiones de las herramientas registradas por el workflow.
#' @return data.frame(tool, versión, path) o NULL.
read_tool_versions <- function(out_dir) {
  f <- file.path(out_dir, "versions.tsv")
  if (!file.exists(f)) return(NULL)
  tryCatch(utils::read.delim(f, stringsAsFactors = FALSE, quote = ""),
           error = function(e) NULL)
}

#' Checksums de las entradas registrados por el workflow.
#' @return data.frame(file, size_bytes, md5) o NULL.
read_input_checksums <- function(out_dir) {
  f <- file.path(out_dir, "checksums.tsv")
  if (!file.exists(f)) return(NULL)
  tryCatch(utils::read.delim(f, stringsAsFactors = FALSE, quote = ""),
           error = function(e) NULL)
}

#' Renderiza un badge HTML con el estado de una ejecución.
#'
#' Usa las clases `.pill` de la hoja de estilos en lugar de colores en línea:
#' antes cada sitio que mostraba un estado lo pintaba con su propia paleta, de
#' modo que el mismo "completado" se veia distinto en la portada, en las
#' metricas y en la tabla de interpretación.
status_badge <- function(status) {
  clase <- switch(status,
    completado = "pill pill-ok",
    error      = "pill pill-bad",
    incompleto = "pill pill-warn",
    "pill pill-neutral"
  )
  tags$span(class = clase, status)
}

#' Infiere si una run es Paired-end o Single-end mirando 02_trimmed_reads/.
#' Devuelve "Paired-end" o "Single-end" (texto). Fallback: "Paired-end".
infer_read_type_from_dir <- function(out_dir) {
  trimmed <- file.path(out_dir, "02_trimmed_reads")
  if (!dir.exists(trimmed)) return("Paired-end")
  files <- list.files(trimmed, pattern = "_trimmed\\.fastq\\.gz$", full.names = FALSE)
  if (!length(files)) return("Paired-end")
  has_r1 <- any(grepl("_R1_trimmed\\.fastq\\.gz$", files))
  has_r2 <- any(grepl("_R2_trimmed\\.fastq\\.gz$", files))
  if (has_r1 && has_r2) return("Paired-end")
  # Si hay *_trimmed.fastq.gz sin R1/R2, es single-end
  has_se <- any(!grepl("_R[12]_trimmed\\.fastq\\.gz$", files))
  if (has_se && !has_r1 && !has_r2) return("Single-end")
  "Paired-end"
}

#' Lee el run_params.tsv que deja workflow.sh en el directorio de salida.
#' Devuelve una lista con lo que haya, o list() si el fichero no existe (las
#' ejecuciones anteriores a su introducción no lo tienen).
#' @param filename fichero clave-valor a leer. El mismo formato lo usan
#'   `run_params.tsv` (parámetros del pipeline) y `deg_params.tsv` (parámetros
#'   de un análisis diferencial persistido), así que comparten lector.
read_run_params_file <- function(out_dir, filename = "run_params.tsv") {
  f <- file.path(out_dir, filename)
  if (!file.exists(f)) return(list())
  df <- tryCatch(
    utils::read.delim(f, header = FALSE, col.names = c("key", "value"),
                      stringsAsFactors = FALSE, quote = ""),
    error = function(e) NULL
  )
  if (is.null(df) || !nrow(df)) return(list())
  stats::setNames(as.list(trimws(df$value)), trimws(df$key))
}

#' Ruta del fichero de anotación de una ejecución guardada, si se conoce.
#' La necesitan tximport (para el mapa transcrito-gen) y la metrica de rRNA.
annotation_file_for_run <- function(out_dir) {
  p <- read_run_params_file(out_dir)
  af <- p$annotation_file %||% ""
  if (nzchar(af) && file.exists(af)) af else NULL
}

#' Infiere parámetros de una run pasada solo a partir del directorio de salida.
#' Útil cuando el usuario abre una carpeta de outputs/ sin contexto de sesión.
infer_result_params <- function(out_dir, workflow_path) {
  saved_params <- read_run_params_file(out_dir)
  tool <- if (dir.exists(file.path(out_dir, "03_alignments", "bowtie2"))) {
    "bowtie2"
  } else if (dir.exists(file.path(out_dir, "03_alignments", "subjunc"))) {
    "subjunc"
  } else if (dir.exists(file.path(out_dir, "03_alignments", "salmon"))) {
    "salmon"
  } else if (dir.exists(file.path(out_dir, "03_alignments", "kallisto"))) {
    "kallisto"
  } else if (isTRUE(saved_params$tool %in% c("bowtie2", "subjunc", "salmon", "kallisto"))) {
    # La carpeta 03_alignments puede no estar: se borra para ahorrar espacio (los
    # BAM son lo más pesado de una ejecución) o la ejecución se copio sin ella.
    # El workflow deja la herramienta escrita en run_params.tsv, así que se
    # respeta antes de darla por desconocida. Sin esto, una ejecución con su
    # matriz de conteos intacta no se podia cargar.
    saved_params$tool
  } else if (file.exists(file.path(out_dir, "04_counts", "count_matrix.tsv"))) {
    # Último recurso: hay matriz por gen pero no consta como se genero. Se lee
    # igual que la de featureCounts, que es un TSV de genes x muestras.
    "bowtie2"
  } else {
    "desconocida"
  }
  analysis <- if (tool %in% c("bowtie2", "subjunc")) "alignment" else "pseudo"
  saved <- saved_params
  counts <- tryCatch(
    load_counts_from_workflow(out_dir, tool, annotation_file = annotation_file_for_run(out_dir)),
    error = function(e) NULL
  )
  stat <- file.info(out_dir)
  dash <- function(x) if (is.null(x) || !nzchar(x %||% "")) "—" else x
  list(
    analysis_type = analysis,
    tool = tool,
    input_dir = dash(saved$input_dir),
    output_dir = out_dir,
    genome_file = dash(saved$genome_file),
    annotation_file = dash(saved$annotation_file),
    n_samples = if (!is.null(counts)) ncol(counts) else "—",
    read_type = infer_read_type_from_dir(out_dir),
    started_at = stat$mtime %||% Sys.time(),
    r_version = paste(R.version$major, R.version$minor, sep = ".")
  )
}

# ── Metricas de coste de la ejecución ───────────────────────────────────────
#
# El workflow mide dos cosas que hacen falta para decidir si un conjunto de
# datos cabe en el equipo que se tiene: cuanto tarda y cuanta memoria pide.
#
# Sobre la memoria: el dato útil es el PICO, no un promedio ni un mínimo. Un
# promedio no dice nada —la mayor parte del tiempo el pipeline esta escribiendo
# a disco— y el mínimo sería el consumo en reposo. Lo que responde a "cuanta
# RAM necesito" es el máximo que llego a ocupar el árbol de procesos.

#' Metricas por paso que escribe el workflow.
#'
#' @return data.frame(paso, segundos, duración, pico_rss_mb) o NULL.
read_run_metrics <- function(out_dir) {
  f <- file.path(out_dir, "metrics.tsv")
  if (!file.exists(f)) return(NULL)
  df <- tryCatch(
    utils::read.delim(f, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(NULL)
  df
}

#' Duración en segundos como texto legible.
fmt_duracion <- function(segundos) {
  s <- suppressWarnings(as.numeric(segundos))
  if (length(s) != 1L || is.na(s) || s < 0) return("—")
  h <- s %/% 3600; m <- (s %% 3600) %/% 60; seg <- round(s %% 60)
  if (h > 0) sprintf("%dh %02dm", h, m)
  else if (m > 0) sprintf("%dm %02ds", m, seg)
  else sprintf("%ds", seg)
}

#' Memoria en MB como texto legible, pasando a GB cuando toca.
fmt_memoria <- function(mb) {
  v <- suppressWarnings(as.numeric(mb))
  if (length(v) != 1L || is.na(v) || v <= 0) return("—")
  if (v >= 1024) sprintf("%.1f GB", v / 1024) else sprintf("%.0f MB", v)
}

# ── Herramientas del pipeline ───────────────────────────────────────────────
#
# El workflow comprueba en su primer paso que las ocho están en el PATH y aborta
# si falta alguna. Hasta ahora la interfaz no lo comprobaba, de modo que dejaba
# lanzar una ejecución condenada a fallar y el motivo quedaba enterrado en el
# log: "Error (código 1)" arriba y la causa real veinte líneas más abajo.
#
# La causa habitual no es que falten instaladas, sino arrancar la aplicación sin
# activar el entorno donde viven. Por eso, cuando no se encuentran en el PATH se
# buscan en los sitios donde conda las deja: poder decir "están instaladas pero
# no en el PATH" ahorra el rato de comprobar si hay que instalar algo.

#' Herramientas que exige cada estrategia, en el mismo orden en que las
#' comprueba workflow.sh.
herramientas_requeridas <- function(analysis_type = "alignment", tool = "bowtie2") {
  comunes <- c("fastqc", "fastp", "multiqc")
  if (identical(analysis_type, "alignment")) {
    # subjunc y subread-buildindex vienen en el mismo paquete que featureCounts,
    # de modo que esta ruta no exige instalar nada nuevo.
    alineador <- if (identical(tool, "subjunc")) c("subjunc", "subread-buildindex")
                 else "bowtie2"
    return(c(comunes, alineador, "samtools", "featureCounts"))
  }
  c(comunes, if (identical(tool, "kallisto")) "kallisto" else "salmon")
}

#' Entornos de conda donde buscar cuando una herramienta no está en el PATH.
entornos_conda_probables <- function() {
  bases <- c(Sys.getenv("CONDA_PREFIX", ""),
             path.expand("~/miniforge3/envs"), path.expand("~/miniconda3/envs"),
             path.expand("~/anaconda3/envs"), path.expand("~/mambaforge/envs"),
             "/opt/miniforge3/envs", "/opt/conda/envs")
  bases <- bases[nzchar(bases)]
  dirs <- unlist(lapply(bases, function(b) {
    if (!dir.exists(b)) return(character(0))
    # CONDA_PREFIX apunta al entorno, no a la carpeta de entornos.
    if (dir.exists(file.path(b, "bin"))) return(b)
    list.dirs(b, recursive = FALSE)
  }))
  unique(dirs[dir.exists(file.path(dirs, "bin"))])
}

#' Estado de las herramientas del pipeline.
#'
#' @return list(faltan, encontradas_fuera, entorno). `faltan` son las que no
#'   están en el PATH; `entorno` es la ruta de un entorno que las contiene
#'   todas, si existe, para poder decir exactamente que hacer.
#' `necesarias` y `entornos` entran como argumentos con valor por defecto para
#' poder probar la lógica sin depender de que herramientas haya instaladas en la
#' maquina donde corren los tests.
comprobar_herramientas <- function(analysis_type = "alignment", tool = "bowtie2",
                                   entornos = entornos_conda_probables(),
                                   necesarias = herramientas_requeridas(analysis_type, tool)) {
  faltan <- necesarias[!nzchar(Sys.which(necesarias))]
  if (!length(faltan)) return(list(faltan = character(0), entorno = NULL))

  # Un entorno solo sirve si tiene TODAS las que faltan: mandar al usuario a uno
  # que resuelve la mitad del problema es peor que no decir nada.
  entorno <- NULL
  for (e in entornos) {
    if (all(file.exists(file.path(e, "bin", faltan)))) { entorno <- e; break }
  }
  list(faltan = faltan, entorno = entorno)
}

#' Mensaje de error para la lista de validación, o NULL si está todo.
mensaje_herramientas <- function(estado) {
  if (!length(estado$faltan)) return(NULL)
  lista <- paste(estado$faltan, collapse = ", ")
  if (!is.null(estado$entorno)) {
    return(paste0(
      "Herramientas no encontradas en el PATH (", lista, "). ",
      "Están instaladas en ", estado$entorno, ", pero la aplicación se ha ",
      "arrancado sin ese entorno activo: cierrala y usa lanzar_app.sh."))
  }
  paste0("Herramientas no encontradas (", lista,
         "). Instalalas con requirements.sh o activa el entorno donde esten.")
}
