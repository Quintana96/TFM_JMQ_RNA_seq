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

#' Infiere el estado de una ejecucion: "completado" / "error" / "incompleto" / "sin log".
#'
#' Prioriza el fichero `exit_status.tsv` que escribe el workflow con un `trap
#' EXIT`, porque es un dato explicito. El respaldo por texto del log queda para
#' ejecuciones anteriores a que ese fichero existiera, y es fragil por dos
#' motivos: depende de una frase concreta en ingles, y clasifica como fallida
#' cualquier ejecucion cuyo log contenga la palabra "Error", aunque venga de un
#' aviso inocuo de una herramienta.
status_from_log <- function(out_dir) {
  st <- read_exit_status(out_dir)
  if (!is.null(st) && !is.null(st$status)) {
    return(switch(st$status, success = "completado", error = "error", "incompleto"))
  }
  log_file <- file.path(out_dir, "workflow_live.log")
  if (!file.exists(log_file)) return("sin log")
  txt <- read_tail_text(log_file, max_bytes = 512000L)
  if (grepl("Analysis completed successfully|Analisis finalizado OK", txt, ignore.case = TRUE))
    return("completado")
  if (grepl("ERROR|Error \\(codigo|fallo en la linea", txt, ignore.case = TRUE))
    return("error")
  "incompleto"
}

#' Versiones de las herramientas registradas por el workflow.
#' @return data.frame(tool, version, path) o NULL.
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

#' Renderiza un badge HTML con color pastel segun el estado
status_badge <- function(status) {
  color <- switch(status,
    completado = "#B8D8BA",
    error = "#F4A6A6",
    incompleto = "#F6D58A",
    "#D7EEF1"
  )
  text_color <- switch(status,
    error = "#5A2323",
    incompleto = "#5C4A16",
    "#20332A"
  )
  tags$span(
    style = paste0(
      "display:inline-block;padding:3px 9px;border-radius:999px;",
      "color:", text_color, ";background:", color, ";font-weight:700;"
    ),
    status
  )
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
#' ejecuciones anteriores a su introduccion no lo tienen).
read_run_params_file <- function(out_dir) {
  f <- file.path(out_dir, "run_params.tsv")
  if (!file.exists(f)) return(list())
  df <- tryCatch(
    utils::read.delim(f, header = FALSE, col.names = c("key", "value"),
                      stringsAsFactors = FALSE, quote = ""),
    error = function(e) NULL
  )
  if (is.null(df) || !nrow(df)) return(list())
  stats::setNames(as.list(trimws(df$value)), trimws(df$key))
}

#' Ruta del fichero de anotacion de una ejecucion guardada, si se conoce.
#' La necesitan tximport (para el mapa transcrito-gen) y la metrica de rRNA.
annotation_file_for_run <- function(out_dir) {
  p <- read_run_params_file(out_dir)
  af <- p$annotation_file %||% ""
  if (nzchar(af) && file.exists(af)) af else NULL
}

#' Infiere parametros de una run pasada solo a partir del directorio de salida.
#' Util cuando el usuario abre una carpeta de outputs/ sin contexto de sesion.
infer_result_params <- function(out_dir, workflow_path) {
  saved_params <- read_run_params_file(out_dir)
  tool <- if (dir.exists(file.path(out_dir, "03_alignments", "bowtie2"))) {
    "bowtie2"
  } else if (dir.exists(file.path(out_dir, "03_alignments", "salmon"))) {
    "salmon"
  } else if (dir.exists(file.path(out_dir, "03_alignments", "kallisto"))) {
    "kallisto"
  } else if (isTRUE(saved_params$tool %in% c("bowtie2", "salmon", "kallisto"))) {
    # La carpeta 03_alignments puede no estar: se borra para ahorrar espacio (los
    # BAM son lo mas pesado de una ejecucion) o la ejecucion se copio sin ella.
    # El workflow deja la herramienta escrita en run_params.tsv, asi que se
    # respeta antes de darla por desconocida. Sin esto, una ejecucion con su
    # matriz de conteos intacta no se podia cargar.
    saved_params$tool
  } else if (file.exists(file.path(out_dir, "04_counts", "count_matrix.tsv"))) {
    # Ultimo recurso: hay matriz por gen pero no consta como se genero. Se lee
    # igual que la de featureCounts, que es un TSV de genes x muestras.
    "bowtie2"
  } else {
    "desconocida"
  }
  analysis <- if (identical(tool, "bowtie2")) "alignment" else "pseudo"
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
