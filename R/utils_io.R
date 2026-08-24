#' utils_io.R
#' Lectura de archivos, manejo de uploads y descubrimiento de runs en outputs/.

#' Lee la cola (tail) de un archivo grande sin cargarlo entero en memoria
read_tail_text <- function(path, max_bytes = 262144L) {
  if (!file.exists(path)) return("")
  size <- file.info(path)$size %||% 0
  if (is.na(size) || size <= 0) return("")

  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  offset <- max(0, size - max_bytes)
  seek(con, where = offset, origin = "start")
  terminal_text(rawToChar(readBin(con, what = "raw", n = size - offset)))
}

#' Copia un archivo subido por Shiny a una ruta persistente dentro del output_dir
prepare_uploaded_input_file <- function(upload, output_dir, prefix = NULL, optional = FALSE) {
  if (is.null(upload) || nrow(upload) == 0) {
    if (optional) return("/dev/null")
    return("")
  }
  if (!nzchar(output_dir)) stop("Output directory required to store uploaded input files.")
  upload_dir <- file.path(output_dir, "uploaded_inputs")
  if (!dir.exists(upload_dir)) dir.create(upload_dir, recursive = TRUE, showWarnings = FALSE)
  dest_name <- basename(upload$name)
  if (!is.null(prefix) && nzchar(prefix)) {
    dest_name <- paste0(prefix, "_", dest_name)
  }
  dest_path <- file.path(upload_dir, dest_name)
  # Decisión: copiar por defecto. Solo saltamos la copia si destino existe Y
  # tanto el tamaño como el hash MD5 coinciden con el origen. Cualquier fallo
  # en la comparación fuerza copia (más seguro que asumir igualdad).
  should_copy <- TRUE
  if (file.exists(dest_path)) {
    same_size <- tryCatch(
      identical(file.info(upload$datapath)$size, file.info(dest_path)$size),
      error = function(e) FALSE
    )
    if (isTRUE(same_size)) {
      same_hash <- tryCatch({
        src_h <- unname(tools::md5sum(upload$datapath))
        dst_h <- unname(tools::md5sum(dest_path))
        !is.na(src_h) && !is.na(dst_h) && identical(src_h, dst_h)
      }, error = function(e) FALSE)
      if (isTRUE(same_hash)) should_copy <- FALSE
    }
  }
  if (should_copy) {
    if (!file.copy(upload$datapath, dest_path, overwrite = TRUE)) {
      stop(sprintf("No se pudo copiar %s a %s", upload$datapath, dest_path))
    }
  }
  dest_path
}

#' Directorio base de outputs (relativo al working dir)
outputs_base_dir <- function() file.path(getwd(), "outputs")

#' Etiqueta segura para una run: <fecha>_<análisis>_<tool>
safe_run_label <- function(analysis_type, tool, time = Sys.time()) {
  analysis <- if (identical(analysis_type, "alignment")) "alineamiento" else "pseudoalineamiento"
  paste(format(time, "%Y%m%d_%H%M%S"), analysis, tool, sep = "_")
}

#' Crea un directorio nuevo dentro de base_dir; si ya existe, le añade _2, _3, ...
create_run_output_dir <- function(base_dir, analysis_type, tool) {
  if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  label <- safe_run_label(analysis_type, tool)
  out <- file.path(base_dir, label)
  i <- 1L
  while (dir.exists(out)) {
    i <- i + 1L
    out <- file.path(base_dir, paste0(label, "_", i))
  }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  normalizePath(out, mustWork = TRUE)
}

#' Lista carpetas de resultados ordenadas por mtime (más reciente primero)
list_result_dirs <- function(base_dir) {
  if (!dir.exists(base_dir)) return(character(0))
  dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[file.info(dirs)$isdir %||% FALSE]
  if (!length(dirs)) return(character(0))
  mt <- file.info(dirs)$mtime
  dirs[order(mt, decreasing = TRUE, na.last = TRUE)]
}

#' Devuelve un setNames(paths, label) listo para selectizeInput
result_choices <- function(base_dir) {
  dirs <- list_result_dirs(base_dir)
  if (!length(dirs)) return(character(0))
  stats <- file.info(dirs)
  labels <- paste0(basename(dirs), "  (", format(stats$mtime, "%Y-%m-%d %H:%M"), ")")
  setNames(dirs, labels)
}

#' Construye data.frame de archivos / tamaños a partir de paths relativos
file_table_for_files <- function(out_dir, files) {
  if (!length(files)) {
    return(data.frame(Archivo = character(), Tamaño = character(),
                      stringsAsFactors = FALSE, check.names = FALSE))
  }
  sz <- file.info(file.path(out_dir, files))$size
  data.frame(Archivo = files, Tamaño = sapply(sz, fmt_bytes),
             stringsAsFactors = FALSE, check.names = FALSE)
}

#' Lista todos los archivos (recursive) de un directorio para mostrar en la app
file_table_for_dir <- function(out_dir) {
  if (!nzchar(out_dir) || !dir.exists(out_dir)) {
    return(data.frame(Archivo = "Directorio no existe.", `Tamaño` = "—",
                      stringsAsFactors = FALSE, check.names = FALSE))
  }
  files <- list.files(out_dir, recursive = TRUE, full.names = FALSE)
  if (!length(files)) {
    return(data.frame(Archivo = "Sin archivos.", `Tamaño` = "—",
                      stringsAsFactors = FALSE, check.names = FALSE))
  }
  file_table_for_files(out_dir, files)
}

#' Wrapper seguro alrededor de read.delim que devuelve NULL en error
read_tsv_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = ""),
    error = function(e) NULL
  )
}

#' Tamaño total (bytes) de un directorio (recursivo)
dir_size <- function(path) {
  files <- list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  if (!length(files)) return(0)
  sum(file.info(files)$size, na.rm = TRUE)
}
