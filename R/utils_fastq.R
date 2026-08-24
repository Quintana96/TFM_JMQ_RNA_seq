#' utils_fastq.R
#' Deteccion y nombrado de muestras FASTQ paired/single end.

#' Genera todas las rutas candidatas para un sample y un set de sufijos
sample_fastq_paths <- function(dir_path, samples, suffixes) {
  file.path(dir_path, as.vector(outer(samples, suffixes, paste0)))
}

#' Quita el sufijo R1 / _1 del nombre base
sample_name_from_r1 <- function(files) {
  sub(FASTQ_R1_PATTERN, "", basename(files), ignore.case = TRUE)
}

#' Quita cualquier sufijo FASTQ (.fastq[.gz], _R1.fastq[.gz], etc.)
sample_name_from_fastq <- function(files) {
  x <- basename(files)
  x <- sub(FASTQ_R1_PATTERN, "", x, ignore.case = TRUE)
  x <- sub(FASTQ_ANY_PATTERN, "", x, ignore.case = TRUE)
  x
}

#' "Paired-end" / "Single-end" en castellano
read_type_label <- function(read_type) {
  if (identical(read_type, "se")) "Single-end" else "Paired-end"
}

#' Detecta muestras (caracter vector) en un directorio segun el read_type
detect_samples <- function(dir_path, read_type = "pe") {
  if (!nzchar(dir_path) || !dir.exists(dir_path)) return(character(0))
  if (identical(read_type, "se")) {
    files <- list.files(
      dir_path,
      pattern = FASTQ_ANY_PATTERN,
      full.names = FALSE, ignore.case = TRUE
    )
    files <- files[!grepl(FASTQ_R2_PATTERN, files, ignore.case = TRUE)]
    unique(sample_name_from_fastq(files))
  } else {
    files <- list.files(
      dir_path,
      pattern = FASTQ_R1_PATTERN,
      full.names = FALSE, ignore.case = TRUE
    )
    unique(sample_name_from_r1(files))
  }
}

#' Vector logico: TRUE si a la muestra le falta el archivo R2
missing_r2 <- function(dir_path, samples) {
  if (length(samples) == 0) return(logical(0))
  vapply(samples, function(s) {
    !any(file.exists(sample_fastq_paths(dir_path, s, FASTQ_R2_SUFFIXES)))
  }, logical(1))
}

#' Suma de tamaños (bytes) de los FASTQ asociados a cada muestra
sample_fastq_sizes <- function(dir_path, samples, read_type = "pe") {
  if (length(samples) == 0 || !dir.exists(dir_path)) return(setNames(numeric(0), samples))
  sapply(samples, function(s) {
    suffixes <- if (identical(read_type, "se")) {
      c(".fastq.gz", ".fastq", FASTQ_R1_SUFFIXES)
    } else {
      c(FASTQ_R1_SUFFIXES, FASTQ_R2_SUFFIXES)
    }
    candidates <- sample_fastq_paths(dir_path, s, suffixes)
    sum(file.info(candidates)$size, na.rm = TRUE)
  }, USE.NAMES = TRUE)
}

#' Devuelve los nombres con caracteres problematicos para shell
bad_sample_chars <- function(names) {
  grep("[^A-Za-z0-9_.-]", names, value = TRUE)
}
