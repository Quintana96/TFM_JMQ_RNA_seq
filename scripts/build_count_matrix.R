#!/usr/bin/env Rscript
#' build_count_matrix.R
#' Escribe la matriz de conteos por GEN de una ejecución de salmon o kallisto.
#'
#' Uso:
#'   Rscript scripts/build_count_matrix.R <output_dir> <tool> <destino.tsv> [anotación]
#'
#' No reimplementa nada: reutiliza `load_counts_from_workflow()`, que es la misma
#' función que usa la aplicación. Así la matriz que queda en disco y la que la
#' app analiza no pueden divergir — construirlas por caminos distintos sería una
#' fuente garantizada de discrepancias.
#'
#' Vive en scripts/ y no en R/ a propósito: Shiny sourcea R/ entero al arrancar,
#' y un fichero con código de nivel superior se ejecutaria al abrir la app.
#'
#' Códigos de salida: 0 matriz escrita; 1 error de uso o de entrada; 2 no se ha
#' podido construir la matriz (el workflow continua, pero sin este artefacto).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Uso: build_count_matrix.R <output_dir> <tool> <destino.tsv> [anotación]\n",
      file = stderr())
  quit(status = 1)
}
output_dir <- args[[1]]
tool       <- args[[2]]
dest       <- args[[3]]
annotation <- if (length(args) >= 4 && nzchar(args[[4]])) args[[4]] else NULL

if (!dir.exists(output_dir)) {
  cat("No existe el directorio de la ejecución: ", output_dir, "\n",
      sep = "", file = stderr())
  quit(status = 1)
}
if (!tool %in% c("salmon", "kallisto")) {
  cat("Herramienta no soportada por este script: ", tool, "\n",
      sep = "", file = stderr())
  quit(status = 1)
}

# La raíz del proyecto es el directorio padre de scripts/.
app_root <- normalizePath(file.path(dirname(
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."),
  mustWork = FALSE)
if (!file.exists(file.path(app_root, "global.R"))) app_root <- normalizePath(".")

suppressWarnings(suppressMessages({
  sys.source(file.path(app_root, "global.R"), envir = globalenv())
  for (f in sort(list.files(file.path(app_root, "R"), pattern = "[.]R$",
                            full.names = TRUE))) {
    sys.source(f, envir = globalenv())
  }
}))

m <- tryCatch(
  load_counts_from_workflow(output_dir, tool, annotation_file = annotation),
  error = function(e) {
    cat("Error al construir la matriz: ", conditionMessage(e), "\n",
        sep = "", file = stderr())
    NULL
  }
)

if (is.null(m) || !nrow(m)) {
  cat("No se ha podido construir la matriz de conteos por gen.\n", file = stderr())
  quit(status = 2)
}

src <- attr(m, "counts_source")
out <- data.frame(gene_id = rownames(m), m, check.names = FALSE,
                  stringsAsFactors = FALSE)
utils::write.table(out, dest, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Matriz escrita: ", dest, " (", nrow(m), " genes x ", ncol(m), " muestras)\n",
    sep = "")
if (!is.null(src)) {
  cat("  Método: ", src$method %||% "—", "\n", sep = "")
  # Si el resumen a gen ha degradado a la via de respaldo hay que decirlo: la
  # matriz existe igualmente, pero no lleva los offsets de longitud efectiva.
  if (!isTRUE(src$ok)) cat("  AVISO: ", src$detail %||% "", "\n", sep = "")
}
quit(status = 0)
