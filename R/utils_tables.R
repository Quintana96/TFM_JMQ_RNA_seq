#' utils_tables.R
#' DataTables y tablas auxiliares (FastQC, alineamiento, conteos, artefactos, log tail).

#' Diccionario de DataTables en español.
#'
#' Se define en local y no por URL: la versión oficial se descarga de un CDN, y
#' una app que puede correr sin conexión no debe quedarse con los textos en
#' ingles (o peor, sin tabla) por no alcanzarlo. Va aquí, en el único wrapper de
#' tablas, para que traducirlas sea un solo punto de cambio.
DT_ES <- list(
  processing = "Procesando...",
  search = "Buscar:",
  lengthMenu = "Mostrar _MENU_ filas",
  info = "Mostrando _START_ a _END_ de _TOTAL_ filas",
  infoEmpty = "Sin filas que mostrar",
  infoFiltered = "(filtrado de _MAX_ filas totales)",
  loadingRecords = "Cargando...",
  zeroRecords = "No se han encontrado resultados",
  emptyTable = "Sin datos disponibles en la tabla",
  paginate = list(first = "Primera", previous = "Anterior",
                  `next` = "Siguiente", last = "Última"),
  aria = list(sortAscending = ": activar para ordenar de forma ascendente",
              sortDescending = ": activar para ordenar de forma descendente")
)

#' Wrapper estandar de DT::datatable con scrollX y dom ftip.
#'
#' Formatea ademas los numeros en castellano. Va aqui, en el unico envoltorio de
#' tablas, por el mismo motivo que el diccionario de arriba: que la convencion
#' numerica de la aplicacion sea un solo punto de cambio y no pueda divergir
#' entre pestanas.
dt_table <- function(data, page_length = 10, filter = "none", decimales = 3) {
  dt <- datatable(
    data,
    filter = filter,
    options = list(pageLength = page_length, scrollX = TRUE, dom = "ftip",
                   language = DT_ES),
    rownames = FALSE
  )
  formatear_numeros_es(dt, data, decimales)
}

#' Aplica a un DataTable la convencion numerica castellana.
#'
#' Coma decimal y punto de millares: "12.345,679". Es la del documento del TFM y
#' la que ya usan las tablas del harness, y las capturas de la aplicacion van
#' dentro de ese documento.
#'
#' Tres tratamientos, segun lo que contenga la columna:
#'
#'   - **P-valores**: tres cifras SIGNIFICATIVAS, no decimales.
#'     `round(4.42e-46, 3)` es 0, y ese 4,42e-46 es el p ajustado de CRISPLD2 en
#'     GSE52778, uno de los genes de control de la memoria. Ese cero no seria un
#'     numero redondeado sino el dato perdido, y dejaria un gen significativo
#'     indistinguible de uno que no lo es.
#'   - **Enteros**: sin decimales, pero con separador de millares. Una libreria
#'     de 12.345.678 lecturas es ilegible escrita del tiron, y escribirle tres
#'     decimales es peor.
#'   - **El resto**: `decimales` decimales.
#'
#' Se formatea en el lado del cliente, con los formateadores de DT, y no
#' redondeando el data.frame. Los formateadores solo actuan cuando el tipo es
#' 'display', asi que la ordenacion y el filtro por rango siguen usando el valor
#' completo: redondeando el dato, dos genes con padj 0,0499 y 0,0501 se
#' ordenarian como iguales. Y las descargas entregan el data.frame intacto, de
#' modo que lo que cambia es la LECTURA y no el resultado.
formatear_numeros_es <- function(dt, data, decimales = 3) {
  # `datatable()` tambien acepta matrices, donde recorrer columnas como lista no
  # funciona. No formatear es preferible a romper la tabla.
  if (!is.data.frame(data) || !ncol(data)) return(dt)
  num <- names(data)[vapply(data, is.numeric, logical(1))]
  if (!length(num)) return(dt)
  p_val <- intersect(num, COLS_P_VALOR)
  resto <- setdiff(num, p_val)
  enteros <- resto[vapply(data[resto], function(x)
    is.integer(x) || all(is.na(x) | x == round(x)), logical(1))]
  continuas <- setdiff(resto, enteros)

  if (length(continuas))
    dt <- formatRound(dt, continuas, digits = decimales, mark = ".", dec.mark = ",")
  if (length(enteros))
    dt <- formatRound(dt, enteros, digits = 0, mark = ".", dec.mark = ",")
  if (length(p_val))
    dt <- formatSignif(dt, p_val, digits = decimales, mark = ".", dec.mark = ",")
  dt
}

#' Tabla FastQC filtrada a columnas relevantes
fastqc_table <- function(out_dir) {
  fq <- fastqc_stats(out_dir)
  if (is.null(fq)) return(data.frame(Mensaje = "No se encontro multiqc_fastqc.txt."))
  qc_cols <- intersect(c(
    "basic_statistics", "per_base_sequence_quality", "per_sequence_quality_scores",
    "per_base_sequence_content", "per_sequence_gc_content", "per_base_n_content",
    "sequence_length_distribution", "sequence_duplication_levels",
    "overrepresented_sequences", "adapter_content"
  ), names(fq))
  keep <- intersect(c("Sample", "Filename", "Total Sequences", "%GC", "avg_sequence_length", qc_cols), names(fq))
  fq[, keep, drop = FALSE]
}

#' Tabla de metricas de alineamiento o cuantificación según la tool
alignment_table <- function(out_dir, tool) {
  if (identical(tool, "salmon")) {
    df <- salmon_stats(out_dir)
    if (is.null(df)) return(message_df("No se encontro multiqc_salmon.txt."))
    keep <- intersect(c("Sample", "salmon_version", "library_types", "num_processed", "num_mapped",
                        "percent_mapped", "frag_length_mean", "frag_length_sd"), names(df))
  } else if (identical(tool, "kallisto")) {
    df <- kallisto_stats(out_dir)
    if (is.null(df)) return(message_df("No se encontro multiqc_kallisto.txt."))
    keep <- intersect(c("Sample", "kallisto_version", "n_targets", "n_bootstraps",
                        "fragments", "pseudoaligned", "percent_pseudoaligned",
                        "n_processed", "n_pseudoaligned", "n_unique",
                        "fragments_processed", "pseudoaligned_reads",
                        "frag_length_mean", "frag_length_sd"), names(df))
    if (!length(keep)) keep <- names(df)
  } else if (identical(tool, "bowtie2")) {
    df <- bowtie2_stats(out_dir)
    if (is.null(df)) return(message_df("No se encontro multiqc_bowtie2.txt."))
    keep <- names(df)
  } else {
    df <- featurecounts_stats(out_dir)
    if (is.null(df)) return(message_df("No se encontraron metricas de alineamiento/cuantificación."))
    keep <- names(df)
  }
  df <- df[, keep, drop = FALSE]
  round_numeric_columns(df)
}

#' Tablas resumen de conteos: librerias por muestra y top 30 genes
counts_tables <- function(out_dir, tool) {
  counts <- tryCatch(load_counts_from_workflow(out_dir, tool, annotation_file = annotation_file_for_run(out_dir)), error = function(e) NULL)
  if (is.null(counts) || !length(counts)) {
    return(list(
      libs = data.frame(Mensaje = "No se pudo cargar una matriz de conteos útil."),
      top = data.frame(Mensaje = "Sin conteos disponibles.")
    ))
  }
  counts <- as.matrix(counts)
  storage.mode(counts) <- "numeric"
  libs <- data.frame(
    Muestra = colnames(counts),
    Lecturas_asignadas = round(colSums(counts, na.rm = TRUE)),
    Genes_detectados = colSums(counts > 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  total <- rowSums(counts, na.rm = TRUE)
  top <- head(data.frame(
    Gen = rownames(counts),
    Conteo_total = round(total),
    Media = round(rowMeans(counts, na.rm = TRUE), 2),
    stringsAsFactors = FALSE
  )[order(total, decreasing = TRUE), ], 30)
  list(libs = libs, top = top)
}

#' Lista de artefactos importantes (MultiQC, log, matriz, HTML fastp/fastqc)
important_artifacts <- function(out_dir) {
  candidates <- c(
    file.path(out_dir, "multiqc_report.html"),
    file.path(out_dir, "workflow_live.log"),
    file.path(out_dir, "04_counts", "count_matrix.tsv"),
    list.files(file.path(out_dir, "02_trimmed_reads"), pattern = "_fastp\\.html$", full.names = TRUE),
    list.files(file.path(out_dir, "01_quality"), pattern = "_fastqc\\.html$", full.names = TRUE)
  )
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) return(data.frame(Mensaje = "No se encontraron informes principales."))
  data.frame(
    Tipo = ifelse(grepl("multiqc_report", candidates), "MultiQC",
           ifelse(grepl("workflow_live\\.log", candidates), "Log",
           ifelse(grepl("count_matrix", candidates), "Matriz de conteos",
           ifelse(grepl("fastp\\.html", candidates), "fastp HTML", "FastQC HTML")))),
    Archivo = basename(candidates),
    Ruta = candidates,
    Tamaño = sapply(file.info(candidates)$size, fmt_bytes),
    stringsAsFactors = FALSE
  )
}

#' Devuelve las últimas n líneas del workflow_live.log de una run
log_tail_text <- function(out_dir, n = 140) {
  f <- file.path(out_dir, "workflow_live.log")
  if (!file.exists(f)) return("No se encontro workflow_live.log.")
  lines <- strsplit(read_tail_text(f), "\n", fixed = TRUE)[[1]]
  paste(tail(lines, n), collapse = "\n")
}
