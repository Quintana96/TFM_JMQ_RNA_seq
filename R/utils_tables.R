#' utils_tables.R
#' DataTables y tablas auxiliares (FastQC, alineamiento, conteos, artefactos, log tail).

#' Wrapper estandar de DT::datatable con scrollX y dom ftip
dt_table <- function(data, page_length = 10, filter = "none") {
  datatable(
    data,
    filter = filter,
    options = list(pageLength = page_length, scrollX = TRUE, dom = "ftip"),
    rownames = FALSE
  )
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

#' Tabla de metricas de alineamiento o cuantificacion segun la tool
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
    if (is.null(df)) return(message_df("No se encontraron metricas de alineamiento/cuantificacion."))
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
      libs = data.frame(Mensaje = "No se pudo cargar una matriz de conteos util."),
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
    Tamano = sapply(file.info(candidates)$size, fmt_bytes),
    stringsAsFactors = FALSE
  )
}

#' Devuelve las ultimas n lineas del workflow_live.log de una run
log_tail_text <- function(out_dir, n = 140) {
  f <- file.path(out_dir, "workflow_live.log")
  if (!file.exists(f)) return("No se encontro workflow_live.log.")
  lines <- strsplit(read_tail_text(f), "\n", fixed = TRUE)[[1]]
  paste(tail(lines, n), collapse = "\n")
}
