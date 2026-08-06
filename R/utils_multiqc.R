#' utils_multiqc.R
#' Lectores de los reports MultiQC y resumenes agregados por run.

#' Carpeta multiqc_data dentro de una run
multiqc_dir <- function(out_dir) file.path(out_dir, "multiqc_data")

#' Lecturas individuales de los TSV de MultiQC
general_stats        <- function(out_dir) read_tsv_safe(file.path(multiqc_dir(out_dir), "multiqc_general_stats.txt"))
fastqc_stats         <- function(out_dir) read_tsv_safe(file.path(multiqc_dir(out_dir), "multiqc_fastqc.txt"))
salmon_stats         <- function(out_dir) read_tsv_safe(file.path(multiqc_dir(out_dir), "multiqc_salmon.txt"))
kallisto_stats       <- function(out_dir) read_tsv_safe(file.path(multiqc_dir(out_dir), "multiqc_kallisto.txt"))
bowtie2_stats        <- function(out_dir) read_tsv_safe(file.path(multiqc_dir(out_dir), "multiqc_bowtie2.txt"))
featurecounts_stats  <- function(out_dir) read_tsv_safe(file.path(multiqc_dir(out_dir), "multiqc_featurecounts.txt"))

#' Redondea numericamente todas las columnas excepto las indicadas en `exclude`
round_numeric_columns <- function(df, exclude = "Sample", digits = 3) {
  for (nm in setdiff(names(df), exclude)) {
    x <- num_or_na(df[[nm]])
    if (any(!is.na(x))) df[[nm]] <- round(x, digits)
  }
  df
}

#' Resumen agregado de una run (estado, conteos, metricas clave) para Tab 3
summarise_result <- function(out_dir, params) {
  gs <- general_stats(out_dir)
  fq <- fastqc_stats(out_dir)
  tool <- params$tool %||% "desconocida"
  status <- status_from_log(out_dir)
  counts <- tryCatch(load_counts_from_workflow(out_dir, tool, annotation_file = annotation_file_for_run(out_dir)), error = function(e) NULL)

  mapped_col <- intersect(c("salmon-percent_mapped", "kallisto-percent_pseudoaligned", "bowtie2-overall_alignment_rate"), names(gs))[1] %||% ""
  mapped_vals <- if (nzchar(mapped_col)) num_or_na(gs[[mapped_col]]) else numeric(0)
  trim_vals <- if (!is.null(gs) && "fastp-pct_surviving" %in% names(gs)) num_or_na(gs[["fastp-pct_surviving"]]) else numeric(0)
  q30_vals <- if (!is.null(gs) && "fastp-after_filtering_q30_rate" %in% names(gs)) num_or_na(gs[["fastp-after_filtering_q30_rate"]]) else numeric(0)
  fastqc_fail <- if (!is.null(fq)) {
    qc_cols <- setdiff(names(fq), c("Sample", "Filename", "File type", "Encoding", "Total Sequences", "Total Bases",
                                    "Sequences flagged as poor quality", "Sequence length", "%GC",
                                    "total_deduplicated_percentage", "avg_sequence_length", "median_sequence_length"))
    sum(as.matrix(fq[, qc_cols, drop = FALSE]) == "fail", na.rm = TRUE)
  } else NA_integer_

  list(
    status = status,
    tool = tool,
    out_dir = out_dir,
    total_size = fmt_bytes(dir_size(out_dir)),
    n_files = nrow(file_table_for_dir(out_dir)),
    n_samples = if (!is.null(counts)) ncol(counts) else params$n_samples %||% "—",
    n_features = if (!is.null(counts)) nrow(counts) else "—",
    mean_mapped = if (length(mapped_vals)) mean(mapped_vals, na.rm = TRUE) else NA_real_,
    mean_trim_survival = if (length(trim_vals)) mean(trim_vals, na.rm = TRUE) else NA_real_,
    mean_q30 = if (length(q30_vals)) mean(q30_vals, na.rm = TRUE) else NA_real_,
    fastqc_fail = fastqc_fail,
    has_multiqc = file.exists(file.path(out_dir, "multiqc_report.html"))
  )
}

#' Tabla resumida de multiqc_general_stats.txt filtrando columnas vacias y duplicados R1/R2
result_general_table <- function(out_dir) {
  gs <- general_stats(out_dir)
  if (is.null(gs)) return(data.frame(Mensaje = "No se encontro multiqc_general_stats.txt."))
  keep <- intersect(c(
    "Sample",
    "salmon-percent_mapped", "salmon-num_mapped", "salmon-library_types",
    "kallisto-percent_pseudoaligned", "kallisto-n_processed", "kallisto-n_pseudoaligned",
    "bowtie2-overall_alignment_rate",
    "fastp-pct_surviving", "fastp-pct_adapter", "fastp-after_filtering_q30_rate",
    "fastp-after_filtering_gc_content",
    "fastqc-total_sequences", "fastqc-percent_gc", "fastqc-percent_duplicates", "fastqc-percent_fails"
  ), names(gs))
  df <- gs[, keep, drop = FALSE]

  # Filtrar filas (muestras) que terminen en _1, _2, _R1 o _R2
  if ("Sample" %in% names(df)) {
    df <- df[!grepl("(_1|_2|_R1|_R2)$", df$Sample), , drop = FALSE]
  }

  # Filtrar columnas donde todos los valores sean NA o vacios
  cols_to_keep <- intersect("Sample", names(df))
  for (nm in setdiff(names(df), "Sample")) {
    x <- num_or_na(df[[nm]])
    if (any(!is.na(x))) {
      cols_to_keep <- c(cols_to_keep, nm)
    }
  }
  df <- df[, cols_to_keep, drop = FALSE]
  round_numeric_columns(df)
}
