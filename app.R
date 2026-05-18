# ============================================================
# app.R  —  RNA-seq Workflow Runner
# ============================================================

library(shiny)
library(bslib)
library(shinyjs)
library(DT)
library(shinyFiles)

# Allow large uploads (approx 10 GB)
options(shiny.maxRequestSize = 10 * 1024^3)
options(sass.cache = file.path(tempdir(), "rnaseq_shiny_sass_cache"))

# ── Paquetes opcionales — deteccion en tiempo de carga ───────────────────────
pkg_ok <- function(p) requireNamespace(p, quietly = TRUE)
HAS_PROCESSX   <- pkg_ok("processx")
HAS_TXIMPORT   <- pkg_ok("tximport")


# ── Operador nulo-coalescente ─────────────────────────────────────────────────
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  HELPERS                                                                 ║
# ╚══════════════════════════════════════════════════════════════════════════╝

# ── Deteccion de muestras FASTQ ───────────────────────────────
FASTQ_R1_SUFFIXES <- c("_1.fastq.gz", "_R1.fastq.gz", "_1.fastq", "_R1.fastq")
FASTQ_R2_SUFFIXES <- c("_2.fastq.gz", "_R2.fastq.gz", "_2.fastq", "_R2.fastq")
FASTQ_R1_PATTERN <- "(_1\\.fastq\\.gz$|_R1\\.fastq\\.gz$|_1\\.fastq$|_R1\\.fastq$)"
FASTQ_R2_PATTERN <- "(_2\\.fastq\\.gz$|_R2\\.fastq\\.gz$|_2\\.fastq$|_R2\\.fastq$)"
FASTQ_ANY_PATTERN <- "(\\.fastq\\.gz$|\\.fastq$)"

sample_fastq_paths <- function(dir_path, samples, suffixes) {
  file.path(dir_path, as.vector(outer(samples, suffixes, paste0)))
}

sample_name_from_r1 <- function(files) {
  sub(FASTQ_R1_PATTERN, "", basename(files), ignore.case = TRUE)
}

sample_name_from_fastq <- function(files) {
  x <- basename(files)
  x <- sub(FASTQ_R1_PATTERN, "", x, ignore.case = TRUE)
  x <- sub(FASTQ_ANY_PATTERN, "", x, ignore.case = TRUE)
  x
}

read_type_label <- function(read_type) {
  if (identical(read_type, "se")) "Single-end" else "Paired-end"
}

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

missing_r2 <- function(dir_path, samples) {
  if (length(samples) == 0) return(logical(0))
  vapply(samples, function(s) {
    !any(file.exists(sample_fastq_paths(dir_path, s, FASTQ_R2_SUFFIXES)))
  }, logical(1))
}

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

fmt_bytes <- function(b) {
  b <- as.numeric(b)
  if (is.na(b))    return("—")
  if (b >= 1e9)    return(sprintf("%.1f GB", b / 1e9))
  if (b >= 1e6)    return(sprintf("%.1f MB", b / 1e6))
  if (b >= 1e3)    return(sprintf("%.1f KB", b / 1e3))
  paste0(b, " B")
}

ts_log <- function(msg) paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)

terminal_text <- function(x) {
  if (length(x) == 0 || is.null(x)) return("")
  x <- paste(x, collapse = "\n")
  x <- gsub("\r\n?", "\n", x)
  Encoding(x) <- "UTF-8"
  x
}

trim_log_text <- function(x, max_chars = 250000L) {
  if (nchar(x, type = "chars", allowNA = FALSE) <= max_chars) return(x)
  paste0(
    "[log truncado en la app; consulta workflow_live.log para el log completo]\n",
    substring(x, nchar(x) - max_chars + 1L)
  )
}

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

#' Formatea segundos en cadena legible (3s, 2m 15s, 1h 4m)
fmt_elapsed <- function(secs) {
  s <- round(as.numeric(secs))
  if (s < 60)   return(sprintf("%ds", s))
  if (s < 3600) return(sprintf("%dm %ds", s %/% 60, s %% 60))
  sprintf("%dh %dm", s %/% 3600, (s %% 3600) %/% 60)
}

#' Valida que el nombre de muestra no tenga caracteres problematicos
bad_sample_chars <- function(names) {
  grep("[^A-Za-z0-9_.-]", names, value = TRUE)
}

#' Carga la matriz de conteos generada por el workflow
load_quant_counts <- function(qfiles, tool) {
  mats <- lapply(qfiles, function(f) {
    x <- tryCatch(read.delim(f, check.names = FALSE), error = function(e) NULL)
    if (is.null(x)) return(NULL)

    if (identical(tool, "salmon") && all(c("Name", "NumReads") %in% names(x))) {
      vals <- x$NumReads
      names(vals) <- x$Name
      return(vals)
    }
    if (identical(tool, "kallisto") && all(c("target_id", "est_counts") %in% names(x))) {
      vals <- x$est_counts
      names(vals) <- x$target_id
      return(vals)
    }
    NULL
  })
  valid <- !vapply(mats, is.null, logical(1))
  if (!any(valid)) return(NULL)

  mats <- mats[valid]
  features <- Reduce(union, lapply(mats, names))
  out <- matrix(0, nrow = length(features), ncol = length(mats),
                dimnames = list(features, names(qfiles)[valid]))
  for (nm in names(mats)) {
    vals <- mats[[nm]]
    out[names(vals), nm] <- vals
  }
  round(out)
}

load_count_matrix_tsv <- function(path) {
  if (!file.exists(path)) return(NULL)
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (!length(lines)) return(NULL)

  header_idx <- grep("^Geneid\\t|^geneid\\t|^Gene\\t|^gene\\t", lines, ignore.case = TRUE)[1]
  if (is.na(header_idx)) header_idx <- 1L

  df <- tryCatch(
    read.delim(
      text = paste(lines[header_idx:length(lines)], collapse = "\n"),
      header = TRUE, row.names = 1, check.names = FALSE, comment.char = ""
    ),
    error = function(e) NULL
  )
  if (is.null(df)) return(NULL)

  names(df) <- sub("\\.bam$", "", basename(names(df)))
  df
}

load_counts_from_workflow <- function(output_dir, tool) {
  if (tool == "bowtie2") {
    f <- file.path(output_dir, "04_counts", "count_matrix.tsv")
    load_count_matrix_tsv(f)
  } else if (HAS_TXIMPORT) {
    aln_dir <- file.path(output_dir, "03_alignments", tool)
    if (!dir.exists(aln_dir)) return(NULL)
    sdirs <- list.dirs(aln_dir, recursive = FALSE, full.names = TRUE)
    if (length(sdirs) == 0) return(NULL)
    qfiles <- if (tool == "salmon")
      file.path(sdirs, "quant.sf")
    else
      file.path(sdirs, "abundance.h5")
    if (!all(file.exists(qfiles)))
      qfiles <- file.path(sdirs, "abundance.tsv")
    valid <- file.exists(qfiles)
    if (!any(valid)) return(NULL)
    qfiles <- qfiles[valid]
    names(qfiles) <- basename(sdirs[valid])
    txi <- tryCatch(
      tximport::tximport(qfiles, type = tool, ignoreTxVersion = TRUE),
      error = function(e) NULL
    )
    if (is.null(txi)) return(load_quant_counts(qfiles, tool))
    round(txi$counts)
  } else if (tool %in% c("salmon", "kallisto")) {
    aln_dir <- file.path(output_dir, "03_alignments", tool)
    if (!dir.exists(aln_dir)) return(NULL)
    sdirs <- list.dirs(aln_dir, recursive = FALSE, full.names = TRUE)
    if (!length(sdirs)) return(NULL)
    qfiles <- if (tool == "salmon") file.path(sdirs, "quant.sf") else file.path(sdirs, "abundance.tsv")
    valid <- file.exists(qfiles)
    if (!any(valid)) return(NULL)
    qfiles <- qfiles[valid]
    names(qfiles) <- basename(sdirs[valid])
    load_quant_counts(qfiles, tool)
  } else NULL
}

#' Copia un archivo subido por Shiny a una ruta persistente
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
  if (!file.exists(dest_path) || file.info(upload$datapath)$size != file.info(dest_path)$size) {
    if (!file.copy(upload$datapath, dest_path, overwrite = TRUE)) {
      stop(sprintf("No se pudo copiar %s a %s", upload$datapath, dest_path))
    }
  }
  dest_path
}

## Reproducibilidad: eliminado según petición del usuario

outputs_base_dir <- function() file.path(getwd(), "outputs")

safe_run_label <- function(analysis_type, tool, time = Sys.time()) {
  analysis <- if (identical(analysis_type, "alignment")) "alineamiento" else "pseudoalineamiento"
  paste(format(time, "%Y%m%d_%H%M%S"), analysis, tool, sep = "_")
}

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

list_result_dirs <- function(base_dir) {
  if (!dir.exists(base_dir)) return(character(0))
  dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[file.info(dirs)$isdir %||% FALSE]
  if (!length(dirs)) return(character(0))
  mt <- file.info(dirs)$mtime
  dirs[order(mt, decreasing = TRUE, na.last = TRUE)]
}

result_choices <- function(base_dir) {
  dirs <- list_result_dirs(base_dir)
  if (!length(dirs)) return(character(0))
  stats <- file.info(dirs)
  labels <- paste0(basename(dirs), "  (", format(stats$mtime, "%Y-%m-%d %H:%M"), ")")
  setNames(dirs, labels)
}

file_table_for_files <- function(out_dir, files) {
  if (!length(files)) return(data.frame(Archivo=character(), Tamano=character(),
                                       stringsAsFactors=FALSE, check.names=FALSE))
  sz <- file.info(file.path(out_dir, files))$size
  data.frame(Archivo = files, Tamano = sapply(sz, fmt_bytes),
             stringsAsFactors = FALSE, check.names = FALSE)
}

file_table_for_dir <- function(out_dir) {
  if (!nzchar(out_dir) || !dir.exists(out_dir)) {
    return(data.frame(Archivo="Directorio no existe.", `Tamano`="—",
                      stringsAsFactors=FALSE, check.names=FALSE))
  }
  files <- list.files(out_dir, recursive = TRUE, full.names = FALSE)
  if (!length(files)) {
    return(data.frame(Archivo="Sin archivos.", `Tamano`="—",
                      stringsAsFactors=FALSE, check.names=FALSE))
  }
  file_table_for_files(out_dir, files)
}

infer_result_params <- function(out_dir, workflow_path) {
  tool <- if (dir.exists(file.path(out_dir, "03_alignments", "bowtie2"))) {
    "bowtie2"
  } else if (dir.exists(file.path(out_dir, "03_alignments", "salmon"))) {
    "salmon"
  } else if (dir.exists(file.path(out_dir, "03_alignments", "kallisto"))) {
    "kallisto"
  } else {
    "desconocida"
  }
  analysis <- if (identical(tool, "bowtie2")) "alignment" else "pseudo"
  counts <- tryCatch(load_counts_from_workflow(out_dir, tool), error = function(e) NULL)
  stat <- file.info(out_dir)
  list(
    analysis_type = analysis,
    tool = tool,
    input_dir = "—",
    output_dir = out_dir,
    genome_file = "—",
    annotation_file = "—",
    n_samples = if (!is.null(counts)) ncol(counts) else "—",
    read_type = "Paired-end",
    started_at = stat$mtime %||% Sys.time(),
    r_version = paste(R.version$major, R.version$minor, sep=".")
  )
}

read_tsv_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = ""),
    error = function(e) NULL
  )
}

num_or_na <- function(x) suppressWarnings(as.numeric(x))

dir_size <- function(path) {
  files <- list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  if (!length(files)) return(0)
  sum(file.info(files)$size, na.rm = TRUE)
}

status_from_log <- function(out_dir) {
  log_file <- file.path(out_dir, "workflow_live.log")
  if (!file.exists(log_file)) return("sin log")
  txt <- read_tail_text(log_file, max_bytes = 512000L)
  if (grepl("Analysis completed successfully|Analisis finalizado OK", txt, ignore.case = TRUE))
    return("completado")
  if (grepl("ERROR|Error \\(codigo|fallo en la linea", txt, ignore.case = TRUE))
    return("error")
  "incompleto"
}

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

multiqc_dir <- function(out_dir) file.path(out_dir, "multiqc_data")

general_stats <- function(out_dir) {
  read_tsv_safe(file.path(multiqc_dir(out_dir), "multiqc_general_stats.txt"))
}

fastqc_stats <- function(out_dir) {
  read_tsv_safe(file.path(multiqc_dir(out_dir), "multiqc_fastqc.txt"))
}

salmon_stats <- function(out_dir) {
  read_tsv_safe(file.path(multiqc_dir(out_dir), "multiqc_salmon.txt"))
}

kallisto_stats <- function(out_dir) {
  read_tsv_safe(file.path(multiqc_dir(out_dir), "multiqc_kallisto.txt"))
}

bowtie2_stats <- function(out_dir) {
  read_tsv_safe(file.path(multiqc_dir(out_dir), "multiqc_bowtie2.txt"))
}

featurecounts_stats <- function(out_dir) {
  read_tsv_safe(file.path(multiqc_dir(out_dir), "multiqc_featurecounts.txt"))
}

summarise_result <- function(out_dir, params) {
  gs <- general_stats(out_dir)
  fq <- fastqc_stats(out_dir)
  tool <- params$tool %||% "desconocida"
  status <- status_from_log(out_dir)
  counts <- tryCatch(load_counts_from_workflow(out_dir, tool), error = function(e) NULL)

  mapped_col <- intersect(c("salmon-percent_mapped", "bowtie2-overall_alignment_rate"), names(gs))[1] %||% ""
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

pct_label <- function(x, digits = 1) {
  if (is.na(x) || !is.finite(x)) return("—")
  sprintf(paste0("%.", digits, "f%%"), x)
}

result_general_table <- function(out_dir) {
  gs <- general_stats(out_dir)
  if (is.null(gs)) return(data.frame(Mensaje = "No se encontro multiqc_general_stats.txt."))
  keep <- intersect(c(
    "Sample",
    "salmon-percent_mapped", "salmon-num_mapped", "salmon-library_types",
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
  
  # Filtrar columnas donde todos los valores sean NA o vacíos
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

round_numeric_columns <- function(df, exclude = "Sample", digits = 3) {
  for (nm in setdiff(names(df), exclude)) {
    x <- num_or_na(df[[nm]])
    if (any(!is.na(x))) df[[nm]] <- round(x, digits)
  }
  df
}

dt_table <- function(data, page_length = 10, filter = "none") {
  datatable(
    data,
    filter = filter,
    options = list(pageLength = page_length, scrollX = TRUE, dom = "ftip"),
    rownames = FALSE
  )
}

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

alignment_table <- function(out_dir, tool) {
  if (identical(tool, "salmon")) {
    df <- salmon_stats(out_dir)
    if (is.null(df)) return(data.frame(Mensaje = "No se encontro multiqc_salmon.txt."))
    keep <- intersect(c("Sample", "salmon_version", "library_types", "num_processed", "num_mapped",
                        "percent_mapped", "frag_length_mean", "frag_length_sd"), names(df))
  } else if (identical(tool, "bowtie2")) {
    df <- bowtie2_stats(out_dir)
    if (is.null(df)) return(data.frame(Mensaje = "No se encontro multiqc_bowtie2.txt."))
    keep <- names(df)
  } else {
    df <- featurecounts_stats(out_dir)
    if (is.null(df)) return(data.frame(Mensaje = "No se encontraron metricas de alineamiento/cuantificacion."))
    keep <- names(df)
  }
  df <- df[, keep, drop = FALSE]
  round_numeric_columns(df)
}

counts_tables <- function(out_dir, tool) {
  counts <- tryCatch(load_counts_from_workflow(out_dir, tool), error = function(e) NULL)
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

log_tail_text <- function(out_dir, n = 140) {
  f <- file.path(out_dir, "workflow_live.log")
  if (!file.exists(f)) return("No se encontro workflow_live.log.")
  lines <- strsplit(read_tail_text(f), "\n", fixed = TRUE)[[1]]
  paste(tail(lines, n), collapse = "\n")
}


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  QC adicional: alineamiento y pseudoalineamiento                         ║
# ╚══════════════════════════════════════════════════════════════════════════╝

# Umbrales configurables para las alertas de QC adicional.
qc_thresholds <- list(
  mapping_rate_warning = 0.70,
  mapping_rate_error = 0.50,
  assigned_rate_warning = 0.60,
  multimapping_rate_warning = 0.20,
  pseudoalignment_rate_warning = 0.70,
  low_reads_fraction_warning = 0.50,
  low_detected_fraction_warning = 0.50,
  near_zero_tpm_fraction_warning = 0.80,
  tpm_distribution_shift_warning = 2
)

rate_fraction <- function(x) {
  x <- num_or_na(x)
  ifelse(!is.na(x) & abs(x) > 1, x / 100, x)
}

find_metric_col <- function(df, patterns) {
  if (is.null(df) || !length(names(df))) return("")
  for (pat in patterns) {
    hit <- grep(pat, names(df), ignore.case = TRUE, value = TRUE)
    if (length(hit)) return(hit[1])
  }
  ""
}

message_df <- function(msg) {
  data.frame(Mensaje = msg, stringsAsFactors = FALSE, check.names = FALSE)
}

has_real_rows <- function(df) {
  !is.null(df) && nrow(df) > 0 && !"Mensaje" %in% names(df)
}

plotly_message <- function(msg) {
  plotly::plot_ly(
    x = 0.5, y = 0.5, type = "scatter", mode = "text",
    text = msg, hoverinfo = "none", textposition = "middle center"
  ) |>
    plotly::layout(
      xaxis = list(visible = FALSE, zeroline = FALSE),
      yaxis = list(visible = FALSE, zeroline = FALSE),
      margin = list(l = 20, r = 20, b = 20, t = 20)
    )
}

sample_column <- function(df) {
  hit <- intersect(c("Sample", "sample", "Muestra", "sample_id"), names(df))
  if (!length(hit) || is.na(hit[1])) "" else hit[1]
}

clean_result_samples <- function(df) {
  sc <- sample_column(df)
  if (!nzchar(sc)) return(df)
  df[!grepl("(_1|_2|_R1|_R2)$", df[[sc]]), , drop = FALSE]
}

featurecounts_assignment_summary <- function(out_dir) {
  fc <- featurecounts_stats(out_dir)
  if (is.null(fc)) return(message_df("No se encontraron metricas de asignacion genica para este analisis."))
  sc <- sample_column(fc)
  if (!nzchar(sc)) return(message_df("No se encontraron metricas de asignacion genica para este analisis."))

  assigned_col <- find_metric_col(fc, c("^Assigned$", "featurecounts.*assigned$"))
  unassigned_cols <- grep("unassigned", names(fc), ignore.case = TRUE, value = TRUE)
  if (!nzchar(assigned_col) && !length(unassigned_cols)) {
    return(message_df("No se encontraron metricas de asignacion genica para este analisis."))
  }

  assigned <- if (nzchar(assigned_col)) num_or_na(fc[[assigned_col]]) else rep(NA_real_, nrow(fc))
  unassigned <- if (length(unassigned_cols)) {
    rowSums(as.data.frame(lapply(fc[, unassigned_cols, drop = FALSE], num_or_na)), na.rm = TRUE)
  } else rep(NA_real_, nrow(fc))
  total <- assigned + unassigned
  data.frame(
    sample_id = fc[[sc]],
    assigned_reads = assigned,
    unassigned_reads = unassigned,
    assigned_rate = ifelse(total > 0, assigned / total, NA_real_),
    stringsAsFactors = FALSE
  )
}

align_qc_summary <- function(out_dir) {
  gs <- clean_result_samples(general_stats(out_dir))
  bt <- bowtie2_stats(out_dir)
  fc_sum <- featurecounts_assignment_summary(out_dir)

  if (is.null(gs) && is.null(bt) && !has_real_rows(fc_sum)) {
    return(message_df("No se encontraron metricas de alineamiento clasico para este analisis."))
  }

  source <- if (!is.null(bt)) bt else gs
  sc <- sample_column(source)
  samples <- if (nzchar(sc)) source[[sc]] else character(0)

  mapping_col <- find_metric_col(source, c("overall_alignment_rate", "percent_mapped", "percent_aligned"))
  mapping_rate <- if (nzchar(mapping_col)) rate_fraction(source[[mapping_col]]) else rep(NA_real_, length(samples))

  total_col <- find_metric_col(source, c("total_reads", "fastqc-total_sequences", "total_sequences", "reads_processed"))
  total_reads <- if (nzchar(total_col)) num_or_na(source[[total_col]]) else rep(NA_real_, length(samples))

  multi_col <- find_metric_col(source, c("multi", "aligned.*>1", "aligned.*multi"))
  multimapping_rate <- if (nzchar(multi_col)) rate_fraction(source[[multi_col]]) else rep(NA_real_, length(samples))

  unique_reads <- ifelse(!is.na(total_reads) & !is.na(mapping_rate),
                         total_reads * pmax(mapping_rate - ifelse(is.na(multimapping_rate), 0, multimapping_rate), 0),
                         NA_real_)
  multimapped_reads <- ifelse(!is.na(total_reads) & !is.na(multimapping_rate),
                              total_reads * multimapping_rate, NA_real_)
  unmapped_reads <- ifelse(!is.na(total_reads) & !is.na(mapping_rate),
                           total_reads * pmax(1 - mapping_rate, 0), NA_real_)

  out <- data.frame(
    sample_id = samples,
    total_reads = total_reads,
    unique_reads = unique_reads,
    multimapped_reads = multimapped_reads,
    unmapped_reads = unmapped_reads,
    mapping_rate = mapping_rate,
    multimapping_rate = multimapping_rate,
    stringsAsFactors = FALSE
  )

  if (has_real_rows(fc_sum)) {
    out <- merge(out, fc_sum, by = "sample_id", all = TRUE)
  } else {
    out$assigned_reads <- NA_real_
    out$unassigned_reads <- NA_real_
    out$assigned_rate <- NA_real_
  }

  if (!nrow(out)) message_df("No se encontraron metricas de alineamiento clasico para este analisis.") else out
}

read_salmon_quant_table <- function(out_dir) {
  qfiles <- list.files(file.path(out_dir, "03_alignments", "salmon"),
                       pattern = "^quant\\.sf$", recursive = TRUE, full.names = TRUE)
  qfiles <- qfiles[grepl("/03_alignments/salmon/[^/]+/quant\\.sf$", qfiles)]
  if (!length(qfiles)) return(NULL)
  pieces <- lapply(qfiles, function(f) {
    x <- read_tsv_safe(f)
    if (is.null(x)) return(NULL)
    x$sample_id <- basename(dirname(f))
    x
  })
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (!length(pieces)) return(NULL)
  do.call(rbind, pieces)
}

read_kallisto_quant_table <- function(out_dir) {
  qfiles <- list.files(file.path(out_dir, "03_alignments", "kallisto"),
                       pattern = "^abundance\\.tsv$", recursive = TRUE, full.names = TRUE)
  qfiles <- qfiles[grepl("/03_alignments/kallisto/[^/]+/abundance\\.tsv$", qfiles)]
  if (!length(qfiles)) return(NULL)
  pieces <- lapply(qfiles, function(f) {
    x <- read_tsv_safe(f)
    if (is.null(x)) return(NULL)
    names(x) <- sub("^target_id$", "Name", names(x))
    names(x) <- sub("^length$", "Length", names(x))
    names(x) <- sub("^eff_length$", "EffectiveLength", names(x))
    names(x) <- sub("^tpm$", "TPM", names(x))
    names(x) <- sub("^est_counts$", "NumReads", names(x))
    x$sample_id <- basename(dirname(f))
    x
  })
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (!length(pieces)) return(NULL)
  do.call(rbind, pieces)
}

infer_rna_type <- function(x) {
  id <- tolower(as.character(x))
  type <- rep("unknown", length(id))
  type[grepl("trna|transfer[-_ ]?rna", id)] <- "tRNA"
  type[grepl("rrna|ribosomal[-_ ]?rna|16s|23s|5s", id)] <- "rRNA"
  type[grepl("tmrna|ssra", id)] <- "tmRNA"
  type[grepl("srna|small[-_ ]?rna", id)] <- "sRNA"
  type[grepl("snrna", id)] <- "snRNA"
  type[grepl("snorna", id)] <- "snoRNA"
  type[grepl("mirna|micro[-_ ]?rna", id)] <- "miRNA"
  type[grepl("lncrna|long[-_ ]?non[-_ ]?coding", id)] <- "lncRNA"
  type[grepl("ncrna|non[-_ ]?coding[-_ ]?rna", id) & type == "unknown"] <- "ncRNA"
  type[grepl("mrna|cds|protein[-_ ]?coding", id) & type == "unknown"] <- "mRNA"
  type
}

pseudo_qc_quant_table <- function(out_dir) {
  df <- read_salmon_quant_table(out_dir)
  if (is.null(df)) df <- read_kallisto_quant_table(out_dir)
  if (is.null(df)) return(message_df("No se encontraron resultados de pseudoalineamiento para este analisis."))
  if ("Name" %in% names(df) && !"Type" %in% names(df)) {
    df$Type <- infer_rna_type(df$Name)
  }
  keep <- intersect(c("Name", "Type", "Length", "EffectiveLength", "TPM", "NumReads", "sample_id"), names(df))
  if (!length(keep)) return(message_df("No se encontraron resultados de pseudoalineamiento para este analisis."))
  df[, keep, drop = FALSE]
}

pseudo_qc_summary <- function(out_dir) {
  salmon <- salmon_stats(out_dir)
  kallisto <- kallisto_stats(out_dir)
  stats <- if (!is.null(salmon)) salmon else kallisto
  quant <- pseudo_qc_quant_table(out_dir)

  if (is.null(stats) && !has_real_rows(quant)) {
    return(message_df("No se encontraron resultados de pseudoalineamiento para este analisis."))
  }

  if (!is.null(stats)) {
    sc <- sample_column(stats)
    samples <- if (nzchar(sc)) stats[[sc]] else character(0)
    rate_col <- find_metric_col(stats, c("percent_mapped", "percent_aligned"))
    processed_col <- find_metric_col(stats, c("num_processed", "total_reads"))
    mapped_col <- find_metric_col(stats, c("num_mapped", "pseudoaligned_reads"))
    out <- data.frame(
      sample_id = samples,
      fragments_processed = if (nzchar(processed_col)) num_or_na(stats[[processed_col]]) else NA_real_,
      pseudoaligned_reads = if (nzchar(mapped_col)) num_or_na(stats[[mapped_col]]) else NA_real_,
      pseudoalignment_rate = if (nzchar(rate_col)) rate_fraction(stats[[rate_col]]) else NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    out <- data.frame(sample_id = unique(quant$sample_id), stringsAsFactors = FALSE)
    out$fragments_processed <- NA_real_
    out$pseudoaligned_reads <- NA_real_
    out$pseudoalignment_rate <- NA_real_
  }

  if (has_real_rows(quant)) {
    q <- quant
    q$TPM <- num_or_na(q$TPM)
    q$NumReads <- num_or_na(q$NumReads)
    detected <- aggregate(q$TPM > 0, list(sample_id = q$sample_id), sum, na.rm = TRUE)
    names(detected)[2] <- "transcripts_detected"
    near_zero <- aggregate(q$TPM <= 0.1, list(sample_id = q$sample_id), mean, na.rm = TRUE)
    names(near_zero)[2] <- "near_zero_tpm_fraction"
    tpm_median <- aggregate(q$TPM, list(sample_id = q$sample_id), median, na.rm = TRUE)
    names(tpm_median)[2] <- "median_tpm"
    out <- Reduce(function(x, y) merge(x, y, by = "sample_id", all = TRUE),
                  list(out, detected, near_zero, tpm_median))
  }

  counts <- tryCatch(load_counts_from_workflow(out_dir, "salmon"), error = function(e) NULL)
  if (is.null(counts)) counts <- tryCatch(load_counts_from_workflow(out_dir, "kallisto"), error = function(e) NULL)
  if (!is.null(counts) && length(counts)) {
    detected_features <- data.frame(
      sample_id = colnames(counts),
      genes_detected = colSums(as.matrix(counts) > 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    out <- merge(out, detected_features, by = "sample_id", all = TRUE)
  }

  if (!nrow(out)) message_df("No se encontraron resultados de pseudoalineamiento para este analisis.") else out
}

align_qc_alerts <- function(out_dir, thresholds = qc_thresholds) {
  x <- align_qc_summary(out_dir)
  if (!has_real_rows(x)) return(x)
  alerts <- list()
  med_reads <- median(x$total_reads, na.rm = TRUE)
  for (i in seq_len(nrow(x))) {
    s <- x$sample_id[i]
    add <- function(level, metric, value, message) {
      alerts[[length(alerts) + 1L]] <<- data.frame(
        sample_id = s, nivel = level, metrica = metric, valor = value,
        alerta = message, stringsAsFactors = FALSE
      )
    }
    if (!is.na(x$mapping_rate[i]) && x$mapping_rate[i] < thresholds$mapping_rate_error)
      add("error", "mapping_rate", round(x$mapping_rate[i], 3), "mapping_rate bajo.")
    else if (!is.na(x$mapping_rate[i]) && x$mapping_rate[i] < thresholds$mapping_rate_warning)
      add("warning", "mapping_rate", round(x$mapping_rate[i], 3), "mapping_rate por debajo del umbral recomendado.")
    if (!is.na(x$multimapping_rate[i]) && x$multimapping_rate[i] > thresholds$multimapping_rate_warning)
      add("warning", "multimapping_rate", round(x$multimapping_rate[i], 3), "multimapping_rate alto.")
    if (!is.na(x$assigned_rate[i]) && x$assigned_rate[i] < thresholds$assigned_rate_warning)
      add("warning", "assigned_rate", round(x$assigned_rate[i], 3), "assigned_rate bajo.")
    if (is.finite(med_reads) && med_reads > 0 && !is.na(x$total_reads[i]) &&
        x$total_reads[i] < med_reads * thresholds$low_reads_fraction_warning)
      add("warning", "total_reads", round(x$total_reads[i]), "Numero de lecturas muy inferior a la mediana del conjunto.")
    if (!is.na(x$total_reads[i]) && is.finite(med_reads) && med_reads > 0 &&
        (x$total_reads[i] < med_reads * 0.5 || x$total_reads[i] > med_reads * 1.5))
      add("info", "total_reads", round(x$total_reads[i]), "Muestra atipica respecto a la mediana del conjunto.")
  }
  if (!length(alerts)) return(message_df("No se detectaron alertas automaticas con los umbrales actuales."))
  do.call(rbind, alerts)
}

pseudo_qc_alerts <- function(out_dir, thresholds = qc_thresholds) {
  x <- pseudo_qc_summary(out_dir)
  if (!has_real_rows(x)) return(x)
  alerts <- list()
  med_processed <- median(x$fragments_processed, na.rm = TRUE)
  med_detected <- median(x$transcripts_detected, na.rm = TRUE)
  med_tpm <- median(x$median_tpm, na.rm = TRUE)
  mad_tpm <- mad(x$median_tpm, na.rm = TRUE)
  for (i in seq_len(nrow(x))) {
    s <- x$sample_id[i]
    add <- function(level, metric, value, message) {
      alerts[[length(alerts) + 1L]] <<- data.frame(
        sample_id = s, nivel = level, metrica = metric, valor = value,
        alerta = message, stringsAsFactors = FALSE
      )
    }
    if (!is.na(x$pseudoalignment_rate[i]) &&
        x$pseudoalignment_rate[i] < thresholds$pseudoalignment_rate_warning)
      add("warning", "pseudoalignment_rate", round(x$pseudoalignment_rate[i], 3), "pseudoalignment_rate bajo.")
    if (is.finite(med_detected) && med_detected > 0 && !is.na(x$transcripts_detected[i]) &&
        x$transcripts_detected[i] < med_detected * thresholds$low_detected_fraction_warning)
      add("warning", "transcripts_detected", x$transcripts_detected[i], "Numero de transcritos detectados anomalamente bajo.")
    if (is.finite(med_processed) && med_processed > 0 && !is.na(x$fragments_processed[i]) &&
        x$fragments_processed[i] < med_processed * thresholds$low_reads_fraction_warning)
      add("warning", "fragments_processed", round(x$fragments_processed[i]), "Bajo numero de fragmentos procesados.")
    if (!is.na(x$near_zero_tpm_fraction[i]) &&
        x$near_zero_tpm_fraction[i] > thresholds$near_zero_tpm_fraction_warning)
      add("info", "near_zero_tpm_fraction", round(x$near_zero_tpm_fraction[i], 3), "Exceso de transcritos con TPM cercano a cero.")
    if (is.finite(mad_tpm) && mad_tpm > 0 && !is.na(x$median_tpm[i]) &&
        abs(x$median_tpm[i] - med_tpm) / mad_tpm > thresholds$tpm_distribution_shift_warning)
      add("warning", "median_tpm", round(x$median_tpm[i], 3), "Distribucion de TPM muy diferente al resto.")
  }
  if (!length(alerts)) return(message_df("No se detectaron alertas automaticas con los umbrales actuales."))
  do.call(rbind, alerts)
}

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  UI                                                                      ║
# ╚══════════════════════════════════════════════════════════════════════════╝

app_theme <- bs_theme(
  version = 5,
  primary = "#BEE8C8",
  secondary = "#D8F1DD",
  success = "#A8DDB8",
  info = "#A8DADC",
  warning = "#F6D58A",
  danger = "#F4A6A6",
  base_font = font_collection("-apple-system", "BlinkMacSystemFont", "Segoe UI", "Roboto", "Arial", "sans-serif"),
  heading_font = font_collection("-apple-system", "BlinkMacSystemFont", "Segoe UI", "Roboto", "Arial", "sans-serif"),
  bg = "#F6FBF7",
  fg = "#20332A"
)

app_css <- HTML("
  :root {
    --pastel-green: #D8F1DD;
    --pastel-green-strong: #A8DDB8;
    --pastel-green-soft: #F1FAF3;
    --pastel-mint: #EAF8EE;
    --pastel-blue: #D7EEF1;
    --pastel-lavender: #E7DDF5;
    --pastel-peach: #F9DCC4;
    --pastel-rose: #FADDE1;
    --pastel-yellow: #FFF1C9;
    --pastel-coral: #FAD1CF;
    --ink: #20332A;
    --muted: #60756A;
    --line: #5F6F66;
    --thin-border: 0.75px solid #5F6F66;
  }

  body,
  .bslib-page-navbar {
    background:
      linear-gradient(180deg, #FAFEFB 0%, #F1FAF3 48%, #FBFEFC 100%);
    color: var(--ink);
  }

  .navbar {
    background: linear-gradient(90deg, #D7EEF1 0%, #E7DDF5 56%, #F9DCC4 100%) !important;
    border: var(--thin-border);
    border-width: 0 0 1px 0;
    box-shadow: 0 8px 24px rgba(93, 137, 108, 0.12);
  }

  .navbar-brand,
  .navbar .nav-link {
    color: var(--ink) !important;
    font-weight: 600;
    letter-spacing: 0;
  }

  .navbar .nav-link.active {
    background: rgba(255, 255, 255, 0.62);
    border: var(--thin-border);
    border-radius: 8px;
  }

  .card {
    border: var(--thin-border);
    border-radius: 8px;
    background: rgba(255, 255, 255, 0.92);
    box-shadow: 0 10px 24px rgba(94, 128, 105, 0.10);
  }

  .card-header {
    background: linear-gradient(90deg, #E7DDF5 0%, #F8FCF9 100%);
    border-bottom: var(--thin-border);
    color: var(--ink);
    font-weight: 600;
  }

  .card:nth-of-type(3n + 1) > .card-header {
    background: linear-gradient(90deg, #D7EEF1 0%, #F8FCF9 100%);
  }

  .card:nth-of-type(3n + 2) > .card-header {
    background: linear-gradient(90deg, #F9DCC4 0%, #F8FCF9 100%);
  }

  .card:nth-of-type(3n) > .card-header {
    background: linear-gradient(90deg, #E7DDF5 0%, #F8FCF9 100%);
  }

  .btn-primary,
  .btn-success {
    background: #A8DDB8;
    border: var(--thin-border);
    color: #173126;
    font-weight: 650;
  }

  .btn-primary:hover,
  .btn-success:hover {
    background: #95D2A8;
    border-color: #5F6F66;
    color: #12261D;
  }

  .btn-secondary,
  .btn-outline-secondary {
    background: #FFFFFF;
    border: var(--thin-border);
    color: #315342;
  }

  .btn-danger {
    background: #F4A6A6;
    border: var(--thin-border);
    color: #4A1F1F;
  }

  .btn-continue-blue {
    background: #A8DADC;
    border: var(--thin-border);
    color: #1F454B;
    font-weight: 600;
    box-shadow: 0 8px 18px rgba(80, 140, 150, 0.18);
  }

  .btn-continue-blue:hover,
  .btn-continue-blue:focus {
    background: #93CDD3;
    border-color: #5F6F66;
    color: #173A40;
  }

  .form-control,
  .form-select,
  .selectize-input {
    border: var(--thin-border) !important;
    border-radius: 8px;
  }

  .btn-default,
  .btn-file,
  .input-group .btn,
  button[id$='_btn'] {
    background: #D7EEF1;
    border: var(--thin-border);
    color: #244A50;
    font-weight: 500;
  }

  .btn-default:hover,
  .btn-file:hover,
  .input-group .btn:hover,
  button[id$='_btn']:hover {
    background: #C6E5EA;
    border-color: #5F6F66;
    color: #1F454B;
  }

  .input-group .form-control {
    background: rgba(255, 255, 255, 0.94);
    color: #20332A;
  }

  .nav-tabs .nav-link.active {
    background: #E7DDF5;
    border-color: #5F6F66 #5F6F66 #E7DDF5;
    color: #20332A;
    font-weight: 600;
  }

  .nav-tabs .nav-link {
    color: #4B3F61;
    font-weight: 500;
  }

  .alert-info {
    background: #D7EEF1;
    border: var(--thin-border);
    color: #244A50;
  }

  .alert-success {
    background: #EAF8EE;
    border: var(--thin-border);
    color: #244B34;
  }

  .alert-warning {
    background: #FFF1C9;
    border: var(--thin-border);
    color: #5C4A16;
  }

  .alert-danger {
    background: #FAD1CF;
    border: var(--thin-border);
    color: #5A2323;
  }

  .metric-card {
    padding: 12px;
    border-radius: 8px;
    min-height: 76px;
    display: flex;
    flex-direction: column;
    justify-content: center;
  }

  .metric-card:nth-child(1) {
    background: linear-gradient(180deg, #FFFFFF 0%, #E7DDF5 100%);
    border: var(--thin-border);
  }

  .metric-card:nth-child(2) {
    background: linear-gradient(180deg, #FFFFFF 0%, #D7EEF1 100%);
    border: var(--thin-border);
  }

  .metric-card:nth-child(3) {
    background: linear-gradient(180deg, #FFFFFF 0%, #F9DCC4 100%);
    border: var(--thin-border);
  }

  .metric-card:nth-child(4) {
    background: linear-gradient(180deg, #FFFFFF 0%, #FADDE1 100%);
    border: var(--thin-border);
  }

  .metric-card-value {
    font-size: 18px;
    font-weight: 650;
    color: #20332A;
  }

  .metric-card-label {
    font-size: 11px;
    color: #60756A;
    margin-top: 6px;
    text-transform: uppercase;
    letter-spacing: .04em;
  }

  table.dataTable thead th {
    background: #E7DDF5;
    color: #20332A;
    font-weight: 600;
    border-bottom: var(--thin-border) !important;
  }

  table.dataTable {
    border: var(--thin-border);
  }

  table.dataTable tbody tr:nth-child(odd) {
    background-color: #F8FCF9 !important;
  }

  table.dataTable tbody tr:nth-child(even) {
    background-color: #F7F1FC !important;
  }

  table.dataTable tbody tr:hover {
    background-color: #FFF1C9 !important;
  }

  pre,
  .shiny-text-output,
  code {
    border-radius: 6px;
    font-family: SFMono-Regular, Menlo, Consolas, 'Liberation Mono', monospace;
  }

  .card-title-download {
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 8px;
    width: 100%;
  }

  .header-download {
    min-width: 30px;
    min-height: 28px;
    padding: 2px 8px;
    line-height: 1;
  }

  .header-download .fa,
  .header-download svg {
    margin: 0;
  }
")

download_header <- function(title, output_id) {
  card_header(
    tags$div(
      class = "card-title-download",
      tags$span(title),
      downloadButton(
        output_id,
        label = NULL,
        icon = icon("download"),
        class = "btn-sm btn-outline-secondary header-download",
        title = paste("Descargar", title)
      )
    )
  )
}

ui <- page_navbar(
  title  = tags$span("RNA-seq Workflow Runner"),
  id     = "main_nav",
  theme  = app_theme,
  header = tagList(useShinyjs(), tags$style(app_css)),

  # ════════════════════════════════════════════════════════════════════════
  # TAB 1  —  CONFIGURACION
  # ════════════════════════════════════════════════════════════════════════
  nav_panel(
    title = "1 · Configuracion",
    value = "tab_config",

        layout_columns(
      col_widths = c(12),
      div(style = "display:flex; justify-content:flex-end; margin-bottom:12px;",
          actionButton("btn_to_processing",
                       tagList("Continuar al procesamiento", icon("arrow-right")),
                       class = "btn-continue-blue btn-lg")
      ),
      tags$div(style = "display:grid; grid-template-columns:repeat(3, minmax(240px, 1fr)); grid-auto-rows:minmax(240px, auto); gap:12px; align-items:stretch;",

        # Card 1: Modo de inicio
        card(
          card_header("Modo de inicio"),
          style = "min-height:240px; display:flex; flex-direction:column; justify-content:space-between; overflow:auto;",
          radioButtons(
            "start_mode", label = NULL,
            choices  = c("Ejecutar workflow completo" = "workflow",
                         "Análisis a partir de matriz de conteos"  = "load"),
            selected = "workflow", inline = TRUE
          ),
          conditionalPanel(
            condition = "input.start_mode === 'load'",
            layout_columns(
              col_widths = c(6, 6),
              div(
                tags$strong("Matriz de conteos (TSV/CSV)"),
                tags$br(),
                tags$small(class="text-muted",
                  "Genes como filas, muestras como columnas. Primera columna = ID de gen."),
                fileInput("upload_counts", label = NULL,
                          accept = c(".tsv", ".csv", ".txt"))
              ),
              div(
                tags$strong("Uso"),
                tags$p(class = "text-muted small",
                       "Carga una matriz externa o usa la pestaña Resultados para revisar ejecuciones previas guardadas automaticamente en outputs/.")
              )
            ),
            div(
              style = "text-align:right; margin-top:8px;",
              actionButton("btn_load_existing", tagList(icon("upload"), " Análisis a partir de matriz de conteos"),
                           class = "btn-success btn-sm")
            )
          )
        ),

        # Card 2: Tipo de análisis
        card(
          card_header("Tipo de analisis"),
          style = "min-height:240px; display:flex; flex-direction:column; justify-content:space-between;",
          radioButtons(
            "analysis_type", label = NULL,
            choices  = c("Alineamiento"      = "alignment",
                         "Pseudoalineamiento" = "pseudo"),
            selected = "alignment"
          ),
          radioButtons(
            "read_type", "Tipo de lectura",
            choices = c("Paired-end" = "pe", "Single-end" = "se"),
            selected = "pe",
            inline = TRUE
          ),
          conditionalPanel(
            "input.analysis_type === 'alignment'",
            tags$p(class = "text-muted small mb-0",
                   icon("circle-info"), " Bowtie2 + featureCounts")
          ),
          conditionalPanel(
            "input.analysis_type === 'pseudo'",
            selectInput("pseudo_tool", "Herramienta",
                        choices  = c("Salmon" = "salmon", "Kallisto" = "kallisto"),
                        selected = "salmon"),
            conditionalPanel(
              "input.pseudo_tool === 'kallisto' && input.read_type === 'se'",
              layout_columns(
                col_widths = c(6, 6),
                numericInput("fragment_length", "Longitud media fragmento", value = 200, min = 1, step = 1),
                numericInput("fragment_sd", "SD fragmento", value = 20, min = 1, step = 1)
              )
            )
          )
        ),

        # Card 3: Rutas y archivos
        card(
          card_header("Rutas y archivos"),
          style = "min-height:240px; display:flex; flex-direction:column; justify-content:space-between; overflow:auto;",
          tags$div(style="display:flex;gap:12px;align-items:flex-start;",
            tags$div(style="flex:1;",
              shinyDirButton("input_dir_btn", "Seleccionar directorio de FASTQs", "Seleccionar..."),
              tags$div(style="margin-top:6px;", textOutput("input_dir_path")),
              uiOutput("genome_label_ui"),
              conditionalPanel(
                "input.analysis_type === 'alignment'",
                fileInput("annotation_file_upload", "Archivo de anotacion GFF/GTF (requerido)",
                          accept = c(".gff", ".gtf", ".gff3", ".gz"), multiple = FALSE)
              ),
              conditionalPanel(
                "input.analysis_type === 'pseudo'",
                tagList(
                  fileInput("annotation_file_pseudo_upload", "Archivo de anotacion GFF/GTF (opcional)",
                            accept = c(".gff", ".gtf", ".gff3", ".gz"), multiple = FALSE),
                  tags$small(class = "text-muted", icon("circle-info"),
                             " Si se omite, se pasa /dev/null al script.")
                )
              ),
              tags$strong("Directorio de salida automatico"),
              tags$div(style="margin-top:6px;", textOutput("output_base_path"))
            )
          )
        ),

        # Card 4: Resumen rapido
        card(
          card_header("Resumen"),
          style = "min-height:240px; display:flex; flex-direction:column; justify-content:space-between;",
          tags$dl(class="row small mb-0",
            tags$dt(class="col-6","Entrada:"), tags$dd(class="col-6", tags$code(class="small", textOutput("input_dir_path", inline = TRUE))),
            tags$dt(class="col-6","Salida (base):"), tags$dd(class="col-6", tags$code(textOutput("output_base_path", inline = TRUE))),
            tags$dt(class="col-6","Workflow:"), tags$dd(class="col-6", tags$code(textOutput("workflow_path_text", inline = TRUE)))
          )
        ),

        # Card 5: Checklist
        card(
          card_header("Checklist"),
          style = "min-height:240px; display:flex; flex-direction:column; justify-content:center; overflow:auto;",
          tags$div(style="flex:1; display:flex; align-items:center; justify-content:center;",
            uiOutput("checklist_ui")
          )
        ),

        # Card 6: Muestras detectadas
        card(
          card_header("Muestras detectadas"),
          style = "min-height:240px; overflow:auto; max-height:420px;",
          uiOutput("sample_preview_ui")
        )
      )
    )
  ), # end tab_config


  # ════════════════════════════════════════════════════════════════════════
  # TAB 2  —  PROCESAMIENTO WORKFLOW
  # ════════════════════════════════════════════════════════════════════════
  nav_panel(
    title = "2 · Procesamiento",
    value = "tab_process",
    uiOutput("tab2_content")
  ),


  # ════════════════════════════════════════════════════════════════════════
  # TAB 3  —  RESULTADOS WORKFLOW
  # ════════════════════════════════════════════════════════════════════════
  nav_panel(
    title = "3 · Resultados",
    value = "tab_results",
    uiOutput("tab3_content")
  )

) # end page_navbar


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SERVER                                                                  ║
# ╚══════════════════════════════════════════════════════════════════════════╝

server <- function(input, output, session) {

  # ── Estado global ─────────────────────────────────────────────────────────
  workflow_path <- normalizePath("workflow.sh", mustWork = FALSE)
  outputs_dir <- outputs_base_dir()
  if (!dir.exists(outputs_dir)) dir.create(outputs_dir, recursive = TRUE, showWarnings = FALSE)
  addResourcePath("saved_outputs", normalizePath(outputs_dir, mustWork = TRUE))

  output$workflow_path_text <- renderText({ workflow_path })

  # Log de ejecucion
  log_text <- reactiveVal(paste0(ts_log("Workflow listo.\n")))

  # Archivos en el directorio de salida
  output_files_rv <- reactiveVal(
    data.frame(Archivo = character(), `Tamano` = character(),
               stringsAsFactors = FALSE, check.names = FALSE)
  )

  # Navegacion entre pestañas
  process_unlocked <- reactiveVal(FALSE)
  analysis_done    <- reactiveVal(FALSE)
  config_snap      <- reactiveVal(list())
  run_params_rv    <- reactiveVal(list())
  pending_output_dir <- reactiveVal("")
  results_refresh <- reactiveVal(Sys.time())

  # Estado del proceso (real-time con processx)          [v3-PROC]
  proc_rv <- reactiveValues(
    proc       = NULL,
    running    = FALSE,
    start_time = NULL,
    checkpoints = character(0),
    cp_idx     = 0,
    samp_stat  = list(),
    cur_sample = NULL,
    n_total    = 0,
    log_file   = NULL,
    log_seen_size = 0,
    last_output_time = NULL,
    last_heartbeat   = NULL
  )

  # Matriz de conteos
  data_rv <- reactiveValues(
    count_matrix = NULL,
    counts_ready = FALSE,
    source       = "workflow"
  )

  # ── Timer de polling (processx)                       [v3-PROC] ────────────
  poll_timer <- reactiveTimer(500)

  progress_tick <- function() {
    if (isTRUE(proc_rv$running)) poll_timer()
  }

  append_run_log <- function(chunk) {
    chunk <- terminal_text(chunk)
    if (!nzchar(chunk)) return(invisible(FALSE))
    if (!is.null(proc_rv$log_file) && nzchar(proc_rv$log_file)) {
      cat(chunk, file = proc_rv$log_file, append = TRUE)
      proc_rv$log_seen_size <- file.info(proc_rv$log_file)$size %||% 0
    }
    log_text(trim_log_text(paste0(log_text(), chunk)))
    proc_rv$last_output_time <- Sys.time()
    invisible(TRUE)
  }

  sync_run_log_file <- function() {
    log_file <- proc_rv$log_file
    if (is.null(log_file) || !nzchar(log_file) || !file.exists(log_file)) return("")
    size <- file.info(log_file)$size %||% 0
    old_size <- proc_rv$log_seen_size %||% 0
    if (is.na(size) || size <= old_size) return("")
    if (old_size < 0 || old_size > size) old_size <- 0

    con <- file(log_file, open = "rb")
    on.exit(close(con), add = TRUE)
    seek(con, where = old_size, origin = "start")
    raw <- readBin(con, what = "raw", n = size - old_size)
    proc_rv$log_seen_size <- size

    chunk <- terminal_text(rawToChar(raw))
    if (nzchar(chunk)) {
      log_text(trim_log_text(paste0(log_text(), chunk)))
      proc_rv$last_output_time <- Sys.time()
    }
    chunk
  }

  update_progress_from_log <- function(chunk) {
    lines <- strsplit(terminal_text(chunk), "\n", fixed = TRUE)[[1]]
    lines <- lines[nzchar(lines)]
    for (line in lines) {
      # cp_idx representa checkpoints completados. El siguiente queda "en curso".
      if (grepl("fastqc", line, ignore.case=TRUE) && proc_rv$cp_idx < 1L)
        proc_rv$cp_idx <- 1L
      if (grepl("Processing sample:", line, ignore.case=TRUE) ||
          grepl("bowtie2|salmon quant|kallisto quant|fastp", line, ignore.case=TRUE)) {
        proc_rv$cp_idx <- max(proc_rv$cp_idx, 2L)
      }
      if (grepl("samtools sort|samtools index|featureCounts|abundance\\.tsv|abundance\\.h5|quant(?:\\.sf)?", line, ignore.case=TRUE))
        proc_rv$cp_idx <- max(proc_rv$cp_idx, 3L)
      if (grepl("Generating count matrix", line, ignore.case=TRUE)) {
        if (!is.null(proc_rv$cur_sample)) {
          ss <- proc_rv$samp_stat
          ss[[proc_rv$cur_sample]] <- "done"
          proc_rv$samp_stat <- ss
          proc_rv$cur_sample <- NULL
        }
        proc_rv$cp_idx <- max(proc_rv$cp_idx, 4L)
      }
      if (grepl("multiqc", line, ignore.case=TRUE))
        proc_rv$cp_idx <- max(proc_rv$cp_idx, 6L)
      if (grepl("Analysis completed|Analysis completed successfully", line, ignore.case=TRUE))
        proc_rv$cp_idx <- length(proc_rv$checkpoints)
      if (grepl("Processing sample:", line, ignore.case=TRUE)) {
        sname <- trimws(sub(".*Processing sample:\\s*", "", line, ignore.case=TRUE))
        cur   <- proc_rv$cur_sample
        if (!is.null(cur)) {
          ss <- proc_rv$samp_stat; ss[[cur]] <- "done"; proc_rv$samp_stat <- ss
          bytes <- proc_rv$sample_sizes[[cur]] %||% 0
          proc_rv$bytes_done <- (proc_rv$bytes_done %||% 0) + bytes
        }
        proc_rv$cur_sample <- sname
        ss <- proc_rv$samp_stat; ss[[sname]] <- "running"; proc_rv$samp_stat <- ss
      }
    }
  }


  # ── Observer de polling ────────────────────────────────────────────────────
  observe({
    req(HAS_PROCESSX, proc_rv$running, !is.null(proc_rv$proc))
    poll_timer()

    proc <- proc_rv$proc
    if (!proc$is_alive()) {
      # Leer salida restante
      remaining <- sync_run_log_file()
      if (nzchar(remaining)) {
        update_progress_from_log(remaining)
      }

      exit_code <- proc$get_exit_status() %||% 0L
      proc_rv$running <- FALSE
      proc_rv$end_time <- Sys.time()
      shinyjs::disable("stop_btn")
      shinyjs::enable("run_btn")
      append_run_log(paste0(ts_log(paste0("Codigo de salida: ", exit_code)), "\n"))

      # Marcar muestra actual como terminada
      if (!is.null(proc_rv$cur_sample)) {
        ss <- proc_rv$samp_stat
        ss[[proc_rv$cur_sample]] <- "done"
        proc_rv$samp_stat <- ss
        bytes <- proc_rv$sample_sizes[[proc_rv$cur_sample]] %||% 0
        proc_rv$bytes_done <- (proc_rv$bytes_done %||% 0) + bytes
        proc_rv$cur_sample <- NULL
      }
      proc_rv$cp_idx <- length(proc_rv$checkpoints)

      if (exit_code == 0) {
        # Auto-cargar matriz de conteos
        p <- run_params_rv()
        counts <- tryCatch(load_counts_from_workflow(p$output_dir, p$tool), error=function(e) NULL)
        if (!is.null(counts)) {
          data_rv$count_matrix <- counts
          data_rv$counts_ready <- TRUE
          append_run_log(paste0(
            ts_log(sprintf("Matriz cargada: %d genes x %d muestras",
                           nrow(counts), ncol(counts))), "\n"))
        }
        # Actualizar lista de archivos
        files <- list.files(p$output_dir, recursive = TRUE, full.names = FALSE)
        if (length(files) > 0) {
          output_files_rv(file_table_for_files(p$output_dir, files))
        }
        append_run_log(paste0(ts_log("=== Analisis finalizado OK ==="), "\n"))
        analysis_done(TRUE)
        results_refresh(Sys.time())
        updateSelectInput(session, "selected_result_dir",
                          choices = result_choices(outputs_dir),
                          selected = p$output_dir)
        showNotification("Workflow finalizado correctamente.", type = "message")
        updateNavbarPage(session, "main_nav", selected = "tab_results")
      } else {
        append_run_log(paste0(ts_log("=== ERROR en el workflow ==="), "\n"))
        showNotification(
          paste0("Error (codigo ", exit_code, "). Revisa el log de ejecucion."),
          type = "error", duration = 12
        )
      }
    } else {
      # Proceso vivo: sincronizar salida redirigida a workflow_live.log.
      new_chunk <- sync_run_log_file()
      if (nzchar(new_chunk)) {
        update_progress_from_log(new_chunk)
      } else if (!is.null(proc_rv$last_output_time)) {
        now <- Sys.time()
        quiet_for <- as.numeric(difftime(now, proc_rv$last_output_time, units = "secs"))
        heartbeat_for <- if (is.null(proc_rv$last_heartbeat)) Inf else
          as.numeric(difftime(now, proc_rv$last_heartbeat, units = "secs"))
        if (quiet_for >= 30 && heartbeat_for >= 30) {
          append_run_log(paste0(
            ts_log(sprintf("Proceso activo sin nueva salida desde hace %s.", fmt_elapsed(quiet_for))),
            "\n"
          ))
          proc_rv$last_heartbeat <- now
        }
      }
    }
  })

  # ── UI dinamica: label genoma/transcriptoma ────────────────────────────────
  output$genome_label_ui <- renderUI({
    if (isTRUE(input$analysis_type == "alignment"))
      fileInput("genome_file_upload", "Genoma de referencia (FASTA)",
                accept = c(".fa", ".fasta", ".fna", ".gz"), multiple = FALSE)
    else
      fileInput("genome_file_upload", "Transcriptoma de referencia (FASTA)",
                accept = c(".fa", ".fasta", ".fna", ".gz"), multiple = FALSE)
  })

  # ── Seleccion de directorios con shinyFiles
  roots <- c(wd = normalizePath(getwd()), home = normalizePath("~"))
  shinyFiles::shinyDirChoose(input, "input_dir_btn", roots = roots, session = session)

  # Selector específico para carpetas de resultados (outputs/)
  roots_results <- c(outputs = normalizePath(outputs_dir, mustWork = FALSE))
  shinyFiles::shinyDirChoose(input, "select_result_dir_btn", roots = roots_results, session = session)

  input_dir_val <- reactive({
    if (!is.null(input$input_dir_btn)) {
      path <- shinyFiles::parseDirPath(roots, input$input_dir_btn)
      as.character(path)
    } else ""
  })

  output_dir_val <- reactive(pending_output_dir() %||% "")

  output$input_dir_path <- renderText({ input_dir_val() })
  output$output_base_path <- renderText({
    paste0(normalizePath(outputs_dir, mustWork = FALSE), "/<fecha_hora>_<tipo>_<herramienta>")
  })

  # ── Deteccion de muestras con debounce (usando shinyFiles)
  input_dir_debounced <- debounce(reactive(input_dir_val()), 600)
  samples_live <- reactive({
    dir_val <- input_dir_debounced()
    detect_samples(dir_val, input$read_type %||% "pe")
  })

  output$sample_preview_ui <- renderUI({
    dir_val <- input_dir_val() %||% ""
    if (!nzchar(dir_val))
      return(div(class="text-muted", icon("folder-open"),
                 " Introduce un directorio de entrada."))
    if (!dir.exists(dir_val))
      return(div(class="text-danger", icon("circle-xmark"), " El directorio no existe."))
    samples <- samples_live()
    if (length(samples) == 0)
      return(div(class="alert alert-warning", icon("triangle-exclamation"),
                 if (isTRUE(input$read_type == "se"))
                   " No se encontraron FASTQ single-end (*.fastq.gz / *.fastq)."
                 else
                   " No se encontraron FASTQ paired-end (*_1.fastq.gz / *_R1.fastq.gz)."))

    # Buenas practicas: advertir nombres con caracteres problematicos
    bad <- bad_sample_chars(samples)
    r2_miss <- if (isTRUE(input$read_type == "se")) rep(FALSE, length(samples)) else missing_r2(dir_val, samples)

    tagList(
      div(class="alert alert-info py-2",
          icon("circle-check"), sprintf(" %d muestra(s) %s detectadas.",
                                        length(samples), read_type_label(input$read_type %||% "pe"))),
      if (length(bad) > 0)
        div(class="alert alert-warning py-2", icon("triangle-exclamation"),
            " Nombres con caracteres especiales (espacios, @, etc.) pueden causar",
            " problemas: ", tags$b(paste(bad, collapse=", "))),
      if (any(r2_miss))
        div(class="alert alert-danger py-2", icon("circle-xmark"),
            " Faltan R2 para: ", tags$b(paste(samples[r2_miss], collapse=", "))),
      tags$table(
        class="table table-sm table-striped",
        tags$thead(tags$tr(tags$th("Muestra"), tags$th("Lectura 1 / FASTQ"), tags$th("R2"))),
        tags$tbody(lapply(seq_along(samples), function(i) {
          tags$tr(
            tags$td(samples[i]),
            tags$td(tags$span(style="color:#315342;font-weight:700;","\u2713")),
            tags$td(if (isTRUE(input$read_type == "se")) tags$span(class="text-muted", "No aplica")
                    else if (!r2_miss[[i]]) tags$span(style="color:#315342;font-weight:700;","\u2713")
                    else tags$span(style="color:#8A2F2F;font-weight:700;","\u2717 falta"))
          )
        }))
      )
    )
  })

  # ── Validaciones Tab 1 ────────────────────────────────────────────────────
  val_errors <- reactive({
    errs <- character(0)
    dir_in  <- input_dir_val() %||% ""
    gf      <- if (!is.null(input$genome_file_upload) && nrow(input$genome_file_upload) > 0)
           input$genome_file_upload$datapath else ""
    aln     <- input$analysis_type %||% "alignment"

    if (!nzchar(dir_in))           errs <- c(errs, "Directorio de FASTQs: campo vacio.")
    else if (!dir.exists(dir_in))  errs <- c(errs, "Directorio de FASTQs: no existe.")
    else {
      samps <- detect_samples(dir_in, input$read_type %||% "pe")
      if (length(samps) == 0)
        errs <- c(errs, if (isTRUE(input$read_type == "se")) "No se encontraron archivos FASTQ single-end." else "No se encontraron archivos FASTQ R1.")
      else if (!isTRUE(input$read_type == "se")) {
        bad <- samps[missing_r2(dir_in, samps)]
        if (length(bad)) errs <- c(errs, paste("Faltan R2 para:", paste(bad, collapse=", ")))
      }
    }
    if (!nzchar(gf))              errs <- c(errs, "FASTA de referencia: campo vacio.")
    else if (!file.exists(gf))    errs <- c(errs, "FASTA de referencia: no existe.")
    if (aln == "alignment") {
      af <- if (!is.null(input$annotation_file_upload) && nrow(input$annotation_file_upload) > 0)
        input$annotation_file_upload$datapath
      else ""
      if (!nzchar(af))            errs <- c(errs, "Anotacion GFF/GTF: campo vacio.")
      else if (!file.exists(af))  errs <- c(errs, "Anotacion GFF/GTF: no existe.")
    }
    if (isTRUE(input$analysis_type == "pseudo") &&
        identical(input$pseudo_tool %||% "salmon", "kallisto") &&
        isTRUE(input$read_type == "se")) {
      fl <- input$fragment_length %||% NA_real_
      fsd <- input$fragment_sd %||% NA_real_
      if (is.na(fl) || fl <= 0)
        errs <- c(errs, "Kallisto single-end: longitud media de fragmento invalida.")
      if (is.na(fsd) || fsd <= 0)
        errs <- c(errs, "Kallisto single-end: desviacion estandar de fragmento invalida.")
    }
    if (!dir.exists(outputs_dir))
      errs <- c(errs, paste0("No se pudo acceder a la carpeta de salidas: ", outputs_dir))
    if (!file.exists(workflow_path))
      errs <- c(errs, paste0("workflow.sh no encontrado en: ", workflow_path))
    errs
  })

  # Checklist de pasos (UI: circulos rojos/verde segun estado)
  checklist_status <- reactive({
    dir_in <- input_dir_val() %||% ""
    dir_ok <- nzchar(dir_in) && dir.exists(dir_in)
    samples_ok <- FALSE
    if (dir_ok) {
      samps <- detect_samples(dir_in, input$read_type %||% "pe")
      samples_ok <- length(samps) > 0 &&
        (isTRUE(input$read_type == "se") || !any(missing_r2(dir_in, samps)))
    }
    genome_ok <- if (!is.null(input$genome_file_upload) && nrow(input$genome_file_upload) > 0) "ok" else "missing"
    # Annotation: required for alignment, optional for pseudo
    if (isTRUE(input$analysis_type == "alignment")) {
      annot_status <- if (!is.null(input$annotation_file_upload) && nrow(input$annotation_file_upload) > 0) "ok" else "missing"
    } else {
      annot_status <- if (!is.null(input$annotation_file_pseudo_upload) && nrow(input$annotation_file_pseudo_upload) > 0) "ok" else "optional"
    }
    valid_ok <- if (length(val_errors()) == 0) "ok" else "missing"
    list(
      "Directorio FASTQ" = if (dir_ok) "ok" else "missing",
      "Muestras detectadas" = if (samples_ok) "ok" else "missing",
      "Genoma/transcriptoma" = genome_ok,
      "Anotación" = annot_status,
      "Validación" = valid_ok
    )
  })

  output$checklist_ui <- renderUI({
    items <- checklist_status()
    tags$div(style="border:1px solid #CFE5D4;padding:10px;border-radius:8px;background:#F8FCF9;min-width:220px;",
      tags$ul(style="list-style:none;margin:0;padding:0;font-size:0.95rem;width:100%;",
        lapply(names(items), function(nm) {
          st <- items[[nm]]
          color <- switch(as.character(st), ok = "#7BBF9A", optional = "#F6D58A", missing = "#F4A6A6", "#D7EEF1")
          tags$li(style="display:flex;align-items:center;justify-content:space-between;gap:8px;padding:4px 0;",
            tags$span(nm),
            tags$span(style=paste0("display:inline-block;width:14px;height:14px;border-radius:50%;background:",color,";"))
          )
        })
      )
    )
  })

  output$validation_ui <- renderUI({
    errs <- val_errors()
    if (length(errs) == 0)
      div(class="alert alert-success mb-0", icon("circle-check"),
          " Todos los campos correctos. Puedes continuar.")
    else
      div(class="alert alert-danger mb-0",
          tags$b(icon("circle-xmark"), " Corrige antes de continuar:"),
          tags$ul(class="mb-0 mt-1", lapply(errs, tags$li)))
  })

  observe({
    if (length(val_errors()) > 0) shinyjs::disable("btn_to_processing")
    else                           shinyjs::enable("btn_to_processing")
  })

  # ── Herramienta y anotacion efectivas ─────────────────────────────────────
  effective_tool <- reactive({
    if (isTRUE(input$analysis_type == "alignment")) "bowtie2"
    else input$pseudo_tool %||% "salmon"
  })

  effective_read_type <- reactive({
    if (isTRUE(input$read_type == "se")) "se" else "pe"
  })

  effective_fragment_length <- reactive({
    input$fragment_length %||% 200
  })

  effective_fragment_sd <- reactive({
    input$fragment_sd %||% 20
  })

  effective_genome_file <- reactive({
    if (!is.null(input$genome_file_upload) && nrow(input$genome_file_upload) > 0)
      input$genome_file_upload$datapath
    else
      ""
  })

  effective_annotation <- reactive({
    if (isTRUE(input$analysis_type == "alignment")) {
      if (!is.null(input$annotation_file_upload) && nrow(input$annotation_file_upload) > 0)
        input$annotation_file_upload$datapath
      else
        ""
    } else {
      if (!is.null(input$annotation_file_pseudo_upload) && nrow(input$annotation_file_pseudo_upload) > 0)
        input$annotation_file_pseudo_upload$datapath
      else "/dev/null"
    }
  })

  # ── Comando del workflow ──────────────────────────────────────────────────
  workflow_cmd <- reactive({
    req(input_dir_val(), output_dir_val())
    sprintf(
      "bash %s --INPUT %s --OUTPUT %s --GENOME_FILE %s --ANNOTATION_FILE %s --ALIGNMENT_TYPE %s --READ_TYPE %s --FRAGMENT_LENGTH %s --FRAGMENT_SD %s",
      shQuote(workflow_path),
      shQuote(input_dir_val()), shQuote(output_dir_val()),
      shQuote(effective_genome_file()), shQuote(effective_annotation()), shQuote(effective_tool()),
      shQuote(effective_read_type()), shQuote(effective_fragment_length()), shQuote(effective_fragment_sd())
    )
  })

  # ── Navegar Tab 1 -> Tab 2 ────────────────────────────────────────────────
  observeEvent(input$btn_to_processing, {
    req(length(val_errors()) == 0)
    process_unlocked(TRUE)
    samps <- detect_samples(input_dir_val(), effective_read_type())
    run_output_dir <- create_run_output_dir(outputs_dir, input$analysis_type, effective_tool())
    pending_output_dir(run_output_dir)
    config_snap(list(
      analysis_type = input$analysis_type, tool = effective_tool(),
      input_dir  = input_dir_val(), output_dir = run_output_dir,
      genome_file = effective_genome_file(), annotation = effective_annotation(),
      n_samples = length(samps), read_type = read_type_label(effective_read_type()),
      fragment_length = effective_fragment_length(), fragment_sd = effective_fragment_sd()
    ))
    log_text(paste(
      ts_log("=== Configuracion validada ==="),
      ts_log(paste("Analisis:", input$analysis_type, "/", effective_tool())),
      ts_log(paste("Muestras:", length(samps))),
      ts_log(paste("Entrada:", input_dir_val())),
      ts_log(paste("Salida:",  run_output_dir)), "",
      sep = "\n"
    ))
    updateNavbarPage(session, "main_nav", selected = "tab_process")
  })

  observeEvent(input$btn_back, {
    updateNavbarPage(session, "main_nav", selected = "tab_config")
  })

  # ── Carga desde resultados previos                    [v3-LOAD] ────────────
  observeEvent(input$btn_load_existing, {
    counts <- NULL
    loaded_dir <- ""
    loaded_tool <- "matriz subida"
    # Option A: upload file
    if (!is.null(input$upload_counts)) {
      counts <- tryCatch(
        read.table(input$upload_counts$datapath, header=TRUE, row.names=1,
                   sep="\t", comment.char="#"),
        error = function(e)
          tryCatch(read.csv(input$upload_counts$datapath, row.names=1),
                   error=function(e2) NULL)
      )
      loaded_dir <- dirname(input$upload_counts$datapath)
    }
    if (is.null(counts)) {
      showNotification("No se pudo cargar la matriz de conteos. Verifica el formato o usa la pestaña Resultados para abrir ejecuciones previas.",
                       type="error"); return()
    }
    counts <- round(counts)
    data_rv$count_matrix <- counts
    data_rv$counts_ready <- TRUE
    data_rv$source       <- "uploaded"

    files <- if (nzchar(loaded_dir) && dir.exists(loaded_dir)) {
      list.files(loaded_dir, recursive=TRUE, full.names=FALSE)
    } else character(0)
    if (length(files)) {
      output_files_rv(file_table_for_files(loaded_dir, files))
    } else {
      upload_size <- if (!is.null(input$upload_counts)) {
        file.info(input$upload_counts$datapath)$size
      } else {
        NA_real_
      }
      output_files_rv(data.frame(
        Archivo="Matriz cargada sin directorio de resultados.",
        `Tamano`=fmt_bytes(upload_size),
        stringsAsFactors=FALSE, check.names=FALSE
      ))
    }
    run_params_rv(list(
      analysis_type = if (identical(loaded_tool, "bowtie2")) "alignment" else "pseudo",
      tool = loaded_tool,
      input_dir = "resultados previos",
      output_dir = loaded_dir,
      genome_file = "—",
      annotation_file = "—",
      n_samples = ncol(counts),
      read_type = "Paired-end",
      started_at = Sys.time(),
      r_version = paste(R.version$major, R.version$minor, sep=".")
    ))
    analysis_done(TRUE)
    showNotification(
      sprintf("Datos cargados: %d genes x %d muestras.",
              nrow(counts), ncol(counts)), type="message", duration=8
    )
    results_refresh(Sys.time())
    updateNavbarPage(session, "main_nav", selected="tab_results")
  })

  # ── Contenido Tab 2 ───────────────────────────────────────────────────────
  output$tab2_content <- renderUI({
    if (!process_unlocked())
      return(div(class="alert alert-info mt-4", icon("arrow-left-long"),
                 " Completa la pestaña 1 antes de continuar."))
    cfg <- config_snap()
    if (length(cfg) == 0) return(NULL)

    total_sz <- {
      cfg_read_type <- if (identical(cfg$read_type, "Single-end")) "se" else "pe"
      cfg_samples <- detect_samples(cfg$input_dir, cfg_read_type)
      sizes <- sample_fastq_sizes(cfg$input_dir, cfg_samples, cfg_read_type)
      if (length(sizes)) fmt_bytes(sum(sizes, na.rm = TRUE)) else "—"
    }

    tagList(
      # ── Fila superior: resumen + progreso ─────────────────────────────
      layout_columns(
        col_widths = c(4, 8),

        # Izquierda: resumen y controles
        card(
          card_header("Resumen del analisis"),
          tags$dl(class="row small mb-1",
            tags$dt(class="col-6","Tipo:"),
            tags$dd(class="col-6", if(cfg$analysis_type=="alignment")
              "Alineamiento (Bowtie2)" else paste0("Pseudoalineamiento (",cfg$tool,")")),
            tags$dt(class="col-6","Muestras:"), tags$dd(class="col-6", cfg$n_samples),
            tags$dt(class="col-6","Lectura:"),  tags$dd(class="col-6", cfg$read_type),
            if (identical(cfg$tool, "kallisto") && identical(cfg$read_type, "Single-end"))
              tagList(
                tags$dt(class="col-6","Fragmento:"),
                tags$dd(class="col-6", paste0(cfg$fragment_length, " ± ", cfg$fragment_sd))
              ),
            tags$dt(class="col-6","Tama\u00f1o est.:"), tags$dd(class="col-6", total_sz),
            tags$dt(class="col-6","Entrada:"),
            tags$dd(class="col-6", tags$code(class="small", cfg$input_dir)),
            tags$dt(class="col-6","Salida:"),
            tags$dd(class="col-6", tags$code(class="small", cfg$output_dir))
          ),
          hr(),
          tags$details(
            tags$summary(tags$small(icon("terminal")," Ver comando")),
            tags$pre(class="small mt-1",
                     style="white-space:pre-wrap;word-break:break-all;background:#f8f9fa;padding:8px;",
                     textOutput("cmd_preview_text", inline=TRUE))
          ),
          hr(),
          div(style="display:flex;gap:8px;align-items:center;",
              actionButton("btn_back", tagList(icon("arrow-left")," Volver"), class="btn-secondary"),
              actionButton("stop_btn", tagList(icon("stop")," Detener"), class="btn-danger", disabled = "disabled"),
              actionButton("run_btn",  tagList(icon("play")," Ejecutar workflow"), class="btn-success btn-lg")
          )
        ),

        # Derecha: progreso detallado                 [v3-PROG]
        card(
          card_header("Progreso en tiempo real"),

          # Metricas de tiempo y muestras
          div(class="d-flex gap-3 mb-3 p-2 rounded",
              style="background:#F7F1FC; border:1px solid #E3D6F2; font-size:.9rem;",
              div(icon("clock"), " Transcurrido: ",
                  tags$b(textOutput("elapsed_text", inline=TRUE))),
              div(icon("hourglass-half"), " Tiempo restante: ",
                  tags$b(textOutput("eta_text", inline=TRUE))),
              div(icon("percent"), " % completado: ",
                  tags$b(textOutput("remaining_pct", inline=TRUE))),
              div(icon("vials"), " Muestras: ",
                  tags$b(textOutput("samp_prog_text", inline=TRUE)))
          ),

          # Lista de checkpoints
          uiOutput("checkpoints_ui"),
          hr(),
          # Tabla de estado por muestra
          uiOutput("sample_status_ui")
        )
      ),

      # ── Fila inferior: log ancho completo            [v3-LAYOUT] ───────
      card(
        card_header(
          "Log de terminal",
          div(class="float-end",
              actionButton("refresh_btn", icon("rotate"),
                           class="btn-sm btn-outline-secondary",
                           title="Actualizar lista de archivos"))
        ),
        tags$small(class="text-muted",
                   "Se muestra stdout/stderr en vivo. Copia completa: workflow_live.log en el directorio de salida."),
        verbatimTextOutput("run_log"),
        max_height = "420px",
        style = "overflow-y:auto; width:100%; background:#ffffff; color:#111827;"
      )
    )
  })

  # ── Outputs de progreso (dependen de poll_timer para actualizarse) ─────────
  output$elapsed_text <- renderText({
    progress_tick()
    if (is.null(proc_rv$start_time)) return("—")
    end_time <- proc_rv$end_time %||% Sys.time()
    fmt_elapsed(as.numeric(difftime(end_time, proc_rv$start_time, units="secs")))
  })

  output$eta_text <- renderText({
    progress_tick()
    if (is.null(proc_rv$start_time) || !proc_rv$running) return("—")
    elapsed <- as.numeric(difftime(Sys.time(), proc_rv$start_time, units="secs"))
    total_bytes <- proc_rv$total_bytes %||% 0
    bytes_done <- proc_rv$bytes_done %||% 0
    if (total_bytes > 0 && bytes_done > 0) {
      speed <- bytes_done / elapsed
      if (speed > 0) {
        remaining <- (total_bytes - bytes_done) / speed
        if (remaining < 0) return("finalizando...")
        return(paste0("~", fmt_elapsed(remaining)))
      }
    }
    n_done  <- sum(vapply(proc_rv$samp_stat, identical, logical(1), "done"))
    n_total <- proc_rv$n_total
    if (n_total == 0) return("—")
    if (n_done == 0) {
      cur <- proc_rv$cur_sample %||% ""
      if (nzchar(cur)) {
        remaining <- elapsed * (n_total - 1)
        if (remaining < 0) return("finalizando...")
        return(paste0("~", fmt_elapsed(remaining)))
      }
      return("calculando...")
    }
    avg <- elapsed / n_done
    remaining <- avg * (n_total - n_done)
    if (remaining < 0) return("finalizando...")
    paste0("~", fmt_elapsed(remaining))
  })

  output$remaining_pct <- renderText({
    progress_tick()
    total_bytes <- proc_rv$total_bytes %||% 0
    bytes_done <- proc_rv$bytes_done %||% 0
    if (total_bytes > 0 && bytes_done > 0) {
      pct_done <- min(100, bytes_done / total_bytes * 100)
      return(paste0(round(pct_done, 1), "%"))
    }
    n_done <- sum(vapply(proc_rv$samp_stat, identical, logical(1), "done"))
    n_total <- proc_rv$n_total
    if (n_total > 0) {
      cur <- proc_rv$cur_sample %||% ""
      progress <- if (nzchar(cur) && n_done < n_total) {
        (n_done + 0.5) / n_total
      } else {
        n_done / n_total
      }
      pct_done <- min(100, max(0, progress * 100))
      return(paste0(round(pct_done, 1), "%"))
    }
    "—"
  })

  output$samp_prog_text <- renderText({
    progress_tick()
    n_done  <- sum(vapply(proc_rv$samp_stat, identical, logical(1), "done"))
    n_total <- proc_rv$n_total
    if (n_total == 0) return("—")
    cur <- proc_rv$cur_sample %||% ""
    if (nzchar(cur)) {
      sprintf("%d / %d completadas (%s en curso)", n_done, n_total, cur)
    } else {
      sprintf("%d / %d completadas", n_done, n_total)
    }
  })

  output$checkpoints_ui <- renderUI({
    progress_tick()
    cps    <- proc_rv$checkpoints
    cp_idx <- proc_rv$cp_idx
    if (length(cps) == 0) return(NULL)
    tags$ul(class="list-unstyled mb-0",
      lapply(seq_along(cps), function(i) {
        if (i <= cp_idx)
          tags$li(style="color:#315342;font-weight:600;", tags$span("\u2713 "), cps[i])
        else if (proc_rv$running && i == cp_idx + 1L)
          tags$li(style="color:#8A6D1C; font-weight:650;",
                  tags$span("\u27f3 "), cps[i])
        else
          tags$li(style="color:#7C9185;", tags$span("\u25f7 "), cps[i])
      })
    )
  })

  output$sample_status_ui <- renderUI({
    progress_tick()
    ss <- proc_rv$samp_stat
    if (length(ss) == 0) return(NULL)
    tags$table(class="table table-sm",
      tags$thead(tags$tr(tags$th("Muestra"), tags$th("Estado"))),
      tags$tbody(lapply(names(ss), function(s) {
        st <- ss[[s]]
        icon_el <- switch(st,
          done    = tags$span(style="color:#315342;font-weight:600;",  "\u2713 Completada"),
          running = tags$span(style="color:#8A6D1C;font-weight:600;", "\u27f3 Procesando..."),
          cancelled = tags$span(style="color:#60756A;", "\u25a0 Cancelada"),
          tags$span(style="color:#7C9185;", "\u25f7 Pendiente")
        )
        tags$tr(tags$td(s), tags$td(icon_el))
      }))
    )
  })

  output$cmd_preview_text <- renderText({
    tryCatch(workflow_cmd(),
             error=function(e) "Faltan campos para generar el comando.")
  })

  # ── Ejecutar workflow ──────────────────────────────────────────────────────
  observeEvent(input$run_btn, {
    errs <- val_errors()
    if (length(errs) > 0) {
      showNotification(paste(errs, collapse="\n"), type="error", duration=8); return()
    }
    if (!file.exists(workflow_path)) {
      showNotification("workflow.sh no encontrado.", type="error"); return()
    }
    if (!nzchar(output_dir_val())) {
      pending_output_dir(create_run_output_dir(outputs_dir, input$analysis_type, effective_tool()))
    }

    genome_path <- tryCatch(
      prepare_uploaded_input_file(input$genome_file_upload, output_dir_val(), prefix = "genome"),
      error = function(e) {
        showNotification(paste0("Error al preparar el archivo FASTA: ", conditionMessage(e)), type = "error")
        NULL
      }
    )
    if (is.null(genome_path)) return()

    annotation_path <- tryCatch({
      if (isTRUE(input$analysis_type == "alignment")) {
        prepare_uploaded_input_file(input$annotation_file_upload, output_dir_val(), prefix = "annotation")
      } else {
        prepare_uploaded_input_file(input$annotation_file_pseudo_upload, output_dir_val(), prefix = "annotation", optional = TRUE)
      }
    }, error = function(e) {
      showNotification(paste0("Error al preparar el archivo de anotación: ", conditionMessage(e)), type = "error")
      NULL
    })
    if (is.null(annotation_path)) return()

    cmd <- sprintf(
      "bash %s --INPUT %s --OUTPUT %s --GENOME_FILE %s --ANNOTATION_FILE %s --ALIGNMENT_TYPE %s --READ_TYPE %s --FRAGMENT_LENGTH %s --FRAGMENT_SD %s",
      shQuote(workflow_path), shQuote(input_dir_val()), shQuote(output_dir_val()),
      shQuote(genome_path), shQuote(annotation_path), shQuote(effective_tool()), shQuote(effective_read_type()),
      shQuote(effective_fragment_length()), shQuote(effective_fragment_sd())
    )

    # Snapshot de parametros para Tab 3
    samps <- detect_samples(input_dir_val(), effective_read_type())
    run_params_rv(list(
      analysis_type = input$analysis_type, tool = effective_tool(),
      input_dir = input_dir_val(), output_dir = output_dir_val(),
      genome_file = genome_path, annotation_file = annotation_path,
      n_samples = length(samps), read_type = read_type_label(effective_read_type()),
      fragment_length = effective_fragment_length(), fragment_sd = effective_fragment_sd(),
      started_at = Sys.time(),
      r_version = paste(R.version$major, R.version$minor, sep=".")
    ))

    analysis_done(FALSE)
    data_rv$counts_ready <- FALSE
    shinyjs::disable("run_btn")
    shinyjs::disable("stop_btn")
    on.exit({ if (!proc_rv$running) shinyjs::enable("run_btn") }, add=TRUE)

    # Checkpoints segun tipo de analisis
    cps <- if (input$analysis_type == "alignment")
      c("Construyendo indice Bowtie2",
        "Control de calidad inicial (FastQC)",
        "Alineamiento de muestras (Bowtie2 + fastp)",
        "Procesando BAM (samtools sort + index)",
        "Conteos por gen (featureCounts)",
        "Control de calidad post-trimming (FastQC)",
        "Informe global (MultiQC)")
    else
      c(paste0("Construyendo indice (", effective_tool(), ")"),
        "Control de calidad inicial (FastQC)",
        paste0("Cuantificacion de muestras (", effective_tool(), " + fastp)"),
        "Importando cuantificaciones",
        "Matriz de conteos",
        "Control de calidad post-trimming (FastQC)",
        "Informe global (MultiQC)")

    # Inicializar proc_rv
    proc_rv$checkpoints <- cps
    proc_rv$cp_idx      <- 0L
    proc_rv$n_total     <- length(samps)
    proc_rv$samp_stat   <- setNames(as.list(rep("pending", length(samps))), samps)
    proc_rv$cur_sample  <- NULL
    proc_rv$sample_sizes <- sample_fastq_sizes(input_dir_val(), samps, effective_read_type())
    proc_rv$total_bytes <- sum(proc_rv$sample_sizes)
    proc_rv$bytes_done <- 0
    proc_rv$start_time  <- Sys.time()
    proc_rv$end_time    <- NULL
    proc_rv$last_output_time <- Sys.time()
    proc_rv$last_heartbeat   <- NULL
    proc_rv$log_file <- file.path(output_dir_val(), "workflow_live.log")
    proc_rv$log_seen_size <- 0
    cat("", file = proc_rv$log_file, append = FALSE)

    append_run_log(paste0(
      ts_log("=== Iniciando analisis ==="), "\n",
      ts_log("Log completo en: "), proc_rv$log_file, "\n",
      ts_log(paste0("Lanzando: ", cmd)), "\n"))

    if (HAS_PROCESSX) {
      # ── Ejecucion no bloqueante con processx          [v3-PROC] ─────────
      proc <- tryCatch(
        processx::process$new("bash", c("-lc", cmd),
                              stdout=proc_rv$log_file, stderr="2>&1",
                              env=c("PATH"=Sys.getenv("PATH"))),
        error=function(e) { showNotification(conditionMessage(e), type="error"); NULL }
      )
      if (!is.null(proc)) {
        proc_rv$proc    <- proc
        proc_rv$running <- TRUE
        shinyjs::enable("stop_btn")
      } else {
        shinyjs::enable("run_btn")
      }
    } else {
      # ── Fallback bloqueante con system2 ────────────────────────────────
      withProgress(message="Ejecutando analisis RNA-seq...", value=0, {
        setProgress(0.05, detail="Preparando comando...")
        setProgress(0.15, detail=cps[1])
        setProgress(0.25, detail=cps[2])
        result <- tryCatch(
          system2("bash", c("-lc", cmd), stdout=TRUE, stderr=TRUE),
          error=function(e) structure(conditionMessage(e), status=1L)
        )
        exit_status <- attr(result, "status") %||% 0L
        append_run_log(paste0(terminal_text(result), "\n",
                              ts_log(paste0("Codigo de salida: ", exit_status)), "\n"))
        setProgress(0.85, detail="Cargando resultados...")

        # Cargar conteos y archivos
        p <- run_params_rv()
        if (exit_status == 0) {
          counts <- tryCatch(load_counts_from_workflow(p$output_dir, p$tool), error=function(e) NULL)
          if (!is.null(counts)) {
            data_rv$count_matrix <- counts
            data_rv$counts_ready <- TRUE
          }
          files <- list.files(p$output_dir, recursive=TRUE, full.names=FALSE)
          if (length(files)) {
            output_files_rv(file_table_for_files(p$output_dir, files))
          }
          append_run_log(paste0(ts_log("=== Analisis finalizado OK ==="), "\n"))
          analysis_done(TRUE)
          results_refresh(Sys.time())
          updateSelectInput(session, "selected_result_dir",
                            choices = result_choices(outputs_dir),
                            selected = p$output_dir)
          proc_rv$cp_idx <- length(proc_rv$checkpoints)
          setProgress(1, detail="Completado.")
          showNotification("Workflow finalizado correctamente.", type="message")
          updateNavbarPage(session, "main_nav", selected="tab_results")
        } else {
          setProgress(1, detail="Finalizado con error.")
          showNotification(paste0("Error (codigo ",exit_status,"). Revisa el log."),
                           type="error", duration=12)
        }
      })
    }
  })

  # ── Boton de detener proceso ───────────────────────────────────────────
  observeEvent(input$stop_btn, {
    req(HAS_PROCESSX, proc_rv$running, !is.null(proc_rv$proc))
    proc <- proc_rv$proc
    if (proc$is_alive()) {
      proc$kill()
      proc$wait(2000)
      proc_rv$running <- FALSE
      proc_rv$end_time <- Sys.time()
      proc_rv$proc <- NULL
      proc_rv$cp_idx <- length(proc_rv$checkpoints)
      proc_rv$cur_sample <- NULL
      ss <- proc_rv$samp_stat
      if (length(ss) > 0L) {
        proc_rv$samp_stat <- lapply(ss, function(x) if (x %in% c("running", "pending")) "cancelled" else x)
      }
      append_run_log(paste0(ts_log("=== Proceso detenido por el usuario ==="), "\n"))
      shinyjs::disable("stop_btn")
      shinyjs::enable("run_btn")
      showNotification("Proceso detenido.", type="warning")
    }
  })

  # ── Log y refresh de archivos ─────────────────────────────────────────────
  output$run_log <- renderText({ log_text() })

  observeEvent(input$refresh_btn, {
    p <- run_params_rv()
    out_dir <- p$output_dir %||% output_dir_val() %||% ""
    if (!nzchar(out_dir) || !dir.exists(out_dir)) {
      output_files_rv(data.frame(Archivo="Directorio no existe.", `Tamano`="—",
                                 stringsAsFactors=FALSE, check.names=FALSE)); return()
    }
    files <- list.files(out_dir, recursive=TRUE, full.names=FALSE)
    if (!length(files)) {
      output_files_rv(data.frame(Archivo="Sin archivos.", `Tamano`="—",
                                 stringsAsFactors=FALSE, check.names=FALSE))
    } else {
      output_files_rv(file_table_for_files(out_dir, files))
    }
  })

  available_result_choices <- reactive({
    results_refresh()
    result_choices(outputs_dir)
  })

  selected_result_dir <- reactive({
    choices <- available_result_choices()
    selected <- input$selected_result_dir %||% ""
    if (nzchar(selected) && selected %in% unname(choices)) return(selected)

    current <- run_params_rv()$output_dir %||% ""
    if (nzchar(current) && current %in% unname(choices)) return(current)

    if (length(choices)) return(unname(choices)[1])

    if (analysis_done() && nzchar(current)) return(current)
    ""
  })

  selected_result_params <- reactive({
    out_dir <- selected_result_dir()
    current <- run_params_rv()
    current_dir <- current$output_dir %||% ""
    if (length(current) && nzchar(current_dir) && nzchar(out_dir) &&
        identical(normalizePath(current_dir, mustWork = FALSE),
                  normalizePath(out_dir, mustWork = FALSE))) {
      return(current)
    }
    if (nzchar(out_dir) && dir.exists(out_dir)) return(infer_result_params(out_dir, workflow_path))
    list()
  })

  selected_result_files <- reactive({
    file_table_for_dir(selected_result_dir())
  })

  selected_result_summary <- reactive({
    p <- selected_result_params()
    out_dir <- selected_result_dir()
    if (!length(p) || !nzchar(out_dir)) return(NULL)
    summarise_result(out_dir, p)
  })

  selected_counts_tables <- reactive({
    p <- selected_result_params()
    out_dir <- selected_result_dir()
    if (!length(p) || !nzchar(out_dir)) {
      return(list(libs = data.frame(Mensaje = "Sin ejecucion seleccionada."),
                  top = data.frame(Mensaje = "Sin ejecucion seleccionada.")))
    }
    counts_tables(out_dir, p$tool %||% "")
  })

  selected_multiqc_href <- reactive({
    out_dir <- selected_result_dir()
    report <- file.path(out_dir, "multiqc_report.html")
    if (!nzchar(out_dir) || !file.exists(report)) return("")

    base <- normalizePath(outputs_dir, winslash = "/", mustWork = TRUE)
    target <- normalizePath(report, winslash = "/", mustWork = TRUE)
    if (!startsWith(target, paste0(base, "/"))) return("")

    rel <- substring(target, nchar(base) + 2L)
    rel_parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
    paste0("saved_outputs/", paste(vapply(rel_parts, URLencode, character(1), reserved = TRUE), collapse = "/"))
  })

  output$result_selector_ui <- renderUI({
    choices <- available_result_choices()
    if (!length(choices)) {
      return(div(class="alert alert-info mb-0", icon("folder-open"),
                 " Aun no hay ejecuciones guardadas en outputs/."))
    }
    tagList(
      div(style="display:flex;gap:10px;align-items:flex-end;width:100%;",
        div(style="flex:1 1 auto;min-width:0;",
          selectizeInput(
            "selected_result_dir",
            "Ejecucion guardada",
            choices = choices,
            selected = selected_result_dir(),
            options = list(placeholder = 'Escribe para buscar...'),
            width = "100%"
          )
        ),
        div(style="flex:0 0 auto;margin-bottom:16px;",
          shinyDirButton("select_result_dir_btn", "Seleccionar carpeta...", "Seleccionar...")
        )
      ),
      tags$div(style="margin-top:6px;", tags$small(class="text-muted", "También puedes escribir para filtrar las ejecuciones guardadas."))
    )
  })

  output$multiqc_open_ui <- renderUI({
    href <- selected_multiqc_href()
    if (!nzchar(href)) {
      return(tags$small(class = "text-muted", "MultiQC no disponible para esta ejecucion."))
    }
    tags$a(
      href = href,
      target = "_blank",
      rel = "noopener noreferrer",
      class = "btn btn-sm btn-primary",
      style = "align-self:flex-start;",
      title = "Abrir el informe MultiQC de esta ejecución en una nueva pestaña",
      tagList(icon("up-right-from-square"), " Abrir MultiQC")
    )
  })

  observeEvent(input$refresh_results_btn, {
    results_refresh(Sys.time())
    choices <- result_choices(outputs_dir)
    updateSelectInput(session, "selected_result_dir",
                      choices = choices,
                      selected = if (length(choices)) unname(choices)[1] else character(0))
  })

  # ── Contenido Tab 3 ───────────────────────────────────────────────────────
  output$tab3_content <- renderUI({
    if (!analysis_done() && !length(available_result_choices()))
      return(div(class="alert alert-info mt-4", icon("clock"),
                 " Ejecuta el workflow en la pestaña 2 para ver los resultados o usa una ejecucion previa guardada en outputs/."))
    p <- selected_result_params()
    if (!length(p)) return(NULL)
    s <- selected_result_summary()
    tagList(
        layout_columns(col_widths = c(6, 6),
        # Left: compact 2x2 boxes placed directly (no card)
        div(style="display:grid; grid-template-columns:repeat(2,1fr); gap:10px;",
          div(class = "metric-card",
            div(class = "metric-card-value",
              if (!is.null(s)) status_badge(s$status) else "—"),
            div(class = "metric-card-label","Estado")
          ),
          div(class = "metric-card",
            div(class = "metric-card-value",
              if (!is.null(s)) paste0((s$n_samples %||% "—"), " / ", (s$n_features %||% "—")) else "—"),
            div(class = "metric-card-label","Muestras / genes")
          ),
          div(class = "metric-card",
            div(class = "metric-card-value",
              if (!is.null(s)) pct_label(s$mean_mapped) else "—"),
            div(class = "metric-card-label","Mapeo medio")
          ),
          div(class = "metric-card",
            div(class = "metric-card-value",
              if (!is.null(s)) (s$n_features %||% "—") else "—"),
            div(class = "metric-card-label","Genes detectados")
          )
        ),
        # Right: vertical selector for results (search + folder chooser)
        card(
          card_header("Resultados"),
          div(style="min-height:160px; display:flex; flex-direction:column; justify-content:center; gap:10px;",
              uiOutput("result_selector_ui"),
              uiOutput("multiqc_open_ui")
          )
        )
      ),
      navset_tab(
        id = "results_tabs",
        nav_panel(
          "Resumen",
          layout_columns(
            col_widths = c(5, 7),
            card(
              card_header("Interpretacion rapida"),
              uiOutput("result_interpretation_ui")
            ),
            card(
              download_header("Estadisticas principales MultiQC", "download_run_stats"),
              DTOutput("run_stats_table")
            )
          )
        ),
        nav_panel(
          "Calidad",
          layout_columns(
            col_widths = c(5, 7),
            card(
              download_header("Estado FastQC", "download_fastqc"),
              DTOutput("fastqc_table")
            ),
            card(
              download_header("Trimming y cuantificacion/alineamiento", "download_alignment"),
              DTOutput("alignment_table")
            )
          )
        ),
        nav_panel(
          "Conteos",
          layout_columns(
            col_widths = c(5, 7),
            card(
              download_header("Librerias por muestra", "download_count_lib"),
              DTOutput("count_lib_table")
            ),
            card(
              download_header("Genes/transcritos mas abundantes", "download_count_top"),
              DTOutput("count_top_table")
            )
          )
        ),
        nav_panel(
          "Informes y archivos",
          card(
            download_header("Informes principales", "download_artifacts"),
            DTOutput("artifact_table")
          ),
          card(
            download_header("Todos los archivos generados", "download_filelist"),
            DTOutput("output_files_table")
          )
        ),
        nav_panel(
          "Log",
          card(
            card_header("Ultimas lineas del workflow"),
            verbatimTextOutput("selected_log_tail")
          )
        )
      ),
      # ── INICIO: Resultados adicionales de alineamiento y pseudoalineamiento ──
      tags$div(
        class = "mt-4",
        tags$h3("Resultados adicionales de alineamiento y pseudoalineamiento"),
        navset_tab(
          id = "method_qc_tabs",
          nav_panel(
            "Alineamiento",
            tags$h4("Control de calidad del alineamiento"),
            layout_columns(
              col_widths = c(6, 6),
              card(download_header("Resumen de alineamiento", "download_align_qc_summary"),
                   DTOutput("align_qc_summary_table")),
              card(download_header("Alertas interpretativas", "download_align_qc_alerts"),
                   DTOutput("align_qc_alerts_table"))
            ),
            layout_columns(
              col_widths = c(6, 6),
              card(download_header("Lecturas unicas, multimapeadas y no alineadas", "download_align_qc_mapping_plot"),
                   plotly::plotlyOutput("align_qc_mapping_plot", height = "320px")),
              card(download_header("Lecturas asignadas y no asignadas", "download_align_qc_assignment_plot"),
                   plotly::plotlyOutput("align_qc_assignment_plot", height = "320px"))
            ),
            layout_columns(
              col_widths = c(6, 6),
              card(download_header("Distribucion exonica, intronica e intergenica", "download_align_qc_region_plot"),
                   plotly::plotlyOutput("align_qc_region_plot", height = "280px")),
              card(download_header("Cobertura 5'-3' del cuerpo genico", "download_align_qc_gene_body_plot"),
                   plotly::plotlyOutput("align_qc_gene_body_plot", height = "280px"))
            )
          ),
          nav_panel(
            "Pseudoalineamiento",
            tags$h4("Control de calidad del pseudoalineamiento"),
            layout_columns(
              col_widths = c(12),
              card(download_header("Resumen de pseudoalineamiento", "download_pseudo_qc_summary"),
                   DTOutput("pseudo_qc_summary_table"))
            ),
            layout_columns(
              col_widths = c(6, 6),
              card(download_header("Pseudoalignment rate por muestra", "download_pseudo_qc_rate_plot"),
                   plotly::plotlyOutput("pseudo_qc_rate_plot", height = "300px")),
              card(download_header("Distribucion de TPM por muestra", "download_pseudo_qc_tpm_plot"),
                   plotly::plotlyOutput("pseudo_qc_tpm_plot", height = "300px"))
            ),
            layout_columns(
              col_widths = c(6, 6),
              card(download_header("Transcritos detectados por muestra", "download_pseudo_qc_detected_plot"),
                   plotly::plotlyOutput("pseudo_qc_detected_plot", height = "300px")),
              card(download_header("TPM vs NumReads", "download_pseudo_qc_scatter_plot"),
                   plotly::plotlyOutput("pseudo_qc_scatter_plot", height = "300px"))
            ),
            layout_columns(
              col_widths = c(12),
              card(download_header("Cuantificacion", "download_pseudo_quant"),
                   DTOutput("pseudo_qc_quant_table"))
            ),
            card(download_header("Alertas interpretativas", "download_pseudo_qc_alerts"),
                 DTOutput("pseudo_qc_alerts_table"))
          )
        )
      ),
      # ── FIN: Resultados adicionales de alineamiento y pseudoalineamiento ──
      ## reproducibilidad eliminada
    )
  })

  output$output_files_table <- renderDT({
    dt_table(selected_result_files(), page_length = 20)
  })

  output$result_interpretation_ui <- renderUI({
    s <- selected_result_summary()
    p <- selected_result_params()
    out_dir <- selected_result_dir()
    if (is.null(s) || !length(p)) return(NULL)

    items <- list(
      tags$li(tags$b("Estado: "), status_badge(s$status)),
      tags$li(tags$b("Carpeta: "), tags$code(out_dir)),
      tags$li(tags$b("Tamano total: "), s$total_size),
      tags$li(tags$b("MultiQC: "), if (isTRUE(s$has_multiqc)) "disponible" else "no encontrado")
    )

    if (!is.na(s$mean_mapped)) {
      items <- c(items, list(tags$li(tags$b("Mapeo medio: "), pct_label(s$mean_mapped),
                                if (s$mean_mapped < 70) tags$span(class="text-danger", "  Revisar: bajo para RNA-seq bacteriano.") else NULL)))
    }
    if (!is.na(s$mean_trim_survival)) {
      items <- c(items, list(tags$li(tags$b("Lecturas retenidas tras fastp: "), pct_label(s$mean_trim_survival),
                                if (s$mean_trim_survival < 80) tags$span(class="text-warning", "  Perdida alta durante trimming.") else NULL)))
    }
    if (!is.na(s$mean_q30)) {
      items <- c(items, list(tags$li(tags$b("Q30 post-trimming: "), pct_label(100 * s$mean_q30),
                                if (s$mean_q30 < 0.85) tags$span(class="text-warning", "  Calidad post-trimming moderada/baja.") else NULL)))
    }
    if (!is.na(s$fastqc_fail)) {
      items <- c(items, list(tags$li(tags$b("Checks FastQC fallidos: "), s$fastqc_fail,
                                if (s$fastqc_fail > 0) tags$span(class="text-warning", "  Revisar pestaña Calidad.") else NULL)))
    }

    tagList(
      tags$ul(class="mb-0", items),
      if (file.exists(file.path(out_dir, "multiqc_report.html")))
        div(class="alert alert-info mt-3 mb-0",
            icon("circle-info"),
            " Abre multiqc_report.html desde la carpeta de resultados para el informe interactivo completo.")
    )
  })

  output$run_stats_table <- renderDT({
    dt_table(result_general_table(selected_result_dir()))
  })

  output$fastqc_table <- renderDT({
    dt_table(fastqc_table(selected_result_dir()), page_length = 12)
  })

  output$alignment_table <- renderDT({
    p <- selected_result_params()
    dt_table(alignment_table(selected_result_dir(), p$tool %||% ""), page_length = 12)
  })

  output$count_lib_table <- renderDT({
    dt_table(selected_counts_tables()$libs, page_length = 12)
  })

  output$count_top_table <- renderDT({
    dt_table(selected_counts_tables()$top, page_length = 15)
  })

  output$artifact_table <- renderDT({
    dt_table(important_artifacts(selected_result_dir()), page_length = 12)
  })

  output$selected_log_tail <- renderText({
    log_tail_text(selected_result_dir())
  })

  ## reproducibilidad eliminada
  csv_download <- function(prefix, data_fun) {
    downloadHandler(
      filename = function() paste0(prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(f) write.csv(data_fun(), f, row.names = FALSE)
    )
  }

  plotly_download <- function(prefix, plot_fun) {
    downloadHandler(
      filename = function() paste0(prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip"),
      content = function(f) {
        if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
          stop("El paquete htmlwidgets es necesario para descargar figuras Plotly.")
        }
        tmp_dir <- tempfile(prefix)
        dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
        html_file <- file.path(tmp_dir, "index.html")
        htmlwidgets::saveWidget(
          plot_fun(),
          file = html_file,
          selfcontained = FALSE,
          libdir = "lib"
        )

        image_files <- character(0)
        image_error <- NULL
        if (requireNamespace("webshot2", quietly = TRUE)) {
          png_file <- file.path(tmp_dir, "plot.png")
          image_error <- tryCatch({
            webshot2::webshot(html_file, png_file, vwidth = 1400, vheight = 900, zoom = 1.4)
            NULL
          }, error = function(e) conditionMessage(e))

          if (file.exists(png_file) && file.info(png_file)$size > 0) {
            image_files <- c(image_files, "plot.png")
          } else {
            jpg_file <- file.path(tmp_dir, "plot.jpg")
            image_error <- tryCatch({
              webshot2::webshot(html_file, jpg_file, vwidth = 1400, vheight = 900, zoom = 1.4)
              NULL
            }, error = function(e) conditionMessage(e))
            if (file.exists(jpg_file) && file.info(jpg_file)$size > 0) {
              image_files <- c(image_files, "plot.jpg")
            }
          }
        } else {
          image_error <- "El paquete webshot2 no esta instalado."
        }

        if (!length(image_files)) {
          writeLines(
            c(
              "No se pudo generar una imagen PNG/JPEG de la figura.",
              "La figura interactiva esta disponible en index.html.",
              "",
              "Para generar imagenes, instala o configura Chrome/Chromium para webshot2/chromote.",
              paste("Detalle:", image_error %||% "sin detalle disponible")
            ),
            file.path(tmp_dir, "README_imagen.txt")
          )
        }

        old_wd <- setwd(tmp_dir)
        on.exit(setwd(old_wd), add = TRUE)
        zip_files <- c("index.html", "lib", image_files)
        if (!length(image_files)) zip_files <- c(zip_files, "README_imagen.txt")
        utils::zip(zipfile = f, files = zip_files, flags = "-r9X")
      }
    )
  }

  output$download_run_stats <- csv_download(
    "multiqc_stats",
    function() result_general_table(selected_result_dir())
  )
  output$download_fastqc <- csv_download(
    "fastqc",
    function() fastqc_table(selected_result_dir())
  )
  output$download_alignment <- csv_download(
    "alignment_metrics",
    function() {
      p <- selected_result_params()
      alignment_table(selected_result_dir(), p$tool %||% "")
    }
  )
  output$download_count_lib <- csv_download(
    "count_libraries",
    function() selected_counts_tables()$libs
  )
  output$download_count_top <- csv_download(
    "top_counts",
    function() selected_counts_tables()$top
  )
  output$download_artifacts <- csv_download(
    "artifacts",
    function() important_artifacts(selected_result_dir())
  )

  output$download_log <- downloadHandler(
    filename=function() paste0("rnaseq_log_",format(Sys.time(),"%Y%m%d_%H%M%S"),".txt"),
    content=function(f) {
      log_file <- file.path(selected_result_dir(), "workflow_live.log")
      if (file.exists(log_file)) file.copy(log_file, f, overwrite = TRUE)
      else writeLines(log_text(), f)
    })
  output$download_filelist <- downloadHandler(
    filename=function() paste0("output_files_",format(Sys.time(),"%Y%m%d_%H%M%S"),".csv"),
    content=function(f) write.csv(selected_result_files(), f, row.names=FALSE))

  # ── INICIO: Server QC adicional de alineamiento y pseudoalineamiento ─────
  output$align_qc_summary_table <- renderDT({
    dt_table(align_qc_summary(selected_result_dir()))
  })

  output$align_qc_alerts_table <- renderDT({
    dt_table(align_qc_alerts(selected_result_dir()))
  })

  output$pseudo_qc_summary_table <- renderDT({
    dt_table(pseudo_qc_summary(selected_result_dir()))
  })

  output$pseudo_qc_quant_table <- renderDT({
    dt_table(pseudo_qc_quant_table(selected_result_dir()), filter = "top")
  })

  output$pseudo_qc_alerts_table <- renderDT({
    dt_table(pseudo_qc_alerts(selected_result_dir()))
  })

  align_qc_mapping_plot_obj <- function() {
    df <- align_qc_summary(selected_result_dir())
    if (!has_real_rows(df)) return(plotly_message(df$Mensaje[1]))
    cols <- c("unique_reads", "multimapped_reads", "unmapped_reads")
    if (!all(cols %in% names(df)) || all(is.na(as.matrix(df[, cols]))))
      return(plotly_message("No se encontraron metricas suficientes para representar lecturas unicas, multimapeadas y no alineadas."))
    df[, cols] <- lapply(df[, cols, drop = FALSE], function(x) ifelse(is.na(x), 0, x))
    plotly::plot_ly(df, x = ~sample_id) |>
      plotly::add_bars(y = ~unique_reads, name = "Unicas", marker = list(color = "#7BBF9A")) |>
      plotly::add_bars(y = ~multimapped_reads, name = "Multimapeadas", marker = list(color = "#F6D58A")) |>
      plotly::add_bars(y = ~unmapped_reads, name = "No alineadas", marker = list(color = "#F4A6A6")) |>
      plotly::layout(
        barmode = "stack",
        xaxis = list(title = "Muestra"),
        yaxis = list(title = "Lecturas"),
        legend = list(orientation = "h", x = 0, y = 1.12),
        margin = list(b = 90)
      )
  }

  align_qc_assignment_plot_obj <- function() {
    df <- align_qc_summary(selected_result_dir())
    if (!has_real_rows(df)) return(plotly_message(df$Mensaje[1]))
    cols <- c("assigned_reads", "unassigned_reads")
    if (!all(cols %in% names(df)) || all(is.na(as.matrix(df[, cols]))))
      return(plotly_message("No se encontraron metricas de asignacion genica para este analisis."))
    df[, cols] <- lapply(df[, cols, drop = FALSE], function(x) ifelse(is.na(x), 0, x))
    plotly::plot_ly(df, x = ~sample_id) |>
      plotly::add_bars(y = ~assigned_reads, name = "Asignadas", marker = list(color = "#8BC9A6")) |>
      plotly::add_bars(y = ~unassigned_reads, name = "No asignadas", marker = list(color = "#D7EEF1")) |>
      plotly::layout(
        barmode = "stack",
        xaxis = list(title = "Muestra"),
        yaxis = list(title = "Lecturas"),
        legend = list(orientation = "h", x = 0, y = 1.12),
        margin = list(b = 90)
      )
  }

  align_qc_region_plot_obj <- function() {
    plotly_message("No se encontraron metricas exonicas, intronicas o intergenicas para este analisis.")
  }

  align_qc_gene_body_plot_obj <- function() {
    plotly_message("No se encontro cobertura 5'-3' del cuerpo genico para este analisis.")
  }

  pseudo_qc_rate_plot_obj <- function() {
    df <- pseudo_qc_summary(selected_result_dir())
    if (!has_real_rows(df)) return(plotly_message(df$Mensaje[1]))
    if (!"pseudoalignment_rate" %in% names(df) || all(is.na(df$pseudoalignment_rate)))
      return(plotly_message("No se encontro pseudoalignment_rate para este analisis."))
    df$pseudoalignment_rate_pct <- 100 * df$pseudoalignment_rate
    plotly::plot_ly(
      df, x = ~sample_id, y = ~pseudoalignment_rate_pct,
      type = "bar", marker = list(color = "#8BC9A6"),
      text = ~paste0(round(pseudoalignment_rate_pct, 2), "%"),
      hovertemplate = "Muestra: %{x}<br>Rate: %{y:.2f}%<extra></extra>"
    ) |>
      plotly::layout(
        xaxis = list(title = "Muestra"),
        yaxis = list(title = "Pseudoalignment rate (%)", range = c(0, 100)),
        shapes = list(list(
          type = "line", xref = "paper", x0 = 0, x1 = 1,
          y0 = 100 * qc_thresholds$pseudoalignment_rate_warning,
          y1 = 100 * qc_thresholds$pseudoalignment_rate_warning,
          line = list(color = "#F4A6A6", dash = "dash")
        )),
        margin = list(b = 90)
      )
  }

  pseudo_qc_tpm_plot_obj <- function() {
    q <- pseudo_qc_quant_table(selected_result_dir())
    if (!has_real_rows(q)) return(plotly_message(q$Mensaje[1]))
    if (!all(c("TPM", "sample_id") %in% names(q)))
      return(plotly_message("No se encontraron valores TPM para este analisis."))
    q$TPM <- num_or_na(q$TPM)
    q$log_tpm <- log10(q$TPM + 1)
    plotly::plot_ly(
      q, x = ~sample_id, y = ~log_tpm,
      type = "box", color = ~sample_id,
      boxpoints = "outliers",
      hovertemplate = "Muestra: %{x}<br>log10(TPM+1): %{y:.3f}<extra></extra>"
    ) |>
      plotly::layout(
        showlegend = FALSE,
        xaxis = list(title = "Muestra"),
        yaxis = list(title = "log10(TPM + 1)"),
        margin = list(b = 90)
      )
  }

  pseudo_qc_detected_plot_obj <- function() {
    df <- pseudo_qc_summary(selected_result_dir())
    if (!has_real_rows(df)) return(plotly_message(df$Mensaje[1]))
    if (!"transcripts_detected" %in% names(df) || all(is.na(df$transcripts_detected)))
      return(plotly_message("No se pudo calcular el numero de transcritos detectados."))
    plotly::plot_ly(
      df, x = ~sample_id, y = ~transcripts_detected,
      type = "bar", marker = list(color = "#7BBF9A"),
      hovertemplate = "Muestra: %{x}<br>Transcritos detectados: %{y}<extra></extra>"
    ) |>
      plotly::layout(
        xaxis = list(title = "Muestra"),
        yaxis = list(title = "Transcritos detectados"),
        margin = list(b = 90)
      )
  }

  pseudo_qc_scatter_plot_obj <- function() {
    q <- pseudo_qc_quant_table(selected_result_dir())
    if (!has_real_rows(q)) return(plotly_message(q$Mensaje[1]))
    if (!all(c("TPM", "NumReads", "sample_id") %in% names(q)))
      return(plotly_message("No se encontraron TPM y NumReads para generar el scatter plot."))
    q$TPM <- num_or_na(q$TPM)
    q$NumReads <- num_or_na(q$NumReads)
    keep <- is.finite(q$TPM) & is.finite(q$NumReads)
    if (!any(keep)) return(plotly_message("No hay valores numericos validos de TPM y NumReads."))
    q <- q[keep, , drop = FALSE]
    q$log_tpm <- log10(q$TPM + 1)
    q$log_reads <- log10(q$NumReads + 1)
    hover_text <- paste0(
      "Muestra: ", q$sample_id,
      "<br>Name: ", q$Name %||% "",
      if ("Type" %in% names(q)) paste0("<br>Type: ", q$Type) else "",
      "<br>TPM: ", round(q$TPM, 3),
      "<br>NumReads: ", round(q$NumReads, 3)
    )
    plotly::plot_ly(
      q, x = ~log_tpm, y = ~log_reads,
      type = "scatter", mode = "markers",
      color = ~sample_id, text = hover_text, hoverinfo = "text",
      marker = list(size = 6, opacity = 0.55)
    ) |>
      plotly::layout(
        xaxis = list(title = "log10(TPM + 1)"),
        yaxis = list(title = "log10(NumReads + 1)"),
        legend = list(orientation = "h", x = 0, y = 1.12)
      )
  }

  output$align_qc_mapping_plot <- plotly::renderPlotly(align_qc_mapping_plot_obj())
  output$align_qc_assignment_plot <- plotly::renderPlotly(align_qc_assignment_plot_obj())
  output$align_qc_region_plot <- plotly::renderPlotly(align_qc_region_plot_obj())
  output$align_qc_gene_body_plot <- plotly::renderPlotly(align_qc_gene_body_plot_obj())
  output$pseudo_qc_rate_plot <- plotly::renderPlotly(pseudo_qc_rate_plot_obj())
  output$pseudo_qc_tpm_plot <- plotly::renderPlotly(pseudo_qc_tpm_plot_obj())
  output$pseudo_qc_detected_plot <- plotly::renderPlotly(pseudo_qc_detected_plot_obj())
  output$pseudo_qc_scatter_plot <- plotly::renderPlotly(pseudo_qc_scatter_plot_obj())

  output$download_align_qc_summary <- downloadHandler(
    filename = function() paste0("align_qc_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(f) write.csv(align_qc_summary(selected_result_dir()), f, row.names = FALSE)
  )

  output$download_align_qc_alerts <- csv_download(
    "align_qc_alerts",
    function() align_qc_alerts(selected_result_dir())
  )

  output$download_pseudo_qc_summary <- downloadHandler(
    filename = function() paste0("pseudo_qc_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(f) write.csv(pseudo_qc_summary(selected_result_dir()), f, row.names = FALSE)
  )

  output$download_pseudo_quant <- downloadHandler(
    filename = function() paste0("pseudo_quantification_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(f) write.csv(pseudo_qc_quant_table(selected_result_dir()), f, row.names = FALSE)
  )

  output$download_pseudo_qc_alerts <- csv_download(
    "pseudo_qc_alerts",
    function() pseudo_qc_alerts(selected_result_dir())
  )

  output$download_align_qc_mapping_plot <- plotly_download(
    "align_qc_mapping",
    align_qc_mapping_plot_obj
  )
  output$download_align_qc_assignment_plot <- plotly_download(
    "align_qc_assignment",
    align_qc_assignment_plot_obj
  )
  output$download_align_qc_region_plot <- plotly_download(
    "align_qc_region",
    align_qc_region_plot_obj
  )
  output$download_align_qc_gene_body_plot <- plotly_download(
    "align_qc_gene_body",
    align_qc_gene_body_plot_obj
  )
  output$download_pseudo_qc_rate_plot <- plotly_download(
    "pseudo_qc_rate",
    pseudo_qc_rate_plot_obj
  )
  output$download_pseudo_qc_tpm_plot <- plotly_download(
    "pseudo_qc_tpm",
    pseudo_qc_tpm_plot_obj
  )
  output$download_pseudo_qc_detected_plot <- plotly_download(
    "pseudo_qc_detected",
    pseudo_qc_detected_plot_obj
  )
  output$download_pseudo_qc_scatter_plot <- plotly_download(
    "pseudo_qc_scatter",
    pseudo_qc_scatter_plot_obj
  )

  output$download_qc_alerts <- downloadHandler(
    filename = function() paste0("qc_alerts_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(f) {
      a <- align_qc_alerts(selected_result_dir())
      p <- pseudo_qc_alerts(selected_result_dir())
      if (has_real_rows(a)) a$bloque <- "alineamiento"
      if (has_real_rows(p)) p$bloque <- "pseudoalineamiento"
      out <- if (has_real_rows(a) && has_real_rows(p)) rbind(a, p)
             else if (has_real_rows(a)) a
             else if (has_real_rows(p)) p
             else message_df("No se detectaron alertas automaticas con los umbrales actuales.")
      write.csv(out, f, row.names = FALSE)
    }
  )
  # ── FIN: Server QC adicional de alineamiento y pseudoalineamiento ───────
  
  # Observador: cuando el usuario elige una carpeta desde el selector de archivos,
  # actualizar el select de ejecuciones guardadas.
  observeEvent(input$select_result_dir_btn, {
    req(input$select_result_dir_btn)
    sel <- tryCatch(shinyFiles::parseDirPath(roots_results, input$select_result_dir_btn), error=function(e) NULL)
    if (!is.null(sel)) {
      sel <- as.character(sel)
      if (nzchar(sel) && dir.exists(sel)) {
        results_refresh(Sys.time())
        updateSelectInput(session, "selected_result_dir",
                          choices = result_choices(outputs_dir), selected = sel)
      }
    }
  })

}

app <- shinyApp(ui, server)
if (sys.nframe() == 0L) {
  shiny::runApp(app, launch.browser = TRUE)
} else {
  app
}
