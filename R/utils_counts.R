#' utils_counts.R
#' Carga de matrices de conteos generadas por el workflow (bowtie2, salmon, kallisto).

#' Combina quant.sf / abundance.tsv en una matriz de conteos (genes x muestras)
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

#' Carga la matriz count_matrix.tsv generada por featureCounts (bowtie2)
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

#' Carga la matriz de conteos generada por el workflow para una tool dada.
#' Si esta disponible, usa tximport para salmon/kallisto.
load_counts_from_workflow <- function(output_dir, tool) {
  if (tool == "bowtie2") {
    f <- file.path(output_dir, "04_counts", "count_matrix.tsv")
    return(load_count_matrix_tsv(f))
  }
  if (HAS_TXIMPORT) {
    aln_dir <- file.path(output_dir, "03_alignments", tool)
    if (!dir.exists(aln_dir)) return(NULL)
    sdirs <- list.dirs(aln_dir, recursive = FALSE, full.names = TRUE)
    if (length(sdirs) == 0) return(NULL)
    qfiles <- if (tool == "salmon") file.path(sdirs, "quant.sf") else file.path(sdirs, "abundance.h5")
    if (!all(file.exists(qfiles))) qfiles <- file.path(sdirs, "abundance.tsv")
    valid <- file.exists(qfiles)
    if (!any(valid)) return(NULL)
    qfiles <- qfiles[valid]
    names(qfiles) <- basename(sdirs[valid])
    txi <- tryCatch(
      tximport::tximport(qfiles, type = tool, ignoreTxVersion = TRUE),
      error = function(e) NULL
    )
    if (is.null(txi)) return(load_quant_counts(qfiles, tool))
    return(round(txi$counts))
  }
  if (tool %in% c("salmon", "kallisto")) {
    aln_dir <- file.path(output_dir, "03_alignments", tool)
    if (!dir.exists(aln_dir)) return(NULL)
    sdirs <- list.dirs(aln_dir, recursive = FALSE, full.names = TRUE)
    if (!length(sdirs)) return(NULL)
    qfiles <- if (tool == "salmon") file.path(sdirs, "quant.sf") else file.path(sdirs, "abundance.tsv")
    valid <- file.exists(qfiles)
    if (!any(valid)) return(NULL)
    qfiles <- qfiles[valid]
    names(qfiles) <- basename(sdirs[valid])
    return(load_quant_counts(qfiles, tool))
  }
  NULL
}
