#' utils_qc.R
#' QC adicional: resumenes y alertas de alineamiento y pseudoalineamiento.

#' Convierte un porcentaje a fraccion (divide entre 100 si |x| > 1)
rate_fraction <- function(x) {
  x <- num_or_na(x)
  ifelse(!is.na(x) & abs(x) > 1, x / 100, x)
}

#' Busca la primera columna que casa con alguno de los patrones (regex)
find_metric_col <- function(df, patterns) {
  if (is.null(df) || !length(names(df))) return("")
  for (pat in patterns) {
    hit <- grep(pat, names(df), ignore.case = TRUE, value = TRUE)
    if (length(hit)) return(hit[1])
  }
  ""
}

#' Devuelve un data.frame con una unica columna Mensaje (para mostrar avisos)
message_df <- function(msg) {
  data.frame(Mensaje = msg, stringsAsFactors = FALSE, check.names = FALSE)
}

#' TRUE si el data.frame tiene filas reales (no es un mensaje placeholder)
has_real_rows <- function(df) {
  !is.null(df) && nrow(df) > 0 && !"Mensaje" %in% names(df)
}

#' Plot plotly vacio con un mensaje centrado
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

#' Detecta la columna que identifica la muestra (Sample/Muestra/sample_id)
sample_column <- function(df) {
  hit <- intersect(c("Sample", "sample", "Muestra", "sample_id"), names(df))
  if (!length(hit) || is.na(hit[1])) "" else hit[1]
}

#' Elimina filas R1/R2 duplicadas a nivel de muestra
clean_result_samples <- function(df) {
  sc <- sample_column(df)
  if (!nzchar(sc)) return(df)
  df[!grepl("(_1|_2|_R1|_R2)$", df[[sc]]), , drop = FALSE]
}

#' Resumen de asignacion genica desde multiqc_featurecounts.txt
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

#' Resumen QC de alineamiento clasico (bowtie2 + featureCounts)
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

  # Buscamos *rate* de multimapping de forma especifica para no capturar
  # columnas de conteo crudo (p.ej. paired_aligned_multi). Si no encontramos
  # un match explicito como rate/porcentaje, dejamos multimapping_rate = NA.
  multi_col <- find_metric_col(source, c(
    "multi_rate", "multimapping_rate", "pct_multi", "percent_multi",
    "aligned.*multi.*rate", "multi.*percent"
  ))
  if (nzchar(multi_col)) {
    raw_multi <- num_or_na(source[[multi_col]])
    # Heuristica conservadora: si el nombre termina en _pct o empieza con percent_
    # o la mediana es > 1, lo tratamos como porcentaje (0-100), si no fraccion.
    is_pct <- grepl("(_pct$|^percent_|percent$|_percent_|pct_)", multi_col, ignore.case = TRUE)
    med <- suppressWarnings(median(raw_multi, na.rm = TRUE))
    if (is.finite(med) && (is_pct || med > 1)) {
      multimapping_rate <- raw_multi / 100
    } else {
      multimapping_rate <- raw_multi
    }
  } else {
    # Si no hay rate explicito, no inventamos uno a partir de conteos: deja NA.
    multimapping_rate <- rep(NA_real_, length(samples))
  }

  # Si conocemos ambos rates podemos calcular unique_reads.
  # Si multimapping_rate es NA, NO etiquetamos como unique_reads (queda NA)
  # y exponemos mapped_reads = total_reads * mapping_rate como columna extra.
  unique_reads <- ifelse(!is.na(total_reads) & !is.na(mapping_rate) & !is.na(multimapping_rate),
                         total_reads * pmax(mapping_rate - multimapping_rate, 0),
                         NA_real_)
  mapped_reads <- ifelse(!is.na(total_reads) & !is.na(mapping_rate),
                         total_reads * mapping_rate, NA_real_)
  multimapped_reads <- ifelse(!is.na(total_reads) & !is.na(multimapping_rate),
                              total_reads * multimapping_rate, NA_real_)
  unmapped_reads <- ifelse(!is.na(total_reads) & !is.na(mapping_rate),
                           total_reads * pmax(1 - mapping_rate, 0), NA_real_)

  out <- data.frame(
    sample_id = samples,
    total_reads = total_reads,
    unique_reads = unique_reads,
    mapped_reads = mapped_reads,
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

#' Lee y concatena todos los quant.sf de un run salmon
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

#' Lee y concatena todos los abundance.tsv de un run kallisto (renombrando columnas a estilo salmon)
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

#' Fraccion de lecturas asignadas a rRNA por muestra.
#'
#' Por que es una metrica de primer nivel y no un detalle (ver
#' docs/REVISION_ESTADISTICA.md, B10): el rRNA supera el 85 % del RNA celular
#' procariota y la depleccion es imperfecta — hasta un 50 % de lecturas
#' residuales de rRNA es un resultado habitual. Con normalizacion por
#' composicion (TMM/RLE), un arrastre de rRNA que varie entre muestras sesga
#' TODOS los factores de tamano, no solo los genes de rRNA. Por eso lo relevante
#' no es el valor absoluto sino la VARIACION entre muestras.
#'
#' @param counts matriz de conteos (genes x muestras)
#' @param rrna_ids identificadores de rRNA de la anotacion. Si es NULL o no casa
#'   con la matriz, se cae a la heuristica de nombres `infer_rna_type()`.
#' @return list(table = data.frame(sample_id, rrna_reads, total_reads, frac),
#'   source, n_rrna_genes, spread, alert)
rrna_fraction_per_sample <- function(counts, rrna_ids = NULL) {
  if (is.null(counts) || !nrow(counts) || !ncol(counts)) return(NULL)
  cm <- as.matrix(counts)
  ids <- rownames(cm)
  if (is.null(ids)) return(NULL)

  is_rrna <- rep(FALSE, nrow(cm))
  source <- "ninguna"
  if (!is.null(rrna_ids) && length(rrna_ids)) {
    is_rrna <- ids %in% rrna_ids
    if (any(is_rrna)) source <- "anotacion"
  }
  if (!any(is_rrna)) {
    is_rrna <- infer_rna_type(ids) == "rRNA"
    if (any(is_rrna)) source <- "heuristica de nombres"
  }
  if (!any(is_rrna)) {
    return(list(table = NULL, source = "ninguna", n_rrna_genes = 0L,
                spread = NA_real_, alert = NULL))
  }

  total <- colSums(cm, na.rm = TRUE)
  rr <- colSums(cm[is_rrna, , drop = FALSE], na.rm = TRUE)
  frac <- ifelse(total > 0, rr / total, NA_real_)
  df <- data.frame(sample_id = colnames(cm), rrna_reads = as.numeric(rr),
                   total_reads = as.numeric(total), frac = as.numeric(frac),
                   stringsAsFactors = FALSE)
  spread <- suppressWarnings(diff(range(df$frac, na.rm = TRUE)))
  alert <- NULL
  if (is.finite(spread) && spread > 0.10) {
    alert <- paste0(
      "El porcentaje de rRNA varia ", round(100 * spread, 1),
      " puntos entre muestras (de ", round(100 * min(df$frac, na.rm = TRUE), 1),
      " % a ", round(100 * max(df$frac, na.rm = TRUE), 1),
      " %). Un arrastre desigual de rRNA sesga los factores de tamano de la ",
      "normalizacion por composicion, asi que afecta a todos los genes. ",
      "Considera excluir el rRNA antes del analisis diferencial.")
  } else if (is.finite(max(df$frac, na.rm = TRUE)) && max(df$frac, na.rm = TRUE) > 0.3) {
    alert <- paste0(
      "Hay muestras con mas del ", round(100 * max(df$frac, na.rm = TRUE), 0),
      " % de lecturas en rRNA. La depleccion no ha sido eficaz, lo que reduce ",
      "la profundidad efectiva sobre el mRNA.")
  }
  list(table = df, source = source, n_rrna_genes = sum(is_rrna),
       spread = spread, alert = alert)
}

#' Heuristica para inferir el tipo de RNA (mRNA, rRNA, tRNA, sRNA, ...)
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

#' Tabla de cuantificacion concatenada (salmon o kallisto) con columnas comunes
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

#' Resumen QC del pseudoalineamiento por muestra
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
    # Mediana excluyendo TPM = 0 (transcritos NO detectados).
    # Util porque median_tpm a menudo es 0 cuando hay muchos features no detectados.
    tpm_median_det <- aggregate(
      q$TPM,
      list(sample_id = q$sample_id),
      function(x) {
        x <- x[!is.na(x) & x > 0]
        if (!length(x)) NA_real_ else median(x)
      }
    )
    names(tpm_median_det)[2] <- "median_tpm_detected"
    out <- Reduce(function(x, y) merge(x, y, by = "sample_id", all = TRUE),
                  list(out, detected, near_zero, tpm_median, tpm_median_det))
  }

  counts <- tryCatch(load_counts_from_workflow(out_dir, "salmon", annotation_file = annotation_file_for_run(out_dir)), error = function(e) NULL)
  if (is.null(counts)) counts <- tryCatch(load_counts_from_workflow(out_dir, "kallisto", annotation_file = annotation_file_for_run(out_dir)), error = function(e) NULL)
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

#' Alertas heuristicas sobre el alineamiento clasico
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
        x$total_reads[i] > med_reads * 1.5)
      add("info", "total_reads", round(x$total_reads[i]), "Numero de lecturas muy superior a la mediana del conjunto.")
  }
  if (!length(alerts)) return(message_df("No se detectaron alertas automaticas con los umbrales actuales."))
  do.call(rbind, alerts)
}

#' Alertas heuristicas sobre el pseudoalineamiento
#' Nota (auditoria): Para el aviso "tpm_distribution_shift_warning" usamos
#' median_tpm_detected cuando esta disponible, porque median_tpm tiende a 0
#' con muchos transcritos no detectados y no discrimina bien shifts reales
#' entre muestras. Si median_tpm_detected no esta (sin tabla quant), caemos
#' de vuelta a median_tpm.
pseudo_qc_alerts <- function(out_dir, thresholds = qc_thresholds) {
  x <- pseudo_qc_summary(out_dir)
  if (!has_real_rows(x)) return(x)
  alerts <- list()
  med_processed <- median(x$fragments_processed, na.rm = TRUE)
  med_detected <- median(x$transcripts_detected, na.rm = TRUE)
  use_det <- "median_tpm_detected" %in% names(x) && any(!is.na(x$median_tpm_detected))
  tpm_vec <- if (use_det) x$median_tpm_detected else x$median_tpm
  med_tpm <- median(tpm_vec, na.rm = TRUE)
  mad_tpm <- mad(tpm_vec, na.rm = TRUE)
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
    tpm_i <- if (use_det) x$median_tpm_detected[i] else x$median_tpm[i]
    metric_name <- if (use_det) "median_tpm_detected" else "median_tpm"
    if (is.finite(mad_tpm) && mad_tpm > 0 && !is.na(tpm_i) &&
        abs(tpm_i - med_tpm) / mad_tpm > thresholds$tpm_distribution_shift_warning)
      add("warning", metric_name, round(tpm_i, 3), "Distribucion de TPM muy diferente al resto.")
  }
  if (!length(alerts)) return(message_df("No se detectaron alertas automaticas con los umbrales actuales."))
  do.call(rbind, alerts)
}
