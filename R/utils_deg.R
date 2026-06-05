#' utils_deg.R
#' Funciones puras (sin Shiny) para analisis de expresion diferencial (DEG).
#' Soporta tres motores: DESeq2, edgeR (glmQLFit + glmQLFTest) y limma-voom.
#' Todas devuelven un data.frame con columnas estandar:
#'   gene, baseMean, log2FC, lfcSE, stat, pvalue, padj

#' Convierte porcentaje o fraccion (>1 o <=1) — no se usa aqui, helper similar.
num_safe <- function(x) suppressWarnings(as.numeric(x))

#' Construye un data.frame "estandar" DEG con NA para columnas que falten
deg_empty_row <- function(gene = character(0)) {
  data.frame(
    gene = gene,
    baseMean = NA_real_,
    log2FC = NA_real_,
    lfcSE = NA_real_,
    stat = NA_real_,
    pvalue = NA_real_,
    padj = NA_real_,
    stringsAsFactors = FALSE
  )[seq_along(gene), , drop = FALSE]
}

#' Valida un samplesheet contra la matriz de conteos.
#' Devuelve list(ok = TRUE/FALSE, errors = character).
validate_samplesheet <- function(df, samples_in_counts) {
  errors <- character(0)
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) {
    return(list(ok = FALSE, errors = "El samplesheet esta vacio."))
  }
  required <- c("sample_id", "condition")
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols)) {
    errors <- c(errors, paste0("Faltan columnas obligatorias: ",
                               paste(missing_cols, collapse = ", "), "."))
  }
  if ("sample_id" %in% names(df)) {
    if (any(duplicated(df$sample_id))) {
      errors <- c(errors, "Hay sample_id duplicados en el samplesheet.")
    }
    if (length(samples_in_counts)) {
      not_in_counts <- setdiff(df$sample_id, samples_in_counts)
      if (length(not_in_counts)) {
        errors <- c(errors, paste0(
          "Muestras del samplesheet no presentes en la matriz: ",
          paste(not_in_counts, collapse = ", "), "."
        ))
      }
    }
  }
  if ("condition" %in% names(df)) {
    lvl <- unique(df$condition[!is.na(df$condition) & nzchar(as.character(df$condition))])
    if (length(lvl) < 2) {
      errors <- c(errors, "La columna condition debe tener al menos 2 niveles distintos.")
    }
  }
  list(ok = !length(errors), errors = errors)
}

#' Reordena las columnas de counts para que coincidan con meta$sample_id.
#' Avisa (mediante atributo "warnings") cuando faltan muestras.
align_counts_to_metadata <- function(counts, meta) {
  if (is.null(counts) || !length(counts)) return(list(counts = NULL, meta = meta, warnings = "Matriz vacia."))
  if (is.null(meta) || !"sample_id" %in% names(meta)) {
    return(list(counts = counts, meta = meta, warnings = "Metadatos sin sample_id."))
  }
  shared <- intersect(meta$sample_id, colnames(counts))
  warns <- character(0)
  if (length(shared) < ncol(counts)) {
    warns <- c(warns, paste0(
      "Muestras en counts no presentes en samplesheet: ",
      paste(setdiff(colnames(counts), shared), collapse = ", ")
    ))
  }
  if (length(shared) < nrow(meta)) {
    warns <- c(warns, paste0(
      "Muestras en samplesheet no presentes en counts: ",
      paste(setdiff(meta$sample_id, shared), collapse = ", ")
    ))
  }
  meta_sub <- meta[match(shared, meta$sample_id), , drop = FALSE]
  counts_sub <- counts[, shared, drop = FALSE]
  list(counts = counts_sub, meta = meta_sub, warnings = warns)
}

#' Filtra genes con poca expresion: al menos `min_count` lecturas en >= `min_samples`.
#' Si min_samples es NULL, usa el tamano del grupo mas pequeno.
prefilter_counts <- function(counts, min_count = 10, min_samples = NULL) {
  if (is.null(counts) || !nrow(counts)) return(counts)
  if (is.null(min_samples) || !is.finite(min_samples) || min_samples < 1) {
    min_samples <- max(2L, floor(ncol(counts) / 2))
  }
  keep <- rowSums(as.matrix(counts) >= min_count) >= min_samples
  counts[keep, , drop = FALSE]
}

#' Construye matriz de diseno con opcional batch
build_design <- function(meta, ref_level = NULL, batch = NULL) {
  meta$condition <- as.factor(meta$condition)
  if (!is.null(ref_level) && ref_level %in% levels(meta$condition)) {
    meta$condition <- stats::relevel(meta$condition, ref = ref_level)
  }
  if (!is.null(batch) && nzchar(batch) && batch %in% names(meta)) {
    meta[[batch]] <- as.factor(meta[[batch]])
    formula_obj <- stats::as.formula(paste0("~ ", batch, " + condition"))
  } else {
    formula_obj <- stats::as.formula("~ condition")
  }
  list(meta = meta, formula = formula_obj)
}

#' Motor DESeq2
run_deg_deseq2 <- function(counts, meta, ref_level = NULL, batch = NULL) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    return(list(table = NULL, error = "DESeq2 no esta instalado."))
  }
  d <- build_design(meta, ref_level, batch)
  out <- tryCatch({
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = round(as.matrix(counts)),
      colData = d$meta,
      design = d$formula
    )
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    res <- DESeq2::results(dds)
    df <- as.data.frame(res)
    df$gene <- rownames(df)
    data.frame(
      gene = df$gene,
      baseMean = df$baseMean,
      log2FC = df$log2FoldChange,
      lfcSE = df$lfcSE,
      stat = df$stat,
      pvalue = df$pvalue,
      padj = df$padj,
      stringsAsFactors = FALSE
    )
  }, error = function(e) e)
  if (inherits(out, "error")) return(list(table = NULL, error = conditionMessage(out)))
  list(table = out, error = NULL)
}

#' Motor edgeR: glmQLFit + glmQLFTest
run_deg_edger <- function(counts, meta, ref_level = NULL, batch = NULL) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    return(list(table = NULL, error = "edgeR no esta instalado."))
  }
  d <- build_design(meta, ref_level, batch)
  out <- tryCatch({
    y <- edgeR::DGEList(counts = round(as.matrix(counts)))
    y <- edgeR::calcNormFactors(y)
    design <- stats::model.matrix(d$formula, data = d$meta)
    y <- edgeR::estimateDisp(y, design)
    fit <- edgeR::glmQLFit(y, design)
    # Coeficiente: el ultimo "condition*" en el diseno
    coef_names <- colnames(design)
    cond_coefs <- grep("^condition", coef_names, value = TRUE)
    coef_test <- if (length(cond_coefs)) cond_coefs[length(cond_coefs)] else coef_names[ncol(design)]
    qlf <- edgeR::glmQLFTest(fit, coef = coef_test)
    tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
    # AveLogCPM ~ baseMean estandarizado en otra escala. Para mantener la
    # interfaz, devolvemos baseMean = 2^AveLogCPM (CPM medio aproximado).
    base_mean <- if ("logCPM" %in% names(tt)) 2 ^ tt$logCPM
                 else if ("AveLogCPM" %in% names(tt)) 2 ^ tt$AveLogCPM
                 else rep(NA_real_, nrow(tt))
    data.frame(
      gene = rownames(tt),
      baseMean = base_mean,
      log2FC = tt$logFC,
      lfcSE = NA_real_,
      stat = tt$F %||% tt$LR %||% rep(NA_real_, nrow(tt)),
      pvalue = tt$PValue,
      padj = tt$FDR,
      stringsAsFactors = FALSE
    )
  }, error = function(e) e)
  if (inherits(out, "error")) return(list(table = NULL, error = conditionMessage(out)))
  list(table = out, error = NULL)
}

#' Motor limma-voom
run_deg_limma <- function(counts, meta, ref_level = NULL, batch = NULL) {
  if (!requireNamespace("limma", quietly = TRUE) ||
      !requireNamespace("edgeR", quietly = TRUE)) {
    return(list(table = NULL, error = "Se requieren limma y edgeR para limma-voom."))
  }
  d <- build_design(meta, ref_level, batch)
  out <- tryCatch({
    y <- edgeR::DGEList(counts = round(as.matrix(counts)))
    y <- edgeR::calcNormFactors(y)
    design <- stats::model.matrix(d$formula, data = d$meta)
    v <- limma::voom(y, design)
    fit <- limma::lmFit(v, design)
    fit <- limma::eBayes(fit)
    coef_names <- colnames(design)
    cond_coefs <- grep("^condition", coef_names, value = TRUE)
    coef_test <- if (length(cond_coefs)) cond_coefs[length(cond_coefs)] else coef_names[ncol(design)]
    tt <- limma::topTable(fit, coef = coef_test, number = Inf, sort.by = "none")
    base_mean <- if ("AveExpr" %in% names(tt)) 2 ^ tt$AveExpr else rep(NA_real_, nrow(tt))
    data.frame(
      gene = rownames(tt),
      baseMean = base_mean,
      log2FC = tt$logFC,
      lfcSE = NA_real_,
      stat = tt$t %||% rep(NA_real_, nrow(tt)),
      pvalue = tt$P.Value,
      padj = tt$adj.P.Val,
      stringsAsFactors = FALSE
    )
  }, error = function(e) e)
  if (inherits(out, "error")) return(list(table = NULL, error = conditionMessage(out)))
  list(table = out, error = NULL)
}

#' Dispatcher: corre el motor solicitado.
#' Devuelve list(table = data.frame_estandar, error = NULL_o_mensaje, method = method).
run_deg <- function(counts, meta, method = c("DESeq2", "edgeR", "limma-voom"),
                    ref_level = NULL, batch = NULL) {
  method <- match.arg(method)
  res <- switch(method,
    "DESeq2"     = run_deg_deseq2(counts, meta, ref_level, batch),
    "edgeR"      = run_deg_edger(counts, meta, ref_level, batch),
    "limma-voom" = run_deg_limma(counts, meta, ref_level, batch)
  )
  res$method <- method
  res
}

#' Aplica filtros (FDR, |log2FC|, baseMean) a una tabla DEG estandar.
apply_deg_filters <- function(deg_df, fdr = 0.05, abs_log2fc = 1, base_mean = 0) {
  if (is.null(deg_df) || !nrow(deg_df)) return(deg_df)
  keep <- rep(TRUE, nrow(deg_df))
  if (!is.null(fdr) && is.finite(fdr)) {
    keep <- keep & !is.na(deg_df$padj) & deg_df$padj <= fdr
  }
  if (!is.null(abs_log2fc) && is.finite(abs_log2fc)) {
    keep <- keep & !is.na(deg_df$log2FC) & abs(deg_df$log2FC) >= abs_log2fc
  }
  if (!is.null(base_mean) && is.finite(base_mean)) {
    keep <- keep & !is.na(deg_df$baseMean) & deg_df$baseMean >= base_mean
  }
  deg_df[keep, , drop = FALSE]
}

#' Transformacion estabilizadora de varianza para visualizacion.
#' Usa DESeq2::vst si nrow >= 1000 (mas rapido), rlog si menos (mas suave).
vst_or_rlog <- function(counts, meta, blind = TRUE) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    # Fallback: log2(counts + 1) normalizado por libreria
    cm <- as.matrix(counts)
    libs <- colSums(cm)
    libs[libs == 0] <- 1
    cpm <- t(t(cm) / libs) * 1e6
    return(log2(cpm + 1))
  }
  out <- tryCatch({
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = round(as.matrix(counts)),
      colData = meta,
      design = ~ 1
    )
    if (nrow(counts) >= 1000) {
      SummarizedExperiment::assay(DESeq2::vst(dds, blind = blind))
    } else {
      SummarizedExperiment::assay(DESeq2::rlog(dds, blind = blind))
    }
  }, error = function(e) {
    cm <- as.matrix(counts)
    libs <- colSums(cm); libs[libs == 0] <- 1
    cpm <- t(t(cm) / libs) * 1e6
    log2(cpm + 1)
  })
  out
}

#' Datos para PCA a partir de una matriz transformada (vst/rlog)
pca_data <- function(vst_mat, meta, ntop = 500) {
  if (is.null(vst_mat) || !nrow(vst_mat) || !ncol(vst_mat)) return(NULL)
  rv <- apply(vst_mat, 1, stats::var, na.rm = TRUE)
  ntop <- min(ntop, sum(is.finite(rv) & rv > 0))
  if (!ntop) return(NULL)
  idx <- order(rv, decreasing = TRUE)[seq_len(ntop)]
  pca <- stats::prcomp(t(vst_mat[idx, , drop = FALSE]), scale. = FALSE)
  varexp <- (pca$sdev ^ 2) / sum(pca$sdev ^ 2)
  df <- data.frame(
    sample_id = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )
  if (!is.null(meta) && "sample_id" %in% names(meta)) {
    df <- merge(df, meta, by = "sample_id", all.x = TRUE, sort = FALSE)
  }
  attr(df, "var_explained") <- varexp[1:min(2, length(varexp))]
  df
}

#' Matriz de distancias entre muestras (euclidiana sobre vst)
sample_distance_matrix <- function(vst_mat) {
  if (is.null(vst_mat) || !ncol(vst_mat)) return(NULL)
  as.matrix(stats::dist(t(vst_mat)))
}

#' Top N genes por varianza (sobre matriz transformada)
top_var_genes <- function(vst_mat, n = 30) {
  if (is.null(vst_mat) || !nrow(vst_mat)) return(character(0))
  rv <- apply(vst_mat, 1, stats::var, na.rm = TRUE)
  n <- min(n, sum(is.finite(rv)))
  if (!n) return(character(0))
  rownames(vst_mat)[order(rv, decreasing = TRUE)[seq_len(n)]]
}
