#' utils_deg.R
#' Funciones puras (sin Shiny) para analisis de expresion diferencial (DEG).
#' Soporta tres motores: DESeq2, edgeR v4 (glmQLFit + glmQLFTest/glmTreat) y
#' limma-voom. Todas devuelven un data.frame con columnas estandar:
#'   gene, baseMean, log2FC, log2FC_shrunk, lfcSE, stat, pvalue, padj
#'
#' Dos criterios estadisticos que atraviesan todo el archivo (ver
#' docs/REVISION_ESTADISTICA.md):
#'
#'   1. El umbral de |log2FC| se aplica DENTRO del test (lfcThreshold en DESeq2,
#'      glmTreat en edgeR, treat en limma) y nunca como filtro posterior sobre
#'      una lista ya corregida por FDR. Filtrar a posteriori reduce la FDR real
#'      en una cantidad desconocida, de modo que la lista deja de tener el
#'      porcentaje de falsos positivos que declara (McCarthy y Smyth, 2009).
#'   2. El nivel de FDR que va a usar quien llama se pasa a los motores, porque
#'      DESeq2 lo necesita en `alpha` para calibrar su filtrado independiente.

#' Convierte porcentaje o fraccion (>1 o <=1) — no se usa aqui, helper similar.
num_safe <- function(x) suppressWarnings(as.numeric(x))

#' Construye un data.frame "estandar" DEG con NA para columnas que falten
deg_empty_row <- function(gene = character(0)) {
  data.frame(
    gene = gene,
    baseMean = NA_real_,
    log2FC = NA_real_,
    log2FC_shrunk = NA_real_,
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

#' Tamano del grupo mas pequeno de `group`. Solo cuando no hay informacion de
#' grupo utilizable cae al fallback historico (la mitad de las muestras).
smallest_group_size <- function(group, n_samples) {
  if (!is.null(group) && length(group) == n_samples) {
    g <- as.character(group)
    g <- g[!is.na(g) & nzchar(g)]
    if (length(g)) return(max(2L, as.integer(min(table(g)))))
  }
  max(2L, floor(n_samples / 2))
}

#' Prefiltrado de genes de baja expresion.
#'
#' `mode = "auto"` (recomendado) delega en `edgeR::filterByExpr()`, que decide
#' en cuantas muestras debe superarse el umbral a partir del tamano del grupo
#' MAS PEQUENO. `mode = "manual"` mantiene el criterio explicito "al menos
#' `min_count` lecturas en >= `min_samples` muestras".
#'
#' Ojo con el fallback de `min_samples`: la version anterior usaba la mitad de
#' las muestras TOTALES, lo que en un diseno desequilibrado (3 control vs 9
#' tratados -> min_samples = 6) descarta precisamente los genes expresados solo
#' en el grupo pequeno. Ahora se deriva del grupo menor via
#' `smallest_group_size()`.
#'
#' Devuelve la matriz filtrada con un atributo "prefilter" que resume la
#' decision tomada (modo, umbrales efectivos, genes antes/despues), para que la
#' interfaz pueda mostrarla.
prefilter_counts <- function(counts, min_count = 10, min_samples = NULL,
                             mode = c("auto", "manual"),
                             design = NULL, group = NULL) {
  mode <- match.arg(mode)
  if (is.null(counts) || !nrow(counts)) return(counts)
  cm <- as.matrix(counts)
  n_before <- nrow(cm)

  if (identical(mode, "auto") && requireNamespace("edgeR", quietly = TRUE)) {
    keep <- tryCatch({
      y <- edgeR::DGEList(counts = round(cm), group = group)
      if (!is.null(design)) edgeR::filterByExpr(y, design = design)
      else                  edgeR::filterByExpr(y, group = group)
    }, error = function(e) NULL)
    if (!is.null(keep) && length(keep) == n_before) {
      out <- cm[keep, , drop = FALSE]
      attr(out, "prefilter") <- list(
        mode = "filterByExpr", min_count = NA_real_, min_samples = NA_real_,
        n_before = n_before, n_after = nrow(out)
      )
      return(out)
    }
  }

  if (is.null(min_samples) || !is.finite(min_samples) || min_samples < 1) {
    min_samples <- smallest_group_size(group, ncol(cm))
  }
  keep <- rowSums(cm >= min_count) >= min_samples
  out <- cm[keep, , drop = FALSE]
  attr(out, "prefilter") <- list(
    mode = "manual", min_count = min_count, min_samples = min_samples,
    n_before = n_before, n_after = nrow(out)
  )
  out
}

#' Construye matriz de diseno con opcional batch.
#' Devuelve tambien los niveles del factor y el de referencia, que hacen falta
#' para etiquetar el contraste realmente testeado.
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
  lvls <- levels(meta$condition)
  list(meta = meta, formula = formula_obj, levels = lvls, ref = lvls[1])
}

#' Coeficiente de `condition` que se va a testear: el ultimo del diseno.
#'
#' Con dos niveles es el unico posible. Con tres o mas es UNA de las
#' comparaciones posibles (ultimo nivel vs referencia), asi que quien llama
#' tiene que etiquetarla explicitamente y avisar al usuario.
pick_condition_coef <- function(coef_names) {
  cond <- grep("^condition", coef_names, value = TRUE)
  if (length(cond)) cond[length(cond)] else coef_names[length(coef_names)]
}

#' Etiqueta legible del contraste realmente testeado.
#' DESeq2 nombra sus coeficientes "condition_trt_vs_ctrl"; `model.matrix`
#' (edgeR/limma) los nombra "conditiontrt". En ambos casos el denominador es el
#' nivel de referencia del factor.
contrast_label <- function(coef_name, ref_level = NULL) {
  if (is.null(coef_name) || !length(coef_name) || !nzchar(coef_name)) return(NA_character_)
  if (grepl("_vs_", coef_name, fixed = TRUE)) {
    parts <- strsplit(sub("^condition_", "", coef_name), "_vs_", fixed = TRUE)[[1]]
    if (length(parts) == 2L) return(paste0(parts[1], " vs ", parts[2]))
  }
  num <- sub("^condition", "", coef_name)
  if (!nzchar(num)) return(coef_name)
  paste0(num, " vs ", ref_level %||% "referencia")
}

#' edgeR >= 4 renombro `calcNormFactors` a `normLibSizes`. Soportamos ambas.
norm_lib_sizes <- function(y) {
  if ("normLibSizes" %in% getNamespaceExports("edgeR")) edgeR::normLibSizes(y)
  else edgeR::calcNormFactors(y)
}

#' TRUE si edgeR es >= 4.0.0. En v4 `glmQLFit()` estima la dispersion NB y la
#' cuasi-dispersion internamente (con devianzas corregidas por sesgo), asi que
#' `estimateDisp()` ya no hace falta; en v3 si.
edger_is_v4 <- function() {
  isTRUE(tryCatch(utils::packageVersion("edgeR") >= "4.0.0", error = function(e) FALSE))
}

#' Elige el mejor tipo de encogido de log2FC disponible.
#' `apeglm` es la recomendacion de la vinieta de DESeq2 (prior de colas anchas:
#' encoge el ruido sin aplastar los efectos grandes reales).
pick_shrink_type <- function() {
  if (requireNamespace("apeglm", quietly = TRUE)) return("apeglm")
  if (requireNamespace("ashr", quietly = TRUE)) return("ashr")
  "normal"
}

#' Motor DESeq2.
#'
#' @param fdr Nivel de FDR objetivo. Se pasa a `results(alpha = fdr)` porque el
#'   filtrado independiente de DESeq2 elige el umbral de expresion que maximiza
#'   el numero de genes significativos A ESE NIVEL. Con el default (0,1) y un
#'   corte real de 0,05 el filtrado se optimiza para el umbral equivocado y se
#'   pierden descubrimientos.
#' @param lfc_threshold Umbral de |log2FC| aplicado dentro del test. 0 = test
#'   clasico contra H0: log2FC = 0.
#' @param shrink Si TRUE, anade `log2FC_shrunk` via `lfcShrink()`. El encogido
#'   cambia el log2FC pero NO los p-valores, asi que se testea con el estimador
#'   de maxima verosimilitud y se visualiza/ordena con el encogido.
run_deg_deseq2 <- function(counts, meta, ref_level = NULL, batch = NULL,
                           fdr = 0.05, lfc_threshold = 0, shrink = TRUE) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    return(list(table = NULL, error = "DESeq2 no esta instalado."))
  }
  d <- build_design(meta, ref_level, batch)
  info <- list(contrast = NA_character_, coef = NA_character_,
               n_levels = length(d$levels), shrink = "ninguno")
  out <- tryCatch({
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = round(as.matrix(counts)),
      colData = d$meta,
      design = d$formula
    )
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    coef_name <- pick_condition_coef(DESeq2::resultsNames(dds))

    res <- DESeq2::results(
      dds, name = coef_name, alpha = fdr,
      lfcThreshold = if (is.finite(lfc_threshold)) lfc_threshold else 0,
      altHypothesis = "greaterAbs"
    )
    df <- as.data.frame(res)

    lfc_shrunk <- rep(NA_real_, nrow(df))
    shrink_used <- "ninguno"
    if (isTRUE(shrink)) {
      type <- pick_shrink_type()
      shr <- tryCatch(
        DESeq2::lfcShrink(dds, coef = coef_name, type = type, quiet = TRUE),
        error = function(e) NULL
      )
      # lfcShrink devuelve las mismas filas y en el mismo orden que results().
      if (!is.null(shr) && nrow(shr) == nrow(df)) {
        lfc_shrunk <- as.data.frame(shr)$log2FoldChange
        shrink_used <- type
      }
    }

    tab <- data.frame(
      gene = rownames(df),
      baseMean = df$baseMean,
      log2FC = df$log2FoldChange,
      log2FC_shrunk = lfc_shrunk,
      lfcSE = df$lfcSE,
      stat = df$stat,
      pvalue = df$pvalue,
      padj = df$padj,
      stringsAsFactors = FALSE
    )
    list(table = tab, coef = coef_name, shrink = shrink_used)
  }, error = function(e) e)

  if (inherits(out, "error")) {
    return(c(list(table = NULL, error = conditionMessage(out)), info))
  }
  info$coef     <- out$coef
  info$contrast <- contrast_label(out$coef, d$ref)
  info$shrink   <- out$shrink
  c(list(table = out$table, error = NULL), info)
}

#' Motor edgeR, pipeline v4.
#'
#' Cambios respecto al pipeline v3 que habia antes:
#'   - `normLibSizes()` en lugar de `calcNormFactors()` (mismo calculo, nombre
#'     nuevo desde v4).
#'   - Sin `estimateDisp()`: en v4 `glmQLFit()` estima la dispersion NB y la
#'     cuasi-dispersion internamente usando devianzas corregidas por sesgo, lo
#'     que corrige la subestimacion historica de las cuasi-dispersiones en
#'     conteos pequenos (Chen et al., NAR 2025) — justo el regimen de un
#'     experimento con pocas replicas.
#'   - `glmTreat()` cuando hay umbral de fold-change, en lugar de filtrar la
#'     tabla despues.
#'
#' `fdr` y `shrink` se aceptan por uniformidad de la interfaz pero no se usan:
#' edgeR no tiene filtrado independiente que calibrar, y el `logFC` que reporta
#' ya lleva el encogido por conteos previos.
run_deg_edger <- function(counts, meta, ref_level = NULL, batch = NULL,
                          fdr = 0.05, lfc_threshold = 0, shrink = TRUE) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    return(list(table = NULL, error = "edgeR no esta instalado."))
  }
  d <- build_design(meta, ref_level, batch)
  info <- list(contrast = NA_character_, coef = NA_character_,
               n_levels = length(d$levels), shrink = "ninguno")
  out <- tryCatch({
    design <- stats::model.matrix(d$formula, data = d$meta)
    y <- edgeR::DGEList(counts = round(as.matrix(counts)), group = d$meta$condition)
    y <- norm_lib_sizes(y)
    fit <- if (edger_is_v4()) {
      edgeR::glmQLFit(y, design)
    } else {
      edgeR::glmQLFit(edgeR::estimateDisp(y, design), design)
    }
    coef_test <- pick_condition_coef(colnames(design))

    test <- if (is.finite(lfc_threshold) && lfc_threshold > 0) {
      edgeR::glmTreat(fit, coef = coef_test, lfc = lfc_threshold)
    } else {
      edgeR::glmQLFTest(fit, coef = coef_test)
    }
    tt <- edgeR::topTags(test, n = Inf, sort.by = "none")$table

    # AveLogCPM ~ baseMean estandarizado en otra escala. Para mantener la
    # interfaz, devolvemos baseMean = 2^AveLogCPM (CPM medio aproximado).
    base_mean <- if (!is.null(tt[["logCPM"]])) 2 ^ tt[["logCPM"]]
                 else if (!is.null(tt[["AveLogCPM"]])) 2 ^ tt[["AveLogCPM"]]
                 else rep(NA_real_, nrow(tt))
    # Ojo: `tt$F` haria partial matching con la columna FDR y colaria el FDR
    # como estadistico cuando glmTreat no devuelve F. Con [[ ]] el match es
    # exacto y devuelve NULL si la columna no existe.
    stat_col <- tt[["F"]] %||% tt[["LR"]] %||% rep(NA_real_, nrow(tt))

    tab <- data.frame(
      gene = rownames(tt),
      baseMean = base_mean,
      log2FC = tt[["logFC"]],
      log2FC_shrunk = NA_real_,
      lfcSE = NA_real_,
      stat = stat_col,
      pvalue = tt[["PValue"]],
      padj = tt[["FDR"]],
      stringsAsFactors = FALSE
    )
    list(table = tab, coef = coef_test)
  }, error = function(e) e)

  if (inherits(out, "error")) {
    return(c(list(table = NULL, error = conditionMessage(out)), info))
  }
  info$coef     <- out$coef
  info$contrast <- contrast_label(out$coef, d$ref)
  c(list(table = out$table, error = NULL), info)
}

#' Motor limma-voom.
#'
#' Dos endurecimientos respecto a la version anterior:
#'   - `voomWithQualityWeights()` en lugar de `voom()`: pondera a la baja las
#'     muestras de peor calidad en vez de descartarlas, lo que importa con pocas
#'     replicas y algun outlier. Si falla, cae a `voom()`.
#'   - `robust = TRUE` en la moderacion empirica-bayesiana, que protege frente a
#'     genes con varianza extrema.
#'
#' No se pasa `trend = TRUE` a proposito: la tendencia media-varianza ya la
#' modelan los pesos de precision de voom, y `trend = TRUE` es la via de
#' limma-trend sobre logCPM. Activar ambas contaria la tendencia dos veces.
run_deg_limma <- function(counts, meta, ref_level = NULL, batch = NULL,
                          fdr = 0.05, lfc_threshold = 0, shrink = TRUE) {
  if (!requireNamespace("limma", quietly = TRUE) ||
      !requireNamespace("edgeR", quietly = TRUE)) {
    return(list(table = NULL, error = "Se requieren limma y edgeR para limma-voom."))
  }
  d <- build_design(meta, ref_level, batch)
  info <- list(contrast = NA_character_, coef = NA_character_,
               n_levels = length(d$levels), shrink = "ninguno")
  out <- tryCatch({
    design <- stats::model.matrix(d$formula, data = d$meta)
    y <- edgeR::DGEList(counts = round(as.matrix(counts)), group = d$meta$condition)
    y <- norm_lib_sizes(y)
    v <- tryCatch(limma::voomWithQualityWeights(y, design),
                  error = function(e) limma::voom(y, design))
    fit <- limma::lmFit(v, design)
    coef_test <- pick_condition_coef(colnames(design))

    if (is.finite(lfc_threshold) && lfc_threshold > 0) {
      # treat() testea H0: |log2FC| <= lfc en lugar de H0: log2FC = 0.
      fit <- limma::treat(fit, lfc = lfc_threshold, robust = TRUE)
      tt <- limma::topTreat(fit, coef = coef_test, number = Inf, sort.by = "none")
    } else {
      fit <- limma::eBayes(fit, robust = TRUE)
      tt <- limma::topTable(fit, coef = coef_test, number = Inf, sort.by = "none")
    }

    base_mean <- if (!is.null(tt[["AveExpr"]])) 2 ^ tt[["AveExpr"]]
                 else rep(NA_real_, nrow(tt))
    # lfcSE recuperable en limma: stdev.unscaled * sqrt(s2.post) para el
    # coeficiente testeado. El indexado es POSICIONAL a proposito: `s2.post` es
    # un vector sin nombres, asi que indexarlo por gen devolveria NA.
    lfc_se <- tryCatch({
      su <- fit$stdev.unscaled
      s2 <- fit$s2.post
      if (!is.null(su) && !is.null(s2) && coef_test %in% colnames(su)) {
        idx <- match(rownames(tt), rownames(su))
        as.numeric(su[idx, coef_test] * sqrt(s2[idx]))
      } else rep(NA_real_, nrow(tt))
    }, error = function(e) rep(NA_real_, nrow(tt)))

    tab <- data.frame(
      gene = rownames(tt),
      baseMean = base_mean,
      log2FC = tt[["logFC"]],
      log2FC_shrunk = NA_real_,
      lfcSE = lfc_se,
      stat = tt[["t"]] %||% rep(NA_real_, nrow(tt)),
      pvalue = tt[["P.Value"]],
      padj = tt[["adj.P.Val"]],
      stringsAsFactors = FALSE
    )
    list(table = tab, coef = coef_test)
  }, error = function(e) e)

  if (inherits(out, "error")) {
    return(c(list(table = NULL, error = conditionMessage(out)), info))
  }
  info$coef     <- out$coef
  info$contrast <- contrast_label(out$coef, d$ref)
  c(list(table = out$table, error = NULL), info)
}

#' Dispatcher: corre el motor solicitado.
#'
#' Devuelve list(table, error, method, contrast, coef, n_levels, shrink, fdr,
#' lfc_threshold). `contrast` es la etiqueta legible de la comparacion que se ha
#' testeado y `n_levels` el numero de niveles de `condition`: con n_levels > 2
#' solo se reporta una de las comparaciones posibles y hay que avisarlo.
run_deg <- function(counts, meta, method = c("DESeq2", "edgeR", "limma-voom"),
                    ref_level = NULL, batch = NULL,
                    fdr = 0.05, lfc_threshold = 0, shrink = TRUE) {
  method <- match.arg(method)
  res <- switch(method,
    "DESeq2"     = run_deg_deseq2(counts, meta, ref_level, batch, fdr, lfc_threshold, shrink),
    "edgeR"      = run_deg_edger(counts, meta, ref_level, batch, fdr, lfc_threshold, shrink),
    "limma-voom" = run_deg_limma(counts, meta, ref_level, batch, fdr, lfc_threshold, shrink)
  )
  res$method        <- method
  res$fdr           <- fdr
  res$lfc_threshold <- lfc_threshold
  res
}

#' Selecciona los genes significativos y aplica los filtros de visualizacion.
#'
#' IMPORTANTE — dos cosas distintas conviven aqui:
#'   - `fdr` es el nivel del test: seleccionar `padj <= fdr` es leer el
#'     resultado, no filtrarlo. Debe ser el MISMO valor que se paso a
#'     `run_deg()`, o el FDR declarado deja de corresponder al calculado.
#'   - `abs_log2fc` y `base_mean` son filtros de VISUALIZACION. Recortan lo que
#'     se muestra y NO llevan garantia estadistica: un corte de fold-change a
#'     posteriori reduce la FDR real en una cantidad desconocida. Para umbralizar
#'     por fold-change con garantia hay que pasar `lfc_threshold` a `run_deg()`,
#'     que lo mete dentro del test. Por eso `abs_log2fc` vale 0 por defecto.
apply_deg_filters <- function(deg_df, fdr = 0.05, abs_log2fc = 0, base_mean = 0) {
  if (is.null(deg_df) || !nrow(deg_df)) return(deg_df)
  keep <- rep(TRUE, nrow(deg_df))
  if (!is.null(fdr) && is.finite(fdr)) {
    keep <- keep & !is.na(deg_df$padj) & deg_df$padj <= fdr
  }
  if (!is.null(abs_log2fc) && is.finite(abs_log2fc) && abs_log2fc > 0) {
    keep <- keep & !is.na(deg_df$log2FC) & abs(deg_df$log2FC) >= abs_log2fc
  }
  if (!is.null(base_mean) && is.finite(base_mean) && base_mean > 0) {
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
