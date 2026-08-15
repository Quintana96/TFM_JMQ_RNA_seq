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

#' Construye matriz de diseno.
#'
#' Tres modos, de menos a mas expresivo:
#'   - `~ condition` (por defecto);
#'   - `~ batch + condition` si se indica `batch`;
#'   - `design_formula` arbitraria, que permite disenos pareados
#'     (`~ subject + condition`), covariables continuas e interacciones. En ese
#'     caso las variables se tipan con `prepare_design_meta()`, que respeta las
#'     numericas en lugar de convertirlas a factor.
#'
#' `condition` siempre acaba como factor releveleado al denominador del
#' contraste, tambien en modo formula libre, porque de ahi sale la etiqueta del
#' contraste y el coeficiente a testear.
#'
#' Devuelve tambien los niveles del factor y el de referencia.
build_design <- function(meta, ref_level = NULL, batch = NULL,
                         design_formula = NULL) {
  if (!is.null(design_formula)) {
    if (is.character(design_formula)) {
      design_formula <- stats::as.formula(design_formula)
    }
    meta <- prepare_design_meta(meta, all.vars(design_formula))
    formula_obj <- design_formula
  } else if (!is.null(batch) && nzchar(batch) && batch %in% names(meta)) {
    meta[[batch]] <- as.factor(meta[[batch]])
    formula_obj <- stats::as.formula(paste0("~ ", batch, " + condition"))
  } else {
    formula_obj <- stats::as.formula("~ condition")
  }
  meta$condition <- as.factor(as.character(meta$condition))
  if (!is.null(ref_level) && ref_level %in% levels(meta$condition)) {
    meta$condition <- stats::relevel(meta$condition, ref = ref_level)
  }
  lvls <- levels(meta$condition)
  list(meta = meta, formula = formula_obj, levels = lvls, ref = lvls[1])
}

#' Resuelve el coeficiente a testear: el explicito si se da y existe, y si no el
#' que corresponde al numerador del contraste.
#'
#' El emparejamiento es tolerante porque los dos estilos de nombre no coinciden
#' en las interacciones: `model.matrix` produce "genotipowt:conditiontrt" y
#' DESeq2 lo renombra a "genotipowt.conditiontrt". `make.names()` traduce entre
#' ambos, asi que un coeficiente pedido en un estilo se encuentra en el otro.
resolve_test_coef <- function(coef_names, test_coef = NULL, num = NULL, ref = NULL) {
  if (!is.null(test_coef) && length(test_coef) && nzchar(test_coef)) {
    if (test_coef %in% coef_names) return(test_coef)
    hit <- coef_names[make.names(coef_names) == make.names(test_coef)]
    if (length(hit)) return(hit[1])
  }
  condition_coef_for(coef_names, num, ref)
}

#' Coeficiente de `condition` que se va a testear.
#'
#' Si se pide un numerador concreto (`num`) se devuelve su coeficiente, que es
#' como se implementa el selector de contraste: el denominador se ha puesto como
#' nivel de referencia del factor, asi que "num vs den" ES un coeficiente del
#' diseno. Eso importa porque `lfcShrink(type = "apeglm")` solo acepta `coef` y
#' no `contrast`, de modo que la via del relevel es la unica que permite
#' contrastes arbitrarios CON encogido.
#'
#' Sin `num`, se devuelve el ultimo coeficiente de `condition`: con dos niveles
#' es el unico posible, con tres o mas es una comparacion arbitraria que quien
#' llama debe etiquetar.
condition_coef_for <- function(coef_names, num = NULL, ref = NULL) {
  cond <- grep("^condition", coef_names, value = TRUE)
  if (!length(cond)) return(coef_names[length(coef_names)])
  if (!is.null(num) && length(num) && nzchar(num)) {
    # DESeq2 nombra "condition_<num>_vs_<ref>"; model.matrix, "condition<num>".
    if (!is.null(ref) && nzchar(ref)) {
      exact <- paste0("condition_", num, "_vs_", ref)
      if (exact %in% cond) return(exact)
    }
    if (paste0("condition", num) %in% cond) return(paste0("condition", num))
    pref <- cond[startsWith(cond, paste0("condition_", num, "_vs_"))]
    if (length(pref)) return(pref[1])
  }
  cond[length(cond)]
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

#' Estructura de metadatos que devuelven todos los motores.
#'
#' Centraliza el contrato: si se anade un campo, se anade una vez y no cinco. El
#' bloque estaba repetido literalmente en cada motor.
deg_engine_info <- function(d) {
  list(contrast = NA_character_, coef = NA_character_,
       n_levels = length(d$levels), shrink = "ninguno",
       padj_method = "BH", disp_data = NULL, cooks = NULL,
       coef_available = character(0))
}

#' Elige el mejor tipo de encogido de log2FC disponible.
#' `apeglm` es la recomendacion de la vinieta de DESeq2 (prior de colas anchas:
#' encoge el ruido sin aplastar los efectos grandes reales).
pick_shrink_type <- function() {
  if (requireNamespace("apeglm", quietly = TRUE)) return("apeglm")
  if (requireNamespace("ashr", quietly = TRUE)) return("ashr")
  "normal"
}

#' Prepara el uso de IHW como funcion de filtrado de DESeq2.
#'
#' IHW pondera cada hipotesis segun una covariable independiente del p-valor bajo
#' la nula (aqui `baseMean`) y gana potencia sobre BH sin perder el control de la
#' FDR. Es una generalizacion del filtrado independiente que DESeq2 ya hace.
#'
#' Necesita S4Vectors ATACHADO, no solo cargado: la ruta interna de IHW llama a
#' `mcols()` sin cualificar, y la app trabaja con prefijos `DESeq2::` sin atachar
#' nada, asi que sin esto falla con "no se pudo encontrar la funcion mcols".
#' Devuelve la funcion o NULL si no se puede usar.
ihw_filter_fun <- function() {
  if (!requireNamespace("IHW", quietly = TRUE)) return(NULL)
  ok <- requireNamespace("S4Vectors", quietly = TRUE) &&
    (("package:S4Vectors" %in% search()) ||
       isTRUE(suppressPackageStartupMessages(
         require("S4Vectors", quietly = TRUE, character.only = TRUE))))
  if (!ok) return(NULL)
  IHW::ihw
}

#' Tabla de dispersiones de un DESeqDataSet ajustado, para dibujar el
#' equivalente de plotDispEsts() con plotly en lugar de graficos base.
deseq_dispersion_data <- function(dds) {
  tryCatch({
    mc <- S4Vectors::mcols(dds)
    needed <- c("baseMean", "dispGeneEst", "dispFit", "dispersion")
    if (!all(needed %in% names(mc))) return(NULL)
    data.frame(
      baseMean    = as.numeric(mc$baseMean),
      dispGeneEst = as.numeric(mc$dispGeneEst),
      dispFit     = as.numeric(mc$dispFit),
      dispersion  = as.numeric(mc$dispersion),
      outlier     = if ("dispOutlier" %in% names(mc)) as.logical(mc$dispOutlier)
                    else rep(FALSE, nrow(mc)),
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)
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
#' @param contrast_num Nivel que actua de numerador. `ref_level` es el
#'   denominador, de modo que el contraste pedido es un coeficiente del diseno.
#' @param use_ihw Si TRUE, sustituye el filtrado independiente por IHW.
run_deg_deseq2 <- function(counts, meta, ref_level = NULL, batch = NULL,
                           fdr = 0.05, lfc_threshold = 0, shrink = TRUE,
                           contrast_num = NULL, use_ihw = FALSE,
                           design_formula = NULL, test_coef = NULL) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    return(list(table = NULL, error = "DESeq2 no esta instalado."))
  }
  d <- build_design(meta, ref_level, batch, design_formula)
  info <- deg_engine_info(d)
  out <- tryCatch({
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = round(as.matrix(counts)),
      colData = d$meta,
      design = d$formula
    )
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    coef_name <- resolve_test_coef(DESeq2::resultsNames(dds), test_coef,
                                   contrast_num, d$ref)

    filter_fun <- if (isTRUE(use_ihw)) ihw_filter_fun() else NULL
    padj_method <- if (!is.null(filter_fun)) "IHW" else "BH"
    res_args <- list(
      dds, name = coef_name, alpha = fdr,
      lfcThreshold = if (is.finite(lfc_threshold)) lfc_threshold else 0,
      altHypothesis = "greaterAbs"
    )
    if (!is.null(filter_fun)) res_args$filterFun <- filter_fun
    res <- do.call(DESeq2::results, res_args)
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
    cooks <- tryCatch({
      ck <- SummarizedExperiment::assays(dds)[["cooks"]]
      if (is.null(ck)) NULL else cooks_sample_summary(ck)
    }, error = function(e) NULL)
    list(table = tab, coef = coef_name, shrink = shrink_used,
         padj_method = padj_method, disp_data = deseq_dispersion_data(dds),
         cooks = cooks, coef_available = DESeq2::resultsNames(dds))
  }, error = function(e) e)

  if (inherits(out, "error")) {
    return(c(list(table = NULL, error = conditionMessage(out)), info))
  }
  info$coef           <- out$coef
  info$contrast       <- contrast_label(out$coef, d$ref)
  info$shrink         <- out$shrink
  info$padj_method    <- out$padj_method
  info$disp_data      <- out$disp_data
  info$cooks          <- out$cooks
  info$coef_available <- out$coef_available %||% character(0)
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
                          fdr = 0.05, lfc_threshold = 0, shrink = TRUE,
                          contrast_num = NULL, use_ihw = FALSE,
                          design_formula = NULL, test_coef = NULL) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    return(list(table = NULL, error = "edgeR no esta instalado."))
  }
  d <- build_design(meta, ref_level, batch, design_formula)
  info <- deg_engine_info(d)
  out <- tryCatch({
    design <- stats::model.matrix(d$formula, data = d$meta)
    y <- edgeR::DGEList(counts = round(as.matrix(counts)), group = d$meta$condition)
    y <- norm_lib_sizes(y)
    fit <- if (edger_is_v4()) {
      edgeR::glmQLFit(y, design)
    } else {
      edgeR::glmQLFit(edgeR::estimateDisp(y, design), design)
    }
    coef_test <- resolve_test_coef(colnames(design), test_coef, contrast_num, d$ref)

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
    list(table = tab, coef = coef_test, coef_available = colnames(design))
  }, error = function(e) e)

  if (inherits(out, "error")) {
    return(c(list(table = NULL, error = conditionMessage(out)), info))
  }
  info$coef           <- out$coef
  info$contrast       <- contrast_label(out$coef, d$ref)
  info$coef_available <- out$coef_available %||% character(0)
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
                          fdr = 0.05, lfc_threshold = 0, shrink = TRUE,
                          contrast_num = NULL, use_ihw = FALSE,
                          design_formula = NULL, test_coef = NULL) {
  if (!requireNamespace("limma", quietly = TRUE) ||
      !requireNamespace("edgeR", quietly = TRUE)) {
    return(list(table = NULL, error = "Se requieren limma y edgeR para limma-voom."))
  }
  d <- build_design(meta, ref_level, batch, design_formula)
  info <- deg_engine_info(d)
  out <- tryCatch({
    design <- stats::model.matrix(d$formula, data = d$meta)
    y <- edgeR::DGEList(counts = round(as.matrix(counts)), group = d$meta$condition)
    y <- norm_lib_sizes(y)
    v <- tryCatch(limma::voomWithQualityWeights(y, design),
                  error = function(e) limma::voom(y, design))
    fit <- limma::lmFit(v, design)
    coef_test <- resolve_test_coef(colnames(design), test_coef, contrast_num, d$ref)

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
    list(table = tab, coef = coef_test, coef_available = colnames(design))
  }, error = function(e) e)

  if (inherits(out, "error")) {
    return(c(list(table = NULL, error = conditionMessage(out)), info))
  }
  info$coef           <- out$coef
  info$contrast       <- contrast_label(out$coef, d$ref)
  info$coef_available <- out$coef_available %||% character(0)
  c(list(table = out$table, error = NULL), info)
}

#' Motor Wilcoxon rank-sum sobre CPM normalizados.
#'
#' Contexto (docs/REVISION_ESTADISTICA.md, B7), que da para discusion y conviene
#' contar completo porque la conclusion ha cambiado tres veces:
#'   1. Li et al. (2022) reportaron que en muestras poblacionales humanas la FDR
#'      real de DESeq2 y edgeR llega a superar el 20 % cuando el objetivo es 5 %,
#'      y recomendaron Wilcoxon para n grande.
#'   2. Hejblum et al. (2024) identificaron un fallo en esa simulacion: los datos
#'      se normalizaban DESPUES de permutar, de modo que ya no cumplian la nula.
#'      Normalizando antes, los tres metodos controlaban bien la FDR.
#'   3. Li et al. (2024) aceptaron el sesgo senalado pero mantienen su conclusion
#'      apoyandose en datos totalmente permutados.
#'
#' Lectura practica: con n pequeno los parametricos son la eleccion correcta y
#' Wilcoxon seria un error por falta de potencia. Con n grande (>= 8-10 por
#' grupo) merece la pena comparar. Por eso este motor existe pero no se sugiere
#' hasta que el tamano muestral lo justifica.
#'
#' CPM normalizados por composicion (TMM cuando edgeR esta disponible).
#'
#' Los motores robustos usaban CPM por tamano de libreria crudo. El benchmark
#' que motiva el motor de Wilcoxon (Li et al., Genome Biology 2022) normaliza
#' con TMM, y sin corregir por composicion una diferencia de composicion entre
#' grupos se convierte en falsos positivos — justo en el regimen de n grande en
#' el que la app recomienda estos motores.
#'
#' @param counts matriz de conteos
#' @return matriz de CPM
normalized_cpm <- function(counts) {
  cm <- as.matrix(counts)
  if (requireNamespace("edgeR", quietly = TRUE)) {
    out <- tryCatch({
      y <- edgeR::DGEList(counts = cm)
      y <- norm_lib_sizes(y)
      edgeR::cpm(y, normalized.lib.sizes = TRUE)
    }, error = function(e) NULL)
    if (!is.null(out)) return(out)
  }
  # Respaldo sin edgeR: CPM por tamano de libreria, peor pero utilizable.
  libs <- colSums(cm); libs[libs == 0] <- 1
  t(t(cm) / libs) * 1e6
}

#' Solo admite dos grupos: es un test de dos muestras, no un modelo, asi que no
#' puede ajustar por batch ni por covariables.
run_deg_wilcoxon <- function(counts, meta, ref_level = NULL, batch = NULL,
                             fdr = 0.05, lfc_threshold = 0, shrink = TRUE,
                             contrast_num = NULL, use_ihw = FALSE,
                             design_formula = NULL, test_coef = NULL) {
  d <- build_design(meta, ref_level, batch)
  info <- deg_engine_info(d)
  num <- if (!is.null(contrast_num) && nzchar(contrast_num %||% "")) contrast_num
         else utils::tail(d$levels, 1)
  den <- d$ref
  if (identical(num, den)) {
    return(c(list(table = NULL, error = "Numerador y denominador coinciden."), info))
  }
  out <- tryCatch({
    g <- as.character(d$meta$condition)
    i_num <- which(g == num); i_den <- which(g == den)
    if (length(i_num) < 2 || length(i_den) < 2) {
      stop("Wilcoxon necesita al menos 2 muestras en cada grupo del contraste.")
    }
    cpm <- normalized_cpm(counts)
    a <- cpm[, i_num, drop = FALSE]; b <- cpm[, i_den, drop = FALSE]
    pv <- vapply(seq_len(nrow(cpm)), function(i) {
      tryCatch(stats::wilcox.test(a[i, ], b[i, ], exact = FALSE)$p.value,
               error = function(e) NA_real_)
    }, numeric(1))
    lfc <- log2((rowMeans(a) + 1) / (rowMeans(b) + 1))
    tab <- data.frame(
      gene = rownames(cpm),
      baseMean = rowMeans(cpm),
      log2FC = lfc,
      log2FC_shrunk = NA_real_,
      lfcSE = NA_real_,
      stat = NA_real_,
      pvalue = pv,
      padj = stats::p.adjust(pv, method = "BH"),
      stringsAsFactors = FALSE
    )
    list(table = tab, coef = paste0("condition", num))
  }, error = function(e) e)
  if (inherits(out, "error")) {
    return(c(list(table = NULL, error = conditionMessage(out)), info))
  }
  info$coef     <- out$coef
  info$contrast <- paste(num, "vs", den)
  c(list(table = out$table, error = NULL), info)
}

#' Motor dearseq: test de componentes de varianza con regresion no parametrica.
#'
#' Controla la FDR sin asumir la distribucion de los conteos (Gauthier et al.,
#' NAR Genomics and Bioinformatics 2020) y admite disenos longitudinales, lo que
#' lo hace el complemento natural de los disenos arbitrarios del item 16.
run_deg_dearseq <- function(counts, meta, ref_level = NULL, batch = NULL,
                            fdr = 0.05, lfc_threshold = 0, shrink = TRUE,
                            contrast_num = NULL, use_ihw = FALSE,
                            design_formula = NULL, test_coef = NULL) {
  d <- build_design(meta, ref_level, batch)
  info <- deg_engine_info(d)
  if (!requireNamespace("dearseq", quietly = TRUE)) {
    return(c(list(table = NULL, error = "dearseq no esta instalado."), info))
  }
  num <- if (!is.null(contrast_num) && nzchar(contrast_num %||% "")) contrast_num
         else utils::tail(d$levels, 1)
  den <- d$ref
  out <- tryCatch({
    keep_s <- as.character(d$meta$condition) %in% c(num, den)
    if (sum(keep_s) < 4) stop("dearseq necesita al menos 4 muestras en el contraste.")
    cm <- round(as.matrix(counts))[, keep_s, drop = FALSE]
    m <- d$meta[keep_s, , drop = FALSE]
    variables2test <- stats::model.matrix(~ condition, data = droplevels(m))[, -1, drop = FALSE]
    covariates <- if (!is.null(batch) && nzchar(batch %||% "") && batch %in% names(m)) {
      stats::model.matrix(stats::as.formula(paste0("~ ", batch)), data = m)[, -1, drop = FALSE]
    } else NULL
    res <- dearseq::dear_seq(
      exprmat = cm, variables2test = variables2test,
      covariates = covariates, which_test = "asymptotic",
      preprocessed = FALSE, parallel_comp = FALSE, progressbar = FALSE
    )
    rt <- res$pvals
    pv <- as.numeric(rt[["rawPval"]])
    names(pv) <- rownames(rt)
    cpm <- normalized_cpm(cm)
    i_num <- which(as.character(m$condition) == num)
    i_den <- which(as.character(m$condition) == den)
    lfc <- log2((rowMeans(cpm[, i_num, drop = FALSE]) + 1) /
                  (rowMeans(cpm[, i_den, drop = FALSE]) + 1))
    genes <- rownames(rt)
    tab <- data.frame(
      gene = genes,
      baseMean = rowMeans(cpm)[genes],
      log2FC = lfc[genes],
      log2FC_shrunk = NA_real_,
      lfcSE = NA_real_,
      stat = NA_real_,
      pvalue = pv[genes],
      padj = stats::p.adjust(pv[genes], method = "BH"),
      stringsAsFactors = FALSE
    )
    list(table = tab, coef = paste0("condition", num))
  }, error = function(e) e)
  if (inherits(out, "error")) {
    return(c(list(table = NULL, error = conditionMessage(out)), info))
  }
  info$coef     <- out$coef
  info$contrast <- paste(num, "vs", den)
  c(list(table = out$table, error = NULL), info)
}

#' Motores disponibles, y cuales se recomiendan para el tamano muestral dado.
#'
#' Con n >= 8 por grupo la controversia sobre el control de la FDR con muestras
#' poblacionales deja de ser academica y conviene comparar con un metodo robusto.
#' Por debajo, los parametricos son la eleccion correcta.
DEG_METHODS_PARAMETRIC <- c("DESeq2", "edgeR", "limma-voom")
DEG_METHODS_ROBUST     <- c("Wilcoxon", "dearseq")

suggest_robust_comparison <- function(meta, min_per_group = 8L) {
  if (is.null(meta) || !"condition" %in% names(meta)) return(NULL)
  g <- as.character(meta$condition)
  g <- g[!is.na(g) & nzchar(g)]
  if (!length(g)) return(NULL)
  sizes <- table(g)
  if (min(sizes) < min_per_group) return(NULL)
  list(min_group = as.integer(min(sizes)),
       message = paste0(
         "Hay ", min(sizes), " muestras en el grupo mas pequeno. Con n >= ",
         min_per_group, " por grupo merece la pena comparar el resultado con un ",
         "metodo robusto (Wilcoxon o dearseq): en muestras poblacionales grandes ",
         "se ha reportado que la FDR real de los metodos parametricos puede ",
         "superar el objetivo declarado."))
}

#' Solapamiento entre dos listas de significativos, para comparar metodos.
deg_method_overlap <- function(tab_a, tab_b, fdr = 0.05,
                               name_a = "A", name_b = "B") {
  sig <- function(t) if (is.null(t)) character(0) else
    t$gene[!is.na(t$padj) & t$padj <= fdr]
  a <- sig(tab_a); b <- sig(tab_b)
  inter <- intersect(a, b)
  un <- union(a, b)
  list(name_a = name_a, name_b = name_b,
       n_a = length(a), n_b = length(b), n_common = length(inter),
       only_a = length(setdiff(a, b)), only_b = length(setdiff(b, a)),
       jaccard = if (length(un)) length(inter) / length(un) else NA_real_,
       genes_common = inter)
}

#' Dispatcher: corre el motor solicitado.
#'
#' Devuelve list(table, error, method, contrast, coef, n_levels, shrink,
#' padj_method, disp_data, cooks, fdr, lfc_threshold). `contrast` es la etiqueta
#' legible de la comparacion testeada.
#'
#' El contraste se especifica como `contrast_num` (numerador) y `ref_level`
#' (denominador). Poner el denominador como nivel de referencia del factor hace
#' que la comparacion pedida sea un coeficiente del diseno, lo que permite usar
#' `lfcShrink(coef = ...)`: `apeglm` no admite `contrast`. La diferencia
#' numerica frente a `results(contrast = ...)` es ruido del ajuste iterativo
#' (del orden de 1e-6 en log2FC, sin cambios en las llamadas de significacion).
run_deg <- function(counts, meta,
                    method = c("DESeq2", "edgeR", "limma-voom", "Wilcoxon", "dearseq"),
                    ref_level = NULL, batch = NULL,
                    fdr = 0.05, lfc_threshold = 0, shrink = TRUE,
                    contrast_num = NULL, use_ihw = FALSE,
                    design_formula = NULL, test_coef = NULL) {
  method <- match.arg(method)
  lvls <- unique(as.character(meta$condition[!is.na(meta$condition) &
                                               nzchar(as.character(meta$condition))]))
  if (!is.null(contrast_num) && length(contrast_num) && nzchar(contrast_num)) {
    if (!contrast_num %in% lvls) {
      return(list(table = NULL, method = method,
                  error = paste0("El numerador '", contrast_num,
                                 "' no es un nivel de condition.")))
    }
    if (!is.null(ref_level) && identical(contrast_num, ref_level)) {
      return(list(table = NULL, method = method,
                  error = "El numerador y el denominador del contraste son el mismo nivel."))
    }
  }
  # Con formula libre se valida ANTES de ajustar, para cambiar un error criptico
  # de DESeq2 en ingles por un diagnostico que dice cual es el problema.
  if (!is.null(design_formula)) {
    v <- validate_design_formula(design_formula, meta)
    if (!isTRUE(v$ok)) {
      return(list(table = NULL, method = method,
                  error = paste(v$errors, collapse = " ")))
    }
    design_formula <- v$formula
  }
  res <- switch(method,
    "DESeq2"     = run_deg_deseq2(counts, meta, ref_level, batch, fdr, lfc_threshold,
                                  shrink, contrast_num, use_ihw, design_formula, test_coef),
    "edgeR"      = run_deg_edger(counts, meta, ref_level, batch, fdr, lfc_threshold,
                                 shrink, contrast_num, use_ihw, design_formula, test_coef),
    "limma-voom" = run_deg_limma(counts, meta, ref_level, batch, fdr, lfc_threshold,
                                 shrink, contrast_num, use_ihw, design_formula, test_coef),
    "Wilcoxon"   = run_deg_wilcoxon(counts, meta, ref_level, batch, fdr, lfc_threshold,
                                    shrink, contrast_num, use_ihw, design_formula, test_coef),
    "dearseq"    = run_deg_dearseq(counts, meta, ref_level, batch, fdr, lfc_threshold,
                                   shrink, contrast_num, use_ihw, design_formula, test_coef)
  )
  res$method        <- method
  res$fdr           <- fdr
  res$lfc_threshold <- lfc_threshold
  # IC del log2FC donde el motor haya dado error estandar (DESeq2 y limma).
  if (!is.null(res$table)) res$table <- add_lfc_confidence_interval(res$table)
  # El diseno REPORTADO tiene que ser el que el motor ajusto de verdad. Wilcoxon
  # y dearseq no consumen `design_formula`: Wilcoxon es un test de dos muestras
  # sin modelo, y dearseq recibe la condicion y, como mucho, el batch. Declarar
  # la formula libre en esos casos hacia que el banner y el informe describieran
  # un ajuste que no habia ocurrido.
  res$design <- if (identical(method, "Wilcoxon")) {
    "sin modelo (test de dos muestras sobre CPM normalizados)"
  } else if (identical(method, "dearseq")) {
    if (!is.null(batch) && nzchar(batch %||% "")) paste0("~ ", batch, " + condition")
    else "~ condition"
  } else if (!is.null(design_formula)) {
    deparse1(design_formula)
  } else if (!is.null(batch) && nzchar(batch %||% "")) {
    paste0("~ ", batch, " + condition")
  } else "~ condition"

  # Y si se pidio una formula libre a un motor que no la usa, se dice: callarlo
  # deja al usuario creyendo que su diseno pareado se ha tenido en cuenta.
  if (!is.null(design_formula) && method %in% DEG_METHODS_ROBUST) {
    res$design_warning <- paste0(
      "El motor ", method, " no admite formulas de diseno arbitrarias: se ha ",
      "ajustado ", res$design, ". Si necesitas el diseno completo, usa DESeq2, ",
      "edgeR o limma-voom.")
  }
  res
}

#' Anade intervalo de confianza de Wald al log2FC, donde haya error estandar.
#'
#' Cierra A9 (docs/REVISION_ESTADISTICA.md): la interfaz declaraba una columna
#' `lfcSE` que solo DESeq2 rellenaba, y sin error estandar "no se pueden dibujar
#' intervalos de confianza". Ahora la rellenan DESeq2 y limma, asi que el
#' intervalo se puede calcular y exportar.
#'
#' Nota de diseno: el IC se anade a la TABLA, no como barras de error en el
#' volcano. Con miles de genes, las barras de error se solapan hasta hacer el
#' grafico ilegible y esconden justo lo que el volcano sirve para ver. El valor
#' del IC esta en la tabla y en los datos exportados, donde se puede leer gen a
#' gen.
#'
#' En edgeR-QL no hay un equivalente directo del SE, asi que sus filas quedan sin
#' intervalo; es coherente con usar `glmTreat`, que responde a "es el efecto
#' mayor que este umbral" sin necesitar el SE.
add_lfc_confidence_interval <- function(deg_df, level = 0.95) {
  if (is.null(deg_df) || !nrow(deg_df)) return(deg_df)
  if (!all(c("log2FC", "lfcSE") %in% names(deg_df))) return(deg_df)
  if (all(is.na(deg_df$lfcSE))) return(deg_df)
  z <- stats::qnorm(1 - (1 - level) / 2)
  deg_df$log2FC_lower <- deg_df$log2FC - z * deg_df$lfcSE
  deg_df$log2FC_upper <- deg_df$log2FC + z * deg_df$lfcSE
  attr(deg_df, "ci_level") <- level
  deg_df
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

#' Lista de genes significativos para el ENRIQUECIMIENTO.
#'
#' Solo aplica el FDR con el que se ajusto el modelo. Deliberadamente NO acepta
#' los filtros de |log2FC| ni de baseMean: esos son de visualizacion y no
#' recortan el universo, de modo que usarlos aqui cambiaria el resultado del ORA
#' al mover un deslizador declarado cosmetico, y ademas dejaria lista y fondo
#' definidos con criterios distintos.
#'
#' @param deg_df tabla de resultados del motor
#' @param fdr FDR objetivo del ajuste
#' @return data.frame con las filas significativas
deg_significant_genes <- function(deg_df, fdr = 0.05) {
  apply_deg_filters(deg_df, fdr = fdr, abs_log2fc = 0, base_mean = 0)
}

#' Universo (fondo) del enriquecimiento: los genes EVALUABLES.
#'
#' Son los que tienen p-valor ajustado. Un gen sin `padj` —descartado por el
#' filtrado independiente, con conteo cero o marcado como outlier de Cook— nunca
#' habria podido entrar en la lista de significativos, asi que contarlo como
#' fondo infla el enriquecimiento (Wijesooriya et al., 2022: el fondo mal
#' definido es el error mas extendido en analisis de sobre-representacion).
#'
#' @param deg_df tabla de resultados del motor
#' @return vector de identificadores de gen
deg_testable_universe <- function(deg_df) {
  if (is.null(deg_df) || !nrow(deg_df)) return(character(0))
  u <- deg_df$gene[!is.na(deg_df$padj)]
  # Un motor que no rellene padj dejaria el universo vacio, lo que romperia el
  # enriquecimiento entero; en ese caso es preferible el fondo completo.
  if (!length(u)) deg_df$gene else u
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
