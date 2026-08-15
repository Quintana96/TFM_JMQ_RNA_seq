#' utils_replicability.R
#' Evaluacion por bootstrap de la replicabilidad de un analisis diferencial.
#'
#' Por que existe (docs/REVISION_ESTADISTICA.md, B8): el estudio de
#' replicabilidad en cohortes pequenas (PLOS Comput Biol) submuestreo 18.000
#' experimentos sobre 18 datasets y encontro que con 3 replicas la
#' replicabilidad es pobre (mediana < 0,5) en casi todos, que el minimo
#' aceptable esta en 5-7, y que con 10 la mayoria supera precision 0,9. Su
#' recomendacion es operativa: estimar por BOOTSTRAP la replicabilidad de tus
#' propios resultados antes de publicarlos, con Spearman > 0,9 como senal de
#' precision alta y < 0,8 como aviso de probables falsos positivos.
#'
#' Ninguna de las apps Shiny comparables ofrece esto.
#'
#' El coste computacional es el obstaculo real: son N reajustes del modelo, por
#' lo que el numero de remuestreos es configurable.

#' Umbrales de interpretacion del estudio de cohortes pequenas.
REPLICABILITY_HIGH <- 0.9
REPLICABILITY_LOW  <- 0.8

#' Remuestrea muestras dentro de cada grupo, conservando el diseno.
#'
#' El remuestreo es ESTRATIFICADO por condicion: un bootstrap ingenuo sobre
#' todas las muestras puede dejar un grupo vacio y hacer imposible el contraste.
#' Se exige que cada nivel conserve al menos 2 muestras distintas, porque con
#' una sola replica no hay dispersion dentro del grupo que estimar.
#' @param batch nombre de la columna de batch si esta en el modelo. Cuando la
#'   hay, la estratificacion es por condicion x batch: si solo se estratifica
#'   por condicion, un remuestreo puede dejar un batch con una sola muestra o
#'   colineal con la condicion, y ese remuestreo falla. Como los fallos solo se
#'   cuentan, la estimacion acaba apoyada en los remuestreos "faciles" y sale
#'   optimista.
bootstrap_sample_indices <- function(meta, min_per_group = 2L, batch = NULL) {
  if (is.null(meta) || !nrow(meta) || !"condition" %in% names(meta)) return(NULL)
  usar_batch <- !is.null(batch) && nzchar(batch %||% "") && batch %in% names(meta)
  estrato <- if (usar_batch) {
    paste(as.character(meta$condition), as.character(meta[[batch]]), sep = "\r")
  } else {
    as.character(meta$condition)
  }
  # Con estratos muy pequenos (frecuente en condicion x batch) exigir 2 muestras
  # distintas por estrato es imposible; el minimo se aplica entonces por
  # CONDICION, que es lo que el contraste necesita de verdad.
  min_estrato <- if (usar_batch) 1L else min_per_group
  idx <- integer(0)
  for (lv in unique(estrato)) {
    pos <- which(estrato == lv)
    if (length(pos) < min_estrato) return(NULL)
    s <- pos
    for (attempt in 1:20) {
      # sample.int y no sample(pos, ...): con un estrato de un solo elemento,
      # sample(pos, 1) interpreta `pos` como 1:pos y devuelve un indice
      # cualquiera en ese rango en lugar de esa muestra. Los estratos de tamano
      # 1 son frecuentes al cruzar condicion x batch.
      s <- pos[sample.int(length(pos), length(pos), replace = TRUE)]
      if (length(unique(s)) >= min(min_estrato, length(pos))) break
    }
    idx <- c(idx, s)
  }
  if (usar_batch) {
    # Comprobacion final a nivel de condicion: cada grupo del contraste tiene
    # que conservar al menos dos muestras distintas.
    g <- as.character(meta$condition)[idx]
    for (lv in unique(as.character(meta$condition))) {
      if (length(unique(idx[g == lv])) < min_per_group) return(NULL)
    }
  }
  sort(idx)
}

#' Concordancia entre dos tablas DEG del mismo experimento.
#'
#' Tres metricas complementarias:
#'   - `spearman`: correlacion de rangos del estadistico sobre los genes comunes.
#'     Es la que el estudio de cohortes pequenas usa como criterio.
#'   - `jaccard_topn`: estabilidad de la lista top-N por p-valor.
#'   - `overlap_sig`: solapamiento de las listas de significativos.
deg_concordance <- function(ref, boot, fdr = 0.05, top_n = 100,
                            rank_by = c("stat", "pvalue")) {
  rank_by <- match.arg(rank_by)
  out <- list(spearman = NA_real_, jaccard_topn = NA_real_, overlap_sig = NA_real_,
              n_common = 0L, n_sig_ref = 0L, n_sig_boot = 0L)
  if (is.null(ref) || is.null(boot) || !nrow(ref) || !nrow(boot)) return(out)
  common <- intersect(ref$gene, boot$gene)
  out$n_common <- length(common)
  if (length(common) < 10) return(out)
  ri <- match(common, ref$gene); bi <- match(common, boot$gene)

  v_ref <- if (rank_by == "stat" && !all(is.na(ref$stat))) ref$stat[ri] else
    sign(ref$log2FC[ri]) * -log10(pmax(ref$pvalue[ri], .Machine$double.xmin))
  v_boot <- if (rank_by == "stat" && !all(is.na(boot$stat))) boot$stat[bi] else
    sign(boot$log2FC[bi]) * -log10(pmax(boot$pvalue[bi], .Machine$double.xmin))
  ok <- is.finite(v_ref) & is.finite(v_boot)
  if (sum(ok) >= 10) {
    out$spearman <- suppressWarnings(
      stats::cor(v_ref[ok], v_boot[ok], method = "spearman"))
  }

  top_of <- function(df, n) {
    d <- df[!is.na(df$pvalue), , drop = FALSE]
    if (!nrow(d)) return(character(0))
    d <- d[order(d$pvalue), , drop = FALSE]
    d$gene[seq_len(min(n, nrow(d)))]
  }
  t_ref <- top_of(ref, top_n); t_boot <- top_of(boot, top_n)
  if (length(t_ref) && length(t_boot)) {
    out$jaccard_topn <- length(intersect(t_ref, t_boot)) /
      length(union(t_ref, t_boot))
  }
  s_ref  <- ref$gene[!is.na(ref$padj) & ref$padj <= fdr]
  s_boot <- boot$gene[!is.na(boot$padj) & boot$padj <= fdr]
  out$n_sig_ref <- length(s_ref); out$n_sig_boot <- length(s_boot)
  if (length(s_ref)) out$overlap_sig <- length(intersect(s_ref, s_boot)) / length(s_ref)
  out
}

#' Interpretacion de la replicabilidad segun los umbrales publicados.
interpret_replicability <- function(median_spearman) {
  if (is.na(median_spearman)) {
    return(list(level = "desconocida", label = "No se ha podido estimar",
                detail = "No hay remuestreos validos suficientes."))
  }
  if (median_spearman >= REPLICABILITY_HIGH) {
    return(list(level = "alta", label = "Replicabilidad alta",
      detail = paste0("La correlacion de Spearman mediana entre remuestreos es ",
                      round(median_spearman, 3), " (>= ", REPLICABILITY_HIGH,
                      "). Los rangos de los genes son estables: la lista aguanta ",
                      "el remuestreo de las muestras.")))
  }
  if (median_spearman < REPLICABILITY_LOW) {
    return(list(level = "baja", label = "Replicabilidad baja",
      detail = paste0("La correlacion de Spearman mediana es ",
                      round(median_spearman, 3), " (< ", REPLICABILITY_LOW,
                      "). Es la senal de aviso que el estudio de cohortes ",
                      "pequenas asocia a probables falsos positivos: quitar una ",
                      "muestra cambia sustancialmente el resultado. Con este ",
                      "tamano muestral, interpreta la lista con cautela.")))
  }
  list(level = "intermedia", label = "Replicabilidad intermedia",
    detail = paste0("La correlacion de Spearman mediana es ",
                    round(median_spearman, 3), ", entre ", REPLICABILITY_LOW,
                    " y ", REPLICABILITY_HIGH, ". Ni estable ni alarmante: los ",
                    "genes del top de la lista son mas fiables que la cola."))
}

#' Panel de replicabilidad por bootstrap.
#'
#' @param counts,meta matriz y metadatos ya alineados y prefiltrados
#' @param n_boot numero de remuestreos (20-50 es abordable)
#' @param progress funcion opcional f(i, n) para informar del avance
#' @return list(reference, per_boot, summary, interpretation, n_ok, n_failed)
bootstrap_replicability <- function(counts, meta, method = "DESeq2",
                                    ref_level = NULL, contrast_num = NULL,
                                    batch = NULL, fdr = 0.05, lfc_threshold = 0,
                                    design_formula = NULL, test_coef = NULL,
                                    n_boot = 20L, top_n = 100L, seed = 1L,
                                    progress = NULL) {
  fail <- function(msg) list(reference = NULL, per_boot = NULL, summary = NULL,
                             interpretation = NULL, n_ok = 0L, n_failed = 0L,
                             error = msg)
  if (is.null(counts) || is.null(meta) || !nrow(meta)) return(fail("Sin datos."))
  n_boot <- max(2L, as.integer(n_boot))

  # Con menos de 3 muestras por grupo el bootstrap es DEGENERADO: el unico
  # remuestreo que conserva 2 muestras distintas es el original, asi que la
  # correlacion sale 1 y se reportaria "replicabilidad alta" cuando en realidad
  # no se ha medido nada. Es justo el tipo de conclusion enganosa que este panel
  # existe para evitar, asi que se rechaza en lugar de devolver un numero bonito.
  g <- as.character(meta$condition)
  g <- g[!is.na(g) & nzchar(g)]
  if (!length(g)) return(fail("El samplesheet no tiene condiciones utilizables."))
  min_group <- min(table(g))
  if (min_group < 3) {
    return(fail(paste0(
      "El grupo mas pequeno tiene ", min_group, " muestras. Con menos de 3 por ",
      "grupo el remuestreo solo puede reproducir las mismas muestras, de modo ",
      "que la correlacion saldria 1 sin haber medido nada. Este diagnostico ",
      "necesita al menos 3 replicas por grupo (la literatura situa el minimo ",
      "aceptable para replicabilidad real en 5-7).")))
  }

  # El analisis de referencia: el mismo que se remuestrea, sin encogido (no
  # aporta nada a la correlacion de rangos y cuesta tiempo x n_boot).
  base_args <- list(counts = counts, meta = meta, method = method,
                    ref_level = ref_level, batch = batch, fdr = fdr,
                    lfc_threshold = lfc_threshold, shrink = FALSE,
                    contrast_num = contrast_num, design_formula = design_formula,
                    test_coef = test_coef)
  ref <- tryCatch(do.call(run_deg, base_args), error = function(e) NULL)
  if (is.null(ref) || is.null(ref$table)) {
    return(fail(paste0("El analisis de referencia fallo: ",
                       ref$error %||% "error desconocido")))
  }

  # with_seed y no set.seed a secas: fijar la semilla sin restaurarla altera el
  # RNG del resto de la sesion, que es justo lo que global.R declara que NO pasa.
  # Los indices se generan de una vez para que el resultado sea reproducible
  # incluso si la ejecucion se paraleliza.
  idx_list <- withr::with_seed(seed, lapply(seq_len(n_boot), function(i)
    bootstrap_sample_indices(meta, batch = batch)))
  valid <- !vapply(idx_list, is.null, logical(1))
  if (!any(valid)) {
    return(fail(paste0("No se pueden generar remuestreos con al menos 2 muestras ",
                       "por grupo. Hacen falta mas replicas.")))
  }
  idx_list <- idx_list[valid]

  one_boot <- function(idx) {
    m <- meta[idx, , drop = FALSE]
    # Los sample_id se renombran porque el remuestreo con reemplazo los repite y
    # DESeq2 exige nombres unicos. Hay que reponer tambien los rownames: al
    # subindexar con repeticiones, R los deja como "8", "8.1", y DESeq2 exige que
    # coincidan exactamente con los colnames de la matriz.
    m$sample_id <- paste0("r", seq_along(idx))
    rownames(m) <- m$sample_id
    cm <- counts[, idx, drop = FALSE]
    colnames(cm) <- m$sample_id
    a <- base_args; a$counts <- cm; a$meta <- m
    r <- tryCatch(do.call(run_deg, a), error = function(e) NULL)
    if (is.null(r) || is.null(r$table)) return(NULL)
    deg_concordance(ref$table, r$table, fdr = fdr, top_n = top_n)
  }

  res <- vector("list", length(idx_list))
  for (i in seq_along(idx_list)) {
    if (is.function(progress)) progress(i, length(idx_list))
    res[[i]] <- one_boot(idx_list[[i]])
  }
  keep <- !vapply(res, is.null, logical(1))
  n_failed <- sum(!keep) + sum(!valid)
  if (!any(keep)) return(fail("Todos los remuestreos han fallado."))

  per_boot <- do.call(rbind, lapply(res[keep], function(x) as.data.frame(x)))
  per_boot$boot <- seq_len(nrow(per_boot))
  qs <- function(v) stats::quantile(v, c(0.25, 0.5, 0.75), na.rm = TRUE)
  summary <- data.frame(
    metrica = c("Spearman del estadistico", "Jaccard del top-N",
                "Solapamiento de significativos"),
    q1 = c(qs(per_boot$spearman)[1], qs(per_boot$jaccard_topn)[1], qs(per_boot$overlap_sig)[1]),
    mediana = c(qs(per_boot$spearman)[2], qs(per_boot$jaccard_topn)[2], qs(per_boot$overlap_sig)[2]),
    q3 = c(qs(per_boot$spearman)[3], qs(per_boot$jaccard_topn)[3], qs(per_boot$overlap_sig)[3]),
    stringsAsFactors = FALSE
  )
  list(reference = ref$table, per_boot = per_boot, summary = summary,
       interpretation = interpret_replicability(summary$mediana[1]),
       n_ok = sum(keep), n_failed = n_failed, top_n = top_n, error = NULL)
}
