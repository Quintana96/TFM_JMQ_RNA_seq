#' utils_diag.R
#' Funciones puras (sin Shiny) para los diagnosticos post-ajuste del analisis
#' diferencial: histograma de p-valores, estimacion de pi0, dispersiones, RLE y
#' distribucion de outliers de Cook.
#'
#' Por que existe este archivo (docs/REVISION_ESTADISTICA.md, B1): Pall et al.
#' (PLoS Biology 2023) revisaron 4.616 datasets de GEO y encontraron que solo el
#' 25 % producia histogramas de p-valores con la forma teoricamente esperada, y
#' que el 37 % declaraba implicitamente que mas de la mitad de los genes
#' cambian (pi0 < 0,5). Son fallos que no se ven en la tabla de resultados: solo
#' aparecen si se miran los diagnosticos, y ninguna de las apps Shiny
#' comparables los muestra.

# ── Desglose de los NA en padj ──────────────────────────────────────────────

#' Clasifica por que un gen no tiene padj.
#'
#' En DESeq2 un `padj = NA` tiene tres causas distintas y muy informativas, que
#' se distinguen por el patron de baseMean/pvalue/padj:
#'   - conteo cero en todas las muestras -> baseMean == 0
#'   - outlier detectado por distancia de Cook -> pvalue tambien es NA
#'   - descartado por el filtrado independiente -> pvalue existe, padj es NA
#' Colapsarlas a "no significativo" oculta informacion: un gen con un outlier
#' extremo puede ser el hallazgo mas interesante, o la senal de que una muestra
#' esta mal.
padj_na_breakdown <- function(deg_df) {
  if (is.null(deg_df) || !nrow(deg_df)) return(NULL)
  bm <- deg_df$baseMean
  pv <- deg_df$pvalue
  pa <- deg_df$padj
  zero     <- !is.na(bm) & bm == 0
  outlier  <- !zero & is.na(pv)
  filtered <- !zero & !is.na(pv) & is.na(pa)
  list(
    n_total    = nrow(deg_df),
    n_zero     = sum(zero, na.rm = TRUE),
    n_outlier  = sum(outlier, na.rm = TRUE),
    n_filtered = sum(filtered, na.rm = TRUE),
    n_tested   = sum(!is.na(pa), na.rm = TRUE)
  )
}

#' Texto de una linea con el desglose, listo para mostrar bajo la tabla.
padj_na_breakdown_text <- function(b) {
  if (is.null(b)) return(NULL)
  paste0(
    fmt_int(b$n_tested), " genes con p-valor ajustado; ",
    fmt_int(b$n_filtered), " descartados por el filtrado independiente; ",
    fmt_int(b$n_outlier), " marcados como outlier (Cook); ",
    fmt_int(b$n_zero), " con conteo cero en todas las muestras."
  )
}

# ── Outliers de Cook por muestra ────────────────────────────────────────────

#' Reparto entre muestras de los maximos de distancia de Cook.
#'
#' Si una sola muestra concentra los outliers, el problema no es de genes sino
#' de esa muestra. Con `n` muestras el reparto esperado por azar es 1/n, asi que
#' se marca dominancia cuando una supera el doble de su cuota y ademas el 40 %.
cooks_sample_summary <- function(cooks_mat) {
  if (is.null(cooks_mat) || !length(cooks_mat)) return(NULL)
  cm <- as.matrix(cooks_mat)
  if (!ncol(cm) || !nrow(cm)) return(NULL)
  wm <- apply(cm, 1, function(r) if (all(is.na(r))) NA_integer_ else which.max(r))
  wm <- wm[!is.na(wm)]
  if (!length(wm)) return(NULL)
  tb <- table(factor(colnames(cm)[wm], levels = colnames(cm)))
  df <- data.frame(
    sample_id = names(tb),
    n_max = as.integer(tb),
    frac = as.numeric(tb) / sum(tb),
    stringsAsFactors = FALSE
  )
  expected <- 1 / ncol(cm)
  worst <- df[which.max(df$frac), , drop = FALSE]
  df$max_cooks <- suppressWarnings(apply(cm, 2, max, na.rm = TRUE))[df$sample_id]
  list(
    table = df,
    expected_frac = expected,
    dominant = if (worst$frac > 2 * expected && worst$frac > 0.4) worst$sample_id else NA_character_
  )
}

# ── pi0 y forma del histograma de p-valores ─────────────────────────────────

#' Estima la proporcion de hipotesis nulas ciertas (pi0).
#'
#' `limma::propTrueNull` es la via por defecto porque no anade dependencias y es
#' estable con pocos p-valores; `qvalue::qvalue` se usa como comprobacion
#' cruzada y devuelve NA cuando hay muy pocos tests.
estimate_pi0 <- function(pvalues) {
  p <- pvalues[!is.na(pvalues) & is.finite(pvalues)]
  out <- list(pi0 = NA_real_, pi0_qvalue = NA_real_, n = length(p), method = NA_character_)
  if (length(p) < 20) return(out)
  out$pi0 <- tryCatch(limma::propTrueNull(p), error = function(e) NA_real_)
  if (!is.na(out$pi0)) out$method <- "limma::propTrueNull"
  out$pi0_qvalue <- tryCatch(qvalue::qvalue(p)$pi0, error = function(e) NA_real_)
  out
}

#' Datos del histograma de p-valores (conteos por bin, en lugar del grafico).
pvalue_hist_data <- function(pvalues, bins = 40) {
  p <- pvalues[!is.na(pvalues) & is.finite(pvalues)]
  if (!length(p)) return(NULL)
  brks <- seq(0, 1, length.out = bins + 1L)
  h <- graphics::hist(p, breaks = brks, plot = FALSE)
  data.frame(
    mid = h$mids, count = h$counts,
    # Suelo esperado bajo la nula: n * pi0 / bins. Se anade fuera con el pi0.
    stringsAsFactors = FALSE
  )
}

#' Veredicto automatico sobre la forma del histograma de p-valores.
#'
#' Traduce a lenguaje llano las formas que describe Pall et al.: pico a la
#' izquierda con suelo uniforme (esperada), suelo hundido con exceso de
#' p-valores altos (conservadora), acumulacion en ambos extremos (bimodal, tipica
#' de supuestos violados o efectos batch no modelados), y pi0 implausiblemente
#' bajo.
diagnose_pvalue_shape <- function(pvalues, pi0 = NA_real_) {
  p <- pvalues[!is.na(pvalues) & is.finite(pvalues)]
  if (length(p) < 50) {
    return(list(verdict = "insuficiente",
                label = "Pocos p-valores para diagnosticar",
                detail = "Hacen falta al menos 50 genes testeados."))
  }
  n <- length(p)
  first <- mean(p <= 0.05)            # densidad relativa en el primer bin
  last  <- mean(p >= 0.95)
  mid   <- mean(p > 0.4 & p < 0.6)
  # Bajo la nula, cada tramo de anchura w contiene ~w de la masa.
  r_first <- first / 0.05
  r_last  <- last / 0.05
  r_mid   <- mid / 0.20

  if (!is.na(pi0) && pi0 < 0.5) {
    return(list(verdict = "sospechoso",
      label = "pi0 implausiblemente bajo",
      detail = paste0("pi0 = ", round(pi0, 3), " implica que mas de la mitad de los genes ",
                      "cambian. Es poco creible en la mayoria de disenos y suele indicar ",
                      "estructura no modelada (batch, muestras mal asignadas).")))
  }
  if (r_last > 1.5 && r_last > r_first) {
    return(list(verdict = "bimodal",
      label = "Exceso de p-valores altos (bimodal)",
      detail = paste0("Hay ", round(r_last, 1), " veces mas p-valores cerca de 1 de lo ",
                      "esperado. Suele venir de supuestos violados, de un efecto batch no ",
                      "modelado o de conteos muy bajos sin prefiltrar.")))
  }
  if (r_mid < 0.6) {
    return(list(verdict = "conservador",
      label = "Distribucion conservadora",
      detail = paste0("El suelo del histograma esta hundido (", round(r_mid, 2),
                      " veces lo esperado): el test esta siendo mas conservador de lo ",
                      "que asume la correccion por FDR.")))
  }
  if (r_first < 1.2) {
    return(list(verdict = "plano",
      label = "Sin senal aparente",
      detail = "No hay pico a la izquierda: pocos o ningun gen diferencialmente expresado."))
  }
  list(verdict = "esperado",
    label = "Forma esperada",
    detail = paste0("Pico a la izquierda (", round(r_first, 1),
                    " veces lo esperado) sobre un suelo aproximadamente uniforme. ",
                    "Es la forma que deberia tener un analisis bien especificado."))
}

#' Distribuciones de referencia para las miniaturas de la interfaz.
#' Simuladas de forma determinista (semilla fija) para que las tres formas que
#' describe la literatura se puedan comparar de un vistazo con la propia.
reference_pvalue_shapes <- function(n = 4000) {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(42)
  on.exit({
    if (!is.null(old)) assign(".Random.seed", old, envir = globalenv())
  }, add = TRUE)
  list(
    "Esperada (pico + suelo uniforme)" = c(stats::runif(n * 0.85),
                                           stats::rbeta(n * 0.15, 0.3, 8)),
    "Conservadora (suelo hundido)"     = stats::rbeta(n, 1.6, 1),
    "Bimodal (supuestos violados)"     = c(stats::rbeta(n * 0.5, 0.4, 3),
                                           stats::rbeta(n * 0.5, 6, 0.7))
  )
}

# ── RLE (Relative Log Expression) ───────────────────────────────────────────

#' Diagnostico de sesgo de longitud en el analisis de sobre-representacion.
#'
#' Por que existe (docs/REVISION_ESTADISTICA.md, B3c): Young, Wakefield, Smyth y
#' Oshlack (2010) mostraron que el ORA estandar da resultados sesgados en RNA-seq
#' porque los transcritos largos tienen mas probabilidad de ser detectados como
#' diferenciales, lo que arrastra a las categorias GO que los contienen.
#'
#' Decision de diseno: en lugar de ofrecer `goseq` a ciegas, se MIDE el sesgo en
#' los datos del usuario. Es la informacion que decide: si en este dataset la
#' probabilidad de ser DE no depende de la longitud, corregir no aporta nada; si
#' depende, el enriquecimiento GO hay que leerlo con cautela. Se calcula tambien
#' la funcion de ponderacion por probabilidad (la curva que dibuja `goseq`), sin
#' anadir dependencia.
#'
#' @param deg_df tabla DEG estandar
#' @param lengths vector de longitudes con nombres = identificadores de gen
#' @param fdr umbral de significacion
#' @param n_bins numero de bins de longitud para la curva
#' @return list(table, spearman, p_value, n_used, odds_ratio, verdict)
length_bias_diagnostic <- function(deg_df, lengths, fdr = 0.05, n_bins = 12L) {
  out <- list(table = NULL, spearman = NA_real_, p_value = NA_real_,
              n_used = 0L, odds_ratio = NA_real_, verdict = NULL)
  if (is.null(deg_df) || !nrow(deg_df) || is.null(lengths) || !length(lengths)) return(out)
  idx <- match(deg_df$gene, names(lengths))
  ok <- !is.na(idx) & !is.na(deg_df$padj)
  if (sum(ok) < 50) return(out)

  len <- as.numeric(lengths[idx[ok]])
  de  <- as.integer(deg_df$padj[ok] <= fdr)
  out$n_used <- sum(ok)

  # Correlacion de rangos entre longitud y evidencia de expresion diferencial.
  # Se usa el rango del p-valor invertido para que "mas alto = mas DE".
  ev <- -log10(pmax(deg_df$pvalue[ok], .Machine$double.xmin))
  ct <- suppressWarnings(tryCatch(
    stats::cor.test(len, ev, method = "spearman", exact = FALSE),
    error = function(e) NULL))
  if (!is.null(ct)) {
    out$spearman <- as.numeric(ct$estimate)
    out$p_value  <- as.numeric(ct$p.value)
  }

  # Curva de proporcion de DE por bin de longitud (la PWF de goseq).
  brks <- unique(stats::quantile(len, probs = seq(0, 1, length.out = n_bins + 1L),
                                 na.rm = TRUE))
  if (length(brks) < 3) return(out)
  bin <- cut(len, breaks = brks, include.lowest = TRUE, labels = FALSE)
  df <- data.frame(
    bin = sort(unique(bin[!is.na(bin)])),
    stringsAsFactors = FALSE
  )
  df$n_genes  <- as.integer(tapply(de, bin, length)[as.character(df$bin)])
  df$n_de     <- as.integer(tapply(de, bin, sum)[as.character(df$bin)])
  df$prop_de  <- df$n_de / df$n_genes
  df$len_median <- as.numeric(tapply(len, bin, stats::median)[as.character(df$bin)])
  out$table <- df

  # Odds ratio entre el tercio mas largo y el mas corto: cuantifica el sesgo en
  # una cifra interpretable.
  q <- stats::quantile(len, c(1/3, 2/3), na.rm = TRUE)
  short <- len <= q[1]; long <- len >= q[2]
  a <- sum(de[long]); b <- sum(!de[long]); c <- sum(de[short]); d <- sum(!de[short])
  if (a > 0 && b > 0 && c > 0 && d > 0) out$odds_ratio <- (a / b) / (c / d)

  out$verdict <- interpret_length_bias(out$spearman, out$p_value, out$odds_ratio)
  out
}

#' Traduce el diagnostico de sesgo de longitud a una recomendacion concreta.
interpret_length_bias <- function(rho, p_value, odds_ratio) {
  if (is.na(rho)) {
    return(list(level = "desconocido", label = "No evaluable",
                detail = "No hay suficientes genes con longitud conocida."))
  }
  strong <- abs(rho) >= 0.15 || (!is.na(odds_ratio) && (odds_ratio >= 1.5 || odds_ratio <= 1/1.5))
  signif <- !is.na(p_value) && p_value < 0.05
  or_txt <- if (is.na(odds_ratio)) "" else paste0(
    " Los genes del tercio mas largo tienen ", round(odds_ratio, 2),
    " veces las probabilidades de salir diferenciales que los del tercio mas corto.")
  if (strong && signif) {
    return(list(level = "aviso", label = "Sesgo de longitud detectado",
      detail = paste0(
        "La correlacion entre longitud del gen y evidencia de expresion ",
        "diferencial es rho = ", round(rho, 3), " (p = ", signif(p_value, 3), ").",
        or_txt, " El analisis de sobre-representacion (GO/KEGG) esta sesgado hacia ",
        "las categorias que contienen genes largos: interpretalo con cautela y ",
        "prioriza GSEA, que usa el ranking completo, o una correccion por longitud ",
        "tipo goseq.")))
  }
  if (signif) {
    return(list(level = "leve", label = "Sesgo de longitud leve",
      detail = paste0(
        "La correlacion es estadisticamente detectable (rho = ", round(rho, 3),
        ", p = ", signif(p_value, 3), ") pero pequena.", or_txt,
        " Con decenas de miles de genes casi cualquier correlacion sale ",
        "significativa; la magnitud importa mas que el p-valor.")))
  }
  list(level = "ok", label = "Sin sesgo de longitud apreciable",
    detail = paste0(
      "No hay relacion apreciable entre longitud del gen y evidencia de expresion ",
      "diferencial (rho = ", round(rho, 3), ", p = ", signif(p_value %||% NA, 3), ").",
      or_txt, " Corregir por longitud no aportaria nada en este dataset."))
}

#' Resumen RLE por muestra: log2CPM menos la mediana del gen.
#'
#' Se devuelven cuantiles y no los valores crudos porque con decenas de miles de
#' genes el grafico no necesita los puntos. Una muestra bien normalizada tiene la
#' mediana cerca de 0 y un IQR estrecho; una mediana desplazada indica fallo de
#' normalizacion o un problema con esa muestra.
rle_summary <- function(counts) {
  if (is.null(counts) || !nrow(counts) || !ncol(counts)) return(NULL)
  cm <- as.matrix(counts)
  libs <- colSums(cm, na.rm = TRUE)
  libs[libs == 0 | is.na(libs)] <- 1
  lcpm <- log2(t(t(cm) / libs) * 1e6 + 1)
  med <- apply(lcpm, 1, stats::median, na.rm = TRUE)
  rle <- lcpm - med
  qs <- apply(rle, 2, stats::quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
              na.rm = TRUE)
  df <- data.frame(
    sample_id = colnames(rle),
    p05 = qs[1, ], q1 = qs[2, ], med = qs[3, ], q3 = qs[4, ], p95 = qs[5, ],
    stringsAsFactors = FALSE
  )
  df$iqr <- df$q3 - df$q1
  # Muestras cuya mediana se aleja de 0 mas que el doble del IQR tipico.
  thr <- max(0.1, 2 * stats::median(df$iqr, na.rm = TRUE))
  df$flag <- abs(df$med) > thr
  attr(df, "median_threshold") <- thr
  df
}
