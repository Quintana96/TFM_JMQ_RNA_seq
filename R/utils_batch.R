#' utils_batch.R
#' Funciones puras (sin Shiny) para tratar la variacion no deseada: variables
#' sustitutas (sva), ajuste de conteos (ComBat-seq) y eliminacion del efecto
#' batch para VISUALIZAR (limma::removeBatchEffect).
#'
#' La distincion que la app debe ensenar explicitamente (docs/REVISION_ESTADISTICA.md,
#' B6), porque es una fuente clasica de error:
#'
#'   - Para TESTEAR, lo correcto es incluir el batch como covariable en el
#'     diseno y dejar los conteos intactos. El modelo estima y descuenta el
#'     efecto sin destruir la estructura de la varianza.
#'   - Para VISUALIZAR (PCA, heatmap), se corrige la matriz, porque un grafico no
#'     puede "incluir una covariable".
#'
#' Corregir la matriz y luego testear sobre ella infla los falsos positivos: el
#' test no sabe que los datos ya han sido ajustados y cuenta como reales unos
#' grados de libertad que ya se gastaron.

#' Estima variables sustitutas con svaseq.
#'
#' Es lo apropiado cuando hay estructura extra en los datos pero no se sabe que
#' la causa: sva la estima a partir de los propios datos y devuelve covariables
#' que se anaden al diseno.
#'
#' @param counts matriz de conteos (genes x muestras)
#' @param meta metadatos alineados con las columnas de counts
#' @param full_formula formula del modelo de interes (p. ej. ~ condition)
#' @param n_sv numero de variables sustitutas; NULL = lo estima sva
#' @param min_residual_df grados de libertad residuales que se reservan. Cada
#'   variable sustituta consume uno, y `num.sv` con el metodo de Leek tiende a
#'   proponer tantas que el diseno se queda sin capacidad de estimar nada: con 12
#'   muestras llego a proponer 9, que sumadas a la condicion dejan 1 g.l.
#'   residual. Se usa el metodo "be" (por permutacion, mas conservador) y se
#'   recorta ademas para respetar este minimo.
#' @return list(sv = matriz muestras x n_sv, n_sv, n_sv_estimated, error)
estimate_surrogate_vars <- function(counts, meta, full_formula = ~ condition,
                                    n_sv = NULL, min_residual_df = 3L) {
  if (!requireNamespace("sva", quietly = TRUE)) {
    return(list(sv = NULL, n_sv = 0L, error = "sva no esta instalado."))
  }
  out <- tryCatch({
    cm <- as.matrix(counts)
    # svaseq trabaja sobre datos no normalizados pero pide quitar los genes sin
    # expresion, porque un gen constante no aporta y desestabiliza la SVD.
    keep <- rowMeans(cm) > 1
    cm <- cm[keep, , drop = FALSE]
    if (!nrow(cm)) stop("No quedan genes con expresion suficiente para sva.")
    mod  <- stats::model.matrix(full_formula, data = meta)
    mod0 <- stats::model.matrix(~ 1, data = meta)
    n_est <- if (is.null(n_sv) || !is.finite(n_sv) || n_sv < 1) {
      tryCatch(sva::num.sv(cm, mod, method = "be"), error = function(e) 1L)
    } else as.integer(n_sv)
    # Tope: no dejar el diseno sin grados de libertad residuales.
    n_max <- ncol(cm) - ncol(mod) - as.integer(min_residual_df)
    if (n_max < 1) {
      stop(paste0("Con ", ncol(cm), " muestras y ", ncol(mod), " coeficientes no ",
                  "hay margen para variables sustitutas."))
    }
    n <- max(1L, min(as.integer(n_est), n_max))
    sv <- sva::svaseq(cm, mod, mod0, n.sv = n)$sv
    sv <- as.matrix(sv)
    colnames(sv) <- paste0("SV", seq_len(ncol(sv)))
    rownames(sv) <- colnames(counts)
    list(sv = sv, n_sv = ncol(sv), n_sv_estimated = as.integer(n_est), error = NULL)
  }, error = function(e) {
    msg <- conditionMessage(e)
    # sva falla con errores de algebra lineal cuando hay pocas muestras. El
    # mensaje crudo ("Lapack routine dgesv: system is exactly singular") no le
    # dice nada a quien usa la app.
    if (grepl("singular|Lapack|dgesv", msg, ignore.case = TRUE)) {
      msg <- paste0("sva no ha podido estimar variables sustitutas con ",
                    ncol(counts), " muestras: el sistema queda singular. ",
                    "Hacen falta mas muestras, o indica el batch a mano si lo ",
                    "conoces. (Detalle: ", msg, ")")
    }
    list(sv = NULL, n_sv = 0L, n_sv_estimated = NA_integer_, error = msg)
  })
  out
}

#' Ajusta los conteos con ComBat-seq.
#'
#' A diferencia de ComBat clasico, usa regresion binomial negativa y devuelve
#' una matriz de conteos ENTERA, compatible con DESeq2/edgeR aguas abajo
#' (Zhang, Parmigiani y Johnson, 2020). Aun asi, el uso recomendado es
#' visualizar: para testear, el batch va en el diseno.
combat_seq_counts <- function(counts, batch, group = NULL) {
  if (!requireNamespace("sva", quietly = TRUE)) {
    return(list(counts = NULL, error = "sva no esta instalado."))
  }
  out <- tryCatch({
    cm <- round(as.matrix(counts))
    b <- as.factor(as.character(batch))
    if (nlevels(b) < 2) stop("La variable de batch necesita al menos 2 niveles.")
    if (length(b) != ncol(cm)) stop("El batch no tiene una entrada por muestra.")
    # ComBat-seq falla con genes de conteo cero en todas las muestras.
    keep <- rowSums(cm) > 0
    adj <- sva::ComBat_seq(cm[keep, , drop = FALSE], batch = b,
                           group = if (!is.null(group)) as.factor(as.character(group)) else NULL)
    full <- cm
    full[keep, ] <- adj
    list(counts = full, error = NULL, n_genes = sum(keep))
  }, error = function(e) list(counts = NULL, error = conditionMessage(e)))
  out
}

#' Elimina el efecto batch de una matriz transformada, SOLO para graficos.
#'
#' `design` preserva los efectos de interes: sin el, removeBatchEffect tambien
#' se llevaria por delante la senal de la condicion.
remove_batch_for_plots <- function(mat, batch, design = NULL, covariates = NULL) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    return(list(mat = mat, error = "limma no esta instalado."))
  }
  out <- tryCatch({
    m <- as.matrix(mat)
    args <- list(x = m)
    if (!is.null(batch)) args$batch <- as.factor(as.character(batch))
    if (!is.null(covariates)) args$covariates <- as.matrix(covariates)
    if (!is.null(design)) args$design <- design
    if (is.null(args$batch) && is.null(args$covariates)) {
      stop("Hace falta una variable de batch o covariables para corregir.")
    }
    list(mat = do.call(limma::removeBatchEffect, args), error = NULL)
  }, error = function(e) list(mat = mat, error = conditionMessage(e)))
  out
}
