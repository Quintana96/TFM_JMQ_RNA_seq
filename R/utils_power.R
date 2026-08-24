#' utils_power.R
#' Calculo de potencia y tamaño muestral a priori para RNA-seq.
#'
#' Por que existe (docs/REVISION_ESTADISTICA.md, B8): el tamaño muestral es el
#' determinante mas fuerte de la calidad del resultado y ningun componente de la
#' app orientaba sobre el. Schurch et al. (2016), con 48 replicas por condicion
#' en levadura, concluyeron que hacen falta >= 6 replicas por condicion para
#' deteccion robusta y >= 12 para capturar la mayoria de los DEG a todos los
#' fold-changes; con 3 replicas, 9 de 11 herramientas recuperaban solo el 20-40 %
#' de los genes que se detectan con 42.
#'
#' LIMITACION QUE HAY QUE MOSTRAR, no esconder: la revision en Briefings in
#' Bioinformatics concluye que ninguna herramienta de calculo de tamaño muestral
#' es fiable cuando se exigen efectos pequeños y potencias altas, porque no se
#' pueden fijar parametros razonables a partir de datos piloto limitados. El
#' resultado es una orientacion, no una garantia.

#' Referencias empiricas de Schurch et al. (2016), para contrastar con el calculo.
POWER_REFERENCE_NOTE <- paste(
  "Referencia empirica (Schurch et al., RNA 2016, con 48 replicas en levadura):",
  ">= 6 replicas por condicion para deteccion robusta, >= 12 para capturar la",
  "mayoria de los DEG a todos los fold-changes. Con 3 replicas se recupera",
  "tipicamente entre el 20 % y el 40 % de lo detectable."
)

#' Estima los parametros de potencia A PARTIR de la matriz cargada.
#'
#' Pedir el coeficiente de variacion y la profundidad "a ojo" es la parte mas
#' fragil del calculo: son justo los valores que el usuario no conoce, y de los
#' que depende todo el resultado. Si hay una matriz de conteos cargada, ambos se
#' pueden medir en lugar de adivinarse:
#'
#'   - `cv` es la raiz cuadrada de la dispersion biologica comun (BCV) que
#'     estima edgeR, que es exactamente la definicion del coeficiente de
#'     variacion biologico que espera RNASeqPower.
#'   - `depth` es la mediana de conteos por gen entre los genes expresados; usar
#'     la media la infla por unos pocos genes muy expresados.
#'
#' @param counts matriz de conteos
#' @param meta samplesheet con `condition` (para el diseño del BCV)
#' @return list(cv, depth, n_por_grupo, origen) o NULL si no se puede estimar
estimate_power_params <- function(counts, meta = NULL) {
  if (is.null(counts) || !nrow(counts) || ncol(counts) < 2) return(NULL)
  cm <- round(as.matrix(counts))
  out <- tryCatch({
    depth <- stats::median(cm[cm > 0], na.rm = TRUE)
    cv <- NA_real_
    if (requireNamespace("edgeR", quietly = TRUE)) {
      grupo <- if (!is.null(meta) && "condition" %in% names(meta) &&
                   length(unique(stats::na.omit(meta$condition))) > 1) {
        as.factor(as.character(meta$condition))
      } else NULL
      y <- edgeR::DGEList(counts = cm, group = grupo)
      keep <- edgeR::filterByExpr(y, group = grupo)
      # El metodo `[` de DGEList no admite `drop`.
      if (sum(keep) > 50) y <- y[keep, ]
      y <- norm_lib_sizes(y)
      design <- if (!is.null(grupo)) {
        stats::model.matrix(~ grupo)
      } else {
        matrix(1, ncol(y), 1)
      }
      y <- edgeR::estimateDisp(y, design, robust = TRUE)
      # BCV = sqrt(dispersion comun); es el CV biologico.
      cv <- sqrt(y$common.dispersion %||% NA_real_)
    }
    n_grp <- if (!is.null(meta) && "condition" %in% names(meta)) {
      tb <- table(as.character(meta$condition))
      if (length(tb)) min(tb) else NA_integer_
    } else NA_integer_
    list(cv = cv, depth = depth, n_por_grupo = n_grp,
         origen = "estimados de la matriz cargada")
  }, error = function(e) NULL)
  out
}

#' Potencia para un tamaño muestral dado, via RNASeqPower.
#'
#' @param n replicas por grupo
#' @param cv coeficiente de variacion biologico (0,1 lineas celulares; 0,4 humano)
#' @param effect fold-change minimo a detectar (escala lineal, p. ej. 2)
#' @param depth profundidad media por gen (conteos)
#' @param alpha nivel de significacion por test
#' @return list(power, error)
power_for_n <- function(n, cv = 0.4, effect = 2, depth = 20, alpha = 0.05) {
  if (!requireNamespace("RNASeqPower", quietly = TRUE)) {
    return(list(power = NA_real_, error = "RNASeqPower no esta instalado."))
  }
  out <- tryCatch({
    p <- RNASeqPower::rnapower(depth = depth, n = n, cv = cv,
                               effect = effect, alpha = alpha)
    list(power = as.numeric(p)[1], error = NULL)
  }, error = function(e) list(power = NA_real_, error = conditionMessage(e)))
  out
}

#' Replicas necesarias para alcanzar una potencia objetivo.
n_for_power <- function(power = 0.8, cv = 0.4, effect = 2, depth = 20,
                       alpha = 0.05) {
  if (!requireNamespace("RNASeqPower", quietly = TRUE)) {
    return(list(n = NA_real_, error = "RNASeqPower no esta instalado."))
  }
  out <- tryCatch({
    v <- RNASeqPower::rnapower(depth = depth, cv = cv, effect = effect,
                               alpha = alpha, power = power)
    list(n = as.numeric(v)[1], error = NULL)
  }, error = function(e) list(n = NA_real_, error = conditionMessage(e)))
  out
}

#' Curva potencia-vs-n, para dibujar.
#' Se calcula punto a punto para que un fallo en un n no tumbe la curva entera.
power_curve <- function(n_range = 2:20, cv = 0.4, effect = 2, depth = 20,
                        alpha = 0.05) {
  pw <- vapply(n_range, function(n) {
    r <- power_for_n(n, cv = cv, effect = effect, depth = depth, alpha = alpha)
    r$power %||% NA_real_
  }, numeric(1))
  df <- data.frame(n = n_range, power = pw, stringsAsFactors = FALSE)
  if (all(is.na(df$power))) return(NULL)
  df
}

#' Interpretacion del resultado, cruzando el calculo con la referencia empirica.
#'
#' El calculo puede decir que 3 replicas bastan para un fold-change de 4 y ser
#' formalmente correcto, mientras la evidencia empirica dice que con 3 replicas
#' se recupera una fraccion pequeña de los DEG reales. Conviene decir las dos
#' cosas.
interpret_power <- function(n, power) {
  if (is.na(power)) {
    return(list(level = "desconocida", label = "No se ha podido calcular",
                detail = "Revisa los parametros introducidos."))
  }
  pct <- round(100 * power, 1)
  if (n < 6) {
    return(list(level = "aviso",
      label = paste0("Potencia estimada ", pct, " % con n = ", n),
      detail = paste("Aunque el calculo salga favorable, con menos de 6 replicas",
                     "por condicion la evidencia empirica es desfavorable.",
                     POWER_REFERENCE_NOTE)))
  }
  if (power >= 0.8) {
    return(list(level = "ok",
      label = paste0("Potencia estimada ", pct, " % con n = ", n),
      detail = paste("Por encima del 80 % habitual y con al menos 6 replicas.",
                     POWER_REFERENCE_NOTE)))
  }
  list(level = "bajo",
    label = paste0("Potencia estimada ", pct, " % con n = ", n),
    detail = paste("Por debajo del 80 % habitual: el experimento detectara solo",
                   "una parte de los genes realmente diferenciales.",
                   POWER_REFERENCE_NOTE))
}
