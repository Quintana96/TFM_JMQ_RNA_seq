#' utils_design.R
#' Funciones puras (sin Shiny) para construir y VALIDAR diseños experimentales
#' arbitrarios a partir de una formula escrita por el usuario.
#'
#' Por que existe (docs/REVISION_ESTADISTICA.md, B5): la app solo generaba
#' `~ condition` o `~ batch + condition`, con todo convertido a factor. Eso deja
#' fuera casos muy comunes:
#'   - diseños pareados (`~ subject + condition`), probablemente el diseño mas
#'     frecuente que no se podia analizar, y el que mas potencia gana al
#'     modelarse bien cuando hay pocas replicas;
#'   - covariables continuas (edad, dosis, tiempo), que con `as.factor()` gastan
#'     grados de libertad absurdamente;
#'   - interacciones (`~ genotipo * tratamiento`), es decir la pregunta "el
#'     efecto del tratamiento depende del genotipo?";
#'   - varias covariables a la vez.
#'
#' La parte importante no es aceptar la formula, es VALIDARLA antes de ajustar:
#' un diseño singular produce un error criptico de DESeq2 en ingles, y una
#' confusion perfecta entre dos covariables produce resultados que parecen
#' correctos y no lo son.

#' Columnas del samplesheet utilizables como terminos del diseño.
design_candidate_vars <- function(meta) {
  if (is.null(meta) || !is.data.frame(meta)) return(character(0))
  setdiff(names(meta), c("sample_id"))
}

#' TRUE si la columna debe tratarse como covariable continua.
#'
#' Criterio: numerica y con mas de `min_levels` valores distintos. Una columna
#' numerica con dos o tres valores (p. ej. dosis 0/1/2 o un batch codificado
#' 1/2) casi siempre se quiere como factor, asi que se deja como factor.
is_continuous_var <- function(x, min_levels = 5L) {
  if (is.factor(x) || is.character(x) || is.logical(x)) return(FALSE)
  v <- suppressWarnings(as.numeric(x))
  if (all(is.na(v))) return(FALSE)
  length(unique(v[!is.na(v)])) >= min_levels
}

#' Prepara `meta` para el diseño: convierte a factor lo que no sea continuo.
#' Devuelve tambien el tipo asignado a cada variable, para poder mostrarlo.
prepare_design_meta <- function(meta, vars = NULL, continuous = NULL) {
  if (is.null(meta)) return(NULL)
  vars <- vars %||% design_candidate_vars(meta)
  types <- character(0)
  for (v in intersect(vars, names(meta))) {
    force_cont <- !is.null(continuous) && v %in% continuous
    if (force_cont || is_continuous_var(meta[[v]])) {
      meta[[v]] <- suppressWarnings(as.numeric(meta[[v]]))
      types[v] <- "continua"
    } else {
      meta[[v]] <- as.factor(as.character(meta[[v]]))
      types[v] <- "factor"
    }
  }
  attr(meta, "var_types") <- types
  meta
}

#' Detecta pares de variables categoricas cuyos efectos no se pueden separar.
#'
#' Dos situaciones distintas, ambas causa de singularidad, y el error de DESeq2
#' ("model matrix is not full rank") no dice cual es el par culpable:
#'
#'   - CONFUSION: cada nivel de una se corresponde con exactamente uno de la otra
#'     (biyeccion). Son la misma variable con dos nombres.
#'   - ANIDAMIENTO: cada nivel de A determina el de B, pero no al reves (p. ej.
#'     cada sujeto pertenece a un solo batch, y cada batch tiene varios sujetos).
#'     Es el caso mas frecuente en la practica y tambien rompe el rango.
#'
#' @return lista de list(vars = c(a, b), kind = "confusion"|"anidamiento", detail)
find_confounded_pairs <- function(meta, vars) {
  vars <- intersect(vars, names(meta))
  cats <- vars[vapply(vars, function(v) !is.numeric(meta[[v]]), logical(1))]
  out <- list()
  if (length(cats) < 2) return(out)
  for (i in seq_len(length(cats) - 1L)) {
    for (j in seq(i + 1L, length(cats))) {
      a <- as.character(meta[[cats[i]]]); b <- as.character(meta[[cats[j]]])
      ok <- !is.na(a) & !is.na(b)
      if (!any(ok)) next
      tb <- table(a[ok], b[ok])
      if (nrow(tb) < 2 || ncol(tb) < 2) next
      a_determines_b <- all(rowSums(tb > 0) <= 1)
      b_determines_a <- all(colSums(tb > 0) <= 1)
      if (a_determines_b && b_determines_a) {
        out[[length(out) + 1L]] <- list(
          vars = c(cats[i], cats[j]), kind = "confusion",
          detail = paste0("'", cats[i], "' y '", cats[j],
                          "' son la misma particion de las muestras"))
      } else if (a_determines_b) {
        out[[length(out) + 1L]] <- list(
          vars = c(cats[i], cats[j]), kind = "anidamiento",
          detail = paste0("cada nivel de '", cats[i], "' cae siempre en el mismo '",
                          cats[j], "'"))
      } else if (b_determines_a) {
        out[[length(out) + 1L]] <- list(
          vars = c(cats[j], cats[i]), kind = "anidamiento",
          detail = paste0("cada nivel de '", cats[j], "' cae siempre en el mismo '",
                          cats[i], "'"))
      }
    }
  }
  out
}

#' Texto con los pares problematicos, para incrustar en un mensaje de error.
confounded_pairs_text <- function(pairs) {
  if (!length(pairs)) return("")
  paste(vapply(pairs, function(p) p$detail, character(1)), collapse = "; ")
}

#' Valida una formula de diseño contra el samplesheet SIN ajustar el modelo.
#'
#' Comprueba, en este orden:
#'   1. que la formula sea sintacticamente valida y solo tenga lado derecho;
#'   2. que todas sus variables existan en el samplesheet;
#'   3. que ninguna tenga NA en las muestras que se van a usar;
#'   4. que `model.matrix()` se pueda construir;
#'   5. que la matriz tenga RANGO COMPLETO (si no, señala los pares
#'      confundidos que lo explican);
#'   6. que queden grados de libertad residuales (n - rango > 0).
#'
#' Devuelve list(ok, errors, warnings, formula, design, meta, coef_names,
#' rank, n_samples, residual_df, var_types).
validate_design_formula <- function(formula_text, meta, continuous = NULL) {
  res <- list(ok = FALSE, errors = character(0), warnings = character(0),
              formula = NULL, design = NULL, meta = NULL,
              coef_names = character(0), rank = NA_integer_,
              n_samples = NA_integer_, residual_df = NA_integer_,
              var_types = character(0))
  txt <- trimws(as.character(formula_text %||% ""))
  if (!nzchar(txt)) {
    res$errors <- "La formula esta vacia."
    return(res)
  }
  if (!startsWith(txt, "~")) txt <- paste0("~ ", txt)
  if (grepl("~[^~]*~", txt)) {
    res$errors <- "La formula solo debe tener lado derecho (empezar por ~)."
    return(res)
  }
  f <- tryCatch(stats::as.formula(txt), error = function(e) NULL)
  if (is.null(f) || length(f) != 2L) {
    res$errors <- paste0("Formula no valida: '", txt, "'.")
    return(res)
  }
  res$formula <- f

  vars <- all.vars(f)
  if (!length(vars)) {
    res$errors <- "La formula no menciona ninguna variable."
    return(res)
  }
  missing <- setdiff(vars, names(meta))
  if (length(missing)) {
    res$errors <- paste0("Estas columnas no estan en el samplesheet: ",
                         paste(missing, collapse = ", "), ".")
    return(res)
  }

  m <- prepare_design_meta(meta, vars, continuous)
  res$var_types <- attr(m, "var_types") %||% character(0)

  has_na <- vapply(vars, function(v) any(is.na(m[[v]])), logical(1))
  if (any(has_na)) {
    res$errors <- paste0("Hay valores vacios en: ",
                         paste(vars[has_na], collapse = ", "),
                         ". Completa el samplesheet o quita esas muestras.")
    return(res)
  }
  # Un factor con un solo nivel no aporta nada y rompe model.matrix.
  one_level <- vapply(vars, function(v) {
    is.factor(m[[v]]) && nlevels(droplevels(m[[v]])) < 2
  }, logical(1))
  if (any(one_level)) {
    res$errors <- paste0("Estas variables tienen un solo nivel: ",
                         paste(vars[one_level], collapse = ", "), ".")
    return(res)
  }

  mm <- tryCatch(stats::model.matrix(f, data = m), error = function(e) e)
  if (inherits(mm, "error")) {
    res$errors <- paste0("model.matrix fallo: ", conditionMessage(mm))
    return(res)
  }
  res$design     <- mm
  res$meta       <- m
  res$coef_names <- colnames(mm)
  res$n_samples  <- nrow(mm)
  res$rank       <- as.integer(qr(mm)$rank)
  res$residual_df <- as.integer(nrow(mm) - res$rank)

  if (res$rank < ncol(mm)) {
    conf <- find_confounded_pairs(m, vars)
    detail <- if (length(conf)) paste0(" Causa probable: ",
                                       confounded_pairs_text(conf), ".") else ""
    res$errors <- paste0(
      "El diseño es singular: la matriz tiene ", ncol(mm), " coeficientes pero ",
      "rango ", res$rank, ", asi que sus efectos no se pueden separar.", detail
    )
    return(res)
  }
  if (res$residual_df <= 0) {
    res$errors <- paste0(
      "No quedan grados de libertad residuales (", res$n_samples,
      " muestras y ", res$rank, " coeficientes). Simplifica el diseño o añade ",
      "muestras."
    )
    return(res)
  }
  if (res$residual_df < 2) {
    res$warnings <- c(res$warnings, paste0(
      "Solo queda ", res$residual_df, " grado de libertad residual: la ",
      "estimacion de la dispersion sera muy inestable."))
  }
  # Asociacion fuerte que no llega a romper el rango: conviene avisar igual.
  conf <- find_confounded_pairs(m, vars)
  if (length(conf)) {
    res$warnings <- c(res$warnings, paste0(
      "Variables fuertemente asociadas entre si: ", confounded_pairs_text(conf),
      ". Sus efectos se estimaran con poca precision."))
  }
  res$ok <- TRUE
  res
}

#' Resumen legible de una validacion, para mostrar en la interfaz.
design_summary_text <- function(v) {
  if (is.null(v)) return(NULL)
  if (!isTRUE(v$ok)) return(paste(v$errors, collapse = " "))
  types <- if (length(v$var_types)) paste0(
    "  ·  ", paste(paste0(names(v$var_types), " (", v$var_types, ")"),
                   collapse = ", ")) else ""
  paste0(v$n_samples, " muestras  ·  ", length(v$coef_names), " coeficientes  ·  ",
         v$residual_df, " g.l. residuales", types)
}

