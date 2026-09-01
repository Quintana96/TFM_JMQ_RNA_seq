#' test-diseno-pareado.R
#' El tipado de las variables del diseño: factor o covariable continua.
#'
#' `is_continuous_var()` decide por heuristica —numerica y con cinco o más
#' valores distintos— y esa heuristica se equivoca justo en el diseño pareado,
#' que es el caso que motivó la formula libre. Un identificador de sujeto
#' codificado 1..6 es numerico y tiene seis valores distintos, así que pasa por
#' covariable continua: el modelo le ajusta UNA PENDIENTE LINEAL en lugar de un
#' bloque por sujeto.
#'
#' El fallo no es ruidoso. No hay error, no hay aviso y sale una tabla
#' plausible; simplemente se ha ajustado otro modelo. Por eso el arreglo no es
#' cambiar la heuristica —seguiria adivinando— sino permitir declararlo, y que
#' lo declarado gane.

#' Fixture pareado: seis sujetos, cada uno medido en las dos condiciones.
#' Seis y no cuatro porque `is_continuous_var()` exige cinco valores distintos
#' para dar una columna por continua: con cuatro sujetos el defecto no aparece
#' y el test pasaria sin probar nada.
meta_pareado <- function(counts) {
  meta <- make_test_meta(counts)
  n_por_grupo <- sum(meta$condition == "ctrl")
  meta$subject <- c(seq_len(n_por_grupo), seq_len(n_por_grupo))
  meta
}

test_that("sin declarar nada, un identificador de sujeto pasa por continuo", {
  counts <- make_test_counts(n_per_group = 6)
  meta <- meta_pareado(counts)

  expect_length(unique(meta$subject), 6)          # el fixture destapa el caso
  expect_true(is_continuous_var(meta$subject))

  m <- prepare_design_meta(meta, c("subject", "condition"))
  expect_identical(unname(attr(m, "var_types")[["subject"]]), "continua")
  expect_true(is.numeric(m$subject))
})

test_that("declarar el tipado gana a la heuristica en los dos sentidos", {
  counts <- make_test_counts(n_per_group = 6)
  meta <- meta_pareado(counts)

  # character(0) = "ninguna es continua". Antes esto no se podia expresar:
  # `continuous` solo AÑADIA, y la heuristica seguia decidiendo.
  m_factor <- prepare_design_meta(meta, c("subject", "condition"), character(0))
  expect_identical(unname(attr(m_factor, "var_types")[["subject"]]), "factor")
  expect_true(is.factor(m_factor$subject))

  # Y al reves: lo declarado continuo lo es aunque la heuristica dijera factor.
  meta$dosis <- rep(c(0, 1), length.out = nrow(meta))   # 2 niveles: factor
  expect_false(is_continuous_var(meta$dosis))
  m_cont <- prepare_design_meta(meta, c("dosis", "condition"), "dosis")
  expect_identical(unname(attr(m_cont, "var_types")[["dosis"]]), "continua")
  expect_true(is.numeric(m_cont$dosis))
})

test_that("el tipado cambia el RANGO de la matriz de diseño, no solo la etiqueta", {
  counts <- make_test_counts(n_per_group = 6)
  meta <- meta_pareado(counts)

  cont <- validate_design_formula("~ subject + condition", meta)             # heuristica
  fact <- validate_design_formula("~ subject + condition", meta, character(0))

  expect_true(cont$ok)
  expect_true(fact$ok)

  # Pendiente lineal: intercepto + subject + condition = 3 columnas.
  expect_identical(cont$rank, 3L)
  # Bloque por sujeto: intercepto + 5 indicadoras + condition = 7 columnas.
  expect_identical(fact$rank, 7L)

  # Y por tanto distintos grados de libertad residuales: es OTRO modelo.
  expect_gt(cont$residual_df, fact$residual_df)
})

test_that("solo se ofrecen como continuas las columnas numericas", {
  counts <- make_test_counts(n_per_group = 6)
  meta <- meta_pareado(counts)
  meta$lote <- rep(c("A", "B"), length.out = nrow(meta))

  cand <- design_numeric_vars(meta)
  expect_true("subject" %in% cand)
  expect_false("lote" %in% cand)        # texto: no puede ser continua
  expect_false("condition" %in% cand)   # es el factor del contraste
  expect_false("sample_id" %in% cand)
})

test_that("run_deg propaga el tipado hasta el ajuste real", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  counts <- make_test_counts(n_per_group = 6)
  meta <- meta_pareado(counts)

  como_pendiente <- run_deg(counts, meta, method = "DESeq2", ref_level = "ctrl",
                            contrast_num = "trt", design_formula = "~ subject + condition",
                            shrink = FALSE)
  como_bloque <- run_deg(counts, meta, method = "DESeq2", ref_level = "ctrl",
                         contrast_num = "trt", design_formula = "~ subject + condition",
                         shrink = FALSE, continuous = character(0))

  expect_null(como_pendiente$error)
  expect_null(como_bloque$error)

  # No basta con que el argumento llegue: tiene que cambiar el ajuste. Si los
  # p-valores fueran identicos, `continuous` no estaria conectado a nada.
  p1 <- como_pendiente$table$pvalue
  p2 <- como_bloque$table$pvalue
  expect_equal(length(p1), length(p2))
  expect_false(isTRUE(all.equal(p1, p2)))
})
