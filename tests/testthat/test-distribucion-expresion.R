#' test-distribución-expresión.R
#' La densidad y el diagrama de cajas de la expresión son diagnósticos, no
#' adornos: se leen ANTES de ajustar el modelo (vinieta de limma-voom, workflow
#' de edgeR) y de ellos depende que el usuario decida prefiltrar o descartar una
#' muestra. Lo que se verifica aquí es que lo que dibujan es comparable entre
#' muestras y que lo que dejan fuera se declara.

test_that("los genes sin conteo en ninguna muestra se excluyen y se cuentan", {
  cm <- make_test_counts(n_genes = 200, n_per_group = 3)
  cm[1:20, ] <- 0L

  lc <- expression_logcpm(cm)
  expect_equal(lc$n_total, 200L)
  expect_equal(lc$n_dropped, 20L)
  expect_equal(lc$n_genes, 180L)
  # Un gen con ceros en ALGUNAS muestras se conserva: ahí el cero es un dato,
  # no ausencia de información.
  cm2 <- cm; cm2[21, 1] <- 0L
  expect_equal(expression_logcpm(cm2)$n_genes, 180L)
})

test_that("todas las muestras comparten rejilla y ancho de banda", {
  # Es la propiedad que hace comparable el gráfico: density() elige por defecto
  # un ancho por vector, de modo que la muestra menos dispersa saldría más
  # picuda por construcción y no porque su distribución lo sea.
  cm <- make_test_counts(n_genes = 300, n_per_group = 3)
  d <- expression_distribution(cm)

  expect_false(is.null(d))
  rejillas <- split(d$density$x, d$density$sample_id)
  expect_true(length(rejillas) == ncol(cm))
  for (r in rejillas[-1]) expect_equal(r, rejillas[[1]])
  expect_true(is.finite(d$bw) && d$bw > 0)
})

test_that("los cuantiles del diagrama de cajas están ordenados", {
  cm <- make_test_counts(n_genes = 300, n_per_group = 3)
  b <- expression_distribution(cm)$box

  expect_equal(nrow(b), ncol(cm))
  expect_true(all(b$p05 <= b$q1))
  expect_true(all(b$q1  <= b$med))
  expect_true(all(b$med <= b$q3))
  expect_true(all(b$q3  <= b$p95))
  expect_true(all(b$iqr >= 0))
})

test_that("la escala normalizada y la cruda no dan lo mismo con librerias dispares", {
  # Se duplica la profundidad de la mitad de las muestras: sin normalizar, sus
  # distribuciones se separan; normalizando por composicion, se acercan. Si las
  # dos vistas coincidieran, el selector de escala no estaria diciendo nada.
  cm <- make_test_counts(n_genes = 400, n_per_group = 3)
  cm[, 4:6] <- cm[, 4:6] * 3L

  crudo <- expression_distribution(cm, normalized = FALSE)$box
  norm  <- expression_distribution(cm, normalized = TRUE)$box

  sep_crudo <- abs(mean(crudo$med[1:3]) - mean(crudo$med[4:6]))
  sep_norm  <- abs(mean(norm$med[1:3])  - mean(norm$med[4:6]))
  expect_lt(sep_norm, sep_crudo)
})

test_that("el grupo se asigna por sample_id y no por posición", {
  # La trampa clasica: el samplesheet llega en otro orden que las columnas de la
  # matriz. Emparejar por posición invierte los grupos en silencio y colorea el
  # gráfico al reves, que es peor que no colorearlo.
  cm <- make_test_counts(n_genes = 200, n_per_group = 3)
  meta <- make_test_meta(cm)
  meta_desordenada <- meta[rev(seq_len(nrow(meta))), , drop = FALSE]

  b <- expression_distribution(cm)$box
  g <- distribution_add_group(b, meta_desordenada, "condition")

  esperado <- meta$condition[match(g$sample_id, meta$sample_id)]
  expect_equal(g$grupo, esperado)
})

test_that("las muestras que no están en el samplesheet quedan sin grupo, no fuera", {
  cm <- make_test_counts(n_genes = 200, n_per_group = 3)
  meta <- make_test_meta(cm)[1:4, , drop = FALSE]

  g <- distribution_add_group(expression_distribution(cm)$box, meta, "condition")
  expect_equal(nrow(g), ncol(cm))
  expect_true(any(g$grupo == "(sin grupo)"))
})

test_that("el pie declara cuantos genes se han dejado fuera", {
  cm <- make_test_counts(n_genes = 200, n_per_group = 3)
  cm[1:15, ] <- 0L
  txt <- distribution_caption(expression_distribution(cm))

  expect_true(grepl("15", txt, fixed = TRUE))
  expect_true(grepl("TMM", txt, fixed = TRUE))
  expect_true(grepl("sin normalizar",
                    distribution_caption(expression_distribution(cm, normalized = FALSE)),
                    fixed = TRUE))
})

test_that("una matriz vacia o degenerada devuelve NULL en lugar de fallar", {
  expect_null(expression_distribution(NULL))
  expect_null(expression_distribution(matrix(0L, nrow = 10, ncol = 3,
                                             dimnames = list(paste0("g", 1:10),
                                                             paste0("s", 1:3)))))
  # Menos de 5 genes: no hay densidad que estimar.
  cm <- make_test_counts(n_genes = 200, n_per_group = 3)[1:3, , drop = FALSE]
  expect_null(expression_distribution(cm))
})

# ── Integración del modulo server ───────────────────────────────────────────
#
# Los reactivos internos del modulo no son accesibles desde testServer (viven en
# el marco de server_tab_deg_diag, no en el de la función server de prueba), así
# que se comprueba lo único que es realmente contrato: lo que sale por los
# outputs.

test_that("el modulo server produce los dos gráficos y respeta el selector de escala", {
  skip_if_not_installed("plotly")
  cm <- make_test_counts(n_genes = 300, n_per_group = 3)
  # Profundidades muy dispares entre grupos: garantiza que normalizar cambie el
  # resultado de forma observable.
  cm[, 4:6] <- cm[, 4:6] * 3L
  meta <- make_test_meta(cm)

  srv <- function(input, output, session) {
    state <- new.env(parent = emptyenv())
    state$deg_rv <- shiny::reactiveValues(counts = cm, meta = meta, results = NULL,
                                          disp_data = NULL, cooks = NULL)
    ctx <- new.env(parent = emptyenv())
    server_tab_deg_diag(input, output, session, state, ctx)
  }

  shiny::testServer(srv, {
    session$setInputs(deg_diag_expr_scale = "norm", deg_condition_col = "condition",
                      deg_diag_pv_subset = "tested", deg_source = "current")
    # Evaluar el output es la prueba: si el gráfico fallara, esto lanzaria.
    cajas_norm <- output$deg_expr_box
    expect_true(nzchar(output$deg_expr_density))
    expect_true(nzchar(cajas_norm))
    expect_true(grepl("TMM", as.character(output$deg_expr_dist_caption$html)))
    # Las muestras se identifican por su nombre y los grupos del samplesheet
    # llegan al gráfico.
    expect_true(grepl("ctrl1", cajas_norm, fixed = TRUE))
    expect_true(grepl("trt1", cajas_norm, fixed = TRUE))

    # El selector cambia lo que se CALCULA, no solo la etiqueta: la figura
    # serializada tiene que ser distinta.
    session$setInputs(deg_diag_expr_scale = "raw")
    cajas_raw <- output$deg_expr_box
    expect_true(nzchar(cajas_raw))
    expect_false(identical(cajas_norm, cajas_raw))
    expect_true(grepl("sin normalizar", as.character(output$deg_expr_dist_caption$html)))
  })
})

test_that("sin conteos cargados los gráficos explican que falta en lugar de fallar", {
  skip_if_not_installed("plotly")
  srv <- function(input, output, session) {
    state <- new.env(parent = emptyenv())
    state$deg_rv <- shiny::reactiveValues(counts = NULL, meta = NULL, results = NULL,
                                          disp_data = NULL, cooks = NULL)
    ctx <- new.env(parent = emptyenv())
    server_tab_deg_diag(input, output, session, state, ctx)
  }
  shiny::testServer(srv, {
    session$setInputs(deg_diag_expr_scale = "norm", deg_condition_col = "condition",
                      deg_diag_pv_subset = "tested", deg_source = "current")
    # El mensaje sustituye al gráfico; lo que no puede pasar es que el output
    # lance y se lleve por delante la pestana entera de diagnósticos.
    expect_true(nzchar(output$deg_expr_density))
    expect_true(nzchar(output$deg_expr_box))
    expect_null(output$deg_expr_dist_caption)
  })
})
