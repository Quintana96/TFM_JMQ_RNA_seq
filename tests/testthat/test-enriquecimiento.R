#' test-enriquecimiento.R
#' La entrada del enriquecimiento —lista de genes y universo— determina el
#' resultado tanto como el propio test. Wijesooriya et al. (2022) documentan que
#' el fondo mal definido es el error mas extendido en analisis de
#' sobre-representacion.

tabla_deg <- function() {
  data.frame(
    gene     = sprintf("g%02d", 1:10),
    baseMean = c(500, 400, 5, 300, 8, 900, 700, 2, 600, 100),
    log2FC   = c(3, 0.2, 4, -2.5, 0.1, 0.4, -3, 5, 0.3, 1.5),
    padj     = c(0.001, 0.01, 0.02, 0.001, NA, 0.03, 0.004, NA, 0.02, 0.6),
    pvalue   = c(1e-5, 1e-3, 1e-3, 1e-5, NA, 2e-3, 1e-4, 0.4, 1e-3, 0.5),
    stringsAsFactors = FALSE
  )
}

test_that("la lista del ORA no depende de los filtros de visualizacion", {
  tab <- tabla_deg()

  # Significativos a FDR 0,05: todos los que tienen padj <= 0,05.
  sig <- deg_significant_genes(tab, fdr = 0.05)
  esperado <- tab$gene[!is.na(tab$padj) & tab$padj <= 0.05]
  expect_setequal(sig$gene, esperado)

  # Los filtros visuales recortan la tabla que se DIBUJA...
  visual <- apply_deg_filters(tab, fdr = 0.05, abs_log2fc = 2, base_mean = 100)
  expect_lt(nrow(visual), nrow(sig))
  # ...pero la lista del enriquecimiento debe seguir siendo la misma.
  expect_setequal(deg_significant_genes(tab, fdr = 0.05)$gene, esperado)
})

test_that("el universo excluye los genes sin padj", {
  tab <- tabla_deg()
  u <- deg_testable_universe(tab)

  # g05 y g08 no tienen padj: descartados por filtrado independiente, conteo
  # cero u outlier de Cook. Nunca podrian haber entrado en la lista.
  expect_false("g05" %in% u)
  expect_false("g08" %in% u)
  expect_equal(length(u), sum(!is.na(tab$padj)))

  # Y el universo debe contener a todos los significativos: una lista con genes
  # fuera del fondo rompe el test hipergeometrico.
  sig <- deg_significant_genes(tab, fdr = 0.05)$gene
  expect_true(all(sig %in% u))
})

test_that("el universo cae al fondo completo si ningun gen tiene padj", {
  tab <- tabla_deg()
  tab$padj <- NA_real_
  # Un motor que no rellene padj dejaria el universo vacio y romperia el
  # enriquecimiento entero; el respaldo es preferible a no poder ejecutarlo.
  expect_equal(length(deg_testable_universe(tab)), nrow(tab))
})

test_that("tabla vacia o nula no rompe el universo", {
  expect_equal(deg_testable_universe(NULL), character(0))
  expect_equal(deg_testable_universe(tabla_deg()[0, ]), character(0))
})
