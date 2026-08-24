#' test-informe.R
#' El informe es el documento que queda cuando el análisis ya no está en
#' pantalla. Se estructura según el contenido mínimo que la guia ACMG/AMP exige
#' a un informe de diagnóstico genomico: que se evaluo, con que bases de datos y
#' pipeline (con sus VERSIONES), los hallazgos, la fecha, un resumen en lenguaje
#' llano y las limitaciones. Un informe negativo se trata igual que uno positivo.

rv_minimo <- function(sig = TRUE, n_per_group = 4) {
  counts <- make_test_counts(n_genes = 200, n_per_group = n_per_group)
  meta <- make_test_meta(counts)
  res <- run_deg(counts, meta, method = "DESeq2", ref_level = "ctrl",
                 contrast_num = "trt", fdr = 0.05, shrink = FALSE)
  tab <- res$table
  if (!sig) tab$padj <- 1
  list(results = tab, meta = meta, method = "DESeq2", design = res$design,
       coef = res$coef, contrast = res$contrast, fdr = 0.05, lfc_threshold = 0,
       shrink = "ninguno", padj_method = "BH", prefilter = NULL,
       run_at = Sys.time(), seeds = list())
}

test_that("el informe incluye resumen en lenguaje llano y fecha", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  h <- suppressMessages(build_deg_report_html(rv_minimo()))
  expect_true(grepl("Resumen", h, fixed = TRUE))
  expect_true(grepl("Análisis realizado el", h, fixed = TRUE))
  expect_true(grepl("Limitaciones", h, fixed = TRUE))
  expect_true(grepl("Reproducibilidad", h, fixed = TRUE))
})

test_that("el informe fija fondo claro para no depender del tema del lector", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  h <- suppressMessages(build_deg_report_html(rv_minimo()))
  # Sin fondo explícito, un navegador en modo oscuro deja texto oscuro sobre
  # negro: el documento se descarga y se lee en cualquier parte.
  expect_true(grepl("background:#FFFFFF", h, fixed = TRUE))
})

test_that("un resultado negativo se documenta en lugar de quedar en blanco", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  h <- suppressMessages(build_deg_report_html(rv_minimo(sig = FALSE)))
  expect_true(grepl("no demuestra ausencia de efecto", h, fixed = TRUE))
  # Se listan igualmente los genes de menor p-valor: hay que documentar que se
  # evaluo, no solo que no salio nada.
  expect_true(grepl("Genes con menor p-valor ajustado", h, fixed = TRUE))
})

test_that("las limitaciones se derivan del análisis, no son texto fijo", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  rv <- rv_minimo(n_per_group = 4)
  lims <- deg_report_limitations(rv, diagnostics = NULL, sig_n = 10L)
  # 4 replicas < 6 recomendadas -> debe avisar
  expect_true(any(grepl("replicas", lims)))
  # Sin bootstrap -> debe avisar
  expect_true(any(grepl("replicabilidad por remuestreo", lims)))

  # Con un resumen a gen degradado, aparece esa limitación y no antes
  expect_false(any(grepl("via recomendada", lims)))
  rv$counts_source <- list(method = "est_counts", ok = FALSE, detail = "Sin anotación.")
  expect_true(any(grepl("via recomendada", deg_report_limitations(rv, NULL, 10L))))
})

test_that("la procedencia y las semillas aparecen cuando existen", {
  skip_if_not(requireNamespace("DESeq2", quietly = TRUE), "DESeq2 no instalado")
  rv <- rv_minimo()
  rv$counts_origin <- list(tipo = "Matriz subida", ruta = "x.tsv",
                           md5 = "d41d8cd98f00b204e9800998ecf8427e")
  rv$seeds <- list(sva = 1L, n_sv = 2L)
  h <- suppressMessages(build_deg_report_html(rv))
  expect_true(grepl("Procedencia de los datos", h, fixed = TRUE))
  expect_true(grepl("d41d8cd98f00b204e9800998ecf8427e", h, fixed = TRUE))
  expect_true(grepl("Semilla de sva", h, fixed = TRUE))
})

test_that("orgdb_source_info devuelve las fechas de las fuentes", {
  skip_if_not(requireNamespace("org.Hs.eg.db", quietly = TRUE), "org.Hs.eg.db no instalado")
  info <- orgdb_source_info("org.Hs.eg.db")
  expect_true("GOSOURCEDATE" %in% names(info))
  # El campo KEGG de los OrgDb está congelado: por eso la app no lo usa y
  # consulta la API. Que aparezca en el informe lo deja documentado.
  expect_true("KEGGSOURCEDATE" %in% names(info))
  expect_identical(orgdb_source_info(NULL)$OrgDb, "no disponible")
})

test_that("app_provenance identifica versión y commit", {
  p <- app_provenance()
  expect_true(all(c("Aplicación", "Commit de git", "Versión de R") %in% names(p)))
  expect_identical(p[["Versión de R"]], as.character(getRversion()))
})
