#' test-ejecución-guardada.R
#' Una ejecución guardada debe poder cargarse aunque ya no conserve la carpeta
#' 03_alignments. Los BAM y los índices son lo más pesado de una ejecución y es
#' habitual borrarlos para ahorrar espacio; la matriz de conteos, que es lo que
#' el análisis diferencial necesita, sobrevive.
#'
#' Antes la herramienta se deducia UNICAMENTE de la existencia de
#' 03_alignments/<tool>/, así que sin esa carpeta la ejecución quedaba como
#' "desconocida" y la matriz no se cargaba, con el mensaje genérico de que no
#' había matriz de conteos.

crear_run <- function(root, con_run_params = TRUE, con_alignments = FALSE,
                      tool = "bowtie2") {
  dir.create(file.path(root, "04_counts"), recursive = TRUE, showWarnings = FALSE)
  m <- data.frame(gene_id = sprintf("g%03d", 1:50))
  for (s in c("s1", "s2", "s3", "s4")) m[[s]] <- sample(10:900, 50, replace = TRUE)
  utils::write.table(m, file.path(root, "04_counts", "count_matrix.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  if (con_run_params) {
    writeLines(c(paste0("tool\t", tool)), file.path(root, "run_params.tsv"))
  }
  if (con_alignments) {
    dir.create(file.path(root, "03_alignments", tool), recursive = TRUE,
               showWarnings = FALSE)
  }
  root
}

test_that("una ejecución con 03_alignments se infiere como siempre", {
  root <- withr::local_tempdir()
  crear_run(root, con_run_params = FALSE, con_alignments = TRUE, tool = "salmon")
  p <- infer_result_params(root, "workflow.sh")
  expect_identical(p$tool, "salmon")
})

test_that("sin 03_alignments se respeta la herramienta de run_params.tsv", {
  root <- withr::local_tempdir()
  crear_run(root, con_run_params = TRUE, con_alignments = FALSE, tool = "bowtie2")
  p <- infer_result_params(root, "workflow.sh")
  expect_identical(p$tool, "bowtie2")

  cm <- load_counts_from_workflow(root, p$tool, annotation_file = NULL)
  expect_false(is.null(cm))
  expect_equal(ncol(cm), 4)
  expect_equal(nrow(cm), 50)
})

test_that("con solo la matriz de conteos la ejecución sigue siendo utilizable", {
  root <- withr::local_tempdir()
  crear_run(root, con_run_params = FALSE, con_alignments = FALSE)
  p <- infer_result_params(root, "workflow.sh")
  # No se puede saber como se genero, pero la matriz por gen se lee igual.
  expect_false(identical(p$tool, "desconocida"))
  cm <- load_counts_from_workflow(root, p$tool, annotation_file = NULL)
  expect_false(is.null(cm))
  expect_equal(ncol(cm), 4)
})

test_that("un directorio sin nada utilizable sigue siendo desconocido", {
  root <- withr::local_tempdir()
  p <- infer_result_params(root, "workflow.sh")
  expect_identical(p$tool, "desconocida")
})
