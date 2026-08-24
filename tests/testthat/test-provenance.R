#' test-provenance.R
#' El pipeline debe dejar constancia de con que se ejecuto: estado de salida,
#' versiones de las herramientas y huella de las entradas. Sin eso es imposible
#' reconstruir que produjo un resultado, y cambiar la versión de una herramienta
#' basta para alterar la lista de genes diferenciales.

escribir_run <- function(root, status = "success", code = 0L) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  writeLines(c(paste0("exit_code\t", code),
               paste0("status\t", status),
               "finished_at\t2026-08-15 12:00:00"),
             file.path(root, "exit_status.tsv"))
  writeLines(c("tool\tversion\tpath",
               "fastqc\tFastQC v0.12.1\t/usr/local/bin/fastqc",
               "salmon\t(no instalado)\t—"),
             file.path(root, "versions.tsv"))
  writeLines(c("file\tsize_bytes\tmd5",
               "/datos/a_1.fastq.gz\t120\tb3e00341f3af5fdcb802500f7959bb12"),
             file.path(root, "checksums.tsv"))
  root
}

test_that("el estado sale del fichero explícito y no de adivinar el log", {
  root <- withr::local_tempdir()
  escribir_run(root, status = "error", code = 1L)
  # El log dice lo contrario a propósito: debe ganar el fichero de estado.
  writeLines("Analysis completed successfully", file.path(root, "workflow_live.log"))

  expect_identical(status_from_log(root), "error")
  st <- read_exit_status(root)
  expect_identical(st$exit_code, "1")
})

test_that("una ejecución correcta se reporta como completada", {
  root <- withr::local_tempdir()
  escribir_run(root, status = "success", code = 0L)
  expect_identical(status_from_log(root), "completado")
})

test_that("sin exit_status.tsv se recurre al log (ejecuciones antiguas)", {
  root <- withr::local_tempdir()
  writeLines("Analysis completed successfully", file.path(root, "workflow_live.log"))
  expect_identical(status_from_log(root), "completado")

  root2 <- withr::local_tempdir()
  writeLines("fallo en la línea 42", file.path(root2, "workflow_live.log"))
  expect_identical(status_from_log(root2), "error")

  expect_identical(status_from_log(withr::local_tempdir()), "sin log")
})

test_that("se leen versiones y checksums cuando existen", {
  root <- withr::local_tempdir()
  escribir_run(root)

  v <- read_tool_versions(root)
  expect_true(is.data.frame(v))
  expect_true(all(c("tool", "version", "path") %in% names(v)))
  expect_true("fastqc" %in% v$tool)

  ck <- read_input_checksums(root)
  expect_true(is.data.frame(ck))
  expect_true(all(c("file", "size_bytes", "md5") %in% names(ck)))
  expect_match(ck$md5[1], "^[0-9a-f]{32}$")
})

test_that("ausencia de ficheros de provenance no rompe nada", {
  root <- withr::local_tempdir()
  expect_null(read_exit_status(root))
  expect_null(read_tool_versions(root))
  expect_null(read_input_checksums(root))
})
