#' test-matriz-conteos.R
#' La matriz de conteos es el entregable principal que queda en disco tras una
#' ejecucion. Para salmon y kallisto se escribia con SOLO la cabecera, asi que
#' quien se llevaba la carpeta obtenia un fichero corrupto presentado como el
#' resultado del pipeline.
#'
#' El test es autocontenido: genera una anotacion y unas cuantificaciones
#' sinteticas, porque data/ esta en .gitignore y no puede asumirse presente.

#' Escribe un GFF minimo con la misma convencion que la anotacion de E. coli:
#' CDS con ID = "cds-<accession>.<version>" y Parent = "gene-<locus_tag>".
#' Esa convencion es justo la que rompia a tximport con ignoreTxVersion = TRUE,
#' porque el corte por el primer punto confunde una version de Ensembl con un
#' accession de GenBank.
write_mini_gff <- function(path, n_genes = 40) {
  lt  <- sprintf("b%04d", seq_len(n_genes))
  acc <- sprintf("AAC%05d.1", seq_len(n_genes))
  gene_rows <- sprintf(
    "chr\tRefSeq\tgene\t%d\t%d\t.\t+\t.\tID=gene-%s;locus_tag=%s",
    seq(1, by = 1000, length.out = n_genes),
    seq(900, by = 1000, length.out = n_genes), lt, lt)
  cds_rows <- sprintf(
    "chr\tRefSeq\tCDS\t%d\t%d\t.\t+\t0\tID=cds-%s;Parent=gene-%s;locus_tag=%s;protein_id=%s",
    seq(1, by = 1000, length.out = n_genes),
    seq(900, by = 1000, length.out = n_genes), acc, lt, lt, acc)
  writeLines(c("##gff-version 3", gene_rows, cds_rows), path)
  list(locus_tags = lt, tx_ids = paste0("cds-", acc))
}

#' Estructura de salida de salmon: 03_alignments/salmon/<muestra>/quant.sf
write_fake_salmon <- function(root, tx_ids, samples = c("m1", "m2", "m3")) {
  withr::with_seed(11, {
    for (s in samples) {
      d <- file.path(root, "03_alignments", "salmon", s)
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
      df <- data.frame(
        Name = tx_ids,
        Length = sample(300:2000, length(tx_ids), replace = TRUE))
      df$EffectiveLength <- df$Length - 20
      df$TPM <- round(stats::runif(nrow(df), 0, 900), 3)
      df$NumReads <- round(stats::runif(nrow(df), 0, 4000), 3)
      utils::write.table(df, file.path(d, "quant.sf"), sep = "\t",
                         quote = FALSE, row.names = FALSE)
    }
  })
  samples
}

test_that("las cuantificaciones de salmon se resumen a una matriz por gen real", {
  skip_if_not(requireNamespace("tximport", quietly = TRUE), "tximport no instalado")
  root <- withr::local_tempdir()
  gff <- file.path(root, "anot.gff")
  ann <- write_mini_gff(gff)
  samples <- write_fake_salmon(root, ann$tx_ids)

  m <- suppressMessages(
    load_counts_from_workflow(root, "salmon", annotation_file = gff))

  expect_false(is.null(m))
  # Una matriz REAL: una fila por gen y una columna por muestra, no una cabecera.
  expect_equal(ncol(m), length(samples))
  expect_gt(nrow(m), 1)
  expect_setequal(colnames(m), samples)
  # El identificador canonico es el locus_tag, igual que featureCounts -g
  # locus_tag, para que alineamiento y pseudoalineamiento sean comparables.
  expect_true(all(rownames(m) %in% ann$locus_tags))

  # Y debe haber ido por la ruta buena, no por el respaldo de est_counts.
  src <- attr(m, "counts_source")
  expect_true(isTRUE(src$ok))
  expect_match(src$method, "tximport", fixed = TRUE)
})

test_that("sin anotacion utilizable la degradacion es explicita, no silenciosa", {
  skip_if_not(requireNamespace("tximport", quietly = TRUE), "tximport no instalado")
  root <- withr::local_tempdir()
  ann <- write_mini_gff(file.path(root, "anot.gff"))
  write_fake_salmon(root, ann$tx_ids)

  # Anotacion que no corresponde al transcriptoma cuantificado.
  otro <- file.path(root, "otra.gff")
  writeLines(c("##gff-version 3",
               "chr\tRefSeq\tCDS\t1\t100\t.\t+\t0\tID=cds-XXX.1;Parent=gene-zzz;locus_tag=zzz"),
             otro)

  m <- suppressMessages(load_counts_from_workflow(root, "salmon", annotation_file = otro))
  src <- attr(m, "counts_source")

  # Se devuelve algo utilizable, pero marcado como degradado y con el motivo.
  expect_false(isTRUE(src$ok))
  expect_true(nzchar(src$detail))
})

test_that("el script de linea de comandos escribe la matriz y falla limpiamente", {
  skip_if_not(requireNamespace("tximport", quietly = TRUE), "tximport no instalado")
  root <- withr::local_tempdir()
  gff <- file.path(root, "anot.gff")
  ann <- write_mini_gff(gff)
  write_fake_salmon(root, ann$tx_ids)
  dest <- file.path(root, "count_matrix.tsv")

  script <- file.path(app_root, "scripts", "build_count_matrix.R")
  skip_if_not(file.exists(script), "scripts/build_count_matrix.R no encontrado")

  st <- system2("Rscript", c(shQuote(script), shQuote(root), "salmon",
                             shQuote(dest), shQuote(gff)),
                stdout = FALSE, stderr = FALSE)
  expect_equal(st, 0L)
  expect_true(file.exists(dest))

  tab <- utils::read.delim(dest, check.names = FALSE)
  expect_gt(nrow(tab), 1)                 # no solo la cabecera
  expect_true("gene_id" %in% names(tab))

  # Herramienta no soportada: codigo de salida de error de uso, sin escribir nada.
  dest2 <- file.path(root, "no.tsv")
  st2 <- system2("Rscript", c(shQuote(script), shQuote(root), "bowtie2",
                              shQuote(dest2), shQuote(gff)),
                 stdout = FALSE, stderr = FALSE)
  expect_equal(st2, 1L)
  expect_false(file.exists(dest2))
})
