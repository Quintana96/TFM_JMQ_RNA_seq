#' test-privacidad.R
#' Los identificadores de muestra de un estudio clínico suelen llevar
#' información identificativa, y en esta app viajan a gráficos, tablas, informes
#' y ficheros persistidos: basta exportar una figura para difundirlos.

test_that("la seudonimizacion mantiene alineados matriz y samplesheet", {
  counts <- matrix(1:12, nrow = 3,
                   dimnames = list(paste0("g", 1:3),
                                   c("TCGA-E2-A15C", "TCGA-AR-A250", "SP4007", "SP3995")))
  meta <- data.frame(sample_id = colnames(counts),
                     condition = c("a", "a", "b", "b"), stringsAsFactors = FALSE)

  r <- pseudonymize_dataset(counts, meta)

  # Si se renombra la matriz sin el samplesheet, el alineamiento se rompe y el
  # análisis falla de formas difíciles de diagnosticar.
  expect_identical(colnames(r$counts), r$meta$sample_id)
  expect_equal(nrow(r$map), 4)
  expect_false(any(grepl("TCGA", colnames(r$counts))))
  # Los datos no se tocan, solo las etiquetas.
  expect_equal(unname(r$counts), unname(counts))
})

test_that("la correspondencia es reversible y estable", {
  ids <- c("m3", "m1", "m2", "m1")
  map <- build_pseudonym_map(ids)
  expect_equal(nrow(map), 3)                      # sin duplicar
  expect_identical(map$original[1], "m3")         # orden de aparicion
  alias <- apply_pseudonyms(ids, map)
  # El mismo identificador recibe siempre el mismo alias.
  expect_identical(alias[2], alias[4])
  # Deshacer la correspondencia devuelve los originales.
  vuelta <- map$original[match(alias, map$alias)]
  expect_identical(vuelta, ids)
})

test_that("los identificadores desconocidos se dejan intactos", {
  map <- build_pseudonym_map(c("a", "b"))
  # Convertirlos en NA en silencio sería peor que dejarlos visibles.
  expect_identical(apply_pseudonyms(c("a", "z"), map), c("S01", "z"))
})

test_that("se detectan columnas potencialmente identificativas", {
  meta <- data.frame(
    sample_id = paste0("s", 1:4),
    condition = c("a", "a", "b", "b"),
    nhc_paciente = c("1", "2", "3", "4"),
    centro = c("H1", "H1", "H2", "H2"),
    stringsAsFactors = FALSE)
  det <- detect_identifying_columns(meta)
  expect_true("nhc_paciente" %in% det)
  # condition y centro se repiten entre muestras: no reidentifican por si solas.
  expect_false("condition" %in% det)
  expect_false("centro" %in% det)
  # sample_id se excluye siempre: es la clave, no un dato de más.
  expect_false("sample_id" %in% det)
})

test_that("una columna única por muestra se marca aunque su nombre no lo delate", {
  meta <- data.frame(
    sample_id = paste0("s", 1:4),
    condition = c("a", "a", "b", "b"),
    comentario = c("caso raro Madrid", "control", "caso Bilbao", "control 2"),
    stringsAsFactors = FALSE)
  expect_true("comentario" %in% detect_identifying_columns(meta))
})

test_that("entradas vacias no rompen nada", {
  expect_equal(nrow(build_pseudonym_map(character(0))), 0)
  expect_equal(detect_identifying_columns(NULL), character(0))
  r <- pseudonymize_dataset(NULL, NULL)
  expect_equal(nrow(r$map), 0)
})
