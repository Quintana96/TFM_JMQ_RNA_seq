# Traducción de identificadores con la anotación (R/utils_annotation.R).
#
# El pipeline llama a `featureCounts -g locus_tag`, de modo que la matriz de
# conteos trae locus tags y ningun OrgDb los conoce. La anotación es la única
# fuente que relaciona BW25113_RS00005 con thrL. Sin traducir, el
# enriquecimiento mapea el 0 % y devuelve "sin términos", que no se distingue
# de la ausencia de señal.

# Un GTF mínimo con los mismos atributos que el de RefSeq, para no depender de
# ningun fichero externo.
gtf_de_prueba <- function() {
  filas <- c(
    # Comentario al estilo de los GTF de RefSeq. Un '##gff-versión 3' haría que
    # rtracklayer avisara de que la versión no cuadra con la extensión .gtf.
    '#!genome-build ASM75055v1',
    'chr\tRefSeq\tgene\t1\t100\t.\t+\t.\tgene_id "LT_001"; locus_tag "LT_001"; old_locus_tag "OLD_1"; gene "aaaA";',
    'chr\tRefSeq\tCDS\t1\t100\t.\t+\t0\tgene_id "LT_001"; locus_tag "LT_001"; gene "aaaA"; protein_id "WP_1";',
    'chr\tRefSeq\tgene\t200\t300\t.\t+\t.\tgene_id "LT_002"; locus_tag "LT_002"; old_locus_tag "OLD_2"; gene "bbbB";',
    'chr\tRefSeq\tCDS\t200\t300\t.\t+\t0\tgene_id "LT_002"; locus_tag "LT_002"; gene "bbbB"; protein_id "WP_2";',
    # Tercer gen que comparte simbolo con el segundo: es el caso de los
    # parologos y los fragmentos anotados por separado.
    'chr\tRefSeq\tgene\t400\t500\t.\t+\t.\tgene_id "LT_003"; locus_tag "LT_003"; gene "bbbB";',
    # Cuarto gen SIN simbolo: no todos lo tienen, y por eso featureCounts no
    # puede agrupar por `gene` (aborta si el atributo falta en algun registro).
    'chr\tRefSeq\tgene\t600\t700\t.\t+\t.\tgene_id "LT_004"; locus_tag "LT_004";'
  )
  f <- tempfile(fileext = ".gtf")
  writeLines(filas, f)
  f
}

test_that("se deduce el atributo al que corresponden unos identificadores", {
  skip_if_not(HAS_RTRACKLAYER, "falta rtracklayer")
  gtf <- gtf_de_prueba()

  d <- detect_annotation_keytype(c("LT_001", "LT_002", "LT_003"), gtf)
  expect_equal(d$rate, 1)
  expect_true(d$attr %in% c("gene_id", "locus_tag"))

  # Con simbolos, el atributo deducido tiene que ser otro.
  d2 <- detect_annotation_keytype(c("aaaA", "bbbB"), gtf)
  expect_equal(d2$attr, "gene")

  # Identificadores que no están en ninguna columna: cobertura cero, no error.
  d3 <- detect_annotation_keytype(c("NADA_1", "NADA_2"), gtf)
  expect_equal(d3$rate, 0)
})

test_that("la traducción respeta los genes sin simbolo y colapsa los repetidos", {
  skip_if_not(HAS_RTRACKLAYER, "falta rtracklayer")
  gtf <- gtf_de_prueba()

  tr <- translate_ids_with_annotation(c("LT_001", "LT_002", "LT_003", "LT_004"),
                                      gtf, to = "gene")
  # LT_004 no tiene simbolo y LT_003 comparte el de LT_002.
  expect_setequal(tr$ids, c("aaaA", "bbbB"))
  expect_equal(tr$mapping$n_input, 4L)
  expect_equal(tr$mapping$n_mapped, 2L)
  expect_equal(tr$mapping$n_colapsados, 1L)

  # `back` tiene que permitir volver a los identificadores de la matriz.
  expect_equal(unname(tr$back[["aaaA"]]), "LT_001")
})

test_that("traducir a lo mismo es la identidad y no pierde genes", {
  skip_if_not(HAS_RTRACKLAYER, "falta rtracklayer")
  gtf <- gtf_de_prueba()
  ids <- c("LT_001", "LT_004")
  tr <- translate_ids_with_annotation(ids, gtf, from = "locus_tag", to = "locus_tag")
  expect_equal(tr$ids, ids)
  expect_equal(tr$mapping$rate, 1)
  # Sin este caso, un usuario que pide traducir a lo que ya tiene perdería
  # LT_004, que no tiene simbolo, sin ninguna razón.
})

test_that("los errores explican la causa en vez de devolver una lista vacia", {
  skip_if_not(HAS_RTRACKLAYER, "falta rtracklayer")
  gtf <- gtf_de_prueba()

  expect_match(translate_ids_with_annotation(character(0), gtf)$error, "vacia")
  expect_match(translate_ids_with_annotation("LT_001", gtf, to = "noexiste")$error,
               "no trae el atributo")
  expect_match(translate_ids_with_annotation(c("X1", "X2"), gtf, to = "gene")$error,
               "Ningun atributo")
  expect_match(translate_ids_with_annotation("LT_001", "/no/existe.gtf")$error,
               "No se pudo leer")
})

test_that("el ranking conserva el valor más extremo al colapsar y queda ordenado", {
  skip_if_not(HAS_RTRACKLAYER, "falta rtracklayer")
  gtf <- gtf_de_prueba()

  # LT_002 y LT_003 caen ambos en bbbB. Debe ganar LT_003 por |valor|.
  rk <- c(LT_001 = 1.5, LT_002 = -0.5, LT_003 = -4.0, LT_004 = 9.9)
  r <- translate_ranking_with_annotation(rk, gtf, to = "gene")

  expect_setequal(names(r$ranked), c("aaaA", "bbbB"))
  expect_equal(unname(r$ranked[["bbbB"]]), -4.0)
  expect_equal(unname(r$ranked[["aaaA"]]), 1.5)
  # LT_004 no tiene simbolo: se pierde aunque sea el valor mayor.
  expect_false("LT_004" %in% names(r$ranked))
  # GSEA recorre el ranking en orden: si la traducción lo desordena, el NES no
  # corresponde a la curva que dibuja el running score.
  expect_false(is.unsorted(rev(r$ranked)))
  expect_equal(r$mapping$n_colapsados, 1L)
})

test_that("solo se ofrecen como destino los atributos que la anotación trae", {
  skip_if_not(HAS_RTRACKLAYER, "falta rtracklayer")
  attrs <- annotation_available_attrs(gtf_de_prueba(),
                                      names(ANNOTATION_TARGET_ATTRS))
  expect_true(all(c("gene", "locus_tag", "old_locus_tag", "protein_id") %in% attrs))
  expect_equal(annotation_available_attrs("/no/existe.gtf"), character(0))
})

test_that("el informe registra la traducción, y también su ausencia", {
  rv <- list(
    results = data.frame(gene = c("a", "b"), log2FC = c(1, -1), padj = c(0.01, 0.02)),
    method = "DESeq2", fdr = 0.05, contrast = "Galio vs Control",
    enrich = list(
      enfoque = "ORA (sobre-representación)", ontologia = "GO: Procesos biológicos",
      keytype = "SYMBOL", n_lista = 1502L, n_universo = 3289L, n_terminos = 15L,
      traduccion = "gene_id -> gene (1502 de 1535, 97.9 %; 7 colapsados)",
      anotacion_traduccion = "/ruta/anotacion.gtf",
      mapeo = list(rate = 0.963, n_mapped = 1446L, n_input = 1502L),
      metrica = NA, simplify = FALSE, direccional = FALSE, simbolos = TRUE))

  h <- build_deg_report_html(rv)
  expect_true(grepl("Traducción de IDs", h, fixed = TRUE))
  expect_true(grepl("gene_id -&gt; gene", h, fixed = TRUE))
  expect_true(grepl("anotacion.gtf", h, fixed = TRUE))

  # Y sin traducir tiene que DECIRLO, no callar: un informe que no menciona la
  # traducción no permite saber si los identificadores del enriquecimiento son
  # los de la matriz o unos traducidos, que es justo lo que hay que reproducir.
  rv$enrich$traduccion <- NA_character_
  rv$enrich$anotacion_traduccion <- NA_character_
  h2 <- build_deg_report_html(rv)
  expect_true(grepl("se usaron los IDs de la matriz", h2, fixed = TRUE))
})
