#' test-reactome.R
#' Reactome anade dos cosas que las otras colecciones no tenian y que son las que
#' pueden fallar en silencio: un catalogo CERRADO de organismos (no hay
#' procariotas) y una TRADUCCION de identificadores a ENTREZID. Un fallo en la
#' primera devuelve una tabla vacia que parece "sin enriquecimiento"; uno en la
#' segunda devuelve un resultado calculado sobre una fraccion de los genes sin
#' que se note.

test_that("el organismo de Reactome se deduce del OrgDb, y para procariotas no existe", {
  expect_equal(reactome_organism_for_orgdb("org.Hs.eg.db"), "human")
  expect_equal(reactome_organism_for_orgdb("org.Mm.eg.db"), "mouse")
  expect_equal(reactome_organism_for_orgdb("org.Sc.sgd.db"), "yeast")

  # E. coli es el organismo del dataset de ejemplo de la app: Reactome no lo
  # cubre, y decirlo es mas util que devolver una tabla vacia.
  expect_null(reactome_organism_for_orgdb("org.EcK12.eg.db"))
  expect_null(reactome_organism_for_orgdb(NULL))
  expect_null(reactome_organism_for_orgdb(""))
})

test_that("todo organismo deducido de un OrgDb esta en el catalogo de Reactome", {
  # Si las dos listas se separan, el selector ofreceria un organismo que
  # enrichPathway rechaza.
  orgdbs <- c("org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "org.Dr.eg.db",
              "org.Dm.eg.db", "org.Ce.eg.db", "org.Sc.sgd.db")
  for (o in orgdbs) {
    org <- reactome_organism_for_orgdb(o)
    expect_true(org %in% REACTOME_ORGANISMOS, info = o)
  }
})

test_that("un organismo fuera del catalogo se rechaza con un mensaje accionable", {
  res <- run_enrichment_reactome(c("TP53", "BRCA1"), organism = "ecoli",
                                 OrgDb = "org.EcK12.eg.db")
  expect_null(res$table)
  expect_true(nzchar(res$error))
  # El mensaje tiene que decir que hacer, no solo que ha fallado: sin
  # ReactomePA instalado dice cual falta; con el instalado, que organismos hay.
  expect_true(grepl("ReactomePA|Reactome no cubre", res$error))
})

test_that("sin genes o sin OrgDb la degradacion es explicita", {
  expect_true(nzchar(run_enrichment_reactome(character(0))$error))

  tr <- translate_to_entrez(c("TP53", "BRCA1"), OrgDb = NULL, keyType = "SYMBOL")
  expect_length(tr$ids, 0)
  expect_true(nzchar(tr$error))
  # La tasa de mapeo viaja igualmente: es lo que distingue "no mapea nada" de
  # "no se ha intentado".
  expect_equal(tr$mapping$n_input, 2L)
  expect_equal(tr$mapping$n_mapped, 0L)
})

test_that("con IDs ya en ENTREZID la traduccion es la identidad y no pierde genes", {
  ids <- c("7157", "672", "324")
  tr <- translate_to_entrez(ids, OrgDb = NULL, keyType = "ENTREZID")

  expect_setequal(tr$ids, ids)
  expect_null(tr$error)
  expect_equal(tr$mapping$rate, 1)
  # `back` permite deshacer la traduccion al mostrar los resultados.
  expect_equal(unname(tr$back[ids]), ids)
})

test_that("el ranking traducido conserva el orden por valor y no duplica genes", {
  rk <- c("7157" = 5.2, "672" = -3.1, "324" = 1.4, "207" = 0.2)
  out <- translate_ranking_to_entrez(rk, OrgDb = NULL, keyType = "ENTREZID")

  expect_false(is.null(out$ranked))
  expect_equal(names(out$ranked), names(sort(rk, decreasing = TRUE)))
  expect_false(any(duplicated(names(out$ranked))))
  expect_null(out$error)
})

test_that("un ranking vacio no rompe el GSEA de Reactome", {
  out <- translate_ranking_to_entrez(NULL, OrgDb = NULL, keyType = "ENTREZID")
  expect_null(out$ranked)
  expect_true(nzchar(out$error))

  res <- run_gsea(NULL, ont = "REACTOME", organism = "human")
  expect_null(res$table)
  expect_true(nzchar(res$error))
})

test_that("el informe recibe una etiqueta legible de la coleccion, no el codigo", {
  expect_equal(enrich_collection_label("REACTOME"), "Reactome")
  expect_equal(enrich_collection_label("BP"), "GO: procesos biologicos")
  expect_equal(enrich_collection_label("GMT"), "Gene sets propios (GMT)")
  expect_equal(enrich_collection_label("KEGG"), "KEGG")
  # Un codigo desconocido se muestra tal cual antes que romper el informe.
  expect_equal(enrich_collection_label("XYZ"), "XYZ")
})

test_that("los genes de una ruta se piden a reactome.db, no a la red", {
  # Con reactome.db instalado devuelve ENTREZIDs; sin el, un error que nombra el
  # paquete que falta. En ninguno de los dos casos se consulta una API en linea,
  # que es la diferencia operativa con KEGG.
  tg <- gsea_term_genes("R-HSA-68886", ont = "REACTOME")   # M Phase
  if (requireNamespace("reactome.db", quietly = TRUE)) {
    # Exigir genes de verdad, no "genes o un error": aceptar las dos cosas
    # convierte el test en uno que no puede fallar, y de hecho no fallo cuando
    # la consulta devolvia vacio por un argumento mal pasado. La consecuencia
    # visible era un running score plano sobre un conjunto que si existe.
    expect_null(tg$error)
    expect_gt(length(tg$genes), 50)
    # Vienen en ENTREZID, que es el espacio en el que corre el GSEA de Reactome.
    expect_true(all(grepl("^[0-9]+$", head(tg$genes, 20))))
  } else {
    expect_true(grepl("reactome.db", tg$error, fixed = TRUE))
  }
  # Una ruta inexistente se declara como tal, no devuelve un conjunto vacio mudo.
  malo <- gsea_term_genes("R-HSA-000000", ont = "REACTOME")
  expect_length(malo$genes, 0)
  expect_true(nzchar(malo$error))
})
