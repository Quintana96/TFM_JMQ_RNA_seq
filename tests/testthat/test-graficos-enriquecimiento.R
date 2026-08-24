# Gráficos de red y distribución del enriquecimiento (R/utils_enrich_plots.R).
#
# Se apoyan en la API de enrichplot, que ha cambiado entre versiones de forma no
# retrocompatible, así que estos tests valen sobre todo como detección temprana:
# si una actualización mueve un argumento, fallan aquí y no en una captura de la
# memoria.

test_that("el catalogo distingue ORA de GSEA", {
  ora  <- structure(list(), class = "enrichResult")
  gsea <- structure(list(), class = "gseaResult")

  expect_equal(enrich_obj_enfoque(ora), "ora")
  expect_equal(enrich_obj_enfoque(gsea), "gsea")
  expect_true(is.na(enrich_obj_enfoque(NULL)))

  # El ridge necesita el NES que solo trae GSEA, y las barras cuentan genes de
  # una lista umbralizada que solo tiene el ORA. Ofrecerlos donde no aplican
  # produce un error en tiempo de dibujo, que es justo lo que se evita.
  expect_true("ridge" %in% enrich_plot_choices(gsea))
  expect_false("ridge" %in% enrich_plot_choices(ora))
  expect_true("barra" %in% enrich_plot_choices(ora))
  expect_false("barra" %in% enrich_plot_choices(gsea))

  # Las tres vistas comunes tienen que estar en los dos enfoques.
  for (k in c("cnet", "emap", "upset")) {
    expect_true(k %in% enrich_plot_choices(ora), info = k)
    expect_true(k %in% enrich_plot_choices(gsea), info = k)
  }
})

test_that("cada tipo tiene una ayuda escrita", {
  for (k in names(ENRICH_PLOT_TIPOS)) {
    expect_true(nzchar(enrich_plot_ayuda(k)), info = k)
  }
  expect_equal(enrich_plot_ayuda("inexistente"), "")
})

test_that("se rechaza dibujar sin objeto o con un tipo que no aplica", {
  expect_error(enrich_make_network_plot(NULL, "emap"), "No hay resultado")
  expect_error(enrich_make_network_plot(structure(list(), class = "enrichResult"),
                                        "quesoyyo"), "Tipo desconocido")
})

test_that("pathview avisa en vez de romperse cuando no hay nada que pintar", {
  r <- enrich_pathview_png("eco00020", NULL)
  expect_null(r$path)
  expect_match(r$error, "log2FC")

  r2 <- enrich_pathview_png("", c(a = 1))
  expect_null(r2$path)
  expect_match(r2$error, "Selecciona")
})

# Los siguientes necesitan un enriquecimiento de verdad: se construye uno
# mínimo con datos sinteticos para no depender ni de la red ni de una ejecución
# previa del pipeline.
test_that("los cinco gráficos se dibujan sobre un enriquecimiento real", {
  skip_if_not(HAS_CLUSTERPROFILER && HAS_ORGECDB && HAS_ENRICHPLOT,
              "faltan clusterProfiler, org.EcK12.eg.db o enrichplot")

  orgdb <- getFromNamespace("org.EcK12.eg.db", "org.EcK12.eg.db")
  universo <- AnnotationDbi::keys(orgdb, keytype = "SYMBOL")
  skip_if(length(universo) < 500, "OrgDb demasiado pequeño")

  # Una lista enriquecida de verdad: los genes de un término GO concreto más
  # ruido. Sin señal, enrichGO devuelve una tabla vacia y no habría que dibujar.
  set.seed(42)
  tca <- tryCatch(
    AnnotationDbi::select(orgdb, keys = "GO:0006099", keytype = "GOALL",
                          columns = "SYMBOL")$SYMBOL,
    error = function(e) character(0)
  )
  skip_if(length(tca) < 10, "no se pudieron recuperar los genes del término")
  lista <- unique(c(tca, sample(universo, 150)))

  ego <- suppressWarnings(clusterProfiler::enrichGO(
    lista, universe = universo, OrgDb = orgdb, keyType = "SYMBOL",
    ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.5))
  skip_if(is.null(ego) || !nrow(as.data.frame(ego)), "sin términos enriquecidos")

  fc <- stats::setNames(stats::rnorm(length(lista)), lista)
  for (k in enrich_plot_choices(ego)) {
    p <- enrich_make_network_plot(ego, k, top_n = 8, fold_change = fc)
    expect_s3_class(p, "ggplot")
  }
})

test_that("showCategory llega a enrichplot como doble, no como entero", {
  # ridgeplot() de enrichplot 1.30 da por inválido un showCategory de tipo
  # "integer": avisa de que "should be a number of pathways" y muere después con
  # "objeto 'selected' no encontrado". El error no menciona el tipo, así que sin
  # este test la causa se vuelve a perder. Y el camino por el que se colaba era
  # nrow(), que devuelve entero.
  skip_if_not(HAS_ENRICHPLOT, "falta enrichplot")

  # Un doble de prueba con su propio método de as.data.frame: mockear el de base
  # provoca recursion infinita, porque lo llama medio R.
  falso <- structure(list(), class = c("gseaResultDePrueba", "gseaResult"))
  registerS3method("as.data.frame", "gseaResultDePrueba",
                   function(x, ...) data.frame(ID = letters[1:20]),
                   envir = environment())

  capturado <- NULL
  local_mocked_bindings(
    ridgeplot = function(x, showCategory, ...) {
      capturado <<- showCategory
      ggplot2::ggplot()
    },
    .package = "enrichplot"
  )

  enrich_make_network_plot(falso, "ridge", top_n = 15L)
  expect_false(is.integer(capturado))
  expect_equal(capturado, 15)

  # Y el recorte al número de términos disponibles tampoco puede reintroducirlo.
  enrich_make_network_plot(falso, "ridge", top_n = 500L)
  expect_false(is.integer(capturado))
  expect_equal(capturado, 20)
})

# ── Traducción de identificadores para KEGG ─────────────────────────────────

test_that("KEGG traduce a ENTREZID cuando se le pide ncbi-geneid", {
  skip_if_not(HAS_CLUSTERPROFILER && HAS_ORGECDB, "faltan clusterProfiler u org.EcK12.eg.db")

  # No se llama a enrichKEGG de verdad —necesita red—: se comprueba que lo que
  # le llega ya está traducido. Es el punto exacto donde fallaba: con simbolos,
  # KEGG mapeaba el 0 % y la app decía "sin términos enriquecidos", que no
  # distingue "no hay señal" de "los IDs no eran los correctos".
  recibido <- NULL
  local_mocked_bindings(
    enrichKEGG = function(gene, universe, organism, keyType, ...) {
      recibido <<- list(gene = gene, universe = universe)
      NULL
    },
    .package = "clusterProfiler"
  )
  genes <- c("thrA", "thrB", "thrC", "dnaK")
  run_enrichment_kegg(genes, universe = c(genes, "lacZ", "recA"),
                      organism = "eco", keyType = "ncbi-geneid",
                      OrgDb = "org.EcK12.eg.db", from_keytype = "SYMBOL")

  expect_false(is.null(recibido))
  # Los ENTREZID de E. coli son numéricos; los simbolos no.
  expect_true(all(grepl("^[0-9]+$", recibido$gene)))
  expect_true(all(grepl("^[0-9]+$", recibido$universe)))
  # El fondo tiene que traducirse también: una lista en un espacio y un universo
  # en otro dan un enriquecimiento sin sentido, no un error.
  expect_gt(length(recibido$universe), length(recibido$gene))
})

test_that("sin from_keytype, KEGG pasa los genes tal cual", {
  skip_if_not(HAS_CLUSTERPROFILER, "falta clusterProfiler")
  recibido <- NULL
  local_mocked_bindings(
    enrichKEGG = function(gene, ...) { recibido <<- gene; NULL },
    .package = "clusterProfiler"
  )
  # keyType "kegg" espera los b-numbers: no hay traducción que hacer.
  run_enrichment_kegg(c("b0001", "b0002"), organism = "eco", keyType = "kegg")
  expect_equal(recibido, c("b0001", "b0002"))
})

test_that("una traducción imposible se explica en vez de devolver una tabla vacia", {
  skip_if_not(HAS_ORGECDB, "falta org.EcK12.eg.db")
  r <- run_enrichment_kegg(c("noexisto1", "noexisto2"), organism = "eco",
                           keyType = "ncbi-geneid", OrgDb = "org.EcK12.eg.db",
                           from_keytype = "SYMBOL")
  expect_null(r$table)
  expect_match(r$error, "traducir")
})
