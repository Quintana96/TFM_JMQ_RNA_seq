#' test-enriquecimiento-avanzado.R
#' Cubre lo que se anadio al enriquecimiento por encima del ORA/GSEA basico:
#' parametros de GSEA expuestos, gene sets propios en GMT, ORA direccional,
#' running score y comparacion ORA/GSEA.
#'
#' Los casos sinteticos comprueban la mecanica (que los parametros llegan al
#' motor y que las tablas se combinan bien) y los casos con el dataset humano
#' GSE52778 (airway) comprueban lo que solo se ve con anotacion real: que el
#' keyType decide el mapeo y que ORA y GSEA no responden lo mismo.

tiene_cp   <- function() requireNamespace("clusterProfiler", quietly = TRUE)
tiene_fgsea <- function() requireNamespace("fgsea", quietly = TRUE)

#' Escribe un .gmt temporal a partir de una lista nombrada de vectores de genes.
gmt_temporal <- function(sets) {
  f <- tempfile(fileext = ".gmt")
  writeLines(vapply(names(sets), function(n)
    paste(c(n, paste0("descripcion de ", n), sets[[n]]), collapse = "\t"),
    character(1)), f)
  f
}

#' Ranking descendente y sin empates sobre `n` genes sinteticos.
ranking_sintetico <- function(n = 2000) {
  stats::setNames(seq(3, -3, length.out = n), sprintf("g%04d", seq_len(n)))
}

# ── Gene sets propios (GMT) ─────────────────────────────────────────────────

test_that("read_gene_sets_gmt devuelve TERM2GENE y el rango de tamaños", {
  skip_if_not(tiene_cp())
  genes <- sprintf("g%04d", 1:500)
  f <- gmt_temporal(list(chico = genes[1:12], grande = genes[100:179]))

  gs <- read_gene_sets_gmt(f)
  expect_null(gs$error)
  expect_equal(gs$n_sets, 2L)
  expect_setequal(names(gs$term2gene), c("term", "gene"))
  # El rango de tamaños es lo que permite ver de un vistazo si el GMT es
  # utilizable: conjuntos de 2 genes o de 5.000 no dan nada interpretable.
  expect_equal(sort(gs$sizes), c(12L, 80L))
  expect_equal(gs$n_genes, 92L)
  expect_true(grepl("2 conjuntos", gmt_summary_text(gs)))
  expect_true(grepl("12 y 80 genes", gmt_summary_text(gs)))
  # read.gmt puede devolver factores segun la version; el resto del codigo hace
  # intersect() contra character y un factor lo romperia en silencio.
  expect_type(gs$term2gene$gene, "character")
})

test_that("un GMT ausente o ilegible se reporta, no rompe", {
  gs <- read_gene_sets_gmt(NULL)
  expect_null(gs$term2gene)
  expect_true(nzchar(gs$error))
  expect_equal(gmt_summary_text(gs), gs$error)

  vacio <- tempfile(fileext = ".gmt")
  writeLines(character(0), vacio)
  skip_if_not(tiene_cp())
  expect_true(nzchar(read_gene_sets_gmt(vacio)$error))
})

test_that("mapping_against_sets mide el solapamiento contra el GMT", {
  t2g <- data.frame(term = rep("a", 3), gene = c("x", "y", "z"),
                    stringsAsFactors = FALSE)
  mp <- mapping_against_sets(c("x", "y", "fuera", "tampoco"), t2g)
  expect_equal(mp$n_input, 4L)
  expect_equal(mp$n_mapped, 2L)
  expect_equal(mp$rate, 0.5)
  # Sin OrgDb no hay keyType: el "diccionario" es el propio fichero.
  expect_equal(mp$source, "GMT")
})

test_that("el ORA sobre gene sets propios corre sin ningun OrgDb", {
  skip_if_not(tiene_cp())
  genes <- sprintf("g%04d", 1:2000)
  gs <- read_gene_sets_gmt(gmt_temporal(list(
    objetivo = genes[1:40],
    ruido    = genes[1000:1059]
  )))
  lista <- genes[1:35]

  res <- run_enrichment_gmt(lista, universe = genes, term2gene = gs$term2gene)
  expect_null(res$error)
  expect_true("objetivo" %in% res$table$ID)
  expect_false("ruido" %in% res$table$ID)
  expect_lt(res$table$p.adjust[res$table$ID == "objetivo"], 1e-10)
  expect_equal(res$mapping$n_mapped, 35L)
})

test_that("sin GMT cargado, ORA y GSEA piden el fichero en vez de fallar", {
  expect_true(grepl("GMT", run_enrichment_gmt(c("a", "b"), term2gene = NULL)$error))
  skip_if_not(tiene_fgsea() && tiene_cp())
  expect_true(grepl("GMT", run_gsea(ranking_sintetico(50), ont = "GMT",
                                    term2gene = NULL)$error))
})

# ── Parametros de GSEA expuestos ────────────────────────────────────────────

test_that("eps = 0 sustituye el p-valor truncado por el exacto", {
  skip_if_not(tiene_cp() && tiene_fgsea())
  rk <- ranking_sintetico()
  # Conjunto repartido por la cabeza del ranking: suficiente para bajar del
  # suelo de 1e-10 sin que el multilevel se vuelva caro.
  gs <- read_gene_sets_gmt(gmt_temporal(list(
    objetivo = names(rk)[seq(1, 300, by = 15)],
    ruido    = names(rk)[seq(2, 1800, by = 30)]
  )))

  truncado <- suppressWarnings(run_gsea(rk, ont = "GMT", term2gene = gs$term2gene,
                                        exponent = 0, pvalueCutoff = 1, eps = 1e-10))
  exacto <- suppressWarnings(run_gsea(rk, ont = "GMT", term2gene = gs$term2gene,
                                      exponent = 0, pvalueCutoff = 1, eps = 0))

  p_trunc <- truncado$table$pvalue[truncado$table$ID == "objetivo"]
  p_exact <- exacto$table$pvalue[exacto$table$ID == "objetivo"]
  # Con el default, todos los conjuntos muy significativos empatan en 1e-10 y su
  # orden relativo se pierde; con eps = 0 se separan.
  expect_equal(p_trunc, 1e-10)
  expect_lt(p_exact, 1e-10)
  # El NES no depende del eps: cambia la estimacion del p-valor, no el estadistico.
  expect_equal(truncado$table$NES[truncado$table$ID == "objetivo"],
               exacto$table$NES[exacto$table$ID == "objetivo"])
})

test_that("pvalueCutoff = 1 distingue 'nada llega a 0,05' de un fallo", {
  skip_if_not(tiene_cp() && tiene_fgsea())
  rk <- ranking_sintetico()
  set.seed(11)
  gs <- read_gene_sets_gmt(gmt_temporal(list(
    azar1 = sample(names(rk), 60),
    azar2 = sample(names(rk), 80)
  )))

  # Conjuntos al azar: ninguno deberia pasar el corte por defecto.
  estricto <- suppressWarnings(run_gsea(rk, ont = "GMT", term2gene = gs$term2gene,
                                        exponent = 0, pvalueCutoff = 0.05))
  todo <- suppressWarnings(run_gsea(rk, ont = "GMT", term2gene = gs$term2gene,
                                    exponent = 0, pvalueCutoff = 1))
  expect_null(estricto$table)
  expect_equal(nrow(todo$table), 2L)
  # Y es la tabla completa la que permite ver que SI se testearon los conjuntos:
  # con el corte puesto, ese caso es indistinguible de un error de anotacion.
  expect_true(all(todo$table$p.adjust > 0.05))
})

test_that("minGSSize y maxGSSize acotan los conjuntos testeados", {
  skip_if_not(tiene_cp() && tiene_fgsea())
  rk <- ranking_sintetico()
  gs <- read_gene_sets_gmt(gmt_temporal(list(
    diminuto = names(rk)[seq(1, 60, by = 5)],
    mediano  = names(rk)[seq(1, 600, by = 10)]
  )))

  ambos <- suppressWarnings(run_gsea(rk, ont = "GMT", term2gene = gs$term2gene,
                                     exponent = 0, pvalueCutoff = 1,
                                     minGSSize = 5, maxGSSize = 500))
  expect_setequal(ambos$table$ID, c("diminuto", "mediano"))

  solo_grande <- suppressWarnings(run_gsea(rk, ont = "GMT", term2gene = gs$term2gene,
                                           exponent = 0, pvalueCutoff = 1,
                                           minGSSize = 30, maxGSSize = 500))
  expect_equal(solo_grande$table$ID, "mediano")

  ninguno <- suppressWarnings(run_gsea(rk, ont = "GMT", term2gene = gs$term2gene,
                                       exponent = 0, pvalueCutoff = 1,
                                       minGSSize = 5, maxGSSize = 8))
  expect_null(ninguno$table)
})

# ── ORA direccional ─────────────────────────────────────────────────────────

test_that("separar por direccion recupera la señal que la lista mezclada diluye", {
  skip_if_not(tiene_cp())
  genes <- sprintf("g%04d", 1:2000)
  set.seed(3)
  # Los conjuntos de relleno no son decorativos: enricher() restringe el fondo a
  # los genes que aparecen en el GMT, asi que con solo dos conjuntos el universo
  # se reduce a sus 80 genes y no queda contraste que medir.
  gs <- read_gene_sets_gmt(gmt_temporal(list(
    inducido  = genes[1:40],
    reprimido = genes[1961:2000],
    relleno1  = sample(genes[200:1800], 120),
    relleno2  = sample(genes[200:1800], 150)
  )))
  # Lista con las dos direcciones mezcladas, como la que recibe el ORA normal.
  sig <- data.frame(
    gene   = c(genes[1:35], genes[1970:2000]),
    log2FC = c(rep(2, 35), rep(-2, 31)),
    stringsAsFactors = FALSE
  )
  runner <- function(g) run_enrichment_gmt(g, universe = genes,
                                           term2gene = gs$term2gene)

  res <- run_ora_directional(sig, runner)
  expect_null(res$error)
  expect_setequal(unique(res$table$Direccion), c("Conjunto", "Al alza", "A la baja"))
  expect_equal(names(res$table)[1], "Direccion")
  expect_equal(res$n_genes[["Al alza"]], 35L)
  expect_equal(res$n_genes[["A la baja"]], 31L)

  p_mezcla <- res$table$p.adjust[res$table$Direccion == "Conjunto" &
                                   res$table$ID == "inducido"]
  p_alza <- res$table$p.adjust[res$table$Direccion == "Al alza" &
                                 res$table$ID == "inducido"]
  # El mismo conjunto es varios ordenes de magnitud mas significativo cuando la
  # lista no arrastra los 31 genes de la direccion contraria: eso es la dilucion.
  expect_lt(p_alza, p_mezcla)
  # Y cada mitad solo ve su conjunto.
  expect_false("reprimido" %in% res$table$ID[res$table$Direccion == "Al alza"])
})

test_that("run_ora_directional aisla las direcciones sin genes o sin terminos", {
  # Runner falso: no hace falta clusterProfiler para comprobar la combinacion.
  runner <- function(g) {
    if (length(g) < 3) return(list(table = NULL, error = "Sin terminos enriquecidos.",
                                   mapping = NULL))
    list(table = data.frame(ID = "T1", Description = "termino", p.adjust = 0.01,
                            Count = length(g), stringsAsFactors = FALSE),
         error = NULL, mapping = list(n_input = length(g), n_mapped = length(g),
                                      rate = 1, keytype = "X", source = "test"))
  }
  sig <- data.frame(gene = c("a", "b", "c", "d"), log2FC = c(1, 2, 3, 4),
                    stringsAsFactors = FALSE)

  res <- run_ora_directional(sig, runner)
  # Sin genes a la baja, esa direccion se reporta como tal en vez de tumbar todo.
  expect_true(grepl("Sin genes", res$errores[["A la baja"]]))
  expect_setequal(unique(res$table$Direccion), c("Conjunto", "Al alza"))
  # El mapeo que se muestra es el del conjunto completo, no el de una mitad.
  expect_equal(res$mapping$n_input, 4L)

  expect_true(grepl("log2FC", run_ora_directional(
    data.frame(gene = c("a", "b")), runner)$error))
  expect_true(nzchar(run_ora_directional(NULL, runner)$error))
})

test_that("el dotplot reparte el top entre direcciones y no colapsa etiquetas", {
  df <- data.frame(
    Direccion   = rep(c("Al alza", "A la baja"), each = 4),
    ID          = paste0("T", 1:8),
    Description = rep(c("ruta A", "ruta B", "ruta C", "ruta D"), 2),
    p.adjust    = c(1e-8, 1e-7, 1e-6, 1e-5, 0.01, 0.02, 0.03, 0.04),
    Count       = 5:12,
    GeneRatio   = "5/50",
    stringsAsFactors = FALSE
  )
  top <- enrichment_dotplot_data(df, top_n = 2)
  # Sin reparto, las cuatro filas mas significativas serian todas "Al alza" y la
  # comparacion entre direcciones desapareceria del grafico.
  expect_equal(as.integer(table(top$Direccion)), c(2L, 2L))
  expect_equal(length(unique(top$plot_label)), nrow(top))
  expect_true(all(grepl("\\[", top$plot_label)))

  # Sin columna Direccion se mantiene el comportamiento anterior.
  simple <- enrichment_dotplot_data(
    df[df$Direccion == "Al alza", setdiff(names(df), "Direccion")], top_n = 3)
  expect_equal(nrow(simple), 3L)
  expect_equal(simple$plot_label, simple$Description)
})

# ── Running score y leading edge ────────────────────────────────────────────

test_that("el running score reproduce el enrichment score de la tabla GSEA", {
  skip_if_not(tiene_cp() && tiene_fgsea())
  rk <- ranking_sintetico()
  conjunto <- names(rk)[seq(1, 400, by = 10)]
  gs <- read_gene_sets_gmt(gmt_temporal(list(objetivo = conjunto)))

  res <- suppressWarnings(run_gsea(rk, ont = "GMT", term2gene = gs$term2gene,
                                   exponent = 0, pvalueCutoff = 1))
  fila <- res$table[res$table$ID == "objetivo", ]

  genes <- gsea_term_genes("objetivo", ont = "GMT", term2gene = gs$term2gene)
  expect_setequal(genes$genes, conjunto)

  # gseaParam = exponent: la curva tiene que ser la del estadistico con el que se
  # calculo el NES, o el pico no coincidiria con el ES de la tabla.
  rs <- gsea_running_score(rk, genes$genes, gseaParam = 0)
  expect_null(rs$error)
  expect_equal(rs$n_hits, length(conjunto))
  expect_gt(nrow(rs$curve), 0)
  expect_equal(rs$es, fila$enrichmentScore, tolerance = 1e-6)
  # Una marca por gen del conjunto presente en el ranking.
  expect_equal(nrow(rs$ticks), length(conjunto))

  le <- leading_edge_genes(fila$core_enrichment)
  expect_gt(length(le), 0)
  # El leading edge es un subconjunto del conjunto, no el conjunto entero.
  expect_true(all(le %in% conjunto))
  expect_lte(length(le), length(conjunto))
})

test_that("el running score avisa cuando el conjunto no toca el ranking", {
  skip_if_not(tiene_fgsea())
  rk <- ranking_sintetico(100)
  rs <- gsea_running_score(rk, c("otro1", "otro2"))
  # Es el sintoma tipico de un keyType equivocado, y hay que decirlo en vez de
  # devolver un grafico vacio.
  expect_true(grepl("keyType", rs$error))
  expect_null(rs$curve)
  expect_true(nzchar(gsea_running_score(NULL, "a")$error))
})

test_that("leading_edge_genes parte la cadena core_enrichment", {
  expect_equal(leading_edge_genes("A/B/C"), c("A", "B", "C"))
  expect_equal(leading_edge_genes(NA_character_), character(0))
  expect_equal(leading_edge_genes(NULL), character(0))
  expect_equal(leading_edge_genes(""), character(0))
})

# ── Comparacion ORA / GSEA ──────────────────────────────────────────────────

test_that("compare_ora_gsea calcula el solapamiento y los terminos solo de GSEA", {
  ora <- data.frame(ID = c("T1", "T2", "T3"),
                    Description = c("a", "b", "c"),
                    p.adjust = c(0.001, 0.01, 0.30),
                    stringsAsFactors = FALSE)
  gsea <- data.frame(ID = c("T2", "T4", "T5"),
                     Description = c("b", "d", "e"),
                     NES = c(2.1, -1.8, 1.5),
                     p.adjust = c(0.02, 0.001, 0.04),
                     stringsAsFactors = FALSE)

  cmp <- compare_ora_gsea(ora, gsea)
  # T3 no es significativo, asi que no entra en la comparacion.
  expect_equal(cmp$n_ora, 2L)
  expect_equal(cmp$n_gsea, 3L)
  expect_equal(cmp$n_comun, 1L)
  expect_equal(cmp$jaccard, 1 / 4)
  expect_setequal(cmp$solo_gsea$ID, c("T4", "T5"))
  expect_equal(cmp$solo_ora$ID, "T1")
  expect_equal(cmp$comunes$ID, "T2")
  expect_equal(cmp$comunes$NES, 2.1)

  # Sin terminos en ninguno de los dos, el Jaccard es indefinido y no 0: 0
  # sugeriria discrepancia total entre dos analisis que no han visto nada.
  vacio <- compare_ora_gsea(NULL, NULL)
  expect_true(is.na(vacio$jaccard))
  expect_equal(vacio$n_comun, 0L)
  expect_null(compare_ora_gsea_table(vacio))
})

test_that("compare_ora_gsea_table etiqueta quien ve cada termino", {
  ora <- data.frame(ID = c("T1", "T2"), Description = c("a", "b"),
                    p.adjust = c(0.001, 0.01), stringsAsFactors = FALSE)
  gsea <- data.frame(ID = c("T2", "T3"), Description = c("b", "c"),
                     NES = c(2, -2), p.adjust = c(0.02, 0.03),
                     stringsAsFactors = FALSE)
  tb <- compare_ora_gsea_table(compare_ora_gsea(ora, gsea))

  expect_equal(nrow(tb), 3L)
  expect_setequal(tb$Visto, c("Solo GSEA", "Ambos", "Solo ORA"))
  expect_equal(tb$padj_GSEA[tb$ID == "T3"], 0.03)
  expect_true(is.na(tb$padj_ORA[tb$ID == "T3"]))
  expect_equal(tb$padj_ORA[tb$ID == "T2"], 0.01)
  expect_equal(tb$padj_GSEA[tb$ID == "T2"], 0.02)
})

# ── Integracion del modulo server ───────────────────────────────────────────

test_that("el modulo server encadena GMT, ORA direccional, GSEA y comparacion", {
  skip_if_not(tiene_cp() && tiene_fgsea())
  genes <- sprintf("g%04d", 1:2000)
  set.seed(5)
  gmt <- gmt_temporal(list(
    inducido  = genes[1:40],
    reprimido = genes[1961:2000],
    relleno1  = sample(genes[200:1800], 120),
    relleno2  = sample(genes[200:1800], 150)
  ))
  deg <- data.frame(
    gene   = genes,
    stat   = seq(6, -6, length.out = 2000),
    log2FC = seq(3, -3, length.out = 2000),
    pvalue = 10^-abs(seq(6, -6, length.out = 2000)),
    padj   = c(rep(0.001, 35), rep(0.5, 1934), rep(0.001, 31)),
    stringsAsFactors = FALSE
  )
  sig <- deg[deg$padj <= 0.05, ]

  srv <- function(input, output, session) {
    state <- new.env(parent = emptyenv())
    state$deg_rv <- shiny::reactiveValues(results = deg, fdr = 0.05, enrich = NULL)
    # El modulo busca aqui la anotacion con la que traducir identificadores.
    state$run_params_rv <- shiny::reactiveVal(list())
    ctx <- new.env(parent = emptyenv())
    ctx$deg_universe <- shiny::reactive(genes)
    ctx$deg_significant <- shiny::reactive(sig)
    server_tab_deg_enrich(input, output, session, state, ctx)
    session$userData$state <- state
  }

  shiny::testServer(srv, {
    session$setInputs(
      deg_ontology = "GMT", deg_enrich_approach = "ora",
      deg_gmt_file = data.frame(name = "sets.gmt", size = 1, type = "",
                                datapath = gmt, stringsAsFactors = FALSE),
      deg_ora_directional = TRUE, deg_enrich_readable = FALSE,
      deg_gsea_min_size = 10, deg_gsea_max_size = 500,
      deg_gsea_pcutoff = 1, deg_gsea_eps_exact = TRUE, deg_gsea_metric = "stat"
    )
    expect_true(grepl("4 conjuntos", as.character(output$deg_gmt_summary$html)))

    session$setInputs(deg_run_enrich_btn = 1)
    e <- session$userData$state$deg_rv$enrich
    expect_gt(e$n_terminos, 0)
    expect_true(e$direccional)
    # Sin OrgDb de por medio: el GMT es el que aporta los conjuntos.
    expect_equal(e$orgdb, "no aplica (GMT)")
    expect_true(grepl("4 conjuntos", e$gmt))

    session$setInputs(deg_enrich_approach = "gsea")
    session$setInputs(deg_run_enrich_btn = 2)
    e2 <- session$userData$state$deg_rv$enrich
    expect_equal(e2$enfoque, "GSEA")
    # Los parametros de la interfaz llegan al motor y quedan registrados para el
    # informe: sin ellos el resultado no seria reproducible.
    expect_equal(e2$gsea_eps, 0)
    expect_equal(e2$gsea_pcutoff, 1)
    expect_equal(e2$gsea_min_size, 10)
    expect_equal(e2$gsea_max_size, 500)
    expect_gt(e2$n_terminos, 0)

    session$setInputs(deg_gsea_term = "inducido")
    expect_true(grepl("Leading edge", as.character(output$deg_gsea_leading$html)))

    session$setInputs(deg_run_enrich_compare_btn = 1)
    resumen <- as.character(output$deg_enrich_compare_summary$html)
    expect_true(grepl("Jaccard", resumen))
    expect_true(grepl("ORA", resumen))
  })
})

# ── Datos humanos reales (GSE52778 airway) ──────────────────────────────────

airway_dir <- "/Users/usuario/Desktop/UEM/TFM/datasets_test/GSE52778_airway"

#' Tabla DEG minima a partir del dataset airway.
#'
#' Se usa una t de Welch sobre logCPM y no un motor completo: aqui lo que se
#' testea es el enriquecimiento, y un ranking real de IDs ENSEMBL humanos ya
#' basta para eso. Ahorra un ajuste de DESeq2 sobre 16.000 genes por test.
airway_deg <- function() {
  cm <- utils::read.delim(file.path(airway_dir, "counts_matrix.tsv"),
                          row.names = 1, check.names = FALSE)
  meta <- utils::read.delim(file.path(airway_dir, "samplesheet.tsv"),
                            stringsAsFactors = FALSE)
  cm <- cm[rowSums(cm >= 10) >= 4, , drop = FALSE]
  lg <- log2(t(t(cm) / colSums(cm)) * 1e6 + 1)
  grp <- meta$condition[match(colnames(lg), meta$sample_id)]
  st <- apply(lg, 1, function(x) {
    r <- try(stats::t.test(x[grp == "trt"], x[grp == "untrt"]), silent = TRUE)
    if (inherits(r, "try-error")) NA_real_ else unname(r$statistic)
  })
  df <- data.frame(gene = rownames(lg), stat = st,
                   log2FC = rowMeans(lg[, grp == "trt", drop = FALSE]) -
                            rowMeans(lg[, grp == "untrt", drop = FALSE]),
                   stringsAsFactors = FALSE)
  df$pvalue <- 2 * stats::pt(-abs(df$stat), df = 6)
  df$padj <- stats::p.adjust(df$pvalue, "BH")
  df[!is.na(df$stat), ]
}

test_that("el keyType decide el mapeo, no el OrgDb", {
  skip_if_not(requireNamespace("org.Hs.eg.db", quietly = TRUE))
  skip_if_not(dir.exists(airway_dir))
  deg <- airway_deg()

  ens <- gene_mapping_rate(deg$gene, "org.Hs.eg.db", "ENSEMBL")
  sym <- gene_mapping_rate(deg$gene, "org.Hs.eg.db", "SYMBOL")
  # Los IDs de la matriz son ENSEMBL: con el keyType correcto mapea la gran
  # mayoria y con SYMBOL, el mismo OrgDb y los mismos genes, no mapea nada. Es
  # el fallo silencioso contra el que avisa la interfaz.
  expect_gt(ens$rate, 0.8)
  expect_lt(sym$rate, 0.01)
  expect_true(grepl("ENSEMBL", mapping_rate_text(ens)))
})

test_that("GSEA sobre datos reales: corte de p, running score y ORA que no lo ve", {
  skip_if_not(tiene_cp() && tiene_fgsea())
  skip_if_not(requireNamespace("org.Hs.eg.db", quietly = TRUE))
  skip_if_not(dir.exists(airway_dir))
  deg <- airway_deg()
  rk <- deg_ranking_metric(deg, "stat")
  expect_gt(rk$n, 10000)

  # Corte a 1: la tabla trae TODOS los conjuntos testeados, significativos o no.
  gsea <- suppressWarnings(run_gsea(rk$ranked, ont = "CC", OrgDb = "org.Hs.eg.db",
                                    keyType = "ENSEMBL", exponent = 0,
                                    pvalueCutoff = 1, minGSSize = 50,
                                    maxGSSize = 300, readable = TRUE))
  expect_null(gsea$error)
  n_sig <- sum(gsea$table$p.adjust <= 0.05)
  expect_gt(nrow(gsea$table), n_sig)
  expect_gt(n_sig, 0)
  # readable = TRUE: core_enrichment sale en simbolos, no en ENSG00000...
  expect_false(grepl("ENSG", gsea$table$core_enrichment[1]))

  # Running score del termino mas significativo, reconstruido desde el OrgDb.
  mejor <- gsea$table[order(gsea$table$p.adjust), ][1, ]
  tg <- gsea_term_genes(mejor$ID, ont = "CC", OrgDb = "org.Hs.eg.db",
                        keyType = "ENSEMBL")
  expect_gt(length(tg$genes), 0)
  rs <- gsea_running_score(rk$ranked, tg$genes, gseaParam = 0)
  expect_null(rs$error)
  # El conjunto que recupera GOALL es el mismo que testeo gseGO (por eso se usa
  # GOALL y no GO: con GO faltarian los genes de los terminos descendientes).
  expect_equal(rs$n_hits, mejor$setSize)
  expect_equal(rs$es, mejor$enrichmentScore, tolerance = 1e-6)

  # El mismo contraste por ORA sobre la lista umbralizada.
  sig <- deg[!is.na(deg$padj) & deg$padj <= 0.05, ]
  expect_gt(nrow(sig), 100)
  ora <- run_enrichment_go(sig$gene, universe = deg$gene, OrgDb = "org.Hs.eg.db",
                           ont = "CC", keyType = "ENSEMBL", readable = TRUE)
  cmp <- compare_ora_gsea(ora$table, gsea$table)
  expect_equal(cmp$n_gsea, n_sig)
  # Los terminos que solo ve GSEA son la razon de tener los dos enfoques: señal
  # coordinada que ningun gen individual lleva mas alla del corte de la lista.
  expect_gt(cmp$n_gsea - cmp$n_comun, 0)
  tb <- compare_ora_gsea_table(cmp)
  expect_true("Solo GSEA" %in% tb$Visto)
  expect_equal(nrow(tb), cmp$n_ora + cmp$n_gsea - cmp$n_comun)
})
