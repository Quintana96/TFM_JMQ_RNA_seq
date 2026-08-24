#' test-umbral-lfc-na.R
#' Regresión: el umbral de |log2FC| del test NO siempre es un número.
#'
#' Swish no recibe umbral de fold-change —trabaja sobre las replicas
#' inferenciales, no sobre la matriz prefiltrada—, así que el ajuste lo registra
#' como `NA_real_` para no atribuirle un test que no hizo. `%||%` no captura ese
#' NA (solo NULL y length 0), de modo que las ramas `if (lfc_thr > 0)` abortaban
#' el render con "valor ausente donde TRUE/FALSE es necesario": se caian el panel
#' de estado de la barra lateral y el volcano (y su descarga) en cuanto se
#' ajustaba con Swish.

test_that("has_lfc_threshold() no se cae con NA, NULL ni vacio", {
  # El caso que rompia: lo que deja Swish en state$deg_rv$lfc_threshold.
  expect_false(has_lfc_threshold(NA_real_))
  expect_false(has_lfc_threshold(NULL))
  expect_false(has_lfc_threshold(numeric(0)))
  expect_false(has_lfc_threshold(NA))
  expect_false(has_lfc_threshold(Inf))
  expect_false(has_lfc_threshold(0))
  expect_false(has_lfc_threshold(-1))
  expect_true(has_lfc_threshold(1))
  expect_true(has_lfc_threshold(0.5))
})

#' Tabla DEG mínima con la forma que devuelven los motores.
deg_tabla_minima <- function(n = 60) {
  data.frame(
    gene          = sprintf("g%03d", seq_len(n)),
    baseMean      = seq(10, 1000, length.out = n),
    log2FC        = seq(-3, 3, length.out = n),
    log2FC_shrunk = NA_real_,
    lfcSE         = NA_real_,
    stat          = seq(-4, 4, length.out = n),
    pvalue        = seq(1e-6, 0.9, length.out = n),
    padj          = c(rep(0.01, 10), rep(0.6, n - 10)),
    stringsAsFactors = FALSE
  )
}


test_that("el volcano sigue dibujando las líneas de umbral cuando SI lo hay", {
  deg <- deg_tabla_minima()

  srv <- function(input, output, session) {
    state <- new.env(parent = emptyenv())
    state$deg_rv <- shiny::reactiveValues(
      results = deg, fdr = 0.05, lfc_threshold = 1.5,
      contrast = "trt vs ctrl", method = "DESeq2",
      vst_mat = NULL, meta = NULL)
    ctx <- new.env(parent = emptyenv())
    ctx$deg_filtered <- shiny::reactive(deg)
    server_tab_deg_results(input, output, session, state, ctx)
  }

  shiny::testServer(srv, {
    # El título declara el umbral: es lo que distingue este caso del anterior y
    # asegura que el arreglo del NA no ha desactivado la rama buena.
    expect_true(grepl("umbral del test", output$deg_volcano_plot, fixed = TRUE))
  })
})

test_that("ninguna rama del umbral compara sin proteger contra NA", {
  # Guarda de estilo: el fallo original fue copiar la comparación a tres sitios y
  # protegerla solo en uno. Si vuelve a aparecer un `if (lfc_thr > 0)` suelto,
  # este test lo caza antes que el usuario.
  ficheros <- list.files(file.path(app_root, "R"), pattern = "[.]R$",
                         full.names = TRUE)
  sueltas <- unlist(lapply(ficheros, function(f) {
    l <- readLines(f, warn = FALSE)
    l <- l[!grepl("^\\s*#", l)]
    hits <- grep("\\blfc_thr\\s*>\\s*0", l, value = TRUE)
    hits[!grepl("has_lfc_threshold|is.finite", hits)]
  }))
  expect_equal(sueltas, character(0))
})
