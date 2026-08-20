#' test-reextraccion.R
#' El FDR objetivo y el umbral del test se recalculan en vivo sobre el ajuste
#' guardado, en lugar de exigir relanzar. La propiedad que hay que sostener es
#' una y no admite matices: **reextraer tiene que dar exactamente lo mismo que
#' reajustar**. Si no, la aplicacion muestra en vivo una tabla que ningun
#' analisis reproduciria, que es peor que obligar a pulsar un boton.
#'
#' El riesgo concreto que estos tests vigilan es el contrario al que resuelven:
#' que alguien "optimice" la reextraccion recortando la tabla ya calculada al
#' nuevo umbral. Eso seria un filtro post hoc y la lista no tendria la FDR que
#' declara (McCarthy y Smyth, 2009).

#' Matriz con senal conocida Y cola de baja expresion.
#'
#' La cola no es decorativa: el filtrado independiente de DESeq2 descarta genes
#' poco expresados, y sin ellos no descarta ninguno, de modo que cambiar `alpha`
#' no cambiaria nada y el test de mas abajo pasaria por el motivo equivocado.
#' Una matriz real siempre tiene esa cola.
datos_con_senal <- function(n_genes = 1500, n_de = 200, seed = 11) {
  withr::with_seed(seed, {
    n_s <- 8
    n_bajos <- floor(n_genes / 2)
    base <- c(stats::rnbinom(n_genes - n_bajos, mu = 200, size = 5) + 10,
              stats::rnbinom(n_bajos, mu = 2, size = 2))
    m <- vapply(seq_len(n_s), function(j) {
      mu <- base
      if (j > 4) mu[seq_len(n_de)] <- mu[seq_len(n_de)] * 2.5
      stats::rnbinom(n_genes, mu = mu, size = 10)
    }, numeric(n_genes))
    # Los genes de senal son los primeros, que estan en la parte expresada.
    rownames(m) <- sprintf("g%05d", seq_len(n_genes))
    colnames(m) <- sprintf("s%d", seq_len(n_s))
    list(counts = m,
         meta = data.frame(sample_id = colnames(m),
                           condition = rep(c("ctrl", "trt"), each = 4),
                           stringsAsFactors = FALSE))
  })
}

ajusta <- function(d, metodo, fdr = 0.05, lfc = 0) {
  run_deg(d$counts, d$meta, method = metodo, ref_level = "ctrl",
          contrast_num = "trt", fdr = fdr, lfc_threshold = lfc)
}

for (motor in c("DESeq2", "edgeR", "limma-voom")) {
  test_that(paste0("reextraer equivale a reajustar (", motor, ")"), {
    skip_if_not_installed(if (identical(motor, "limma-voom")) "limma" else motor)
    d <- datos_con_senal()

    ajuste <- ajusta(d, motor, fdr = 0.05, lfc = 0)
    skip_if(is.null(ajuste$table), "el motor no ha producido tabla en este entorno")

    # Referencia: reajuste completo con los parametros nuevos.
    referencia <- ajusta(d, motor, fdr = 0.01, lfc = 1)
    # Candidato: reextraccion desde el ajuste anterior.
    reex <- deg_reextract(ajuste$fit, fdr = 0.01, lfc_threshold = 1)

    expect_null(reex$error)
    expect_equal(nrow(reex$table), nrow(referencia$table))
    cols <- intersect(names(referencia$table), names(reex$table))
    expect_equal(reex$table[, cols], referencia$table[, cols], tolerance = 1e-8)
  })
}

test_that("cambiar el FDR cambia que genes son EVALUABLES, no solo el corte", {
  skip_if_not_installed("DESeq2")
  d <- datos_con_senal()
  ajuste <- ajusta(d, "DESeq2", fdr = 0.05)
  skip_if(is.null(ajuste$table))

  a05 <- deg_reextract(ajuste$fit, fdr = 0.05)$table
  a01 <- deg_reextract(ajuste$fit, fdr = 0.01)$table

  # Esta es la diferencia entre reextraer y recolorear: el filtrado independiente
  # de DESeq2 elige el umbral de expresion que maximiza los significativos AL
  # NIVEL PEDIDO, asi que con alpha distinto el conjunto de genes con padj es
  # otro. Un simple recorte de la tabla no podria producir esto.
  expect_false(identical(is.na(a05$padj), is.na(a01$padj)))
  # Y los p-valores sin ajustar no cambian: el modelo es el mismo.
  expect_equal(a05$pvalue, a01$pvalue)
})

test_that("el umbral del test entra en el modelo, no recorta la tabla despues", {
  skip_if_not_installed("DESeq2")
  d <- datos_con_senal()
  ajuste <- ajusta(d, "DESeq2", fdr = 0.05, lfc = 0)
  skip_if(is.null(ajuste$table))

  sin <- deg_reextract(ajuste$fit, fdr = 0.05, lfc_threshold = 0)$table
  con <- deg_reextract(ajuste$fit, fdr = 0.05, lfc_threshold = 1)$table

  # Mismo numero de filas: no se ha filtrado nada, se ha testeado otra hipotesis.
  expect_equal(nrow(sin), nrow(con))
  # Los p-valores cambian porque cambia H0. Si solo se hubiera recortado la
  # tabla, serian identicos: es la comprobacion que distingue las dos cosas.
  comunes <- !is.na(sin$pvalue) & !is.na(con$pvalue)
  expect_false(isTRUE(all.equal(sin$pvalue[comunes], con$pvalue[comunes])))
})

test_that("los motores sin objeto reutilizable devuelven su tabla intacta", {
  d <- datos_con_senal(n_genes = 400, n_de = 60)
  ajuste <- run_deg(d$counts, d$meta, method = "Wilcoxon", ref_level = "ctrl",
                    contrast_num = "trt", fdr = 0.05)
  skip_if(is.null(ajuste$table))

  expect_equal(ajuste$fit$engine, "estatico")
  reex <- deg_reextract(ajuste$fit, fdr = 0.01)
  expect_null(reex$error)
  # Su padj es la correccion BH de los p-valores del test y no depende del nivel
  # objetivo: cambiar el FDR cambia donde se corta, no lo que se calcula.
  expect_equal(reex$table$padj, ajuste$table$padj)
})

test_that("Wilcoxon y dearseq no declaran un umbral de fold-change que no aplican", {
  d <- datos_con_senal(n_genes = 400, n_de = 60)
  ajuste <- run_deg(d$counts, d$meta, method = "Wilcoxon", ref_level = "ctrl",
                    contrast_num = "trt", fdr = 0.05, lfc_threshold = 1)
  skip_if(is.null(ajuste$table))

  # El motor acepta el argumento por uniformidad de la interfaz pero no lo usa.
  # Declararlo hacia que el banner y el informe afirmaran "H0: |log2FC| <= 1
  # dentro del test" sobre un ajuste que testeo H0: log2FC = 0.
  expect_true(is.na(ajuste$lfc_threshold))
  expect_false(has_lfc_threshold(ajuste$lfc_threshold))
})

test_that("el modo de outliers que cambia el ajuste se rechaza con un mensaje", {
  skip_if_not_installed("DESeq2")
  d <- datos_con_senal(n_genes = 800, n_de = 100)
  ajuste <- run_deg(d$counts, d$meta, method = "DESeq2", ref_level = "ctrl",
                    contrast_num = "trt", fdr = 0.05, outliers = "na")
  skip_if(is.null(ajuste$table))

  # "keep" solo desactiva cooksCutoff en results(): es lectura, se puede.
  expect_null(deg_reextract(ajuste$fit, outliers = "keep")$error)
  # "refit" rebaja minReplicatesForReplace, que es argumento de DESeq(): no.
  malo <- deg_reextract(ajuste$fit, outliers = "refit")
  expect_null(malo$table)
  expect_true(grepl("relanza", malo$error, ignore.case = TRUE))
})

test_that("sin ajuste guardado la reextraccion lo dice en lugar de fallar", {
  r <- deg_reextract(NULL, fdr = 0.05)
  expect_null(r$table)
  expect_true(nzchar(r$error))

  r2 <- deg_reextract(list(engine = "Swish"), fdr = 0.05)
  expect_null(r2$table)
  expect_true(nzchar(r2$error))
})

test_that("el encogido se reutiliza y no se recalcula en cada reextraccion", {
  skip_if_not_installed("DESeq2")
  skip_if_not(isTRUE(HAS_APEGLM) || isTRUE(HAS_ASHR))
  d <- datos_con_senal()
  ajuste <- ajusta(d, "DESeq2", fdr = 0.05)
  skip_if(is.null(ajuste$table))
  skip_if(identical(ajuste$shrink, "ninguno"))

  reex <- deg_reextract(ajuste$fit, fdr = 0.01)
  # El encogido depende del COEFICIENTE, no del nivel de significacion, asi que
  # la columna tiene que ser la misma. Recalcularlo costaria mas que el propio
  # ajuste (6,4 s frente a 5,1 s medidos sobre 20.000 genes).
  expect_equal(reex$table$log2FC_shrunk, ajuste$table$log2FC_shrunk)
  expect_equal(reex$shrink, ajuste$shrink)
})

test_that("las tres guardas de la reextraccion deciden bien", {
  p  <- list(fdr = 0.05, lfc_threshold = 0, use_ihw = FALSE, outliers = "na")
  p2 <- list(fdr = 0.01, lfc_threshold = 0, use_ihw = FALSE, outliers = "na")
  fit <- list(engine = "estatico", table = data.frame())

  # Cambio real y ajuste al dia: se reextrae.
  expect_true(deg_reextract_needed(fit, p2, p, stale = FALSE))
  # Sin cambio: no se toca nada (ni se escribe en el registro de auditoria).
  expect_false(deg_reextract_needed(fit, p, p, stale = FALSE))
  # Sin ajuste guardado: nada que reutilizar.
  expect_false(deg_reextract_needed(NULL, p2, p, stale = FALSE))
  # Ajuste desactualizado: abstenerse. Reextraer aqui daria una tabla que no
  # corresponde a ningun modelo y taparia el aviso de "relanza".
  expect_false(deg_reextract_needed(fit, p2, p, stale = TRUE))
  # Primera extraccion (aun no hay parametros previos).
  expect_true(deg_reextract_needed(fit, p, NULL, stale = FALSE))
})

# ── Integracion en el grafo reactivo ────────────────────────────────────────
#
# Lo anterior comprueba que la reextraccion calcula bien. Esto comprueba lo que
# el usuario ve: que mover el FDR actualiza la tabla sin repetir el ajuste, y
# que cuando el ajuste SI se ha quedado obsoleto la aplicacion lo dice en vez de
# recalcular sobre un modelo que ya no corresponde.

# Nota sobre los avisos: `testServer` vuelca todas las salidas en cada
# `setInputs`, de modo que dibuja el volcano y el MA. Con un fixture que tiene
# cola de baja expresion, plotly avisa "Ignoring N observations" por los genes
# sin padj. Es un aviso del dibujo, no del codigo bajo prueba, y no se silencia
# a proposito: taparlo con suppressWarnings() esconderia tambien los que si
# importan.

servidor_deg_con_ajuste <- function(d, ajuste, tmp, firma = NULL) {
  function(input, output, session) {
    state <- create_app_state(session)
    state$outputs_dir <- tmp
    server_tab_deg(input, output, session, state)
    # Se inyecta el resultado del ajuste como si ya se hubiera pulsado el boton:
    # el objetivo del test es la reextraccion, no volver a probar el ajuste.
    state$deg_rv$counts         <- d$counts
    state$deg_rv$meta           <- d$meta
    state$deg_rv$method         <- "DESeq2"
    state$deg_rv$results        <- ajuste$table
    state$deg_rv$fit            <- ajuste$fit
    state$deg_rv$extract_params <- ajuste$fit$extract
    state$deg_rv$fdr            <- 0.05
    state$deg_rv$lfc_threshold  <- 0
    state$deg_rv$run_at         <- as.POSIXct("2026-08-20 10:00:00")
    state$deg_rv$fit_signature  <- firma
    session$userData$state <- state
  }
}

inputs_base <- list(deg_fdr_target = 0.05, deg_lfc_threshold = 0,
                    deg_use_ihw = FALSE, deg_outliers = "na")

test_that("mover el FDR recalcula la tabla sin repetir el ajuste", {
  skip_if_not_installed("DESeq2")
  d <- datos_con_senal()
  ajuste <- ajusta(d, "DESeq2", fdr = 0.05, lfc = 0)
  skip_if(is.null(ajuste$table))

  tmp <- withr::local_tempdir()
  withr::local_dir(tmp)

  shiny::testServer(servidor_deg_con_ajuste(d, ajuste, tmp), {
    do.call(session$setInputs, inputs_base)
    st <- session$userData$state
    n0        <- sum(st$deg_rv$results$padj <= 0.05, na.rm = TRUE)
    ajuste_en <- st$deg_rv$run_at
    fit_antes <- st$deg_rv$fit

    session$setInputs(deg_fdr_target = 0.01)
    session$elapse(400)   # supera el debounce de 300 ms

    expect_equal(st$deg_rv$fdr, 0.01)
    n1 <- sum(st$deg_rv$results$padj <= 0.01, na.rm = TRUE)
    expect_lt(n1, n0)
    # La marca de la reextraccion existe...
    expect_false(is.null(st$deg_rv$reextracted_at))
    # ...y el ajuste NO se ha repetido: misma hora y el mismo objeto.
    expect_equal(st$deg_rv$run_at, ajuste_en)
    expect_identical(st$deg_rv$fit, fit_antes)
  })
})

test_that("el umbral del test tambien es un control en vivo", {
  skip_if_not_installed("DESeq2")
  d <- datos_con_senal()
  ajuste <- ajusta(d, "DESeq2", fdr = 0.05, lfc = 0)
  skip_if(is.null(ajuste$table))

  tmp <- withr::local_tempdir()
  withr::local_dir(tmp)

  shiny::testServer(servidor_deg_con_ajuste(d, ajuste, tmp), {
    do.call(session$setInputs, inputs_base)
    st <- session$userData$state
    p_antes <- st$deg_rv$results$pvalue

    session$setInputs(deg_lfc_threshold = 1)
    session$elapse(400)

    expect_equal(st$deg_rv$lfc_threshold, 1)
    # Los p-valores CAMBIAN porque cambia la hipotesis nula. Si solo se hubiera
    # recortado la tabla serian los mismos: es la diferencia entre reextraer y
    # filtrar a posteriori.
    expect_false(isTRUE(all.equal(p_antes, st$deg_rv$results$pvalue)))
  })
})

test_that("con el ajuste desactualizado se avisa y no se reextrae", {
  skip_if_not_installed("DESeq2")
  d <- datos_con_senal()
  ajuste <- ajusta(d, "DESeq2", fdr = 0.05, lfc = 0)
  skip_if(is.null(ajuste$table))

  tmp <- withr::local_tempdir()
  withr::local_dir(tmp)

  # Firma que no puede coincidir con la que produce la interfaz: simula haber
  # cambiado un parametro que define el ajuste (motor, diseno, prefiltrado...).
  firma_vieja <- list(metodo = "otro-motor-cualquiera")

  shiny::testServer(servidor_deg_con_ajuste(d, ajuste, tmp, firma = firma_vieja), {
    do.call(session$setInputs, inputs_base)
    st <- session$userData$state
    tabla_antes <- st$deg_rv$results

    # El aviso aparece...
    expect_false(is.null(output$deg_stale_warning))
    expect_true(grepl("otro ajuste", as.character(output$deg_stale_warning$html)))

    # ...y mover el FDR no toca la tabla: reextraer de un ajuste que ya no
    # corresponde daria un resultado de aspecto normal y sin modelo detras.
    session$setInputs(deg_fdr_target = 0.01)
    session$elapse(400)
    expect_identical(st$deg_rv$results, tabla_antes)
    expect_null(st$deg_rv$reextracted_at)
  })
})
