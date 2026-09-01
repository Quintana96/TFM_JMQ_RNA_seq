#' test-tablas-decimales.R
#' Las tablas numericas de la pestana 4 se muestran a tres decimales.
#'
#' El caso que obliga a que esto tenga test es el p-valor. `round(4.42e-46, 3)`
#' es 0, y ese 4,42e-46 es el p ajustado de CRISPLD2 en GSE52778, uno de los
#' genes de control que la memoria reporta. Redondearlo a decimales fijos no lo
#' redondea: lo borra, y lo deja indistinguible de un gen no significativo.
#'
#' La segunda propiedad que se fija aqui es que el formato es de LECTURA. El
#' dato que viaja al widget —y por tanto el que se ordena, se filtra y se
#' descarga— conserva toda su precision.

#' Que formateador de DT le toca a cada columna, leido del propio widget.
#' DT 0.34 los guarda en `columnDefs[[i]]$render`, no en `rowCallback`.
formato_por_columna <- function(dt, df) {
  out <- setNames(rep("ninguno", length(df)), names(df))
  for (e in dt$x$options$columnDefs) {
    if (is.null(e$render)) next
    tipo <- if (grepl("formatSignif", e$render)) "signif"
            else if (grepl("formatRound", e$render)) "round" else "otro"
    for (t in e$targets) out[[t + 1L]] <- tipo   # targets es 0-based
  }
  out
}

tabla_deg <- function() {
  data.frame(
    gene     = c("CRISPLD2", "DUSP1", "SPARCL1"),
    baseMean = c(3379.9876543, 12345.6789, 88.1),
    log2FC   = c(2.6314159, 2.9512345, 4.5678901),
    lfcSE    = c(0.1234567, 0.2, 0.31),
    pvalue   = c(1.2345e-50, 3.3e-12, 0.049912),
    padj     = c(4.42e-46, 7.7e-9, 0.0501),
    Count    = c(18L, 42L, 7L),
    stringsAsFactors = FALSE)
}

test_that("los p-valores no se redondean a decimales, van a cifras significativas", {
  df <- tabla_deg()
  f <- formato_por_columna(dt_table_num(df), df)

  expect_identical(unname(f[["pvalue"]]), "signif")
  expect_identical(unname(f[["padj"]]),   "signif")

  # La comprobacion que da sentido al test: con tres decimales, el padj de
  # CRISPLD2 seria cero.
  expect_identical(round(df$padj[1], 3), 0)
})

test_that("las magnitudes continuas si van a tres decimales", {
  df <- tabla_deg()
  f <- formato_por_columna(dt_table_num(df), df)

  for (nm in c("baseMean", "log2FC", "lfcSE")) {
    expect_identical(unname(f[[nm]]), "round", info = nm)
  }

  # Y con tres, no con los cuatro que habia antes.
  cd <- Filter(function(e) !is.null(e$render) && grepl("formatRound", e$render),
               dt_table_num(df)$x$options$columnDefs)
  expect_true(length(cd) > 0)
  expect_true(grepl("formatRound(data, 3,", as.character(cd[[1]]$render), fixed = TRUE))
})

test_that("ni el texto ni los enteros se tocan", {
  df <- tabla_deg()
  f <- formato_por_columna(dt_table_num(df), df)

  expect_identical(unname(f[["gene"]]), "ninguno")
  # Un recuento de genes escrito "18.000" se leeria como dieciocho mil.
  expect_identical(unname(f[["Count"]]), "ninguno")
})

test_that("el formato es de lectura: el dato conserva toda su precision", {
  df <- tabla_deg()
  dt <- dt_table_num(df)

  # Es lo que hace que ordenar por padj distinga 0,0499 de 0,0501, y lo que
  # permite que la descarga entregue el resultado y no la vista.
  expect_equal(dt$x$data$padj, df$padj)
  expect_equal(dt$x$data$log2FC, df$log2FC)

  # El formateador solo actua en 'display'.
  cd <- Filter(function(e) !is.null(e$render), dt$x$options$columnDefs)
  expect_true(all(vapply(cd, function(e)
    grepl("type !== 'display'", as.character(e$render), fixed = TRUE), logical(1))))
})

test_that("una tabla sin columnas numericas no rompe", {
  expect_silent(dt_table_num(message_df("Sin resultados DEG.")))
})


# ── El mismo criterio, en el informe descargable ────────────────────────────
#
# El informe HTML no usa DT: escribe las celdas a mano. Antes usaba
# `signif(v, 4)` para todo, asi que la tabla de la interfaz y la del informe
# mostraban el mismo gen con distinto numero de cifras.

test_that("fmt_celda_num da tres decimales a las magnitudes continuas", {
  expect_identical(fmt_celda_num(2.6314159, "log2FC"), "2.631")
  expect_identical(fmt_celda_num(-4.5, "log2FC"), "-4.500")
  expect_identical(fmt_celda_num(3379.9876543, "baseMean"), "3379.988")
  expect_identical(fmt_celda_num(0.1234567, "lfcSE"), "0.123")
})

test_that("fmt_celda_num conserva los p-valores pequeños", {
  # El caso que motiva todo: con tres decimales seria "0.000".
  expect_identical(fmt_celda_num(4.42e-46, "padj"), "4.42e-46")
  expect_identical(fmt_celda_num(1.2345e-50, "pvalue"), "1.23e-50")
  # Y uno cerca del umbral sigue siendo legible como tal.
  expect_identical(fmt_celda_num(0.049912, "padj"), "0.0499")
})

test_that("fmt_celda_num no inventa decimales en enteros ni toca el texto", {
  expect_identical(fmt_celda_num(18L, "Count"), "18")
  expect_identical(fmt_celda_num("CRISPLD2", "gene"), "CRISPLD2")
  expect_true(is.na(fmt_celda_num(NA_real_, "log2FC")))
})

test_that("la interfaz y el informe coinciden en que columnas son p-valores", {
  # Un solo diccionario: si se anade un motor con otro nombre de columna, se
  # anade en un sitio y las dos tablas lo respetan.
  expect_true(all(c("pvalue", "padj", "p.adjust", "qvalue") %in% COLS_P_VALOR))
  df <- tabla_deg()
  f <- formato_por_columna(dt_table_num(df), df)
  for (nm in intersect(names(df), COLS_P_VALOR)) {
    expect_identical(unname(f[[nm]]), "signif", info = nm)
    expect_true(grepl("e-", fmt_celda_num(1e-30, nm), fixed = TRUE), info = nm)
  }
})
