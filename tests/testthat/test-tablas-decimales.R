#' test-tablas-decimales.R
#' Convencion numerica de las tablas: coma decimal, punto de millares y tres
#' decimales.
#'
#' Es la del documento del TFM y la que ya usan las tablas del harness, y las
#' capturas de la aplicacion van dentro de ese documento.
#'
#' El caso que obliga a que esto tenga test es el p-valor. `round(4.42e-46, 3)`
#' es 0, y ese 4,42e-46 es el p ajustado de CRISPLD2 en GSE52778, uno de los
#' genes de control que la memoria reporta. Redondearlo a decimales fijos no lo
#' redondea: lo borra, y lo deja indistinguible de un gen no significativo.
#'
#' La segunda propiedad que se fija aqui es que el formato es de LECTURA. El
#' dato que viaja al widget —y por tanto el que se ordena, se filtra y se
#' descarga— conserva toda su precision.

#' Que formateador le toca a cada columna, leido del propio widget.
#' DT 0.34 los guarda en `columnDefs[[i]]$render`, no en `rowCallback`.
#' Devuelve "round:N", "signif:N" o "ninguno".
formato_por_columna <- function(dt, df) {
  out <- setNames(rep("ninguno", length(df)), names(df))
  for (e in dt$x$options$columnDefs) {
    if (is.null(e$render)) next
    txt <- as.character(e$render)
    m <- regmatches(txt, regexec("DTWidget\\.format(Round|Signif)\\(data, ([0-9]+)", txt))[[1]]
    if (!length(m)) next
    etiqueta <- paste0(tolower(m[2]), ":", m[3])
    for (t in e$targets) out[[t + 1L]] <- etiqueta
  }
  out
}

#' Las marcas (millares, decimal) que el widget pasa al formateador.
marcas <- function(dt) {
  for (e in dt$x$options$columnDefs) {
    if (is.null(e$render)) next
    m <- regmatches(as.character(e$render),
                    regexec('data, [0-9]+, [0-9]+, "([^"]*)", "([^"]*)"',
                            as.character(e$render)))[[1]]
    if (length(m)) return(list(millares = m[2], decimal = m[3]))
  }
  NULL
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

test_that("las tablas usan coma decimal y punto de millares", {
  mk <- marcas(dt_table(tabla_deg()))
  expect_identical(mk$millares, ".")
  expect_identical(mk$decimal, ",")
})

test_that("los p-valores no se redondean a decimales, van a cifras significativas", {
  f <- formato_por_columna(dt_table(tabla_deg()), tabla_deg())
  expect_identical(unname(f[["pvalue"]]), "signif:3")
  expect_identical(unname(f[["padj"]]),   "signif:3")

  # La comprobacion que da sentido al test: con tres decimales, el padj de
  # CRISPLD2 seria cero.
  expect_identical(round(tabla_deg()$padj[1], 3), 0)
})

test_that("las magnitudes continuas van a tres decimales", {
  f <- formato_por_columna(dt_table(tabla_deg()), tabla_deg())
  for (nm in c("baseMean", "log2FC", "lfcSE")) {
    expect_identical(unname(f[[nm]]), "round:3", info = nm)
  }
})

test_that("los enteros llevan millares pero no decimales", {
  f <- formato_por_columna(dt_table(tabla_deg()), tabla_deg())
  # "18,000" para un recuento de 18 genes se leeria como dieciocho mil.
  expect_identical(unname(f[["Count"]]), "round:0")
})

test_that("el texto no se toca", {
  f <- formato_por_columna(dt_table(tabla_deg()), tabla_deg())
  expect_identical(unname(f[["gene"]]), "ninguno")
})

test_that("el formato es de lectura: el dato conserva toda su precision", {
  df <- tabla_deg()
  dt <- dt_table(df)

  # Es lo que hace que ordenar por padj distinga 0,0499 de 0,0501, y lo que
  # permite que la descarga entregue el resultado y no la vista.
  expect_equal(dt$x$data$padj, df$padj)
  expect_equal(dt$x$data$log2FC, df$log2FC)

  cd <- Filter(function(e) !is.null(e$render), dt$x$options$columnDefs)
  expect_true(all(vapply(cd, function(e)
    grepl("type !== 'display'", as.character(e$render), fixed = TRUE), logical(1))))
})

test_that("una tabla sin columnas numericas no rompe", {
  expect_silent(dt_table(message_df("Sin resultados DEG.")))
})

test_that("la convencion la aplica el envoltorio, asi que alcanza a toda la app", {
  # No solo a la pestana 4: cualquier tabla que pase por dt_table() la hereda,
  # que es el motivo de que el formato viva en el envoltorio y no en cada
  # llamada. Aqui, una tabla con la forma de las de calidad.
  qc <- data.frame(Sample = c("A", "B"),
                   `Total Sequences` = c(12345678L, 987654L),
                   percent_mapped = c(96.4312, 91.007),
                   check.names = FALSE, stringsAsFactors = FALSE)
  f <- formato_por_columna(dt_table(qc), qc)
  expect_identical(unname(f[["Total Sequences"]]), "round:0")
  expect_identical(unname(f[["percent_mapped"]]), "round:3")
  expect_identical(unname(f[["Sample"]]), "ninguno")
})


# ── El mismo criterio, en el informe descargable ────────────────────────────
#
# El informe HTML no usa DT: escribe las celdas a mano. Antes usaba
# `signif(v, 4)` para todo, asi que la tabla de la interfaz y la del informe
# mostraban el mismo gen con distinto numero de cifras y distinta convencion.

test_that("fmt_celda_num da tres decimales con coma a las magnitudes continuas", {
  expect_identical(fmt_celda_num(2.6314159, "log2FC"), "2,631")
  expect_identical(fmt_celda_num(-4.5, "log2FC"), "-4,500")
  expect_identical(fmt_celda_num(3379.9876543, "baseMean"), "3.379,988")
  expect_identical(fmt_celda_num(12345.6789, "baseMean"), "12.345,679")
})

test_that("fmt_celda_num conserva los p-valores pequeños", {
  # El caso que motiva todo: con tres decimales seria "0,000".
  expect_identical(fmt_celda_num(4.42e-46, "padj"), "4,42e-46")
  expect_identical(fmt_celda_num(1.2345e-50, "pvalue"), "1,23e-50")
  # Y uno cerca del umbral sigue siendo legible como tal.
  expect_identical(fmt_celda_num(0.049912, "padj"), "0,0499")
  expect_identical(fmt_celda_num(0.0501, "padj"), "0,0501")
})

test_that("fmt_celda_num pone millares a los enteros y no toca el texto", {
  expect_identical(fmt_celda_num(12345678, "Lecturas_asignadas"), "12.345.678")
  expect_identical(fmt_celda_num(18L, "Count"), "18")
  expect_identical(fmt_celda_num("CRISPLD2", "gene"), "CRISPLD2")
  expect_true(is.na(fmt_celda_num(NA_real_, "log2FC")))
})

test_that("formatear en castellano no dispara el aviso de big.mark ambiguo", {
  # `big.mark = "."` con el decimal.mark por defecto, que tambien es ".", hace
  # avisar a R. `fmt_int()` ya lo resolvia; el formateador de celdas lo reusa en
  # vez de repetir el arreglo.
  expect_silent(fmt_celda_num(12345678, "Lecturas_asignadas"))
  expect_silent(fmt_celda_num(12345.6789, "baseMean"))
})

test_that("la interfaz y el informe coinciden en que columnas son p-valores", {
  # Un solo diccionario: si se anade un motor con otro nombre de columna, se
  # anade en un sitio y las dos tablas lo respetan.
  expect_true(all(c("pvalue", "padj", "p.adjust", "qvalue") %in% COLS_P_VALOR))
  df <- tabla_deg()
  f <- formato_por_columna(dt_table(df), df)
  for (nm in intersect(names(df), COLS_P_VALOR)) {
    expect_identical(unname(f[[nm]]), "signif:3", info = nm)
    expect_true(grepl("e-", fmt_celda_num(1e-30, nm), fixed = TRUE), info = nm)
  }
})


test_that("el diccionario reconoce los p-valores con sufijo", {
  # No es teorico: la tabla de comparacion ORA/GSEA nombra sus columnas
  # `padj_ORA` y `padj_GSEA`. Con la lista de nombres exactos quedaban fuera y
  # en la interfaz se mostraban como "0,000", que es justo el fallo que el
  # diccionario existe para evitar. Se vio ejecutando la aplicacion.
  expect_true(es_col_p_valor("padj_GSEA"))
  expect_true(es_col_p_valor("padj_ORA"))
  expect_identical(fmt_celda_num(1.55e-05, "padj_GSEA"), "1,55e-05")

  # Y no se lleva por delante columnas que solo se le parecen.
  expect_false(es_col_p_valor("NES"))
  expect_false(es_col_p_valor("setSize"))
  expect_false(es_col_p_valor("enrichmentScore"))
  expect_identical(fmt_celda_num(3.409, "NES"), "3,409")

  df <- data.frame(ID = "R-HSA-445355", padj_ORA = 0.001, padj_GSEA = 1.55e-05,
                   NES = 3.409, stringsAsFactors = FALSE)
  f <- formato_por_columna(dt_table(df), df)
  expect_identical(unname(f[["padj_GSEA"]]), "signif:3")
  expect_identical(unname(f[["NES"]]), "round:3")
})
