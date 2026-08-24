# Metricas de coste de la ejecucion (tiempo y memoria).
#
# El pico de memoria es el dato que responde a "cuanta RAM necesito"; ni el
# promedio ni el minimo dicen nada util. Estos tests fijan el formateo y la
# lectura, que es donde un cambio silencioso pasa desapercibido: un numero mal
# escalado sigue siendo un numero y nadie lo nota.

test_that("las duraciones se escriben en la unidad que toca", {
  expect_equal(fmt_duracion(0), "0s")
  expect_equal(fmt_duracion(45), "45s")
  expect_equal(fmt_duracion(90), "1m 30s")
  expect_equal(fmt_duracion(290), "4m 50s")
  expect_equal(fmt_duracion(3661), "1h 01m")
  # Ausencia y basura dan un guion, no un error ni un cero enganoso.
  expect_equal(fmt_duracion(NA), "—")
  expect_equal(fmt_duracion(NULL), "—")
  expect_equal(fmt_duracion("no soy un numero"), "—")
  expect_equal(fmt_duracion(-5), "—")
})

test_that("la memoria pasa a GB cuando el numero en MB deja de leerse", {
  expect_equal(fmt_memoria(90), "90 MB")
  expect_equal(fmt_memoria(1023), "1023 MB")
  expect_equal(fmt_memoria(1024), "1.0 GB")
  expect_equal(fmt_memoria(2444), "2.4 GB")
  expect_equal(fmt_memoria(NA), "—")
  expect_equal(fmt_memoria(0), "—")
})

test_that("read_run_metrics devuelve NULL sin fichero y la tabla cuando lo hay", {
  expect_null(read_run_metrics(tempfile()))

  d <- withr::local_tempdir()
  writeLines(c("paso\tsegundos\tduracion\tpico_rss_mb",
               "preparacion\t17\t17s\t90",
               "qc_inicial\t38\t38s\t2444",
               "TOTAL\t290\t4m 50s\t2444"),
             file.path(d, "metrics.tsv"))
  m <- read_run_metrics(d)
  expect_equal(nrow(m), 3L)
  expect_equal(names(m), c("paso", "segundos", "duracion", "pico_rss_mb"))
  expect_equal(m$pico_rss_mb[m$paso == "TOTAL"], 2444)
})

test_that("el coste llega al resumen de la ejecucion, tambien si fallo", {
  d <- withr::local_tempdir()
  # exit_status.tsv es lo que el workflow escribe en su trap EXIT: lleva el
  # coste tambien cuando la ejecucion termino mal, que es cuando hace falta
  # para diagnosticar por que.
  writeLines(c("exit_code\t1", "status\terror",
               "finished_at\t2026-08-24 13:32:00",
               "duration_seconds\t290", "peak_rss_mb\t2444"),
             file.path(d, "exit_status.tsv"))
  st <- read_exit_status(d)
  expect_equal(st$status, "error")
  expect_equal(fmt_duracion(st$duration_seconds), "4m 50s")
  expect_equal(fmt_memoria(st$peak_rss_mb), "2.4 GB")
})

test_that("una ejecucion antigua sin metricas no rompe nada", {
  d <- withr::local_tempdir()
  writeLines(c("exit_code\t0", "status\tsuccess",
               "finished_at\t2026-08-01 10:00:00"),
             file.path(d, "exit_status.tsv"))
  st <- read_exit_status(d)
  # Las claves nuevas no existen: el formateador tiene que dar un guion y no
  # un error, o la pestana de resultados se cae con cualquier ejecucion
  # anterior a este cambio.
  expect_equal(fmt_duracion(st$duration_seconds), "—")
  expect_equal(fmt_memoria(st$peak_rss_mb), "—")
  expect_null(read_run_metrics(d))
})
