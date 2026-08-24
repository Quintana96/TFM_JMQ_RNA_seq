# Comprobacion de las herramientas del pipeline (R/utils_status.R).
#
# El workflow ya las verifica en su primer paso y aborta si falta alguna, pero
# para entonces la ejecucion ya esta lanzada: arriba se lee "Error (codigo 1)" y
# la causa aparece veinte lineas mas abajo en el log. Comprobarlo en la
# interfaz convierte eso en una frase antes de empezar.

test_that("cada estrategia exige sus herramientas y ninguna de mas", {
  comunes <- c("fastqc", "fastp", "multiqc")

  aln <- herramientas_requeridas("alignment", "bowtie2")
  expect_true(all(comunes %in% aln))
  expect_true(all(c("bowtie2", "samtools", "featureCounts") %in% aln))
  # Pedir un cuantificador que no se va a usar bloquearia el analisis por una
  # herramienta que no hace falta.
  expect_false(any(c("salmon", "kallisto") %in% aln))

  sal <- herramientas_requeridas("pseudo", "salmon")
  expect_true("salmon" %in% sal)
  expect_false("kallisto" %in% sal)
  expect_false(any(c("bowtie2", "samtools", "featureCounts") %in% sal))

  kal <- herramientas_requeridas("pseudo", "kallisto")
  expect_true("kallisto" %in% kal)
  expect_false("salmon" %in% kal)
})

test_that("con todo en el PATH no se reporta nada", {
  # Herramientas que existen en cualquier sistema, para no depender de que el
  # entorno del pipeline este activo al correr los tests.
  e <- comprobar_herramientas(necesarias = c("ls", "cat"))
  expect_length(e$faltan, 0)
  expect_null(mensaje_herramientas(e))
})

test_that("se nombra lo que falta, no un 'faltan herramientas' generico", {
  e <- comprobar_herramientas(entornos = character(0),
                              necesarias = c("ls", "no_existe_esta_herramienta"))
  expect_equal(e$faltan, "no_existe_esta_herramienta")
  expect_null(e$entorno)
  msg <- mensaje_herramientas(e)
  expect_match(msg, "no_existe_esta_herramienta", fixed = TRUE)
  # Sin entorno donde esten, el consejo es instalarlas.
  expect_match(msg, "requirements.sh", fixed = TRUE)
})

test_that("si estan instaladas pero fuera del PATH, se dice donde", {
  # Es el caso real: la aplicacion arrancada sin activar el entorno de conda.
  # Decir "no encontradas" a secas manda a instalar algo que ya esta instalado.
  d <- withr::local_tempdir()
  bin <- file.path(d, "bin"); dir.create(bin)
  file.create(file.path(bin, "no_existe_esta_herramienta"))

  e <- comprobar_herramientas(entornos = d,
                              necesarias = "no_existe_esta_herramienta")
  expect_equal(e$entorno, d)
  msg <- mensaje_herramientas(e)
  expect_match(msg, d, fixed = TRUE)
  expect_match(msg, "lanzar_app.sh", fixed = TRUE)
})

test_that("un entorno que solo resuelve la mitad no se ofrece", {
  # Mandar al usuario a un entorno que arregla parte del problema le hace
  # arrancar otra vez para volver a fallar, ahora por otra herramienta.
  d <- withr::local_tempdir()
  bin <- file.path(d, "bin"); dir.create(bin)
  file.create(file.path(bin, "falta_una"))

  e <- comprobar_herramientas(entornos = d,
                              necesarias = c("falta_una", "falta_otra"))
  expect_setequal(e$faltan, c("falta_una", "falta_otra"))
  expect_null(e$entorno)
})

test_that("solo se consideran entornos con carpeta bin", {
  d <- withr::local_tempdir()
  dir.create(file.path(d, "sin_bin"))
  con_bin <- file.path(d, "con_bin"); dir.create(file.path(con_bin, "bin"), recursive = TRUE)
  encontrados <- entornos_conda_probables()
  # No se puede fijar el contenido real de la maquina, pero si la invariante.
  expect_true(all(dir.exists(file.path(encontrados, "bin"))))
})
