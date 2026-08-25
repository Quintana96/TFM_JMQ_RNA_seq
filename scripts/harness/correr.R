#!/usr/bin/env Rscript
#' scripts/harness/correr.R
#' Punto de entrada del harness de validación.
#'
#'   Rscript scripts/harness/correr.R GSE273773
#'   Rscript scripts/harness/correr.R --todos
#'   Rscript scripts/harness/correr.R GSE273773 --etapas deg,metricas
#'   Rscript scripts/harness/correr.R GSE273773 --forzar
#'   Rscript scripts/harness/correr.R --tablas
#'
#' Etapas disponibles: pipeline, deg, enriquecimiento, metricas.
#' Sin `--etapas` se ejecutan todas, en ese orden.
#'
#' Cada etapa se salta si ya está hecha, de modo que una ejecución interrumpida
#' se reanuda sin repetir lo caro. `--forzar` las rehace.

setwd(local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  d <- if (length(f)) dirname(dirname(dirname(normalizePath(sub("^--file=", "", f[1]))))) else getwd()
  if (file.exists(file.path(d, "global.R"))) d else getwd()
}))

args <- commandArgs(trailingOnly = TRUE)
tomar <- function(bandera, defecto = NULL) {
  i <- match(bandera, args)
  if (is.na(i) || i == length(args)) defecto else args[i + 1]
}
forzar  <- "--forzar" %in% args
todos   <- "--todos" %in% args
solo_tablas <- "--tablas" %in% args
hilos   <- as.integer(tomar("--hilos", "6"))
etapas  <- strsplit(tomar("--etapas", "pipeline,deg,enriquecimiento,metricas"), ",")[[1]]
ids     <- setdiff(args[!grepl("^--", args)], c(tomar("--etapas"), tomar("--hilos")))

source("scripts/harness/datasets.R")
source("scripts/harness/harness.R")
source("scripts/harness/tablas.R")

if (solo_tablas) {
  msg("Generando las tablas de defensa")
  generar_tablas(if (length(ids)) ids else datasets_disponibles())
  quit(save = "no", status = 0)
}

if (todos) ids <- datasets_disponibles()
if (!length(ids)) {
  cat("Uso: Rscript scripts/harness/correr.R <dataset> [--etapas ...] [--forzar] [--hilos N]\n")
  cat("     Rscript scripts/harness/correr.R --todos\n")
  cat("     Rscript scripts/harness/correr.R --tablas\n\n")
  cat("Datasets configurados:", paste(datasets_disponibles(), collapse = ", "), "\n")
  quit(save = "no", status = 1)
}

msg("Cargando SARA")
cargar_app()

for (id in ids) {
  d <- dataset_config(id)
  cat("\n", strrep("=", 74), "\n", d$id, " — ", d$organismo, "\n",
      d$descripcion, "\n", strrep("=", 74), "\n", sep = "")

  if ("pipeline" %in% etapas && length(d$pipelines)) {
    msg("Etapa 1/4: pipeline (", paste(d$pipelines, collapse = ", "), ")")
    for (pl in d$pipelines) etapa_pipeline(d, pl, hilos = hilos, forzar = forzar)
  }

  degs <- list()
  if ("deg" %in% etapas) {
    msg("Etapa 2/4: expresión diferencial (", paste(d$motores, collapse = ", "), ")")
    degs <- etapa_deg(d, forzar = forzar)
    msg("  ", length(degs), " ajustes disponibles")
  } else {
    f <- list.files(file.path(dir_dataset(d), "deg"), pattern = "\\.rds$", full.names = TRUE)
    degs <- stats::setNames(lapply(f, readRDS), sub("\\.rds$", "", basename(f)))
  }

  if ("enriquecimiento" %in% etapas) {
    msg("Etapa 3/4: enriquecimiento funcional")
    etapa_enriquecimiento(d, degs, forzar = forzar)
  }

  if ("metricas" %in% etapas) {
    msg("Etapa 4/4: métricas de validación")
    salida <- dir_etapa(d, "metricas")

    # Contra lo publicado, un resultado por entrada
    vs <- list()
    for (nm in names(degs)) {
      m <- metricas_vs_publicado(d, degs[[nm]])
      if (!is.null(m)) vs[[nm]] <- m
    }
    if (length(vs)) {
      saveRDS(vs, file.path(salida, "vs_publicado.rds"))
      msg("  contra lo publicado: ", length(vs), " comparaciones")
    } else msg("  sin tabla publicada con la que comparar")

    # Entre motores y entre rutas
    pares <- list()
    nms <- names(degs)
    for (i in seq_along(nms)) for (j in seq_along(nms)) if (i < j) {
      pi <- strsplit(nms[i], "__")[[1]]; pj <- strsplit(nms[j], "__")[[1]]
      # Solo se comparan pares que aíslan UNA variable: o el motor con la misma
      # ruta, o la ruta con el mismo motor. Comparar los dos a la vez no dice
      # de qué viene la diferencia.
      tipo <- if (pi[1] == pj[1]) "motor" else if (pi[2] == pj[2]) "ruta" else NA
      if (is.na(tipo)) next
      m <- metricas_entre(degs[[nms[i]]], degs[[nms[j]]])
      if (!is.null(m)) pares[[paste0(tipo, ": ", nms[i], " vs ", nms[j])]] <- m
    }
    if (length(pares)) {
      saveRDS(pares, file.path(salida, "pares.rds"))
      msg("  acuerdo interno: ", length(pares), " parejas")
    }

    etapa_permutacion(d, forzar = forzar)
    des <- etapa_descomposicion(d, forzar = forzar)
    if (!is.null(des)) msg("  descomposicion: Pearson ", round(des$pearson, 4),
                           " | recuperados ", round(100*des$recuperados, 1), " %")
  }
}

msg("Generando las tablas de defensa")
generar_tablas(ids)
msg("Hecho. Resultados en ", RAIZ_VALIDACION, "/")
