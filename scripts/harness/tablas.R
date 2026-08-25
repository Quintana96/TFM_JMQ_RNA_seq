#' scripts/harness/tablas.R
#' Genera las tablas para defender el TFM.
#'
#' Cada tabla responde a una pregunta que un tribunal puede hacer, y está
#' pensada para copiarse tal cual a la memoria. Se emiten en Markdown (para
#' pegar) y en TSV (para reprocesar).
#'
#' Las siete preguntas:
#'   T1  ¿Con qué datos has trabajado?
#'   T2  ¿Qué coste tiene ejecutar esto?
#'   T3  ¿Coincides con lo publicado?
#'   T4  ¿Coinciden los motores entre sí?
#'   T5  ¿Coinciden las rutas del pipeline entre sí?
#'   T6  ¿Cómo sabes que no estás inventando señal?
#'   T7  ¿De dónde viene el desacuerdo que queda?
#'   T8  ¿Con qué versiones exactas?

fmt <- function(x, d = 4) if (is.null(x) || length(x) == 0 || is.na(x)) "—" else formatC(x, format = "f", digits = d)
pct <- function(x, d = 1) if (is.null(x) || length(x) == 0 || is.na(x)) "—" else paste0(formatC(100*x, format = "f", digits = d), " %")
ent <- function(x) if (is.null(x) || length(x) == 0 || is.na(x)) "—" else format(x, big.mark = ".", decimal.mark = ",")

#' Escribe una tabla en Markdown y en TSV.
escribir <- function(df, id, titulo, nota = NULL, dir_salida = RAIZ_VALIDACION) {
  dir.create(file.path(dir_salida, "tablas"), recursive = TRUE, showWarnings = FALSE)
  base <- file.path(dir_salida, "tablas", id)
  utils::write.table(df, paste0(base, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  md <- c(paste0("### ", id, ". ", titulo), "",
          paste0("| ", paste(names(df), collapse = " | "), " |"),
          paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
          apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
  if (!is.null(nota)) md <- c(md, "", paste0("> ", nota))
  writeLines(md, paste0(base, ".md"))
  invisible(df)
}

leer_rds <- function(d, etapa, f) {
  p <- file.path(RAIZ_VALIDACION, d$id, etapa, f)
  if (file.exists(p)) readRDS(p) else NULL
}

generar_tablas <- function(ids) {
  cfg <- lapply(ids, function(i) tryCatch(dataset_config(i), error = function(e) NULL))
  cfg <- Filter(Negate(is.null), cfg)
  if (!length(cfg)) { message("Sin datasets que tabular"); return(invisible(NULL)) }

  # ── T1: los datos ─────────────────────────────────────────────────────────
  t1 <- do.call(rbind, lapply(cfg, function(d) data.frame(
    Dataset = d$id, Organismo = d$organismo, Comparación = d$descripcion,
    n = length(d$muestras),
    Diseño = if (!is.null(d$batch)) paste0("pareado por ", d$batch) else "independiente",
    Lectura = if (identical(d$read_type, "se")) "single-end" else "paired-end",
    Rutas = if (length(d$pipelines)) paste(d$pipelines, collapse = ", ") else "solo conteos",
    Artículo = d$articulo$cita %||% "—", stringsAsFactors = FALSE)))
  escribir(t1, "T1", "Conjuntos de datos utilizados",
    "Las rutas de alineamiento (bowtie2, subjunc) solo se aplican donde son legítimas: bowtie2 no reconoce uniones exón-exón, así que en eucariotas con intrones se usa subjunc o pseudoalineamiento.")

  # ── T2: coste de ejecución ────────────────────────────────────────────────
  f2 <- list()
  for (d in cfg) for (pl in d$pipelines) {
    dd <- dir_pipeline(d, pl)
    st <- if (dir.exists(dd)) read_exit_status(dd) else NULL
    mt <- if (dir.exists(dd)) read_run_metrics(dd) else NULL
    if (is.null(st)) next
    cm <- tryCatch(matriz_de(d, pl), error = function(e) NULL)
    paso_pico <- if (!is.null(mt) && nrow(mt) > 1) {
      m2 <- mt[mt$paso != "TOTAL", ]; m2$paso[which.max(m2$pico_rss_mb)] } else "—"
    f2[[length(f2)+1]] <- data.frame(
      Dataset = d$id, Ruta = pl, Estado = st$status %||% "?",
      Duración = fmt_duracion(st$duration_seconds),
      `RAM máxima` = fmt_memoria(st$peak_rss_mb),
      `Paso que más pide` = paso_pico,
      Genes = if (is.null(cm)) "—" else ent(nrow(cm)),
      check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (length(f2)) escribir(do.call(rbind, f2), "T2", "Coste de ejecución del pipeline",
    "La memoria es el máximo del árbol de procesos, muestreado una vez por segundo: un pico más corto puede no quedar registrado.")

  # ── T3: concordancia con lo publicado ─────────────────────────────────────
  f3 <- list()
  for (d in cfg) {
    vs <- leer_rds(d, "metricas", "vs_publicado.rds"); if (is.null(vs)) next
    for (nm in names(vs)) { v <- vs[[nm]]; p <- strsplit(nm, "__")[[1]]
      f3[[length(f3)+1]] <- data.frame(
        Dataset = d$id, Ruta = p[1], Motor = p[2],
        `Genes comparados` = ent(v$genes_comunes),
        `DEG propios` = ent(v$deg_propios), `DEG publicados` = ent(v$deg_publicados),
        Recuperados = pct(v$recuperados), Jaccard = fmt(v$jaccard, 3),
        `Pearson log2FC` = fmt(v$pearson), `Mismo signo` = pct(v$mismo_signo),
        `|dif| mediana` = fmt(v$dif_mediana, 3),
        check.names = FALSE, stringsAsFactors = FALSE) }
  }
  if (length(f3)) escribir(do.call(rbind, f3), "T3", "Concordancia con los resultados publicados",
    "La correlación se calcula sobre TODOS los genes comunes y no solo sobre los significativos: restringirla a los diferenciales la infla por selección.")

  # ── T4 y T5: acuerdo interno ──────────────────────────────────────────────
  f4 <- list(); f5 <- list()
  for (d in cfg) {
    pares <- leer_rds(d, "metricas", "pares.rds"); if (is.null(pares)) next
    for (nm in names(pares)) { v <- pares[[nm]]
      tipo <- sub(":.*", "", nm); etiqueta <- trimws(sub("^[^:]*:", "", nm))
      fila <- data.frame(Dataset = d$id, Comparación = etiqueta,
        `Genes comunes` = ent(v$genes_comunes),
        `DEG A` = ent(v$deg_A), `DEG B` = ent(v$deg_B),
        Coincidentes = ent(v$coincidentes), Jaccard = fmt(v$jaccard, 3),
        `Pearson log2FC` = fmt(v$pearson), `Mismo signo` = pct(v$mismo_signo),
        check.names = FALSE, stringsAsFactors = FALSE)
      if (tipo == "motor") f4[[length(f4)+1]] <- fila else f5[[length(f5)+1]] <- fila }
  }
  if (length(f4)) escribir(do.call(rbind, f4), "T4", "Acuerdo entre los motores estadísticos",
    "Sobre la MISMA matriz de conteos, de modo que la única variable es el motor. Es la única tabla que valida la aplicación y no el trabajo ajeno.")
  if (length(f5)) escribir(do.call(rbind, f5), "T5", "Acuerdo entre las rutas del pipeline",
    "Con el MISMO motor, de modo que la única variable es cómo se asignan las lecturas.")

  # ── T6: control de permutación ────────────────────────────────────────────
  f6 <- list()
  for (d in cfg) {
    pm <- leer_rds(d, "metricas", "permutacion.rds"); if (is.null(pm)) next
    f6[[length(f6)+1]] <- data.frame(Dataset = d$id,
      `DEG con etiquetas reales` = ent(pm$deg_real),
      `Reetiquetados probados` = ent(pm$n),
      `DEG mediana` = ent(pm$mediana), `DEG máximo` = ent(pm$maximo),
      `Reetiquetados con 0 DEG` = paste0(pm$ceros, " de ", pm$n),
      check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (length(f6)) escribir(do.call(rbind, f6), "T6", "Control de permutación",
    "Barajando las etiquetas de grupo. Las demás tablas miden concordancia con lo publicado; esta mide que el pipeline no invente señal, que es lo que ninguna otra demuestra.")

  # ── T7: descomposición del error ──────────────────────────────────────────
  f7 <- list()
  for (d in cfg) {
    de <- leer_rds(d, "metricas", "descomposicion.rds")
    vs <- leer_rds(d, "metricas", "vs_publicado.rds")
    if (is.null(de) || is.null(vs)) next
    propio <- vs[[grep("DESeq2$", names(vs))[1]]]
    f7[[length(f7)+1]] <- data.frame(Dataset = d$id,
      `Entrada` = "conteos de los autores", `Pearson log2FC` = fmt(de$pearson),
      Recuperados = pct(de$recuperados), `|dif| mediana` = fmt(de$dif_mediana, 3),
      check.names = FALSE, stringsAsFactors = FALSE)
    if (!is.null(propio)) f7[[length(f7)+1]] <- data.frame(Dataset = d$id,
      `Entrada` = "conteos propios", `Pearson log2FC` = fmt(propio$pearson),
      Recuperados = pct(propio$recuperados), `|dif| mediana` = fmt(propio$dif_mediana, 3),
      check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (length(f7)) escribir(do.call(rbind, f7), "T7", "Descomposición del error",
    "El mismo análisis sobre los conteos de los autores frente a los propios. Si la primera fila es casi perfecta, el desacuerdo de la segunda es atribuible al alineamiento y al conteo, no a la estadística.")

  # ── T8: procedencia ───────────────────────────────────────────────────────
  f8 <- list()
  for (d in cfg) for (pl in d$pipelines) {
    dd <- dir_pipeline(d, pl)
    tv <- if (dir.exists(dd)) read_tool_versions(dd) else NULL
    if (is.null(tv) || !nrow(tv)) next
    inst <- tv[!grepl("no instalado", tv$version), , drop = FALSE]
    f8[[length(f8)+1]] <- data.frame(Dataset = d$id, Ruta = pl,
      Herramientas = paste(paste0(inst$tool, " ", sub("^[^0-9]*", "", inst$version)),
                           collapse = "; "),
      check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (length(f8)) escribir(do.call(rbind, f8), "T8", "Versiones de las herramientas",
    "Leídas de versions.tsv, que el propio workflow escribe invocando cada herramienta: no son las declaradas, son las que realmente se ejecutaron.")

  # Índice
  tablas <- sort(list.files(file.path(RAIZ_VALIDACION, "tablas"), pattern = "\\.md$"))
  writeLines(c("# Tablas para la defensa del TFM", "",
               paste0("Generadas el ", format(Sys.time(), "%d de %B de %Y a las %H:%M"), "."),
               "", "Cada tabla está en Markdown, para pegar en la memoria, y en TSV.", "",
               unlist(lapply(tablas, function(f)
                 readLines(file.path(RAIZ_VALIDACION, "tablas", f))))),
             file.path(RAIZ_VALIDACION, "tablas", "TODAS.md"))
  message("Tablas en ", file.path(RAIZ_VALIDACION, "tablas"), ": ",
          paste(sub("\\.md$", "", tablas), collapse = ", "))
  invisible(TRUE)
}
