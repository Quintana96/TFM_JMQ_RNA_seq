#' utils_report.R
#' Informe reproducible y exportación del script R equivalente.
#'
#' Por qué existe (docs/REVISION_ESTADISTICA.md, B11): una aplicación gráfica
#' agrava el problema de trazabilidad que Pall et al. asocian al sesgo del campo,
#' porque al no haber script no queda registro de que se hizo. Dos salidas de
#' bajo coste y alto valor: un informe con todos los parámetros y los
#' diagnósticos, y el script R equivalente al análisis ejecutado — que convierte
#' la app en herramienta de aprendizaje y permite auditar el resultado fuera de
#' ella.
#'
#' El informe se genera como HTML autocontenido escrito a mano, NO via
#' rmarkdown::render(): pandoc no está disponible en todos los entornos donde
#' corre la app (no lo está en esta maquina), y una descarga que falla la mitad
#' de las veces no sirve. El script R exportado si es texto plano, sin
#' dependencias.

#' Escapa texto para insertarlo en HTML.
html_escape <- function(x) {
  x <- as.character(x %||% "")
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub('"', "&quot;", x, fixed = TRUE)
}

#' Tabla HTML de dos columnas a partir de una lista clave-valor.
html_kv_table <- function(x) {
  if (!length(x)) return("")
  rows <- vapply(names(x), function(k) {
    paste0("<tr><th>", html_escape(k), "</th><td>",
           html_escape(paste(x[[k]], collapse = ", ")), "</td></tr>")
  }, character(1))
  paste0("<table class=\"kv\">", paste(rows, collapse = ""), "</table>")
}

#' Parámetros del análisis DEG en forma de lista ordenada, para el informe y el
#' script. Es la fuente única de verdad de "que se ejecuto".
deg_run_parameters <- function(rv) {
  if (is.null(rv) || is.null(rv$results)) return(NULL)
  pf <- rv$prefilter
  list(
    "Fecha del análisis"       = format(rv$run_at %||% Sys.time(), "%Y-%m-%d %H:%M:%S"),
    # El ajuste y la lectura pueden ser de momentos distintos: el FDR y el
    # umbral del test se reextraen del mismo ajuste sin repetirlo. Cuando ha
    # pasado, la fecha del análisis por si sola no describe la tabla que se está
    # informando, así que se declaran las dos.
    "Última reextraccion"      = if (is.null(rv$reextracted_at)) "—"
                                 else format(rv$reextracted_at, "%Y-%m-%d %H:%M:%S"),
    "Motor"                    = rv$method %||% "—",
    "Diseño"                   = rv$design %||% "~ condition",
    "Contraste"                = rv$contrast %||% "—",
    "Coeficiente testeado"     = rv$coef %||% "—",
    "FDR objetivo (alpha)"     = rv$fdr %||% 0.05,
    "Corrección múltiple"      = rv$padj_method %||% "BH",
    "Umbral |log2FC| del test" = rv$lfc_threshold %||% 0,
    "Encogido de log2FC"       = rv$shrink %||% "ninguno",
    # El modo de outliers cambia que genes tienen padj, así que es un parámetro
    # del test y no una preferencia de visualización.
    "Outliers de Cook"         = switch(rv$outliers %||% "na",
      na = "excluidos del test (por defecto)",
      refit = "valor atípico sustituido y regenado al test",
      keep = "filtro desactivado (se tratan como biologia real)",
      rv$outliers %||% "—"),
    "Prefiltrado"              = if (is.null(pf)) "—" else
      paste0(pf$mode, ": ", pf$n_before, " -> ", pf$n_after, " genes"),
    "Muestras"                 = if (!is.null(rv$meta)) nrow(rv$meta) else NA,
    "Genes analizados"         = if (!is.null(rv$results)) nrow(rv$results) else NA,
    "Corrección solo gráfica"  = rv$viz_note %||% "ninguna"
  )
}

#' Identificación de la versión de la aplicación que genero el informe.
#'
#' Un informe que no dice con que versión del software se produjo no es
#' auditable: es el mismo requisito que la guia ACMG impone a un informe de
#' diagnóstico genomico ("pipeline y guias usadas con versiones").
app_provenance <- function() {
  commit <- tryCatch(
    suppressWarnings(system2("git", c("rev-parse", "--short", "HEAD"),
                             stdout = TRUE, stderr = FALSE)),
    error = function(e) character(0))
  sucio <- tryCatch(
    length(suppressWarnings(system2("git", c("status", "--porcelain"),
                                    stdout = TRUE, stderr = FALSE))) > 0,
    error = function(e) FALSE)
  list(
    "Aplicación" = "SARA (Shiny App for RNA-seq Analysis)",
    "Commit de git" = if (length(commit))
      paste0(commit[1], if (isTRUE(sucio)) " (con cambios sin commitear)" else "")
      else "no disponible",
    "Versión de R" = as.character(getRversion()),
    "Plataforma" = R.version$platform
  )
}

#' Resumen del análisis en lenguaje llano.
#'
#' La guia ACMG exige que un informe incluya un parrafo dirigido a alguien sin
#' formacion especializada. Aquí cumple además una función práctica: obliga a
#' que el informe diga en una frase que se comparo y que salio, sin que haya que
#' interpretar una tabla.
deg_plain_summary <- function(rv, sig_n, up, down) {
  fdr <- rv$fdr %||% 0.05
  ct <- rv$contrast %||% "las dos condiciones comparadas"
  n_tot <- nrow(rv$results)
  base <- paste0(
    "Se han comparado los perfiles de expresión de ", ct,
    " sobre ", fmt_int(n_tot), " genes analizados. ")
  if (sig_n == 0) {
    return(paste0(base,
      "No se ha encontrado ningun gen con diferencias estadisticamente ",
      "significativas al nivel exigido (FDR <= ", fdr, "). Esto NO significa ",
      "que no existan diferencias: puede que el experimento no tenga suficientes ",
      "muestras para detectarlas. Revisa el apartado de limitaciones."))
  }
  paste0(base,
    fmt_int(sig_n), " genes muestran diferencias estadisticamente significativas ",
    "(FDR <= ", fdr, "): ", fmt_int(up), " con mayor expresión y ", fmt_int(down),
    " con menor expresión en el primer grupo del contraste respecto al segundo. ",
    "La cifra de FDR indica la proporción esperada de falsos positivos dentro de ",
    "esa lista; con FDR = ", fdr, ", aproximadamente ",
    fmt_int(round(sig_n * fdr)), " de los ", fmt_int(sig_n),
    " genes podrían serlo por azar.")
}

#' Limitaciones detectadas en ESTE análisis concreto.
#'
#' No es un texto fijo: se deriva del estado real del ajuste. Un informe que no
#' declara sus limitaciones invita a sobreinterpretarlo, y la guia ACMG las exige
#' explicitamente, también (sobre todo) cuando el resultado es negativo.
deg_report_limitations <- function(rv, diagnostics = NULL, sig_n = NA_integer_) {
  lim <- character(0)
  n_meta <- if (!is.null(rv$meta)) nrow(rv$meta) else NA_integer_
  grupos <- if (!is.null(rv$meta) && "condition" %in% names(rv$meta))
    table(as.character(rv$meta$condition)) else NULL
  min_grp <- if (!is.null(grupos) && length(grupos)) min(grupos) else NA_integer_

  if (!is.na(min_grp) && min_grp < 6) {
    lim <- c(lim, paste0(
      "El grupo más pequeño tiene ", min_grp, " replicas. Schurch et al. (RNA 2016) ",
      "recomiendan al menos 6 por condición para una detección robusta, y 12 para ",
      "capturar la mayoria de los genes diferenciales. Con menos replicas la ",
      "potencia es limitada y la lista será menos reproducible."))
  }
  if (!is.null(rv$counts_source) && !isTRUE(rv$counts_source$ok)) {
    lim <- c(lim, paste0(
      "La matriz de conteos NO se resumio a gen por la via recomendada: ",
      rv$counts_source$detail %||% "", " Los resultados pueden estar sesgados en ",
      "genes con varias isoformas."))
  }
  if (!is.null(rv$viz_note) && nzchar(rv$viz_note %||% "")) {
    lim <- c(lim, paste0(
      "Los gráficos llevan una corrección que el test NO lleva: ", rv$viz_note,
      " Es lo correcto (el modelo trata la variación como covariable), pero ",
      "significa que gráficos y p-valores describen matrices distintas."))
  }
  e <- rv$enrich
  if (!is.null(e) && !is.null(e$mapeo) && !is.na(e$mapeo$rate %||% NA) &&
      e$mapeo$rate < 0.5) {
    lim <- c(lim, paste0(
      "En el enriquecimiento solo se mapeo el ", round(100 * e$mapeo$rate, 1),
      " % de los genes contra la anotación. Por debajo de la mitad, el resultado ",
      "funcional no es interpretable."))
  }
  if (identical(rv$method, "Swish")) {
    lim <- c(lim, paste0(
      "El análisis se ha hecho a nivel de TRANSCRITO, no de gen: los ",
      "identificadores no son comparables con los del resto de motores."))
  }
  if (is.null(diagnostics) || is.null(diagnostics$replicability)) {
    lim <- c(lim, paste0(
      "No se ha estimado la replicabilidad por remuestreo. Sin ella no hay ",
      "medida de cuanto dependeria esta lista de las muestras concretas ",
      "analizadas."))
  }
  if (!is.na(sig_n) && sig_n == 0) {
    lim <- c(lim, paste0(
      "Un resultado sin genes significativos no demuestra ausencia de efecto. ",
      "Antes de concluir nada, comprueba la potencia del diseño y el histograma ",
      "de p-valores: si es plano, o no hay señal o el modelo no la captura."))
  }
  lim
}

#' Informe HTML autocontenido del análisis.
#'
#' Estructura inspirada en el contenido mínimo que la guia ACMG/AMP (Richards et
#' al., Genetics in Medicine 2015) exige a un informe de diagnóstico genomico:
#' metricas de calidad, que se evaluo, bases de datos y pipeline CON SUS
#' VERSIONES, hallazgos con su metodología, fecha, resumen en lenguaje llano y
#' limitaciones. Un informe negativo se trata igual que uno positivo.
#'
#' @param rv `state$deg_rv`
#' @param diagnostics list(pi0, verdict, na_breakdown, cooks_dominant, replicability)
build_deg_report_html <- function(rv, diagnostics = NULL) {
  params <- deg_run_parameters(rv)
  if (is.null(params)) return(NULL)

  sig_n <- sum(!is.na(rv$results$padj) & rv$results$padj <= (rv$fdr %||% 0.05))
  up <- sum(!is.na(rv$results$padj) & rv$results$padj <= (rv$fdr %||% 0.05) &
              !is.na(rv$results$log2FC) & rv$results$log2FC > 0)
  down <- sig_n - up

  # ── Procedencia de los datos ─────────────────────────────────────────────
  org <- rv$counts_origin
  cs  <- rv$counts_source
  proc <- list()
  if (!is.null(org)) {
    proc[["Origen de la matriz"]] <- org$tipo %||% "—"
    proc[["Ruta / fichero"]]      <- org$ruta %||% "—"
    if (!is.null(org$md5) && !is.na(org$md5)) proc[["md5 del fichero"]] <- org$md5
    if (!is.null(org$detalle)) proc[["Detalle"]] <- org$detalle
  }
  proc[["Resumen a gen"]] <- if (is.null(cs)) "—" else
    paste0(cs$method %||% "—", if (isTRUE(cs$ok)) "" else "  [DEGRADADO]")
  if (!is.null(cs) && !isTRUE(cs$ok) && nzchar(cs$detail %||% ""))
    proc[["Motivo de la degradación"]] <- cs$detail
  proc_html <- paste0("<h2>Procedencia de los datos</h2>", html_kv_table(proc))

  # Versiones de las herramientas del pipeline y checksums de las entradas,
  # leidos de la ejecución de origen: cierra el ciclo entre lo que hizo el
  # pipeline y lo que hizo el análisis estadístico.
  rd <- rv$run_dir
  tv <- if (!is.null(rd) && nzchar(rd %||% "") && dir.exists(rd)) read_tool_versions(rd) else NULL
  ck <- if (!is.null(rd) && nzchar(rd %||% "") && dir.exists(rd)) read_input_checksums(rd) else NULL
  pipeline_html <- if (!is.null(tv) && nrow(tv)) {
    inst <- tv[!grepl("no instalado", tv$version), , drop = FALSE]
    paste0("<h3>Herramientas del pipeline</h3>",
           html_kv_table(stats::setNames(as.list(inst$version), inst$tool)))
  } else ""
  # Coste de la ejecución del pipeline. Es parte de la procedencia: dice en que
  # equipo y en cuanto tiempo se obtuvo la matriz, que es lo que permite juzgar
  # si el análisis es repetible con los medios de quien lo lee.
  mt <- if (!is.null(rd) && nzchar(rd %||% "") && dir.exists(rd)) read_run_metrics(rd) else NULL
  metricas_html <- if (!is.null(mt) && nrow(mt)) {
    filas <- vapply(seq_len(nrow(mt)), function(i) paste0(
      "<tr><td>", html_escape(mt$paso[i]), "</td><td>",
      html_escape(mt$duracion[i]), "</td><td>",
      html_escape(fmt_memoria(mt$pico_rss_mb[i])), "</td></tr>"), character(1))
    paste0("<h3>Coste de la ejecución del pipeline</h3>",
           "<table><thead><tr><th>Paso</th><th>Duracion</th>",
           "<th>Memoria máxima</th></tr></thead><tbody>",
           paste(filas, collapse = ""), "</tbody></table>",
           "<p class=\"nota\">La memoria es el máximo del árbol de procesos, ",
           "muestreado una vez por segundo: un pico más corto que eso puede no ",
           "quedar registrado.</p>")
  } else ""

  checks_html <- if (!is.null(ck) && nrow(ck)) {
    paste0("<h3>Huella de los ficheros de entrada</h3>",
           "<table><thead><tr><th>Fichero</th><th>Bytes</th><th>md5</th></tr></thead><tbody>",
           paste(vapply(seq_len(nrow(ck)), function(i) paste0(
             "<tr><td>", html_escape(basename(ck$file[i])), "</td><td>",
             html_escape(ck$size_bytes[i]), "</td><td><code>",
             html_escape(ck$md5[i]), "</code></td></tr>"), character(1)), collapse = ""),
           "</tbody></table>")
  } else ""

  # ── Enriquecimiento funcional ────────────────────────────────────────────
  e <- rv$enrich
  enrich_html <- if (is.null(e)) "" else {
    campos <- list(
      "Enfoque"            = e$enfoque %||% "—",
      "Ontologia / base"   = e$ontologia %||% "—",
      "keyType"            = e$keytype %||% "—",
      # La traducción de identificadores cambia que genes entran en el test, no
      # como se presentan: sin ella en el informe, el resultado no se puede
      # reproducir a partir de la matriz de conteos.
      "Traducción de IDs"  = if (is.na(e$traduccion %||% NA)) "no (se usaron los IDs de la matriz)"
                             else e$traduccion,
      # Los simbolos cambian entre versiones de la anotación, así que la ruta
      # del fichero es parte del resultado igual que lo es la versión del OrgDb.
      "Anotación usada"    = if (is.na(e$anotacion_traduccion %||% NA)) "—"
                             else basename(e$anotacion_traduccion),
      "Genes en la lista"  = if (is.na(e$n_lista %||% NA)) "— (GSEA usa el ranking completo)"
                             else fmt_int(e$n_lista),
      "Universo (fondo)"   = fmt_int(e$n_universo %||% 0),
      "Tasa de mapeo"      = if (is.null(e$mapeo) || is.na(e$mapeo$rate %||% NA)) "—"
                             else paste0(round(100 * e$mapeo$rate, 1), " % (",
                                         e$mapeo$n_mapped, "/", e$mapeo$n_input, ")"),
      "Términos obtenidos" = fmt_int(e$n_terminos %||% 0)
    )
    if (!is.na(e$metrica %||% NA)) {
      campos[["Metrica del ranking"]] <- e$metrica
      campos[["Ponderación (exponent)"]] <- e$exponent
      # Estos cuatro cambian que conjuntos se testean y cuales se devuelven, de
      # modo que sin ellos el resultado no es reproducible.
      campos[["Tamaño de conjunto"]] <- paste0(e$gsea_min_size %||% "—", " - ",
                                               e$gsea_max_size %||% "—", " genes")
      campos[["Corte de p ajustado"]] <- e$gsea_pcutoff %||% "—"
      campos[["Precisión del p-valor (eps)"]] <-
        if (identical(e$gsea_eps %||% NA_real_, 0)) "exacto (eps = 0)"
        else paste0("truncado en ", e$gsea_eps %||% "1e-10")
    }
    if (!is.na(e$direccional %||% NA)) campos[["ORA por dirección"]] <- e$direccional
    if (!is.na(e$gmt %||% NA)) campos[["Gene sets propios"]] <- e$gmt
    if (!is.null(e$simbolos)) campos[["IDs a simbolos"]] <- e$simbolos
    if (!is.na(e$simplify %||% NA)) campos[["Colapsar redundantes"]] <- e$simplify
    if (!is.na(e$organismo_kegg %||% NA)) campos[["Organismo KEGG"]] <- e$organismo_kegg
    if (!is.na(e$organismo_reactome %||% NA))
      campos[["Organismo Reactome"]] <- e$organismo_reactome
    if (!is.na(e$error %||% NA)) campos[["Resultado"]] <- e$error
    paste0("<h2>Enriquecimiento funcional</h2>", html_kv_table(campos),
           "<h3>Versión de la anotación</h3>", html_kv_table(e$anotacion %||% list()),
           "<p class=\"nota\">Los resultados de enriquecimiento dependen de la ",
           "versión de la anotación: usar anotaciones desactualizadas altera las ",
           "rutas que salen enriquecidas (Wadi et al., Nature Methods 2016). El ",
           "campo KEGG de los OrgDb está congelado desde 2011 y NO se usa: la ",
           "anotación KEGG se consulta en línea.</p>")
  }

  # ── Semillas y parámetros estocasticos ───────────────────────────────────
  sd <- rv$seeds %||% list()
  seeds_html <- {
    campos <- list()
    if (!is.null(sd$sva))   campos[["Semilla de sva (num.sv)"]] <- sd$sva
    if (!is.null(sd$n_sv))  campos[["Variables sustitutas estimadas"]] <- sd$n_sv
    if (!is.null(sd$swish)) campos[["Semilla de Swish"]] <- sd$swish
    if (!is.null(sd$swish_nperms)) campos[["Permutaciones de Swish"]] <- sd$swish_nperms
    if (!is.null(diagnostics$replicability))
      campos[["Semilla del bootstrap"]] <- ANALYSIS_SEED
    # No se afirma que el análisis sea determinista si no se ha registrado nada:
    # DESeq2 con IHW y el enriquecimiento GSEA tienen componentes estocasticos
    # que no siempre quedan anotados aquí, y declarar determinismo sin haberlo
    # comprobado es peor que no decir nada.
    if (!is.null(rv$enrich)) campos[["Semilla de GSEA"]] <- ANALYSIS_SEED
    if (identical(rv$padj_method %||% "BH", "IHW"))
      campos[["Semilla de IHW"]] <- "1 (valor por defecto de IHW::ihw)"
    if (!length(campos)) campos[["Componentes aleatorios registrados"]] <-
      paste("ninguno. El ajuste de este motor es determinista, pero esta línea",
            "solo cubre lo que la aplicación registra: revisa las versiones de",
            "los paquetes si necesitas reproducir el resultado bit a bit.")
    paste0("<h3>Semillas y componentes aleatorios</h3>", html_kv_table(campos))
  }

  # ── Limitaciones ─────────────────────────────────────────────────────────
  lims <- deg_report_limitations(rv, diagnostics, sig_n)
  lims_html <- paste0(
    "<h2>Limitaciones</h2>",
    if (!length(lims)) "<p>No se han detectado limitaciones destacables.</p>" else
      paste0("<ul>", paste0("<li>", vapply(lims, html_escape, character(1)),
                            "</li>", collapse = ""), "</ul>"))

  diag_html <- ""
  if (!is.null(diagnostics)) {
    bits <- list()
    if (!is.null(diagnostics$verdict)) {
      # Se acepta tanto la lista de diagnose_pvalue_shape() como una cadena
      # suelta: el informe no puede caerse por la forma de un diagnóstico.
      v <- diagnostics$verdict
      bits[["Forma del histograma de p-valores"]] <- if (is.list(v))
        paste0(v$label %||% "—", " — ", v$detail %||% "") else as.character(v)
    }
    if (!is.null(diagnostics$pi0) && !is.na(diagnostics$pi0)) {
      bits[["pi0 (proporción de nulas ciertas)"]] <- round(diagnostics$pi0, 4)
    }
    if (!is.null(diagnostics$na_breakdown)) {
      bits[["Genes sin p-valor ajustado"]] <-
        padj_na_breakdown_text(diagnostics$na_breakdown)
    }
    if (!is.null(diagnostics$cooks_dominant) && !is.na(diagnostics$cooks_dominant)) {
      bits[["Muestra que concentra los outliers de Cook"]] <- diagnostics$cooks_dominant
    }
    # El objeto de replicabilidad puede venir de un intento fallido (sin summary),
    # así que se comprueba antes de leerlo en lugar de asumir que está completo.
    r <- diagnostics$replicability
    if (!is.null(r)) {
      bits[["Replicabilidad (bootstrap)"]] <- if (!is.null(r$summary) &&
                                                  is.numeric(r$summary$mediana)) {
        paste0(r$interpretation$label, " — Spearman mediana ",
               round(r$summary$mediana[1], 3), " sobre ", r$n_ok, " remuestreos")
      } else {
        paste0("No estimada: ", r$error %||% "motivo desconocido")
      }
    }
    if (length(bits)) {
      diag_html <- paste0("<h2>Diagnosticos</h2>", html_kv_table(bits))
    }
  }

  top <- rv$results[!is.na(rv$results$padj), , drop = FALSE]
  top <- top[order(top$padj), , drop = FALSE]
  top <- utils::head(top, 25)
  cols <- intersect(c("gene", "baseMean", "log2FC", "log2FC_shrunk", "pvalue", "padj"),
                    names(top))
  top_rows <- if (nrow(top)) paste(vapply(seq_len(nrow(top)), function(i) {
    paste0("<tr>", paste(vapply(cols, function(c) {
      v <- top[[c]][i]
      paste0("<td>", html_escape(if (is.numeric(v)) signif(v, 4) else v), "</td>")
    }, character(1)), collapse = ""), "</tr>")
  }, character(1)), collapse = "") else ""

  si <- utils::capture.output(utils::sessionInfo())
  pkgs <- c("DESeq2", "edgeR", "limma", "apeglm", "ashr", "IHW", "sva",
            "clusterProfiler", "fgsea", "tximport")
  pkg_versions <- stats::setNames(lapply(pkgs, function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p))
    else "no instalado"
  }), pkgs)

  paste0(
'<!DOCTYPE html><html lang="es"><head><meta charset="utf-8">',
'<title>Informe de expresión diferencial</title><style>',
# El fondo se fija EXPLICITAMENTE: sin el, un navegador en modo oscuro deja el
# informe con texto oscuro sobre fondo negro e ilegible. Es un documento
# descargable que se leera en cualquier parte, así que no puede depender del
# tema de quien lo abra.
'body{font-family:system-ui,-apple-system,Segoe UI,Roboto,sans-serif;max-width:900px;',
'margin:2rem auto;padding:0 1rem;color:#20332A;background:#FFFFFF;line-height:1.5;',
'-webkit-print-color-adjust:exact;print-color-adjust:exact}',
'@media (prefers-color-scheme:dark){body{color:#20332A;background:#FFFFFF}}',
'h1{border-bottom:3px solid #7BBF9A;padding-bottom:.3rem}',
'h2{margin-top:2rem;color:#244B34}',
'table{border-collapse:collapse;width:100%;margin:.5rem 0;font-size:.92rem}',
'th,td{border:1px solid #D9E2DC;padding:.35rem .5rem;text-align:left}',
'table.kv th{width:34%;background:#F2F7F4}',
'thead th{background:#E8F1EC}',
'.metric{display:inline-block;background:#F2F7F4;border:1px solid #D9E2DC;',
'border-radius:6px;padding:.5rem .9rem;margin:.25rem .4rem .25rem 0}',
'.metric b{display:block;font-size:1.3rem;color:#244B34}',
'pre{background:#F7F9F8;border:1px solid #E3EAE6;padding:.6rem;overflow-x:auto;',
'font-size:.8rem}',
'h3{margin-top:1.2rem;color:#315342;font-size:1rem}',
'.resumen{background:#F2F7F4;border-left:4px solid #7BBF9A;padding:.8rem 1rem;',
'margin:1rem 0;border-radius:0 6px 6px 0}',
'.nota{font-size:.85rem;color:#60756A;margin-top:.4rem}',
'ul li{margin-bottom:.4rem}',
'code{font-size:.85em}',
'footer{margin-top:2.5rem;font-size:.82rem;color:#60756A;border-top:1px solid #D9E2DC;',
'padding-top:.6rem}</style></head><body>',
'<h1>Informe de expresión diferencial</h1>',
'<p class="nota">Análisis realizado el ',
html_escape(format(rv$run_at %||% Sys.time(), "%Y-%m-%d %H:%M:%S")), '.</p>',
'<div class="resumen"><b>Resumen</b><br>',
html_escape(deg_plain_summary(rv, sig_n, up, down)), '</div>',
'<p><span class="metric"><b>', fmt_int(nrow(rv$results)), '</b>genes analizados</span>',
'<span class="metric"><b>', fmt_int(sig_n), '</b>significativos a FDR &le; ',
html_escape(rv$fdr %||% 0.05), '</span>',
'<span class="metric"><b>', fmt_int(up), ' / ', fmt_int(down),
'</b>up / down</span></p>',
'<h2>Parámetros del análisis</h2>', html_kv_table(params),
proc_html, pipeline_html, metricas_html, checks_html,
diag_html,
enrich_html,
lims_html,
'<h2>Metadatos de las muestras</h2>',
if (!is.null(rv$meta)) paste0(
  '<table><thead><tr>',
  paste0('<th>', vapply(names(rv$meta), html_escape, character(1)), '</th>', collapse = ''),
  '</tr></thead><tbody>',
  paste(vapply(seq_len(nrow(rv$meta)), function(i) paste0('<tr>',
    paste0('<td>', vapply(rv$meta[i, ], function(v) html_escape(as.character(v)), character(1)),
           '</td>', collapse = ''), '</tr>'), character(1)), collapse = ''),
  '</tbody></table>') else '<p>Sin metadatos.</p>',
if (sig_n > 0) '<h2>Top 25 genes por p-valor ajustado</h2>' else
  paste0('<h2>Genes con menor p-valor ajustado</h2>',
         '<p class="nota">Ninguno alcanza el umbral de significacion; se listan ',
         'igualmente los de menor p-valor ajustado, porque un informe negativo ',
         'debe documentar que se evaluo y con que resultado.</p>'),
'<table><thead><tr>',
paste0('<th>', vapply(cols, html_escape, character(1)), '</th>', collapse = ''),
'</tr></thead><tbody>', top_rows, '</tbody></table>',
'<h2>Reproducibilidad</h2>',
html_kv_table(app_provenance()),
seeds_html,
'<h3>Versiones de los paquetes</h3>', html_kv_table(pkg_versions),
'<h3>sessionInfo()</h3><pre>', html_escape(paste(si, collapse = "\n")), '</pre>',
'<footer>Generado por SARA el ',
html_escape(format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
'. Este informe recoge los parámetros con los que se ejecuto el análisis; ',
'el script R equivalente se puede descargar aparte para reproducirlo fuera de la app.',
'</footer></body></html>')
}

#' Script R equivalente al análisis ejecutado.
#'
#' No es una transcripción de la app: es el código mínimo que reproduce el mismo
#' resultado usando directamente los paquetes de Bioconductor, para que se pueda
#' auditar y correr sin la app.
build_deg_r_script <- function(rv) {
  if (is.null(rv) || is.null(rv$results)) return(NULL)
  method <- rv$method %||% "DESeq2"
  design <- rv$design %||% "~ condition"
  # La formula como CÓDIGO. Para Wilcoxon la etiqueta legible es prosa ("sin
  # modelo..."), y interpolarla dentro de model.matrix() producía un script que
  # no parseaba. Cuando no hay modelo, el prefiltrado usa la via `group=`.
  design_code <- rv$design_code %||%
    (if (identical(rv$method %||% "", "Wilcoxon")) NULL else design)
  # Diseño con el que se PREFILTRO, que no es `design_code`. La app construye la
  # matriz de `filterByExpr` con `build_design(meta, ref, batch)`, es decir sin
  # formula libre y —sobre todo— sin las variables sustitutas, porque estás se
  # estiman DESPUÉS, sobre la matriz ya prefiltrada. Emitir `design_code` aquí
  # producía un `model.matrix(~ ... + SV1)` colocado ANTES del bloque que crea
  # SV1, de modo que el script fallaba con "object 'SV1' not found" en cuanto el
  # análisis usaba sva con el prefiltrado automático (el modo por defecto).
  prefilter_design_code <- if (identical(method, "Wilcoxon")) NULL
    else if (!is.null(rv$batch) && nzchar(rv$batch %||% ""))
      paste0("~ ", rv$batch, " + condition")
    else "~ condition"
  fdr <- rv$fdr %||% 0.05
  lfc <- rv$lfc_threshold %||% 0
  coef <- rv$coef %||% ""
  # Denominador REAL del contraste ajustado. Antes se tomaba el primer nivel en
  # orden ALFABETICO, que solo coincide con el denominador elegido por
  # casualidad: en cuanto el usuario contrastaba, por ejemplo, "ctrl vs trt", el
  # script releveleaba a "ctrl" y luego pedia un coeficiente que no existe, de
  # modo que fallaba al ejecutarse. El contraste es parte del ajuste, no algo
  # que se pueda deducir de los datos a posteriori.
  ref <- rv$ref_level %||%
    # Respaldo para estados guardados antes de que se registrara ref_level:
    # el contraste se muestra como "num vs den".
    (if (!is.null(rv$contrast) && grepl(" vs ", rv$contrast, fixed = TRUE))
       sub("^.* vs ", "", rv$contrast) else NULL) %||%
    (if (!is.null(rv$meta) && "condition" %in% names(rv$meta))
       levels(as.factor(as.character(rv$meta$condition)))[1] else "")
  num <- rv$contrast_num %||%
    (if (!is.null(rv$contrast) && grepl(" vs ", rv$contrast, fixed = TRUE))
       sub(" vs .*$", "", rv$contrast) else "")
  seeds <- rv$seeds %||% list()
  pf <- rv$prefilter

  header <- c(
    "# Script equivalente al análisis ejecutado en SARA",
    paste0("# Generado el ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("# Contraste: ", rv$contrast %||% "—",
           "   |   Diseño: ", design,
           "   |   Motor: ", method),
    "#",
    "# Entradas que hay que proporcionar:",
    if (identical(method, "Swish"))
      "#   quant_dir   directorio 03_alignments de la ejecución (salmon/kallisto)"
    else "#   counts.tsv  matriz de conteos (genes x muestras)",
    "#   meta.tsv    samplesheet con sample_id, condition y las covariables del diseño",
    "",
    "# Semilla de todo lo estocastico del análisis. Sin fijarla, los resultados",
    "# que dependen de permutaciones no son reproducibles.",
    paste0("set.seed(", seeds$sva %||% seeds$swish %||% 1L, ")"),
    "",
    if (identical(method, "Swish")) character(0) else
      'counts <- as.matrix(read.delim("counts.tsv", row.names = 1, check.names = FALSE))',
    'meta   <- read.delim("meta.tsv", stringsAsFactors = FALSE)',
    # Con Swish no hay objeto `counts`: el orden se fija por sample_id y las
    # cuantificaciones se reordenan después contra el.
    if (identical(method, "Swish"))
      'meta   <- meta[order(meta$sample_id), , drop = FALSE]'
    else 'meta   <- meta[match(colnames(counts), meta$sample_id), , drop = FALSE]',
    paste0('# El denominador del contraste es el nivel de referencia del factor.'),
    paste0('meta$condition <- relevel(factor(meta$condition), ref = "', ref, '")'),
    ""
  )

  # Variables sustitutas: la formula del diseño las menciona (SV1, SV2...) pero
  # no están en meta.tsv, así que hay que volver a estimarlas o el script no
  # correria. Se reproducen con la misma semilla y el mismo número.
  #
  # El modelo de interes que se pasa a svaseq es `design_base`, el diseño SIN las
  # SV: pasarle `design_code`, que ya las incluye, sería circular. El respaldo a
  # "~ condition" solo cubre estados guardados antes de que se registrara ese
  # campo; usarlo cuando el ajuste llevaba batch o formula libre estimaba las SV
  # con otro modelo y el script no reproducia el resultado.
  sva_block <- if (!is.null(seeds$n_sv) && seeds$n_sv > 0) c(
    "# El diseño incluye variables sustitutas estimadas con sva. Para reproducir",
    "# el ajuste hay que volver a estimarlas: no viajan en el samplesheet.",
    paste0("mod  <- model.matrix(", rv$design_base %||% "~ condition", ", data = meta)"),
    "mod0 <- model.matrix(~ 1, data = meta)",
    "cm_sv <- as.matrix(counts)[rowMeans(as.matrix(counts)) > 1, , drop = FALSE]",
    paste0("sv <- sva::svaseq(cm_sv, mod, mod0, n.sv = ", seeds$n_sv, ")$sv"),
    "colnames(sv) <- paste0(\"SV\", seq_len(ncol(sv)))",
    "meta <- cbind(meta, as.data.frame(sv))",
    ""
  ) else character(0)

  # Swish no parte de la matriz de conteos, así que el prefiltrado no aplica:
  # emitirlo produciría un script que referencia un objeto `counts` inexistente.
  prefilter <- if (identical(method, "Swish")) character(0) else
    if (!is.null(pf) && identical(pf$mode, "filterByExpr")) c(
    "# Prefiltrado: filterByExpr usa el tamaño del grupo más pequeño",
    if (!is.null(prefilter_design_code))
      paste0("design <- model.matrix(", prefilter_design_code, ", data = meta)")
    else character(0),
    "y <- edgeR::DGEList(counts = round(counts), group = meta$condition)",
    if (!is.null(prefilter_design_code)) "keep <- edgeR::filterByExpr(y, design = design)"
    else "keep <- edgeR::filterByExpr(y, group = meta$condition)",
    "counts <- counts[keep, , drop = FALSE]",
    ""
  ) else if (!is.null(pf)) c(
    paste0("# Prefiltrado manual: >= ", pf$min_count, " conteos en >= ",
           pf$min_samples, " muestras"),
    paste0("keep <- rowSums(counts >= ", pf$min_count, ") >= ", pf$min_samples),
    "counts <- counts[keep, , drop = FALSE]",
    ""
  ) else character(0)

  body <- switch(method,
    "DESeq2" = c(
      "library(DESeq2)",
      # IHW necesita S4Vectors ATACHADO, no solo cargado: sin esto falla con
      # "no se pudo encontrar la función mcols". Es el mismo tropiezo que tuvo
      # la propia app, así que el script no puede omitirlo.
      if (identical(rv$padj_method, "IHW"))
        "library(IHW); library(S4Vectors)  # IHW necesita S4Vectors atachado"
      else character(0),
      paste0("dds <- DESeqDataSetFromMatrix(round(counts), meta, ", design_code, ")"),
      if (identical(rv$outliers %||% "na", "refit"))
        "dds <- DESeq(dds, minReplicatesForReplace = 3)  # outliers: sustituir"
      else "dds <- DESeq(dds)",
      paste0("# alpha = FDR objetivo: calibra el filtrado independiente"),
      paste0("res <- results(dds, name = \"", coef, "\", alpha = ", fdr,
             if (lfc > 0) paste0(", lfcThreshold = ", lfc,
                                 ", altHypothesis = \"greaterAbs\"") else "",
             if (identical(rv$padj_method, "IHW")) ", filterFun = IHW::ihw" else "",
             # cooksCutoff = FALSE trata los valores extremos como biologia real
             if (identical(rv$outliers %||% "na", "keep")) ", cooksCutoff = FALSE" else "",
             ")"),
      if (!identical(rv$shrink %||% "ninguno", "ninguno")) c(
        "# El encogido cambia el log2FC pero no los p-valores",
        paste0("shr <- lfcShrink(dds, coef = \"", coef, "\", type = \"",
               rv$shrink, "\")")) else character(0),
      "res <- as.data.frame(res)"
    ),
    "edgeR" = c(
      "library(edgeR)",
      paste0("design <- model.matrix(", design_code, ", data = meta)"),
      "y <- DGEList(counts = round(counts), group = meta$condition)",
      "y <- normLibSizes(y)                 # edgeR v4; antes calcNormFactors",
      "fit <- glmQLFit(y, design)           # v4 estima las dispersiones aquí",
      if (lfc > 0)
        paste0("test <- glmTreat(fit, coef = \"", coef, "\", lfc = ", lfc, ")")
      else paste0("test <- glmQLFTest(fit, coef = \"", coef, "\")"),
      "res <- topTags(test, n = Inf, sort.by = \"none\")$table"
    ),
    "limma-voom" = c(
      "library(limma); library(edgeR)",
      paste0("design <- model.matrix(", design_code, ", data = meta)"),
      "y <- normLibSizes(DGEList(counts = round(counts), group = meta$condition))",
      "v <- voomWithQualityWeights(y, design)",
      "fit <- lmFit(v, design)",
      if (lfc > 0) c(
        paste0("fit <- treat(fit, lfc = ", lfc, ", robust = TRUE)"),
        paste0("res <- topTreat(fit, coef = \"", coef, "\", number = Inf, sort.by = \"none\")"))
      else c(
        "fit <- eBayes(fit, robust = TRUE)",
        paste0("res <- topTable(fit, coef = \"", coef, "\", number = Inf, sort.by = \"none\")"))
    ),
    "Wilcoxon" = c(
      "# Wilcoxon rank-sum sobre CPM: sin modelo, no ajusta covariables.",
      "# La normalización es TMM y no CPM por tamaño de libreria: sin corregir",
      "# por composicion, una diferencia de composicion entre grupos se",
      "# convierte en falsos positivos (Li et al., Genome Biology 2022).",
      "y <- edgeR::normLibSizes(edgeR::DGEList(counts = round(counts)))",
      "cpm <- edgeR::cpm(y, normalized.lib.sizes = TRUE)",
      "g <- as.character(meta$condition)",
      paste0("num <- \"", sub(" vs .*$", "", rv$contrast %||% ""), "\"; ",
             "den <- \"", sub("^.* vs ", "", rv$contrast %||% ""), "\""),
      "pv <- apply(cpm, 1, function(x) wilcox.test(x[g == num], x[g == den],",
      "                                            exact = FALSE)$p.value)",
      "res <- data.frame(pvalue = pv, padj = p.adjust(pv, method = \"BH\"))"
    ),
    "dearseq" = c(
      "library(dearseq)",
      "# El motor subsetea a los dos niveles del contraste y pasa el batch como",
      "# covariable; reproducirlo sin eso daría otro resultado.",
      paste0("keep_s <- as.character(meta$condition) %in% c(\"", num, "\", \"", ref, "\")"),
      "m <- droplevels(meta[keep_s, , drop = FALSE])",
      "cm <- round(counts)[, keep_s, drop = FALSE]",
      "v2t <- model.matrix(~ condition, data = m)[, -1, drop = FALSE]",
      if (!is.null(rv$batch) && nzchar(rv$batch %||% ""))
        paste0("cov <- model.matrix(~ ", rv$batch, ", data = m)[, -1, drop = FALSE]")
      else "cov <- NULL",
      "res <- dear_seq(exprmat = cm, variables2test = v2t, covariates = cov,",
      "                which_test = \"asymptotic\", preprocessed = FALSE)$pvals"
    ),
    "Swish" = c(
      "library(fishpond); library(tximport); library(SummarizedExperiment)",
      "# Swish no parte de la matriz de conteos: la incertidumbre de la",
      "# cuantificación vive en las replicas inferenciales de salmon/kallisto.",
      "# El análisis queda a nivel de TRANSCRITO (txOut = TRUE).",
      paste0('qfiles <- list.files("quant_dir", pattern = "',
             if (identical(rv$quant_tool %||% "", "kallisto")) "abundance.h5"
             else "quant.sf",
             '", recursive = TRUE, full.names = TRUE)'),
      "names(qfiles) <- basename(dirname(qfiles))",
      "qfiles <- qfiles[match(meta$sample_id, names(qfiles))]",
      paste0('txi <- tximport(qfiles, type = "', rv$quant_tool %||% "salmon",
             '", txOut = TRUE, dropInfReps = FALSE)'),
      "se <- SummarizedExperiment(",
      "  assays = c(list(counts = txi$counts, abundance = txi$abundance,",
      "                  length = txi$length), txi$infReps),",
      "  colData = S4Vectors::DataFrame(meta))",
      paste0('colData(se)$condition <- factor(as.character(meta$condition),',
             ' levels = c("', ref, '", "', num, '"))'),
      "y <- scaleInfReps(se)",
      "y <- labelKeep(y)",
      "y <- y[mcols(y)$keep, ]",
      paste0("y <- swish(y, x = \"condition\", nperms = ",
             seeds$swish_nperms %||% 30L, ")"),
      "res <- as.data.frame(mcols(y))",
      "# Swish reporta qvalue, no padj; se renombra para el resumen de abajo.",
      "res$padj <- res$qvalue"
    ),
    "# Motor no reconocido"
  )

  footer <- c(
    "",
    paste0("# Genes significativos a FDR <= ", fdr),
    paste0("sum(res$padj <= ", fdr, ", na.rm = TRUE)"),
    "",
    "sessionInfo()"
  )
  paste(c(header, prefilter, sva_block, body, footer), collapse = "\n")
}
