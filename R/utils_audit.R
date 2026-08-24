#' utils_audit.R
#' Persistencia de los analisis diferenciales y registro de auditoria.
#'
#' Por que existe: hasta ahora el informe y el script solo existian si el usuario
#' pulsaba descargar. Una ejecucion del pipeline dejaba su rastro en disco
#' (log, run_params.tsv), pero los analisis diferenciales hechos sobre ella no
#' dejaban ninguno: no habia forma de saber cuantos se lanzaron, con que
#' parametros, ni cual produjo la figura que acabo en la memoria.
#'
#' El patron es el del cuaderno de recogida de datos electronico de las normas
#' de buena practica clinica: un registro append-only de quien hizo que, cuando
#' y con que parametros. No pretende ser inviolable —es un fichero de texto—,
#' pero convierte "creo que use FDR 0,05" en un dato consultable.

#' Ruta del registro de auditoria (uno por instalacion, en outputs/).
audit_log_path <- function(outputs_dir = NULL) {
  base <- outputs_dir %||% outputs_base_dir()
  file.path(base, "audit_log.tsv")
}

#' Añade una linea al registro de auditoria.
#'
#' Nunca falla hacia fuera: un problema al escribir el registro no puede tumbar
#' un analisis que ya se ha calculado correctamente.
#'
#' @param accion etiqueta corta de la accion ("deg_run", "deg_export"...)
#' @param detalles lista clave-valor; se serializa como "k=v; k=v"
#' @param outputs_dir directorio base; por defecto el de la app
#' @return TRUE si se escribio, FALSE si no
append_audit_log <- function(accion, detalles = list(), outputs_dir = NULL) {
  out <- tryCatch({
    f <- audit_log_path(outputs_dir)
    dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
    nuevo <- !file.exists(f)
    txt <- paste(vapply(names(detalles), function(k) {
      v <- detalles[[k]]
      v <- if (is.null(v) || !length(v)) "—" else paste(as.character(v), collapse = ",")
      # Ni tabuladores ni saltos: cada evento tiene que caber en una linea para
      # que el fichero siga siendo un TSV legible.
      paste0(k, "=", gsub("[\t\r\n]+", " ", v))
    }, character(1)), collapse = "; ")
    linea <- paste(
      format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
      tryCatch(unname(Sys.info()[["user"]]), error = function(e) "—"),
      accion, txt, sep = "\t")
    con <- file(f, open = "a", encoding = "UTF-8")
    on.exit(close(con), add = TRUE)
    if (nuevo) writeLines("timestamp\tusuario\taccion\tdetalles", con)
    writeLines(linea, con)
    TRUE
  }, error = function(e) FALSE)
  invisible(out)
}

#' Lee el registro de auditoria.
read_audit_log <- function(outputs_dir = NULL, n_max = 500L) {
  f <- audit_log_path(outputs_dir)
  if (!file.exists(f)) return(NULL)
  df <- tryCatch(utils::read.delim(f, stringsAsFactors = FALSE, quote = ""),
                 error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(NULL)
  utils::tail(df, n_max)
}

#' Persiste un analisis diferencial completo en disco.
#'
#' Escribe en `<destino>/05_deg/<timestamp>/`:
#'   - `deg_params.tsv`  parametros del ajuste, en el mismo formato clave-valor
#'                       que run_params.tsv del pipeline
#'   - `resultados.tsv`  la tabla COMPLETA, no solo los significativos
#'   - `samplesheet.tsv` los metadatos tal como se usaron (incluidas las SV)
#'   - `informe.html`    el mismo informe que se descarga desde la interfaz
#'   - `analisis.R`      el script equivalente
#'
#' @param rv `state$deg_rv`
#' @param diagnostics diagnosticos para el informe (puede ser NULL)
#' @param base_dir directorio de la ejecucion de origen; si no hay (matriz
#'   subida), se usa `outputs/analisis_sueltos/`
#' @return ruta del directorio escrito, o NULL si no se pudo
persist_deg_analysis <- function(rv, diagnostics = NULL, base_dir = NULL,
                                 outputs_dir = NULL) {
  if (is.null(rv) || is.null(rv$results)) return(NULL)
  out <- tryCatch({
    raiz <- if (!is.null(base_dir) && nzchar(base_dir %||% "") && dir.exists(base_dir)) {
      base_dir
    } else {
      file.path(outputs_dir %||% outputs_base_dir(), "analisis_sueltos")
    }
    stamp <- format(rv$run_at %||% Sys.time(), "%Y%m%d_%H%M%S")
    dest <- file.path(raiz, "05_deg", stamp)
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)

    # Parametros en el mismo formato clave-valor que usa el pipeline, para que
    # se puedan leer con la misma funcion.
    par <- deg_run_parameters(rv)
    extra <- list(
      numerador = rv$contrast_num, denominador = rv$ref_level,
      formula_libre = rv$design_formula, coef_elegido = rv$test_coef,
      semilla_sva = rv$seeds$sva, n_variables_sustitutas = rv$seeds$n_sv,
      semilla_swish = rv$seeds$swish, permutaciones_swish = rv$seeds$swish_nperms,
      origen_matriz = rv$counts_origin$tipo, ruta_matriz = rv$counts_origin$ruta,
      md5_matriz = rv$counts_origin$md5,
      resumen_a_gen = rv$counts_source$method,
      resumen_a_gen_ok = rv$counts_source$ok
    )
    lineas <- c(
      vapply(names(par), function(k)
        paste0(k, "\t", paste(as.character(par[[k]]), collapse = ",")), character(1)),
      vapply(names(extra), function(k) {
        v <- extra[[k]]
        paste0(k, "\t", if (is.null(v) || !length(v)) "—" else
               paste(as.character(v), collapse = ","))
      }, character(1)))
    writeLines(unname(lineas), file.path(dest, "deg_params.tsv"))

    # La tabla completa, no solo los significativos: recortarla impediria
    # recalcular otros umbrales o rehacer el enriquecimiento mas tarde.
    utils::write.table(rv$results, file.path(dest, "resultados.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
    if (!is.null(rv$meta)) {
      utils::write.table(rv$meta, file.path(dest, "samplesheet.tsv"),
                         sep = "\t", quote = FALSE, row.names = FALSE)
    }
    h <- tryCatch(build_deg_report_html(rv, diagnostics), error = function(e) NULL)
    if (!is.null(h)) writeLines(h, file.path(dest, "informe.html"))
    s <- tryCatch(build_deg_r_script(rv), error = function(e) NULL)
    if (!is.null(s)) writeLines(s, file.path(dest, "analisis.R"))
    dest
  }, error = function(e) NULL)
  out
}
