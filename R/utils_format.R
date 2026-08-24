#' utils_format.R
#' Funciones de formato de texto, bytes, tiempos y logs.

#' Formatea un número de bytes en una cadena legible
fmt_bytes <- function(b) {
  b <- as.numeric(b)
  if (is.na(b))    return("—")
  if (b >= 1e9)    return(sprintf("%.1f GB", b / 1e9))
  if (b >= 1e6)    return(sprintf("%.1f MB", b / 1e6))
  if (b >= 1e3)    return(sprintf("%.1f KB", b / 1e3))
  paste0(b, " B")
}

#' Formatea segundos en cadena legible (3s, 2m 15s, 1h 4m)
fmt_elapsed <- function(secs) {
  s <- round(as.numeric(secs))
  if (s < 60)   return(sprintf("%ds", s))
  if (s < 3600) return(sprintf("%dm %ds", s %/% 60, s %% 60))
  sprintf("%dh %dm", s %/% 3600, (s %% 3600) %/% 60)
}

#' Entero con separador de miles ("18.402").
#' `decimal.mark` se pasa explicitamente porque con big.mark = "." y el
#' decimal.mark por defecto también ".", R avisa de la ambiguedad.
fmt_int <- function(x) {
  if (length(x) == 0 || is.null(x) || any(is.na(x))) return("—")
  formatC(as.integer(x), format = "d", big.mark = ".", decimal.mark = ",")
}

#' Etiqueta de porcentaje "12.3%" o "—"
pct_label <- function(x, digits = 1) {
  if (is.na(x) || !is.finite(x)) return("—")
  sprintf(paste0("%.", digits, "f%%"), x)
}

#' Prefija un mensaje de log con timestamp [HH:MM:SS]
ts_log <- function(msg) paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)

#' Normaliza un fragmento de texto de terminal a UTF-8 con saltos \n
terminal_text <- function(x) {
  if (length(x) == 0 || is.null(x)) return("")
  x <- paste(x, collapse = "\n")
  x <- gsub("\r\n?", "\n", x)
  Encoding(x) <- "UTF-8"
  x
}

#' Recorta el log de la app si supera max_chars (el log completo vive en disco)
trim_log_text <- function(x, max_chars = 250000L) {
  if (nchar(x, type = "chars", allowNA = FALSE) <= max_chars) return(x)
  paste0(
    "[log truncado en la app; consulta workflow_live.log para el log completo]\n",
    substring(x, nchar(x) - max_chars + 1L)
  )
}

#' Conversión segura a numérico (suprime warnings y devuelve NA)
num_or_na <- function(x) suppressWarnings(as.numeric(x))
