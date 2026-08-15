#' utils_privacy.R
#' Seudonimizacion de identificadores de muestra.
#'
#' Por que existe: los identificadores de muestra de un estudio clinico suelen
#' llevar informacion identificativa (codigo de paciente, fecha, centro), y en
#' esta aplicacion viajan a todos los graficos, tablas, informes y ficheros
#' persistidos. Basta con exportar una figura para difundirlos.
#'
#' El matiz que conviene tener claro, y que la interfaz debe decir: esto es
#' SEUDONIMIZACION, no anonimizacion. Existe una tabla de correspondencia que
#' permite deshacerla, asi que los datos siguen siendo datos personales a
#' efectos del RGPD; lo que se consigue es que la informacion identificativa no
#' viaje en los entregables. La tabla de correspondencia se descarga aparte y a
#' proposito, para que la decision de exportarla sea explicita.
#'
#' Y un limite que no se puede ocultar: los datos de EXPRESION son en si mismos
#' reidentificables. Schadt, Woo y Hao (Nature Genetics 2012) demostraron que se
#' pueden inferir genotipos individuales a partir de niveles de expresion. Es
#' decir, renombrar las columnas reduce la exposicion accidental, pero no
#' convierte una matriz de expresion humana en un dato anonimo.

#' Construye una tabla de correspondencia entre identificadores reales y alias.
#'
#' @param ids identificadores originales
#' @param prefijo prefijo del alias
#' @return data.frame(original, alias) en el orden de aparicion
build_pseudonym_map <- function(ids, prefijo = "S") {
  ids <- as.character(ids)
  unicos <- unique(ids[!is.na(ids) & nzchar(ids)])
  if (!length(unicos)) {
    return(data.frame(original = character(0), alias = character(0),
                      stringsAsFactors = FALSE))
  }
  ancho <- max(2L, nchar(as.character(length(unicos))))
  data.frame(
    original = unicos,
    alias = sprintf(paste0(prefijo, "%0", ancho, "d"), seq_along(unicos)),
    stringsAsFactors = FALSE
  )
}

#' Aplica una tabla de correspondencia a un vector de identificadores.
#' Los que no esten en el mapa se dejan intactos: es preferible un identificador
#' sin seudonimizar y visible a uno convertido en NA sin avisar.
apply_pseudonyms <- function(ids, map) {
  if (is.null(map) || !nrow(map)) return(as.character(ids))
  idx <- match(as.character(ids), map$original)
  out <- as.character(ids)
  out[!is.na(idx)] <- map$alias[idx[!is.na(idx)]]
  out
}

#' Seudonimiza a la vez la matriz de conteos y el samplesheet.
#'
#' Se hace en un solo sitio porque las columnas de la matriz y la columna
#' `sample_id` del samplesheet tienen que seguir casando: renombrar una sin la
#' otra rompe el alineamiento y el analisis fallaria de formas dificiles de
#' diagnosticar.
#'
#' @param counts matriz de conteos (genes x muestras)
#' @param meta samplesheet con `sample_id`
#' @param prefijo prefijo del alias
#' @return list(counts, meta, map)
pseudonymize_dataset <- function(counts, meta, prefijo = "S") {
  ids <- if (!is.null(counts)) colnames(counts) else
         if (!is.null(meta) && "sample_id" %in% names(meta)) meta$sample_id else character(0)
  map <- build_pseudonym_map(ids, prefijo)
  if (!nrow(map)) return(list(counts = counts, meta = meta, map = map))
  if (!is.null(counts)) colnames(counts) <- apply_pseudonyms(colnames(counts), map)
  if (!is.null(meta) && "sample_id" %in% names(meta)) {
    meta$sample_id <- apply_pseudonyms(meta$sample_id, map)
    if (!is.null(rownames(meta))) rownames(meta) <- meta$sample_id
  }
  list(counts = counts, meta = meta, map = map)
}

#' Columnas del samplesheet que parecen identificativas.
#'
#' Sirve para avisar, no para borrar nada: la decision es del usuario. Detecta
#' por nombre de columna los campos que el RGPD considera de riesgo en un
#' contexto clinico y los que, por ser unicos por muestra, pueden reidentificar
#' aunque no lo parezcan.
#'
#' @return character con los nombres de columna sospechosos
detect_identifying_columns <- function(meta) {
  if (is.null(meta) || !ncol(meta)) return(character(0))
  patrones <- paste0("(?i)(paciente|patient|nombre|name|apellido|surname|",
                     "nhc|historia|record|dni|nif|ssn|email|correo|telefono|",
                     "phone|direccion|address|fecha_nac|birth|dob|iniciales|",
                     "initials|barcode|codigo_paciente)")
  sospechosas <- grep(patrones, names(meta), value = TRUE, perl = TRUE)
  # Columnas de texto con un valor distinto por muestra: por si solas permiten
  # reidentificar, aunque su nombre no lo delate.
  casi_unicas <- names(meta)[vapply(names(meta), function(nm) {
    v <- meta[[nm]]
    if (!is.character(v) && !is.factor(v)) return(FALSE)
    n <- length(unique(as.character(v)))
    n == nrow(meta) && nrow(meta) > 2
  }, logical(1))]
  unique(c(sospechosas, setdiff(casi_unicas, "sample_id")))
}
