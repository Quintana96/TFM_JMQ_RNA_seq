# Auditoria + Modulo DEG — Changelog

Fecha: 2026-05-19
Contexto: TFM JMQ — refactor modular de la app Shiny RNA-seq.

---

## PARTE 1 — Bugs auditados

### Bug 1 — `infer_result_params()` hardcodea `read_type = "Paired-end"`
- **Estado:** Arreglado.
- **Archivo:** `R/utils_status.R`.
- **Que hago:** Nueva funcion auxiliar `infer_read_type_from_dir(out_dir)` que
  inspecciona `out_dir/02_trimmed_reads/`:
  - Si encuentra `*_R1_trimmed.fastq.gz` y `*_R2_trimmed.fastq.gz` → `"Paired-end"`.
  - Si solo hay `*_trimmed.fastq.gz` sin R1/R2 → `"Single-end"`.
  - Fallback: `"Paired-end"`.
- `infer_result_params()` ahora llama a `infer_read_type_from_dir()` en
  vez de devolver el literal.

### Bug 2 — Regex `multi_col` demasiado amplio en `align_qc_summary()`
- **Estado:** Arreglado.
- **Archivo:** `R/utils_qc.R`.
- **Que hago:** Patrones especificos:
  `c("multi_rate", "multimapping_rate", "pct_multi", "percent_multi",
     "aligned.*multi.*rate", "multi.*percent")`.
  Si no hay match, `multimapping_rate <- NA_real_` (no se inventa a partir de
  conteos crudos). Ademas, heuristica conservadora: si el nombre contiene
  `_pct`, `percent`, etc. o la mediana > 1, dividimos por 100; en caso
  contrario asumimos fraccion.

### Bug 3 — `unique_reads` mal etiquetado cuando `multimapping_rate` es NA
- **Estado:** Arreglado.
- **Archivo:** `R/utils_qc.R` + `R/server_tab_results.R`.
- **Que hago:**
  - Si `multimapping_rate` es NA → `unique_reads <- NA_real_`.
  - Nueva columna `mapped_reads = total_reads * mapping_rate`.
  - El plot `align_qc_mapping_plot_obj()` ahora usa `mapped_reads`/`unmapped_reads`
    cuando no hay descomposicion multi (con anotacion explicita "Sin tasa de
    multimapping disponible").

### Bug 4 — `alignment_table()` sin rama `kallisto`
- **Estado:** Arreglado.
- **Archivo:** `R/utils_tables.R`.
- **Que hago:** Anadido bloque `else if (identical(tool, "kallisto"))` que
  usa `kallisto_stats(out_dir)` y un `keep` con las columnas tipicas de
  MultiQC para kallisto. Si la lectura es NULL, devuelve `message_df()`
  con mensaje claro (he homogeneizado los mensajes con `message_df` en
  todas las ramas).

### Bug 5 — `prepare_uploaded_input_file()` salta copia solo por tamano
- **Estado:** Arreglado.
- **Archivo:** `R/utils_io.R`.
- **Que hago:** Si destino existe **y** mismo tamano **y** mismo `tools::md5sum`
  → no copia. Cualquier error en la comparacion fuerza copia por defecto
  (lado seguro).

### Bug 6 — `median_tpm` incluye ceros en `pseudo_qc_summary()`
- **Estado:** Arreglado.
- **Archivo:** `R/utils_qc.R`.
- **Que hago:** Anadida columna `median_tpm_detected = median(q$TPM[q$TPM > 0])`
  por muestra, junto a la `median_tpm` original. Asi el usuario ve ambas.

### Bug 7 — `pseudo_qc_alerts()` con umbral `tpm_distribution_shift_warning`
- **Estado:** Arreglado.
- **Decision:** Uso `median_tpm_detected` si esta disponible, fallback a
  `median_tpm`. Documentado tambien con un bloque `#'` en el codigo. Razon:
  `median_tpm` colapsa a 0 cuando hay muchos features no detectados y no
  discrimina shifts reales entre muestras.

### Bug 8 — Patron `_fastqc\\.html$` en `important_artifacts()`
- **Estado:** Falso positivo. No requiere cambio.
- **Justificacion:** El codigo actual ya hace
  `list.files(file.path(out_dir, "01_quality"), pattern = "_fastqc\\.html$", full.names = TRUE)`,
  que es exactamente la ruta donde el workflow.sh pone los HTMLs. Confirmado
  por revision del codigo. No se toca.

### Bug 9 — `bindCache(rx, cache_key())` en `server_tab_results.R`
- **Estado:** Arreglado.
- **Archivo:** `R/server_tab_results.R`.
- **Que hago:** Sustituido por la forma idiomatica:
  `bindCache(rx, selected_result_dir(), state$results_refresh())`.
  Aplicado en `selected_result_summary`, `selected_counts_tables` y
  `selected_general_table`. Se elimina la reactive intermedia `cache_key()`.

---

## PARTE 2 — Modulo DEG (Tab 4)

### Archivos nuevos
- `R/utils_deg.R` — funciones puras DEG (validate, align, prefilter, motores
  DESeq2 / edgeR / limma-voom, dispatcher, filtros, vst/rlog, PCA, distancias,
  top-var).
- `R/utils_enrich.R` — funciones puras de enriquecimiento (GO via clusterProfiler,
  KEGG, ordenacion para dotplot).
- `R/ui_tab_deg.R` — UI del nav_panel "4. Expresion diferencial":
  5 cards (Datos, Metadatos, Diseno, Analisis, Filtros) + zona de resultados con
  navset (Tabla, Volcano, MA, PCA, Heatmap top-N, Distancia, Enriquecimiento).
- `R/server_tab_deg.R` — server del Tab 4 con manejo de fuentes (current/saved/upload),
  editor inline de samplesheet (DT editable), autocompletado a partir de muestras
  detectadas, dispatch a motores, cache en `state$deg_rv`, filtros aplicados al
  vuelo, enriquecimiento bajo demanda.

### Archivos modificados
- `R/state.R` — anadido `state$deg_rv <- reactiveValues(...)`.
- `ui.R` — anadido `nav_panel("4 · Expresion diferencial", value = "tab_deg", uiOutput("tab_deg_content"))`.
- `server.R` — `server_tab_deg(input, output, session, state)` al final.
- `global.R` — flags `HAS_DESEQ2`, `HAS_EDGER`, `HAS_LIMMA`, `HAS_CLUSTERPROFILER`,
  `HAS_ENRICHPLOT`, `HAS_PHEATMAP`, `HAS_ORGECDB`.
- `requirements.sh` — anadidos `edgeR`, `limma` a `required_bioc` y
  `org.EcK12.eg.db` como opcional (avisa pero no falla si no esta).

### Decisiones de implementacion
- **Pureza:** ningun objeto Shiny entra en `utils_deg.R` ni `utils_enrich.R`.
  Solo `data.frame` y matrices.
- **Manejo de errores:** cada `run_deg_*` devuelve `list(table = NULL, error = msg)`.
  El server muestra `showNotification(type="error")` y no asigna a `state$deg_rv`.
- **Cache:** los filtros (FDR / |log2FC| / baseMean) no recorren los motores;
  reactivos derivados (`deg_filtered`) reusan la `state$deg_rv$results`
  cacheada.
- **vst/rlog:** `vst_or_rlog()` elige `vst` si `nrow >= 1000`, `rlog` si menos.
  Si DESeq2 no esta, fallback a log2(CPM + 1).
- **Plotly + pheatmap:** volcano, MA, PCA con plotly. Heatmap top-N y matriz
  de distancias con `pheatmap` (renderPlot) — si pheatmap no esta, `stats::heatmap`.
