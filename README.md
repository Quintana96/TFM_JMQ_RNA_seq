# RNA-seq Workflow Runner — refactor modular

Refactorizacion de la aplicacion Shiny original (`app.R`, 3122 lineas) en una
estructura multi-archivo lista para ejecutarse con `shiny::runApp()`. Mantiene
la funcionalidad end-to-end, los mismos IDs de inputs/outputs y la misma
paleta pastel.

## Como ejecutar

Desde el directorio que contiene `workflow.sh` y la carpeta `rnaseq_app/`:

```r
shiny::runApp("rnaseq_app")
```

O directamente:

```sh
Rscript rnaseq_app/app.R
```

Shiny actualiza automaticamente, en este orden:

1. `global.R` (librerias, opciones, constantes, tema, CSS).
2. Todos los `R/*.R` (helpers de utilidades + helpers de UI + modulos server).
3. `ui.R` y `server.R`.
4. `app.R` arranca `shinyApp(ui, server)`.

## Estructura

```
rnaseq_app/
├── app.R                # Entrypoint minimo
├── global.R             # Librerias, opciones, constantes FASTQ_*, qc_thresholds, app_theme, app_css, %||%
├── ui.R                 # page_navbar() con las 3 pestanas
├── server.R             # crea state y delega en server_tab_*
├── README.md            # Este archivo
└── R/                   # Auto-sourceado por Shiny
    ├── state.R                  # create_app_state(): reactiveVal/reactiveValues compartidos
    ├── utils_format.R           # fmt_bytes, fmt_elapsed, pct_label, ts_log, terminal_text, trim_log_text, num_or_na
    ├── utils_fastq.R            # sample_fastq_paths, detect_samples, missing_r2, sample_fastq_sizes, etc.
    ├── utils_io.R               # read_tail_text, prepare_uploaded_input_file, outputs_base_dir, list_result_dirs, file_table_*, dir_size, read_tsv_safe
    ├── utils_counts.R           # load_quant_counts, load_count_matrix_tsv, load_counts_from_workflow
    ├── utils_multiqc.R          # multiqc_dir, *_stats, summarise_result, result_general_table, round_numeric_columns
    ├── utils_status.R           # status_from_log, status_badge, infer_result_params
    ├── utils_tables.R           # dt_table, fastqc_table, alignment_table, counts_tables, important_artifacts, log_tail_text
    ├── utils_qc.R               # rate_fraction, find_metric_col, message_df, has_real_rows, plotly_message,
    │                            # featurecounts_assignment_summary, align_qc_summary, pseudo_qc_summary, alerts, ...
    ├── ui_helpers.R             # download_header, csv_download, plotly_download
    ├── ui_tab_config.R          # ui_tab_config() — Tab 1 (6 cards)
    ├── ui_tab_processing.R      # ui_tab_processing_content(cfg, total_sz), ui_tab_processing_locked()
    ├── ui_tab_results.R         # ui_tab_results_content(s), ui_tab_results_empty()
    ├── server_tab_config.R      # server_tab_config(input, output, session, state)
    ├── server_tab_processing.R  # server_tab_processing(input, output, session, state)
    └── server_tab_results.R     # server_tab_results(input, output, session, state)
```

## Decisiones de diseno

- **No se usan `moduleServer()` / `NS()`**. Las tres pestanas comparten mucho
  estado (logs, configuracion validada, parametros de ejecucion, lista de
  resultados...) y los IDs de inputs/outputs ya estaban documentados en la
  app original. En lugar de renombrarlos todos con namespaces, se pasa un
  objeto `state` con todos los reactiveVal/reactiveValues a las funciones
  server de cada pestana, manteniendo el namespace global de Shiny.
- **`state` se crea en `create_app_state(session)`** y expone tambien rutas
  (`workflow_path`, `outputs_dir`, `roots`, `roots_results`) y un slot
  `state$shared` con los reactivos derivados de Tab 1 que necesitan los otros
  tabs (`effective_tool`, `workflow_cmd`, `samples_eff`, etc.).

## Mejoras aplicadas

### Modularizacion
- Separacion clara entre utilidades puras (sin Shiny) en `utils_*.R` y
  modulos Shiny (`server_tab_*.R`, `ui_tab_*.R`).
- `R/` se sourcea automaticamente; no hay `source()` manual.
- Cabeceras `# ╔════╗` ASCII reemplazadas por cabeceras concisas tipo
  roxygen (`#' @section ...`, `#' Descripcion breve`).
- Docstrings ligeras en cada funcion no trivial.

### Optimizaciones de rendimiento
- **Cache de muestras detectadas**: un unico reactivo `samples_eff()` que
  llama a `detect_samples(input_dir_debounced(), input$read_type)` se reusa
  en `val_errors()`, `checklist_status()` y en el observer de
  `btn_to_processing`. Antes se llamaba 3 veces por cambio de input.
- **`bindCache` opcional** en `selected_result_summary`,
  `selected_counts_tables` y `selected_general_table`, con clave string
  `dir|timestamp`. Si la version de Shiny no expone `bindCache`, el reactivo
  cae al modo normal sin cache (chequeo `exists("bindCache", ...)`).
- **`debounce` 600 ms** sobre `input_dir_val()` mantenido.
- **Lectura tail-only del log** (`read_tail_text`, `sync_run_log_file`)
  mantenida; el log completo se persiste en disco
  (`workflow_live.log`).
- **Refresh de archivos centralizado**: el observer de `refresh_btn` ahora
  llama a `file_table_for_dir(out_dir)` (helper unico) en lugar de duplicar
  la logica `list.files(recursive=TRUE)`.

### Limpieza de codigo
- Eliminados comentarios `## reproducibilidad eliminada` y marcadores
  obsoletos `[v3-PROC]`, `[v3-PROG]`, `[v3-LAYOUT]`, `[v3-LOAD]`.
- Unificadas comillas dobles.
- Unificadas las dos ramas casi identicas del observer de `run_btn`
  (processx vs system2) extrayendo la post-ejecucion en
  `finalize_run(exit_code, output_dir)` — antes ambas ramas duplicaban la
  carga de matriz, el showNotification y el `updateNavbarPage`.
- Helper interno `log_line(msg)` que envuelve
  `append_run_log(paste0(ts_log(msg), "\\n"))`.

### Mejoras UI/UX
- En **Tab 2** bloqueada (sin configuracion validada) se anade un boton
  "Ir a configuracion" (`btn_goto_config`) ademas del aviso original.
- En **Tab 3** sin resultados, el aviso vacio incluye un boton
  "Refrescar lista de outputs/" que dispara `refresh_results_btn`
  (anteriormente solo se rendeaba dentro del selector).
- Tema `bs_theme` y paleta pastel intactos.

## Compatibilidad y limitaciones

- El IDs de Shiny se preservan al 100%: las 19 entradas (`input$*`) y las 63
  salidas (`output$*`) coinciden 1:1 con el original. Solo se anade el ID
  nuevo `input$btn_goto_config` para la mejora UX descrita.
- La firma del comando bash y los nombres de los flags pasados a
  `workflow.sh` se mantienen sin cambios.
- `bindCache` requiere Shiny >= 1.6. Se detecta en tiempo de carga; si no
  esta disponible, los reactivos siguen funcionando sin cache (sin afectar
  la correccion, solo el rendimiento de reentradas a Tab 3).

## TODO / notas

- El bloque `pkg_ok` se mantiene en `global.R` por compatibilidad. Si la app
  se distribuye via renv o un Dockerfile fijo, esa deteccion en runtime
  podria sustituirse por dependencias explicitas.
- Las dos funciones `align_qc_region_plot_obj` y
  `align_qc_gene_body_plot_obj` solo devuelven un mensaje placeholder
  (igual que el original); si en el futuro hay metricas de region o gene
  body, basta con poblar esas funciones en `server_tab_results.R`.
