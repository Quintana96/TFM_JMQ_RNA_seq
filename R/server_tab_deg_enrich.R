#' server_tab_deg_enrich.R
#' Enriquecimiento funcional y GSEA de la pestana 4 (items 13 y 14, B3).
#'
#' Parte del modulo de la pestana 4, separado de server_tab_deg.R por tamano: el
#' fichero original llego a 1.608 lineas en una unica funcion. NO se usa
#' moduleServer(): igual que el resto de la aplicacion, se conservan los IDs
#' Shiny originales y el estado se pasa explicitamente.
#'
#' `ctx` es el contexto compartido del modulo (un environment, como `state`, para
#' que las asignaciones sean visibles entre las partes). Lo crea server_tab_deg()
#' y contiene los reactivos que varias partes necesitan.

server_tab_deg_enrich <- function(input, output, session, state, ctx) {
  # ── Enriquecimiento ────────────────────────────────────────────────────────
  enrich_rv <- reactiveVal(NULL)
  enrich_mapping_rv <- reactiveVal(NULL)
  # Contexto del ultimo GSEA (ranking y parametros). Lo necesita el running
  # score: la tabla de resultados no contiene el ranking, y recalcularlo al
  # vuelo daria una curva que no corresponde al NES mostrado si el usuario ha
  # tocado la metrica entretanto.
  gsea_ctx_rv <- reactiveVal(NULL)
  compare_rv <- reactiveVal(NULL)
  # El objeto S4 del enriquecimiento (enrichResult / gseaResult). La tabla no
  # basta para los graficos de red: no conserva ni la pertenencia gen-termino ni
  # el ranking. Ver R/utils_enrich_plots.R.
  enrich_obj_rv <- reactiveVal(NULL)
  # Resultado de la traduccion de identificadores del ultimo calculo, para
  # poder decir cuantos genes se perdieron por el camino.
  enrich_traduccion_rv <- reactiveVal(NULL)

  # OrgDb elegido en la interfaz. Estaba cableado a org.EcK12.eg.db, de modo que
  # con datos de otro organismo el enriquecimiento GO no podia funcionar aunque
  # su paquete estuviera instalado: el mapeo salia del 0 % y la app avisaba,
  # correctamente, de algo que el usuario no tenia forma de arreglar.
  deg_orgdb <- reactive({
    sel <- input$deg_orgdb %||% ""
    if (nzchar(sel) && sel %in% ORGDBS_DISPONIBLES) return(sel)
    if (length(ORGDBS_DISPONIBLES)) return(ORGDBS_DISPONIBLES[1])
    NULL
  })

  # keyType seleccionable, poblado con los keytypes reales del OrgDb. Antes
  # estaba fijado a "SYMBOL" en el codigo, lo que con locus tags de E. coli
  # (b0001) daba un mapeo casi nulo sin ningun aviso.
  observe({
    kt <- orgdb_keytypes(deg_orgdb())
    if (!length(kt)) kt <- c("SYMBOL", "ENTREZID", "ENSEMBL", "ALIAS", "REFSEQ")
    preferred <- if ("SYMBOL" %in% kt) "SYMBOL" else kt[1]
    updateSelectInput(session, "deg_go_keytype", choices = kt,
                      selected = isolate(input$deg_go_keytype) %||% preferred)
  })

  # ── Anotacion con la que traducir identificadores ──────────────────────────
  # La matriz de conteos viene con los identificadores que el pipeline produjo
  # (locus tags, si se lanzo con `featureCounts -g locus_tag`), y ningun OrgDb
  # los conoce. La anotacion es lo unico que sabe a que simbolo corresponde cada
  # uno. Se busca donde este: primero la de la ejecucion elegida, y si no la de
  # la sesion o la que se haya subido en la pestana de configuracion.
  deg_annotation_path <- reactive({
    src <- input$deg_source %||% "current"
    if (identical(src, "saved")) {
      sel <- input$selected_deg_run_dir %||% ""
      if (nzchar(sel) && dir.exists(sel)) {
        af <- annotation_file_for_run(sel)
        if (!is.null(af)) return(af)
      }
    }
    # `state` es un environment, no una clase: sus campos pueden faltar. Y la
    # anotacion es opcional por contrato —esta funcion devuelve NULL cuando no
    # hay ninguna—, asi que un estado incompleto tiene que dar NULL y no un
    # error que tumbe el resto de la pestana.
    if (is.function(state$run_params_rv)) {
      af <- state$run_params_rv()$annotation_file %||% ""
      if (nzchar(af) && file.exists(af)) return(af)
    }
    up <- input$annotation_file_upload
    if (!is.null(up) && nrow(up)) return(up$datapath[1])
    NULL
  })

  # Atributos que la anotacion trae realmente: no tiene sentido ofrecer
  # `protein_id` como destino si el GTF no lo lleva.
  deg_annotation_attrs <- reactive({
    path <- deg_annotation_path()
    if (is.null(path)) return(character(0))
    annotation_available_attrs(path, names(ANNOTATION_TARGET_ATTRS))
  })

  observe({
    # Depende tambien de la casilla, y no por capricho: mientras el
    # conditionalPanel esta oculto, updateSelectInput() llega a un elemento que
    # el cliente aun no ha terminado de inicializar y las opciones se pierden.
    # El sintoma es un selector vacio justo despues de marcar la casilla. Al
    # depender de ella, las opciones se reenvian en el momento en que se muestra.
    input$deg_translate_ids
    attrs <- deg_annotation_attrs()
    if (!length(attrs)) {
      updateSelectInput(session, "deg_translate_to", choices = character(0))
      return()
    }
    etiquetas <- ANNOTATION_TARGET_ATTRS[attrs]
    ch <- stats::setNames(attrs, ifelse(is.na(etiquetas), attrs, etiquetas))
    sel <- isolate(input$deg_translate_to)
    updateSelectInput(session, "deg_translate_to", choices = ch,
                      selected = if (!is.null(sel) && sel %in% attrs) sel
                                 else if ("gene" %in% attrs) "gene" else attrs[1])
  })

  # Diagnostico ANTES de calcular: dice si la traduccion hace falta y si va a
  # funcionar. Sin esto, el unico sintoma de unos identificadores equivocados es
  # un enriquecimiento vacio, que no se distingue de la ausencia de senal.
  deg_translate_diag <- reactive({
    path <- deg_annotation_path()
    tab <- state$deg_rv$results
    if (is.null(path) || is.null(tab) || !"gene" %in% names(tab)) return(NULL)
    det <- detect_annotation_keytype(tab$gene, path)
    if (is.null(det)) return(NULL)
    list(path = path, det = det)
  })

  output$deg_translate_estado <- renderUI({
    path <- deg_annotation_path()
    if (is.null(path)) {
      return(div(class = "small text-muted",
                 paste("Sin anotacion disponible. Se toma de la ejecucion elegida,",
                       "o de la que subas en el paso 1.")))
    }
    d <- deg_translate_diag()
    if (is.null(d)) return(div(class = "small text-muted", "Anotacion cargada."))
    if (d$det$rate < 0.5) {
      return(div(class = "alert alert-warning py-1 px-2 small mb-0",
                 icon("triangle-exclamation"), " ",
                 sprintf(paste("La anotacion no reconoce estos identificadores",
                               "(el mejor atributo, '%s', cubre el %.0f %%).",
                               "Comprueba que es la misma anotacion con la que se",
                               "contaron los genes."),
                         d$det$attr, 100 * d$det$rate)))
    }
    div(class = "small text-muted",
        sprintf("Los identificadores de la matriz son '%s' (%.0f %% de cobertura en la anotacion).",
                d$det$attr, 100 * d$det$rate))
  })

  # ── Gene sets propios (GMT) ────────────────────────────────────────────────
  gmt_rv <- reactive({
    f <- input$deg_gmt_file
    if (is.null(f) || !nrow(f)) return(NULL)
    read_gene_sets_gmt(f$datapath[1])
  })

  output$deg_gmt_summary <- renderUI({
    gs <- gmt_rv()
    if (is.null(gs)) {
      return(tags$small(class = "text-muted", "Ningun fichero cargado."))
    }
    if (!is.null(gs$error)) {
      return(div(class = "alert alert-warning py-1 px-2 small mb-0",
                 icon("triangle-exclamation"), " ", gs$error))
    }
    div(class = "small text-muted", tags$b("Conjuntos cargados: "),
        gmt_summary_text(gs))
  })

  # El organismo de Reactome se deduce del OrgDb ya elegido. Declararlo dos veces
  # es una fuente segura de incoherencia: mapear identificadores con el OrgDb de
  # raton y pedir rutas humanas devuelve una tabla plausible y equivocada. Se
  # preselecciona, no se impone: sigue siendo editable.
  observeEvent(deg_orgdb(), ignoreNULL = FALSE, {
    org <- reactome_organism_for_orgdb(deg_orgdb())
    if (is.null(org)) return()
    updateSelectInput(session, "deg_reactome_organism", selected = org)
  })

  # ── Parametros de la interfaz ──────────────────────────────────────────────
  # Un unico punto de lectura para que el boton de calcular y el de comparar no
  # puedan divergir en los parametros que usan.
  num_or <- function(x, default) {
    v <- suppressWarnings(as.numeric(x %||% NA))
    if (length(v) != 1L || is.na(v)) default else v
  }

  enrich_inputs <- function() {
    ont <- input$deg_ontology %||% "BP"
    org_code <- trimws(input$deg_kegg_organism %||% "eco")
    if (!nzchar(org_code)) org_code <- "eco"
    # Reactome tiene su propio catalogo de organismos, con nombres distintos a
    # los codigos de tres letras de KEGG. Se resuelve aqui, en el unico punto de
    # lectura de parametros, para que `org_code` signifique siempre "el organismo
    # de la coleccion seleccionada" y ni el ORA ni el GSEA puedan divergir.
    if (identical(ont, "REACTOME")) {
      org_code <- input$deg_reactome_organism %||% "human"
    }
    gs <- gmt_rv()
    list(
      ont         = ont,
      approach    = input$deg_enrich_approach %||% "ora",
      org_code    = org_code,
      keytype     = if (identical(ont, "KEGG")) input$deg_kegg_keytype %||% "kegg"
                    else input$deg_go_keytype %||% "SYMBOL",
      # El selector de KEGG declara lo que KEGG espera recibir; este declara lo
      # que son los identificadores de la matriz. Son cosas distintas y hasta
      # ahora solo se leia el primero, de modo que no habia forma de decirle a la
      # aplicacion que los genes venian en simbolos.
      from_keytype = input$deg_go_keytype %||% "SYMBOL",
      # setReadable necesita OrgDb: no aplica a KEGG ni a conjuntos propios.
      readable    = isTRUE(input$deg_enrich_readable) && !ont %in% c("KEGG", "GMT"),
      simplify    = isTRUE(input$deg_go_simplify),
      directional = isTRUE(input$deg_ora_directional),
      metric      = input$deg_gsea_metric %||% "stat",
      min_size    = max(2, round(num_or(input$deg_gsea_min_size, 10))),
      max_size    = max(3, round(num_or(input$deg_gsea_max_size, 500))),
      pcut        = min(1, max(1e-4, num_or(input$deg_gsea_pcutoff, 0.05))),
      # eps = 0 pide a fgsea el p-valor exacto del multilevel en lugar de
      # truncarlo en 1e-10, donde todos los conjuntos muy significativos empatan.
      eps         = if (isTRUE(input$deg_gsea_eps_exact)) 0 else 1e-10,
      # Traduccion de identificadores con la anotacion. Va aqui, en el unico
      # punto de lectura, para que el ORA, el GSEA y la comparacion entre ambos
      # no puedan usar espacios de identificadores distintos.
      traducir    = isTRUE(input$deg_translate_ids) && !is.null(deg_annotation_path()),
      annot_path  = deg_annotation_path(),
      translate_to = input$deg_translate_to %||% "gene",
      universe    = ctx$deg_universe(),
      term2gene   = if (is.null(gs)) NULL else gs$term2gene,
      gmt         = gs
    )
  }

  # ── Ejecucion ──────────────────────────────────────────────────────────────
  #' Corre el ORA de la ontologia seleccionada sobre una lista de genes.
  ora_runner <- function(p) {
    function(genes) {
      if (identical(p$ont, "KEGG")) {
        run_enrichment_kegg(genes, universe = p$universe, organism = p$org_code,
                            keyType = p$keytype, OrgDb = deg_orgdb(),
                            from_keytype = p$from_keytype)
      } else if (identical(p$ont, "GMT")) {
        run_enrichment_gmt(genes, universe = p$universe, term2gene = p$term2gene,
                           minGSSize = p$min_size, maxGSSize = p$max_size)
      } else if (identical(p$ont, "REACTOME")) {
        run_enrichment_reactome(genes, universe = p$universe, OrgDb = deg_orgdb(),
                                keyType = p$keytype, organism = p$org_code,
                                minGSSize = p$min_size, maxGSSize = p$max_size,
                                readable = p$readable)
      } else {
        run_enrichment_go(genes, universe = p$universe, OrgDb = deg_orgdb(),
                          ont = p$ont, keyType = p$keytype,
                          simplify_terms = p$simplify, readable = p$readable)
      }
    }
  }

  #' Ejecuta un enfoque (ORA o GSEA) con los parametros dados.
  #' Devuelve el list(table, error, mapping) de utils_enrich mas, en GSEA, el
  #' ranking usado.
  enrich_compute <- function(p, approach) {
    if (identical(p$ont, "GMT") && is.null(p$term2gene)) {
      return(list(table = NULL, mapping = NULL,
                  error = "Carga un fichero GMT con los conjuntos de genes."))
    }
    if (identical(approach, "gsea")) {
      # GSEA parte del ranking completo, sin umbralizar: por eso usa la tabla
      # entera y no la lista filtrada.
      rk <- deg_ranking_metric(state$deg_rv$results, p$metric)
      if (is.null(rk$ranked)) {
        return(list(table = NULL, mapping = NULL, ranking = rk,
                    error = paste0("No se pudo construir el ranking: ",
                                   rk$error %||% "—")))
      }
      # El ranking se traduce entero antes de ordenar conjuntos: traducir
      # despues obligaria a reordenar y el NES dejaria de corresponder a la
      # curva que dibuja el running score.
      if (isTRUE(p$traducir)) {
        trk <- translate_ranking_with_annotation(rk$ranked, p$annot_path,
                                                 to = p$translate_to)
        if (is.null(trk$ranked)) {
          return(list(table = NULL, mapping = trk$mapping, ranking = rk,
                      error = paste0("No se pudieron traducir los identificadores: ",
                                     trk$error %||% "—")))
        }
        rk$ranked <- trk$ranked
        rk$traduccion <- trk
      }
      res <- run_gsea(rk$ranked, ont = p$ont, OrgDb = deg_orgdb(),
                      organism = p$org_code, keyType = p$keytype, exponent = 0,
                      pvalueCutoff = p$pcut, minGSSize = p$min_size,
                      maxGSSize = p$max_size, eps = p$eps,
                      term2gene = p$term2gene, readable = p$readable)
      res$ranking <- rk
      return(res)
    }
    # La lista del ORA sale del FDR con el que se AJUSTO el modelo, no de los
    # filtros de la tarjeta 6: esos son de visualizacion y no recortan el
    # universo, asi que dejarlos entrar aqui cambiaba el enriquecimiento al
    # mover un deslizador declarado cosmetico.
    df <- ctx$deg_significant()
    if (is.null(df) || !nrow(df)) {
      return(list(table = NULL, mapping = NULL,
                  error = paste0("No hay genes significativos a FDR <= ",
                                 state$deg_rv$fdr %||% 0.05, ".")))
    }
    # La lista Y el universo se traducen con las mismas reglas: un fondo en un
    # espacio de identificadores y una lista en otro dan un enriquecimiento sin
    # sentido, no un error.
    traduccion <- NULL
    if (isTRUE(p$traducir)) {
      tr <- translate_ids_with_annotation(unique(as.character(df$gene)),
                                          p$annot_path, to = p$translate_to)
      if (!length(tr$ids)) {
        return(list(table = NULL, mapping = tr$mapping,
                    error = paste0("No se pudieron traducir los identificadores: ",
                                   tr$error %||% "—")))
      }
      traduccion <- tr
      tu <- translate_ids_with_annotation(p$universe, p$annot_path,
                                          from = tr$from, to = p$translate_to)
      p$universe <- if (length(tu$ids)) tu$ids else NULL
      # El ORA direccional parte la tabla en tres, asi que la traduccion tiene
      # que ir en la columna, no en el vector: se reetiqueta df$gene.
      mapa <- annotation_id_map(p$annot_path, tr$from, p$translate_to)
      nuevos <- unname(mapa[as.character(df$gene)])
      df <- df[!is.na(nuevos) & nzchar(nuevos), , drop = FALSE]
      nuevos <- nuevos[!is.na(nuevos) & nzchar(nuevos)]
      df$gene <- nuevos
      df <- df[!duplicated(df$gene), , drop = FALSE]
    }

    runner <- ora_runner(p)
    res <- if (isTRUE(p$directional)) run_ora_directional(df, runner)
           else runner(unique(as.character(df$gene)))
    res$n_lista <- nrow(df)
    # La tasa que importa comunicar es la de la traduccion cuando la hay: es la
    # que explica por que la lista encogio.
    if (!is.null(traduccion)) res$traduccion <- traduccion
    res
  }

  observeEvent(input$deg_run_enrich_btn, {
    req(state$deg_rv$results)
    # El enriquecimiento tarda entre segundos y un minuto largo (mas con
    # simplify(), que calcula similitud semantica entre todos los terminos, con
    # el ORA direccional, que lo corre tres veces, y con eps = 0, que sustituye
    # el p-valor truncado por el exacto). Sin feedback la app parece colgada y el
    # usuario vuelve a pulsar, encolando ejecuciones.
    shinyjs::disable("deg_run_enrich_btn")
    on.exit(shinyjs::enable("deg_run_enrich_btn"), add = TRUE)

    p <- enrich_inputs()
    approach <- p$approach

    if (identical(approach, "gsea")) {
      rk_check <- deg_ranking_metric(state$deg_rv$results, p$metric)
      if (!is.null(rk_check$ranked) && !is.na(rk_check$ties_frac) &&
          rk_check$ties_frac > 0.01) {
        showNotification(
          paste0("El ", round(100 * rk_check$ties_frac, 1), " % de los genes comparte ",
                 "valor en el ranking. GSEA no resuelve los empates, asi que su ",
                 "orden interno es arbitrario; considera usar la metrica 'stat'."),
          type = "warning", duration = 16
        )
      }
    }

    res <- withProgress(
      message = if (identical(approach, "gsea")) "Calculando GSEA..."
                else "Calculando enriquecimiento...",
      value = 0.3, {
        setProgress(value = 0.6,
                    detail = if (identical(approach, "gsea")) "permutaciones"
                             else "test hipergeometrico")
        enrich_compute(p, approach)
      })

    enrich_mapping_rv(res$mapping)
    enrich_traduccion_rv(res$traduccion %||% res$ranking$traduccion)
    # Un mapeo bajo hace el resultado no interpretable, asi que se avisa aunque
    # el enriquecimiento haya devuelto terminos.
    mp <- res$mapping
    if (!is.null(mp) && !is.na(mp$rate %||% NA) && mp$rate < 0.5) {
      showNotification(
        paste0("Atencion: solo se ha mapeado el ", round(100 * mp$rate, 1),
               " % de los genes (", mp$n_mapped, "/", mp$n_input,
               "). Revisa el keyType: los IDs de featureCounts suelen ser locus ",
               "tags, no simbolos."),
        type = "warning", duration = 16
      )
    }

    # Se registra el enriquecimiento AUNQUE no haya dado terminos: un resultado
    # negativo tambien forma parte del analisis, y sin sus parametros (ontologia,
    # keyType, universo, tasa de mapeo) no es interpretable ni reproducible.
    state$deg_rv$enrich <- list(
      enfoque    = if (identical(approach, "gsea")) "GSEA" else "ORA (sobre-representacion)",
      ontologia  = enrich_collection_label(p$ont),
      organismo_kegg = if (identical(p$ont, "KEGG")) p$org_code else NA_character_,
      organismo_reactome = if (identical(p$ont, "REACTOME")) p$org_code else NA_character_,
      orgdb      = if (identical(p$ont, "GMT")) "no aplica (GMT)" else deg_orgdb() %||% "—",
      keytype    = if (identical(p$ont, "GMT")) "IDs del fichero GMT" else p$keytype,
      metrica    = if (identical(approach, "gsea")) p$metric else NA_character_,
      exponent   = if (identical(approach, "gsea")) 0 else NA_real_,
      simplify   = if (identical(approach, "gsea")) NA else p$simplify,
      # Los parametros de GSEA cambian que conjuntos se testean y cuales se
      # devuelven, asi que sin ellos el resultado no es reproducible.
      gsea_min_size = if (identical(approach, "gsea")) p$min_size else NA_integer_,
      gsea_max_size = if (identical(approach, "gsea")) p$max_size else NA_integer_,
      gsea_pcutoff  = if (identical(approach, "gsea")) p$pcut else NA_real_,
      gsea_eps      = if (identical(approach, "gsea")) p$eps else NA_real_,
      direccional = if (identical(approach, "gsea")) NA else p$directional,
      simbolos   = p$readable,
      gmt        = if (identical(p$ont, "GMT") && !is.null(p$gmt))
                     paste0(p$gmt$n_sets, " conjuntos (", min(p$gmt$sizes), "-",
                            max(p$gmt$sizes), " genes)")
                   else NA_character_,
      n_lista    = if (identical(approach, "gsea")) NA_integer_
                   else res$n_lista %||% nrow(ctx$deg_significant() %||% data.frame()),
      n_universo = length(p$universe %||% character(0)),
      # La traduccion cambia que genes entran en el test, asi que forma parte de
      # los parametros del analisis, no de su presentacion.
      traduccion = if (!isTRUE(p$traducir)) NA_character_ else {
        tr <- res$traduccion %||% res$ranking$traduccion
        if (is.null(tr)) "solicitada, sin resultado"
        else sprintf("%s -> %s (%d de %d, %.1f %%; %d colapsados)",
                     tr$mapping$keytype %||% "?", p$translate_to,
                     tr$mapping$n_mapped %||% 0L, tr$mapping$n_input %||% 0L,
                     100 * (tr$mapping$rate %||% 0), tr$mapping$n_colapsados %||% 0L)
      },
      anotacion_traduccion = if (isTRUE(p$traducir)) p$annot_path %||% NA_character_
                             else NA_character_,
      mapeo      = res$mapping,
      n_terminos = if (is.null(res$table)) 0L else nrow(res$table),
      error      = res$error %||% NA_character_,
      # Fechas de las fuentes del OrgDb: los resultados de enriquecimiento
      # cambian con la version de la anotacion (Wadi et al., Nat Methods 2016),
      # asi que sin esto el resultado no es reproducible.
      anotacion  = orgdb_source_info(deg_orgdb())
    )

    # El contexto de GSEA se guarda aunque no haya terminos: es lo que permite
    # que el running score sepa sobre que ranking se calculo todo.
    if (identical(approach, "gsea")) {
      # `ranked_used` solo lo devuelven las colecciones que traducen los IDs
      # (hoy Reactome, que trabaja en ENTREZID). El running score tiene que
      # dibujarse sobre ese mismo ranking o las posiciones no corresponden al NES.
      gsea_ctx_rv(list(ranked = res$ranked_used %||% res$ranking$ranked, ont = p$ont,
                       keytype = p$keytype, orgdb = deg_orgdb(),
                       organism = p$org_code, term2gene = p$term2gene,
                       exponent = 0))
    } else {
      gsea_ctx_rv(NULL)
    }

    if (is.null(res$table)) {
      # Con GSEA, "sin conjuntos enriquecidos" puede ser simplemente que ninguno
      # baja de 0,05: se sugiere el corte a 1 para poder distinguirlo de un fallo.
      extra <- if (identical(approach, "gsea") && p$pcut < 1)
        " Pon el corte de p ajustado a 1 para ver todos los conjuntos testeados." else ""
      showNotification(paste0("Enriquecimiento sin resultados: ", res$error %||% "—", extra),
                       type = "warning", duration = 10)
      enrich_rv(NULL); enrich_obj_rv(NULL); return()
    }
    if (isTRUE(p$directional) && length(res$errores %||% character(0))) {
      showNotification(
        paste0("Direcciones sin terminos: ",
               paste0(names(res$errores), " (", res$errores, ")", collapse = "; ")),
        type = "message", duration = 10)
    }
    enrich_rv(res$table)
    enrich_obj_rv(res$obj)
  })

  output$deg_enrich_traduccion <- renderUI({
    tr <- enrich_traduccion_rv()
    if (is.null(tr) || is.null(tr$mapping)) return(NULL)
    m <- tr$mapping
    perdidos <- (m$n_input %||% 0L) - (m$n_mapped %||% 0L)
    div(class = "small text-muted py-1 px-2 mb-1",
        tags$b("Traduccion: "),
        sprintf("%s -> %s, %d de %d genes (%.1f %%)",
                m$keytype %||% "?", m$source %||% "?", m$n_mapped %||% 0L,
                m$n_input %||% 0L, 100 * (m$rate %||% 0)),
        if (perdidos > 0)
          sprintf(". Se quedaron fuera %d: %d sin correspondencia y %d colapsados en un identificador ya presente.",
                  perdidos, perdidos - (m$n_colapsados %||% 0L), m$n_colapsados %||% 0L)
        else NULL)
  })

  output$deg_enrich_mapping <- renderUI({
    mp <- enrich_mapping_rv()
    txt <- mapping_rate_text(mp)
    if (is.null(txt)) return(NULL)
    low <- !is.null(mp) && !is.na(mp$rate %||% NA) && mp$rate < 0.5
    div(class = paste("small py-1 px-2 mb-1 rounded",
                      if (low) "alert alert-warning mb-2" else "text-muted"),
        if (low) tagList(icon("triangle-exclamation"), " ") else NULL,
        tags$b("Tasa de mapeo: "), txt)
  })

  #' El dotplot sirve a dos tipos de tabla: el ORA trae GeneRatio/Count y se
  #' ordena por significacion; GSEA trae NES/setSize, y ahi lo informativo es el
  #' signo del NES (si el conjunto sube o baja), no solo el p-valor.
  make_enrich_dotplot <- function(df) {
    if (is.null(df) || !nrow(df))
      return(plotly_message("Pulsa 'Calcular' para correr el enriquecimiento."))
    top <- enrichment_dotplot_data(df, top_n = 15)
    if (is.null(top) || !nrow(top))
      return(plotly_message("Sin terminos enriquecidos."))

    if ("NES" %in% names(top)) {
      top$size_val <- top$setSize %||% rep(10, nrow(top))
      return(plotly::plot_ly(
        top,
        x = ~NES,
        y = ~stats::reorder(plot_label, NES),
        type = "scatter", mode = "markers",
        size = ~size_val,
        marker = list(color = ~ifelse(NES > 0, "#7BBF9A", "#F4A6A6"),
                      sizemode = "area"),
        text = ~paste0("Conjunto: ", Description,
                       "<br>NES: ", round(NES, 3),
                       "<br>genes en el conjunto: ", size_val,
                       "<br>p.adjust: ", signif(p.adjust, 3)),
        hoverinfo = "text"
      ) |>
        plotly::layout(
          xaxis = list(title = "NES (score de enriquecimiento normalizado)"),
          yaxis = list(title = "", automargin = TRUE),
          margin = list(l = 220),
          shapes = list(list(type = "line", x0 = 0, x1 = 0, yref = "paper",
                             y0 = 0, y1 = 1,
                             line = list(dash = "dot", color = "#60756A")))
        ))
    }

    # Con ORA direccional el color codifica la direccion: es la lectura que se
    # busca al separar las listas, y sin ella las dos mitades se confunden.
    dir_col <- if ("Direccion" %in% names(top)) {
      unname(c("Conjunto" = "#A8DADC", "Al alza" = "#7BBF9A",
               "A la baja" = "#F4A6A6")[top$Direccion])
    } else NULL

    p <- plotly::plot_ly(
      top,
      x = ~ -log10(p.adjust),
      y = ~stats::reorder(plot_label, -log10(p.adjust)),
      type = "scatter", mode = "markers",
      size = ~Count,
      marker = if (is.null(dir_col))
        list(color = ~ -log10(p.adjust),
             colorscale = list(c(0, "#A8DADC"), c(1, "#7BBF9A")),
             sizemode = "area")
      else list(color = dir_col, sizemode = "area"),
      text = ~paste0("Termino: ", Description,
                     if ("Direccion" %in% names(top)) paste0("<br>Direccion: ", Direccion) else "",
                     "<br>Count: ", Count,
                     "<br>p.adjust: ", signif(p.adjust, 3)),
      hoverinfo = "text"
    )
    p |>
      plotly::layout(
        xaxis = list(title = "-log10(p.adjust)"),
        yaxis = list(title = "", automargin = TRUE),
        margin = list(l = 220)
      )
  }

  output$deg_enrich_dotplot <- plotly::renderPlotly({
    make_enrich_dotplot(enrich_rv())
  })

  output$deg_enrich_table <- renderDT({
    df <- enrich_rv()
    if (is.null(df) || !nrow(df)) return(dt_table(message_df("Sin enriquecimiento calculado.")))
    df_r <- df
    for (nm in intersect(c("pvalue", "p.adjust", "qvalue"), names(df_r)))
      df_r[[nm]] <- signif(df_r[[nm]], 3)
    for (nm in intersect(c("enrichmentScore", "NES"), names(df_r)))
      df_r[[nm]] <- round(df_r[[nm]], 3)
    dt_table(df_r, page_length = 15, filter = "top")
  })

  output$download_enrich_table <- csv_download(
    "enriquecimiento",
    function() {
      df <- enrich_rv()
      if (is.null(df)) return(message_df("Sin enriquecimiento calculado."))
      df
    }
  )


  # ── Mapas y redes de terminos ──────────────────────────────────────────────
  # Completan el esquema del pipeline: el dotplot resume, pero no deja ver que
  # los terminos comparten genes entre si. Son ggplot2, no plotly (ver la nota
  # de cabecera de R/utils_enrich_plots.R).

  # log2FC con el que colorear los genes de la red. Sale de la tabla DEG, asi
  # que sus nombres son los identificadores de la matriz de conteos; si el
  # enriquecimiento se corrio en otro espacio (keyType distinto, o readable =
  # TRUE traduciendo a simbolos) no casaran, y en vez de un grafico mudo se
  # avisa de cuantos genes han quedado sin color.
  enrich_fold_change <- reactive({
    tab <- state$deg_rv$results
    if (is.null(tab) || !all(c("gene", "log2FC") %in% names(tab))) return(NULL)
    v <- stats::setNames(tab$log2FC, as.character(tab$gene))
    v <- v[!is.na(v) & nzchar(names(v))]
    # Si el enriquecimiento se corrio traduciendo, el objeto guarda los genes ya
    # traducidos: un vector de log2FC con los identificadores de la matriz no
    # casaria con ninguno y la red saldria en gris.
    p <- enrich_inputs()
    if (isTRUE(p$traducir)) {
      tv <- translate_ranking_with_annotation(v, p$annot_path, to = p$translate_to)
      if (!is.null(tv$ranked)) v <- tv$ranked
    }
    v
  })

  observe({
    obj <- enrich_obj_rv()
    ch <- enrich_plot_choices(obj)
    if (!length(ch)) {
      updateSelectInput(session, "deg_enrich_plot_tipo", choices = character(0))
      return()
    }
    actual <- isolate(input$deg_enrich_plot_tipo)
    updateSelectInput(session, "deg_enrich_plot_tipo", choices = ch,
                      selected = if (!is.null(actual) && actual %in% ch) actual else ch[[1]])
  })

  output$deg_enrich_plot_ayuda <- renderUI({
    tipo <- input$deg_enrich_plot_tipo
    if (is.null(tipo) || !nzchar(tipo %||% "")) return(NULL)
    txt <- enrich_plot_ayuda(tipo)
    if (!nzchar(txt)) return(NULL)
    tags$p(class = "small text-muted mb-2", txt)
  })

  #' Construye el grafico una sola vez para la vista y para la descarga: si se
  #' calculase dos veces, el emapplot —que reordena por similitud— podria salir
  #' con una disposicion distinta en el PNG que en pantalla.
  enrich_netplot_rv <- reactive({
    obj <- enrich_obj_rv()
    tipo <- input$deg_enrich_plot_tipo
    if (is.null(obj) || is.null(tipo) || !nzchar(tipo %||% "")) return(NULL)
    n <- max(2, round(as.numeric(input$deg_enrich_plot_n %||% 15)))
    fc <- if (identical(tipo, "cnet")) enrich_fold_change() else NULL
    tryCatch(
      list(plot = enrich_make_network_plot(obj, tipo, top_n = n, fold_change = fc,
                                           etiquetar_genes = isTRUE(input$deg_cnet_genes)),
           error = NULL),
      error = function(e) list(plot = NULL, error = conditionMessage(e))
    )
  })

  output$deg_enrich_netplot <- renderPlot({
    r <- enrich_netplot_rv()
    validate(need(!is.null(r), "Corre un enriquecimiento para ver los mapas."))
    validate(need(is.null(r$error), r$error))
    r$plot
  }, res = 110)

  output$deg_enrich_netplot_aviso <- renderUI({
    if (!identical(input$deg_enrich_plot_tipo, "cnet")) return(NULL)
    obj <- enrich_obj_rv(); fc <- enrich_fold_change()
    if (is.null(obj) || is.null(fc)) return(NULL)
    df <- tryCatch(as.data.frame(obj), error = function(e) NULL)
    col <- if (is.null(df)) NULL else if ("geneID" %in% names(df)) "geneID"
           else if ("core_enrichment" %in% names(df)) "core_enrichment" else NULL
    if (is.null(col)) return(NULL)
    genes <- unique(unlist(strsplit(as.character(df[[col]]), "/")))
    if (!length(genes)) return(NULL)
    tasa <- mean(genes %in% names(fc))
    if (tasa >= 0.5) return(NULL)
    div(class = "alert alert-warning py-1 px-2 small mb-2",
        icon("triangle-exclamation"), " ",
        sprintf(paste("Solo el %.0f %% de los genes del enriquecimiento tiene log2FC",
                      "asignable, asi que la mayoria quedan sin color. Ocurre cuando",
                      "el enriquecimiento traduce los identificadores (keyType o",
                      "'mostrar simbolos') y la tabla DEG usa los de la matriz."),
                100 * tasa))
  })

  output$download_enrich_netplot <- downloadHandler(
    filename = function() {
      paste0("enriquecimiento_", input$deg_enrich_plot_tipo %||% "grafico", "_",
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      r <- enrich_netplot_rv()
      if (is.null(r) || is.null(r$plot)) stop(r$error %||% "No hay grafico que descargar.")
      # 300 ppp: es lo que necesita una figura impresa, y el renderPlot de
      # pantalla se queda muy por debajo.
      ggplot2::ggsave(file, r$plot, width = 11, height = 7.5, dpi = 300, bg = "white")
    }
  )

  # ── Diagrama de ruta KEGG ──────────────────────────────────────────────────
  # Solo tiene sentido sobre un enriquecimiento KEGG: pathview necesita el
  # identificador de ruta (eco00020) para descargar el diagrama oficial.
  kegg_rutas <- reactive({
    df <- enrich_rv()
    if (is.null(df) || !nrow(df) || !"ID" %in% names(df)) return(character(0))
    ids <- as.character(df$ID)
    ok <- grepl("^[a-z]{3,4}[0-9]{5}$", ids)
    if (!any(ok)) return(character(0))
    stats::setNames(ids[ok], paste0(df$Description[ok], "  (", ids[ok], ")"))
  })

  observe({
    ch <- kegg_rutas()
    updateSelectInput(session, "deg_kegg_pathway", choices = ch,
                      selected = if (length(ch)) ch[[1]] else NULL)
  })

  output$deg_kegg_pathview_estado <- renderUI({
    if (!length(kegg_rutas())) {
      return(div(class = "alert alert-light border py-2 px-2 small mb-2",
                 paste("El diagrama necesita un enriquecimiento KEGG. Elige KEGG como",
                       "coleccion y vuelve a calcular; con GO o Reactome no hay",
                       "identificador de ruta que descargar.")))
    }
    NULL
  })

  kegg_pathview_rv <- eventReactive(input$deg_kegg_pathview_btn, {
    pid <- input$deg_kegg_pathway
    req(pid)
    fc <- enrich_fold_change()
    p <- enrich_inputs()
    if (is.null(fc) || !length(fc)) {
      return(list(path = NULL, error = "No hay log2FC que pintar sobre la ruta."))
    }
    # El diagrama se colorea con los log2FC, que vienen con los identificadores
    # de la matriz. pathview no los traduce: hay que darselos ya en el espacio
    # que declara `gene.idtype`, o pinta un diagrama vacio con "no ID can be
    # mapped". Es el mismo paso que hace el enriquecimiento KEGG.
    idtype <- if (identical(p$keytype, "kegg")) "KEGG" else "ENTREZ"
    if (identical(idtype, "ENTREZ") && !identical(p$from_keytype, "ENTREZID")) {
      tr <- translate_to_entrez(names(fc), deg_orgdb(), keyType = p$from_keytype)
      if (!length(tr$ids)) {
        return(list(path = NULL, error = paste0(
          "No se pudo traducir ningun gen de '", p$from_keytype, "' a ENTREZID: ",
          tr$error %||% "sin coincidencias.")))
      }
      # tr$back va de ENTREZID al identificador original: se invierte para
      # reetiquetar el vector de log2FC sin perder el orden ni los valores.
      originales <- unname(tr$back[tr$ids])
      fc <- stats::setNames(fc[originales], tr$ids)
      fc <- fc[!is.na(fc)]
    }
    withProgress(message = "Descargando el diagrama de KEGG...", value = 0.4, {
      enrich_pathview_png(pid, fc, species = p$org_code, gene_idtype = idtype)
    })
  })

  output$deg_kegg_pathview <- renderImage({
    r <- kegg_pathview_rv()
    validate(need(is.null(r$error), r$error))
    list(src = r$path, contentType = "image/png",
         width = "100%", alt = "Diagrama de la ruta KEGG con los log2FC")
  }, deleteFile = FALSE)

  # ── Running score / leading edge ───────────────────────────────────────────
  # Solo tiene sentido sobre una tabla de GSEA (la que trae NES).
  observe({
    df <- enrich_rv()
    if (is.null(df) || !nrow(df) || !"NES" %in% names(df)) {
      updateSelectInput(session, "deg_gsea_term", choices = character(0))
      return()
    }
    ord <- order(df$p.adjust, na.last = TRUE)
    ids <- as.character(df$ID[ord])
    etiquetas <- paste0(df$Description[ord], "  (NES ", round(df$NES[ord], 2),
                        ", p.adj ", signif(df$p.adjust[ord], 2), ")")
    updateSelectInput(session, "deg_gsea_term",
                      choices = stats::setNames(ids, etiquetas),
                      selected = isolate(input$deg_gsea_term) %||% ids[1])
  })

  gsea_runscore_rv <- reactive({
    gctx <- gsea_ctx_rv()
    term <- input$deg_gsea_term
    if (is.null(gctx) || is.null(term) || !nzchar(term %||% "")) return(NULL)
    tg <- gsea_term_genes(term, ont = gctx$ont, OrgDb = gctx$orgdb,
                          keyType = gctx$keytype, organism = gctx$organism,
                          term2gene = gctx$term2gene)
    if (!length(tg$genes)) {
      return(list(error = tg$error %||% "No se han podido recuperar los genes del conjunto."))
    }
    # gseaParam = exponent: la curva debe ser la del estadistico con el que se
    # calculo el NES, no la ponderada por defecto de fgsea.
    gsea_running_score(gctx$ranked, tg$genes, gseaParam = gctx$exponent %||% 0)
  })

  output$deg_gsea_runscore <- plotly::renderPlotly({
    if (is.null(gsea_ctx_rv()))
      return(plotly_message("Corre un GSEA para ver el running score."))
    rs <- gsea_runscore_rv()
    if (is.null(rs)) return(plotly_message("Selecciona un conjunto."))
    if (!is.null(rs$error)) return(plotly_message(rs$error))
    cur <- rs$curve
    base <- min(cur$ES, 0) - 0.08 * max(abs(range(cur$ES)), 0.1)
    plotly::plot_ly() |>
      plotly::add_lines(x = cur$rank, y = cur$ES, line = list(color = "#3F7D5C"),
                        name = "Running score", hoverinfo = "x+y") |>
      # Cada marca es un gen del conjunto en su posicion del ranking: es lo que
      # deja ver si el conjunto se concentra en un extremo o esta repartido.
      plotly::add_markers(x = rs$ticks$rank, y = rep(base, nrow(rs$ticks)),
                          marker = list(symbol = "line-ns-open", color = "#60756A",
                                        size = 8),
                          name = "Genes del conjunto", hoverinfo = "x") |>
      plotly::layout(
        xaxis = list(title = "Posicion en el ranking"),
        yaxis = list(title = "Enrichment score acumulado"),
        showlegend = FALSE,
        shapes = list(list(type = "line", x0 = 0, x1 = rs$n_ranked, y0 = 0, y1 = 0,
                           line = list(dash = "dot", color = "#B0BDB5")))
      )
  })

  output$deg_gsea_leading <- renderUI({
    df <- enrich_rv()
    term <- input$deg_gsea_term
    if (is.null(df) || !"core_enrichment" %in% names(df) ||
        is.null(term) || !nzchar(term %||% "")) return(NULL)
    fila <- df[!is.na(df$ID) & df$ID == term, , drop = FALSE]
    if (!nrow(fila)) return(NULL)
    genes <- leading_edge_genes(fila$core_enrichment[1])
    rs <- gsea_runscore_rv()
    tagList(
      div(class = "small text-muted mt-2",
          tags$b("Conjunto: "), fila$Description[1],
          " — NES ", round(fila$NES[1], 3),
          ", p ajustado ", signif(fila$p.adjust[1], 3),
          ", ", fila$setSize[1], " genes",
          if (!is.null(rs) && is.null(rs$error))
            paste0(" (", rs$n_hits, " presentes en el ranking)") else ""),
      div(class = "small mt-1",
          tags$b(paste0("Leading edge (", length(genes), " genes): ")),
          tags$span(class = "text-muted",
                    "el subconjunto que va desde el extremo del ranking hasta el pico ",
                    "de la curva y que sostiene el enriquecimiento.")),
      div(class = "small mt-1 p-2 border rounded",
          style = "max-height:150px;overflow-y:auto;font-family:monospace;",
          paste(genes, collapse = ", "))
    )
  })

  # ── ORA frente a GSEA ──────────────────────────────────────────────────────
  observeEvent(input$deg_run_enrich_compare_btn, {
    req(state$deg_rv$results)
    if (!isTRUE(HAS_FGSEA)) {
      showNotification("fgsea no esta instalado: no se puede correr el GSEA.",
                       type = "error"); return()
    }
    shinyjs::disable("deg_run_enrich_compare_btn")
    on.exit(shinyjs::enable("deg_run_enrich_compare_btn"), add = TRUE)

    p <- enrich_inputs()
    # La comparacion exige que las dos ramas vean lo mismo: misma ontologia,
    # mismo universo y misma lista. El ORA direccional partiria la lista en tres
    # y dejaria de ser comparable con un unico GSEA.
    p$directional <- FALSE

    res <- withProgress(message = "Comparando ORA y GSEA...", value = 0.2, {
      ora <- enrich_compute(p, "ora")
      setProgress(value = 0.6, detail = "GSEA (permutaciones)")
      gsea <- enrich_compute(p, "gsea")
      list(ora = ora, gsea = gsea)
    })

    if (is.null(res$ora$table) && is.null(res$gsea$table)) {
      compare_rv(NULL)
      showNotification(paste0("Ninguno de los dos enfoques devolvio terminos. ORA: ",
                              res$ora$error %||% "—", " | GSEA: ",
                              res$gsea$error %||% "—"),
                       type = "warning", duration = 12)
      return()
    }
    cmp <- compare_ora_gsea(res$ora$table, res$gsea$table, padj_cutoff = 0.05)
    cmp$ont <- p$ont
    cmp$error_ora <- res$ora$error
    cmp$error_gsea <- res$gsea$error
    compare_rv(cmp)
  })

  output$deg_enrich_compare_summary <- renderUI({
    cmp <- compare_rv()
    if (is.null(cmp)) {
      return(div(class = "small text-muted",
                 "Pulsa 'Comparar ORA y GSEA' para correr los dos enfoques."))
    }
    jac <- if (is.na(cmp$jaccard)) "—" else paste0(round(100 * cmp$jaccard, 1), " %")
    tagList(
      div(class = "alert alert-secondary py-2 px-2 small mb-2",
          tags$b("Terminos significativos (p ajustado <= 0,05): "),
          "ORA ", fmt_int(cmp$n_ora), " | GSEA ", fmt_int(cmp$n_gsea),
          " | comunes ", fmt_int(cmp$n_comun),
          " | indice de Jaccard ", jac, ".",
          tags$br(),
          if (cmp$n_gsea > cmp$n_comun)
            paste0("GSEA ve ", cmp$n_gsea - cmp$n_comun, " terminos que el ORA no ve: ",
                   "conjuntos desplazados en el ranking cuyos genes no llegan al corte ",
                   "de la lista.")
          else "GSEA no anade ningun termino sobre los del ORA."),
      if (!is.na(cmp$error_ora %||% NA))
        div(class = "small text-muted", tags$b("ORA: "), cmp$error_ora) else NULL,
      if (!is.na(cmp$error_gsea %||% NA))
        div(class = "small text-muted", tags$b("GSEA: "), cmp$error_gsea) else NULL
    )
  })

  output$deg_enrich_compare_plot <- plotly::renderPlotly({
    cmp <- compare_rv()
    if (is.null(cmp)) return(plotly_message("Sin comparacion calculada."))
    d <- data.frame(
      grupo = factor(c("Solo ORA", "Comunes", "Solo GSEA"),
                     levels = c("Solo ORA", "Comunes", "Solo GSEA")),
      n = c(cmp$n_ora - cmp$n_comun, cmp$n_comun, cmp$n_gsea - cmp$n_comun)
    )
    plotly::plot_ly(d, x = ~grupo, y = ~n, type = "bar",
                    marker = list(color = c("#A8DADC", "#7BBF9A", "#F4C68A")),
                    text = ~n, textposition = "outside", hoverinfo = "x+y") |>
      plotly::layout(xaxis = list(title = ""),
                     yaxis = list(title = "Terminos significativos"),
                     showlegend = FALSE)
  })

  compare_table_rv <- reactive({
    cmp <- compare_rv()
    if (is.null(cmp)) return(NULL)
    tb <- compare_ora_gsea_table(cmp)
    if (is.null(tb)) return(NULL)
    for (nm in c("padj_ORA", "padj_GSEA")) tb[[nm]] <- signif(tb[[nm]], 3)
    tb$NES <- round(tb$NES, 3)
    tb
  })

  output$deg_enrich_compare_table <- renderDT({
    tb <- compare_table_rv()
    if (is.null(tb)) return(dt_table(message_df("Sin comparacion calculada.")))
    dt_table(tb, page_length = 15, filter = "top")
  })

  output$download_enrich_compare <- csv_download(
    "comparacion_ora_gsea",
    function() {
      tb <- compare_table_rv()
      if (is.null(tb)) return(message_df("Sin comparacion calculada."))
      tb
    }
  )

  invisible(NULL)
}
