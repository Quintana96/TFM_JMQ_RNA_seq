#!/bin/bash

# -E (errtrace) hace que el trap ERR se herede en funciones y subshells. Sin el,
# un fallo dentro de run_cmd o del pipe bowtie2|samtools terminaba el script
# (correcto) pero sin el mensaje que dice en que línea, que es justo lo que se
# necesita para diagnosticarlo.
set -Eeuo pipefail
shopt -s nullglob

# ============================================================================
# E. coli RNA-seq Pipeline
# ============================================================================
# Requirements:
# Alignment options:
#   - bowtie2: Full alignment (requires samtools, featureCounts)
#   - salmon: Pseudo-alignment (faster, produces direct counts)
#   - kallisto: Pseudo-alignment (faster, produces direct counts)
# Quality control:
#   - fastqc: Quality control
#   - fastp: Trimming, filtering and quality control
#   - multiqc: Report generation
#
# Installation via conda/bioconda:
# conda create -n pipeline_ecoli -c bioconda bowtie2 samtools fastqc fastp subread multiqc salmon kallisto
# conda activate pipeline_ecoli
# ============================================================================

# Usage: workflow.sh --INPUT <path> --OUTPUT <path> --GENOME_FILE <file> --ANNOTATION_FILE <file> [--ALIGNMENT_TYPE bowtie2|salmon|kallisto] [--READ_TYPE pe|se] [--FRAGMENT_LENGTH <mean>] [--FRAGMENT_SD <sd>] [--INFERENTIAL_REPS <n>] [--FEATURE_TYPE gene|CDS] [--FEATURE_ATTR <attr>]

# Default values
ALIGNMENT_TYPE="bowtie2"
READ_TYPE="pe"
FRAGMENT_LENGTH=200
FRAGMENT_SD=20
# Replicas inferenciales de la cuantificación (salmon: --numGibbsSamples,
# kallisto: -b). Las necesita Swish para propagar la incertidumbre de asignación
# de lecturas entre transcritos que comparten secuencia. 0 las desactiva.
INFERENTIAL_REPS=20
# Tipo de feature que cuenta featureCounts (-t) y atributo que agrupa (-g).
# El default se mantiene en gene/locus_tag para no cambiar el comportamiento
# existente. Las buenas prácticas para procariotas (Genome Biology 2021) hacen el
# análisis diferencial a nivel de CDS; con FEATURE_TYPE=CDS se obtiene eso, a
# cambio de excluir los genes no codificantes (rRNA, tRNA) del recuento.
FEATURE_TYPE="gene"
FEATURE_ATTR="locus_tag"
# Orientación de la libreria para featureCounts (-s): 0 sin orientar, 1 directa,
# 2 inversa (el caso de los protocolos dUTP, que son la mayoria hoy). "auto"
# la infiere de los propios datos.
STRANDEDNESS="auto"
# Hilos. Antes estaba cableado a 8, lo que sobresuscribe maquinas más pequeñas
# y no quedaba registrado en ninguna parte.
THREADS_ARG=""

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --INPUT) INPUT="$2"; shift 2 ;;
        --OUTPUT) OUTPUT="$2"; shift 2 ;;
        --GENOME_FILE) GENOME_FILE="$2"; shift 2 ;;
        --ANNOTATION_FILE) ANNOTATION_FILE="$2"; shift 2 ;;
        --ALIGNMENT_TYPE) ALIGNMENT_TYPE="$2"; shift 2 ;;
        --READ_TYPE) READ_TYPE="$2"; shift 2 ;;
        --FRAGMENT_LENGTH) FRAGMENT_LENGTH="$2"; shift 2 ;;
        --FRAGMENT_SD) FRAGMENT_SD="$2"; shift 2 ;;
        --INFERENTIAL_REPS) INFERENTIAL_REPS="$2"; shift 2 ;;
        --FEATURE_TYPE) FEATURE_TYPE="$2"; shift 2 ;;
        --FEATURE_ATTR) FEATURE_ATTR="$2"; shift 2 ;;
        --STRANDEDNESS) STRANDEDNESS="$2"; shift 2 ;;
        --THREADS) THREADS_ARG="$2"; shift 2 ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
done

# Validate required arguments
if [[ -z "${INPUT:-}" || -z "${OUTPUT:-}" || -z "${GENOME_FILE:-}" || -z "${ANNOTATION_FILE:-}" ]]; then
    echo "Error: Missing required arguments."
    echo "Usage: pipeline_ecoli.sh --INPUT <path> --OUTPUT <path> --GENOME_FILE <file> --ANNOTATION_FILE <file> [--ALIGNMENT_TYPE bowtie2|salmon|kallisto] [--READ_TYPE pe|se] [--FRAGMENT_LENGTH <mean>] [--FRAGMENT_SD <sd>] [--INFERENTIAL_REPS <n>] [--FEATURE_TYPE gene|CDS] [--FEATURE_ATTR <attr>]"
    exit 1
fi

# Validate alignment type
if [[ ! "$ALIGNMENT_TYPE" =~ ^(bowtie2|subjunc|salmon|kallisto)$ ]]; then
    echo "Error: ALIGNMENT_TYPE no válido. Opciones: bowtie2, subjunc, salmon, kallisto"
    exit 1
fi

if [[ ! "$READ_TYPE" =~ ^(pe|se)$ ]]; then
    echo "Error: Invalid READ_TYPE. Choose: pe or se"
    exit 1
fi

if [[ ! "$FRAGMENT_LENGTH" =~ ^[0-9]+([.][0-9]+)?$ || ! "$FRAGMENT_SD" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo "Error: FRAGMENT_LENGTH and FRAGMENT_SD must be positive numeric values"
    exit 1
fi

if [[ ! "$INFERENTIAL_REPS" =~ ^[0-9]+$ ]]; then
    echo "Error: INFERENTIAL_REPS must be a non-negative integer"
    exit 1
fi

if [[ -z "$FEATURE_TYPE" || -z "$FEATURE_ATTR" ]]; then
    echo "Error: FEATURE_TYPE and FEATURE_ATTR must not be empty"
    exit 1
fi

# Número de hilos: el indicado, o los nucleos disponibles menos uno para no
# dejar la maquina sin capacidad de respuesta.
if [[ -n "$THREADS_ARG" ]]; then
    THREADS="$THREADS_ARG"
else
    detected=$( (command -v nproc >/dev/null 2>&1 && nproc) || sysctl -n hw.ncpu 2>/dev/null || echo 4 )
    THREADS=$(( detected > 1 ? detected - 1 : 1 ))
fi
if [[ ! "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]]; then
    echo "Error: THREADS debe ser un entero positivo"; exit 1
fi
if [[ ! "$STRANDEDNESS" =~ ^(auto|0|1|2)$ ]]; then
    echo "Error: STRANDEDNESS debe ser auto, 0, 1 o 2"; exit 1
fi

# Directorio del propio script, para localizar scripts/ auxiliares con
# independencia de desde donde se invoque el workflow.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

log() {
    printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"
}

run_cmd() {
    log "+ $*"
    "$@"
}

# Helpers portables: macOS trae `md5` y `stat -f`, Linux `md5sum` y `stat -c`.
md5_of() {
    [[ -f "$1" ]] || { echo "—"; return 0; }
    if command -v md5sum >/dev/null 2>&1; then
        md5sum "$1" | awk '{print $1}'
    elif command -v md5 >/dev/null 2>&1; then
        md5 -q "$1"
    else
        echo "—"
    fi
}

file_size() {
    [[ -e "$1" ]] || { echo 0; return 0; }
    stat -f%z "$1" 2>/dev/null || stat -c%s "$1" 2>/dev/null || echo 0
}


# ── Metricas de ejecución: tiempo y memoria ────────────────────────────────
#
# Responden a la pregunta práctica de "cuanto tarda y cuanta RAM hace falta",
# que es lo que hay que saber para decidir si esto corre en un portatil o
# necesita un servidor.
#
# La memoria se mide con un MUESTREADOR en segundo plano y no con
# `/usr/bin/time`, por dos razones. La primera es de portabilidad: la versión
# BSD de macOS y la GNU de Linux difieren en las opciones (-l frente a -v), en
# las unidades (bytes frente a kilobytes) y en si aceptan -o para escribir a
# fichero, de modo que habría que mantener dos caminos. La segunda es que
# `time` mide un proceso, y aquí interesa el árbol entero: cuando bowtie2
# escribe a samtools por una tuberia los dos están vivos a la vez, y lo que hay
# que reservar es lo que ocupan juntos.
#
# Limitación, que conviene declarar en lugar de aparentar precisión: al
# muestrear una vez por segundo, un pico más corto que eso puede pasar
# desapercibido. Para herramientas que tardan minutos la aproximación es buena;
# para medir con exactitud haría falta cgroups o un `time` específico de la
# plataforma.
RUN_START_EPOCH=$(date +%s)
SAMPLES_FILE=""      # epoch <TAB> rss_kb del árbol de procesos
STEPS_FILE=""        # paso <TAB> inicio <TAB> fin
METRICS_FILE=""      # resumen por paso, ya agregado
SAMPLER_PID=""
RSS_INTERVAL="${RSS_INTERVAL:-1}"

# Memoria residente sumada de un proceso y todos sus descendientes, en KB.
# `ps -eo pid=,ppid=,rss=` se comporta igual en macOS y en Linux, y el RSS
# viene en KB en las dos.
rss_tree_kb() {
    ps -eo pid=,ppid=,rss= 2>/dev/null | awk -v raiz="$1" '
        { padre[$1] = $2; rss[$1] = $3 }
        END {
            descendiente[raiz] = 1
            # El árbol es poco profundo, así que basta con repetir hasta que
            # deje de crecer en vez de ordenarlo topologicamente.
            do {
                nuevo = 0
                for (p in padre)
                    if (!(p in descendiente) && (padre[p] in descendiente)) {
                        descendiente[p] = 1; nuevo = 1
                    }
            } while (nuevo)
            total = 0
            for (p in descendiente) total += rss[p]
            printf "%d\n", total
        }'
}

start_sampler() {
    [[ -n "$SAMPLES_FILE" ]] || return 0
    local raiz=$$
    (
        # El muestreador muere solo si el script desaparece: sin esto, un fallo
        # duro dejaría un proceso huerfano muestreando para siempre.
        while kill -0 "$raiz" 2>/dev/null; do
            printf '%s\t%s\n' "$(date +%s)" "$(rss_tree_kb "$raiz")" >> "$SAMPLES_FILE"
            sleep "$RSS_INTERVAL"
        done
    ) &
    SAMPLER_PID=$!
}

stop_sampler() {
    [[ -n "$SAMPLER_PID" ]] || return 0
    kill "$SAMPLER_PID" 2>/dev/null || true
    wait "$SAMPLER_PID" 2>/dev/null || true
    SAMPLER_PID=""
}

# Marca de paso. Cada llamada cierra el paso anterior y abre uno nuevo, de modo
# que no hay que emparejar inicio y fin a mano en un script de 700 líneas.
CURRENT_STEP=""
CURRENT_STEP_START=0
step_close() {
    [[ -n "$CURRENT_STEP" && -n "$STEPS_FILE" ]] || return 0
    printf '%s\t%s\t%s\n' "$CURRENT_STEP" "$CURRENT_STEP_START" "$(date +%s)" \
        >> "$STEPS_FILE"
    CURRENT_STEP=""
}
step() {
    step_close
    CURRENT_STEP="$1"; shift
    CURRENT_STEP_START=$(date +%s)
    log "$*"
}

# Cruza los pasos con las muestras de memoria y escribe metrics.tsv.
write_metrics() {
    [[ -n "$METRICS_FILE" && -s "$STEPS_FILE" ]] || return 0
    local total=$(( $(date +%s) - RUN_START_EPOCH ))
    {
        printf 'paso\tsegundos\tduracion\tpico_rss_mb\n'
        awk -F'\t' -v muestras="$SAMPLES_FILE" '
            function humano(s,   h, m) {
                h = int(s / 3600); m = int((s % 3600) / 60)
                if (h > 0) return sprintf("%dh %02dm %02ds", h, m, s % 60)
                if (m > 0) return sprintf("%dm %02ds", m, s % 60)
                return sprintf("%ds", s)
            }
            BEGIN {
                n = 0
                while ((getline linea < muestras) > 0) {
                    split(linea, c, "\t")
                    t[n] = c[1] + 0; r[n] = c[2] + 0; n++
                }
            }
            {
                paso = $1; ini = $2 + 0; fin = $3 + 0
                pico = 0
                # Ventana abierta por la izquierda, (inicio, fin).
                #
                # Con los dos extremos incluidos la muestra de la frontera se
                # contaba dos veces, y además el instante en que un paso empieza
                # es justo aquel en el que los procesos del anterior TODAVIA no
                # han terminado de cerrarse. Medido: al arrancar el conteo, la
                # muestra de ese segundo daba 2,4 GB porque FastQC seguía vivo;
                # un segundo después, 82 MB. Atribuir ese pico al conteo es
                # falso, así que la ventana empieza después del instante de
                # arranque. El total de la ejecución no se ve afectado: ese se
                # calcula sobre todas las muestras.
                for (i = 0; i < n; i++)
                    if (t[i] > ini && t[i] < fin && r[i] > pico) pico = r[i]
                printf "%s\t%d\t%s\t%.0f\n", paso, fin - ini, humano(fin - ini), pico / 1024
            }' "$STEPS_FILE"
        # Fila final con el total, para no tener que sumar a ojo.
        awk -F'\t' -v total="$total" -v muestras="$SAMPLES_FILE" '
            function humano(s,   h, m) {
                h = int(s / 3600); m = int((s % 3600) / 60)
                if (h > 0) return sprintf("%dh %02dm %02ds", h, m, s % 60)
                if (m > 0) return sprintf("%dm %02ds", m, s % 60)
                return sprintf("%ds", s)
            }
            BEGIN {
                pico = 0
                while ((getline linea < muestras) > 0) {
                    split(linea, c, "\t")
                    if (c[2] + 0 > pico) pico = c[2] + 0
                }
                printf "TOTAL\t%d\t%s\t%.0f\n", total, humano(total), pico / 1024
            }'
    } > "$METRICS_FILE"
}

trap 'log "ERROR: fallo en la línea ${BASH_LINENO[0]:-$LINENO}${FUNCNAME[0]:+ (función ${FUNCNAME[0]})}. Revisa el comando anterior."' ERR

# Estado de salida como FICHERO, no solo como texto en el log.
#
# La app deducia si una ejecución había terminado bien buscando la frase
# "Analysis completed successfully" en el log y clasificaba como error cualquier
# tail que contuviera "Error". Eso es fragil por partida doble: cambiar el texto
# del log rompe la detección, y el aviso de una herramienta que mencione "Error"
# marca como fallida una ejecución correcta. Con esto queda un dato explícito.
RUN_STATUS_FILE=""
# Señal que interrumpió la ejecución, si la hubo. Sin esto, matar el script
# dejaba escrito "success": bash ejecuta el trap EXIT también al recibir
# SIGTERM, y `$?` en ese instante puede valer 0. El resultado era una ejecución
# a medias registrada como buena, que el harness después SALTABA por creerla
# hecha. Medido: un subjunc interrumpido escribió exit_code 0, status success,
# 1.232 s y 2.750 MB de pico sin haber generado ninguna matriz de conteos.
RUN_SIGNAL=""
trap 'RUN_SIGNAL=SIGTERM' TERM
trap 'RUN_SIGNAL=SIGINT'  INT
trap 'RUN_SIGNAL=SIGHUP'  HUP

write_exit_status() {
    local code=$1
    if [[ -n "$RUN_SIGNAL" ]]; then
        code=143
        [[ "$RUN_SIGNAL" == "SIGINT" ]] && code=130
        [[ "$RUN_SIGNAL" == "SIGHUP" ]] && code=129
    fi
    # Las metricas se cierran aquí y no en el camino feliz: una ejecución que
    # falla a mitad también consumio tiempo y memoria, y saber cuanto es justo
    # lo que hace falta para diagnosticarla.
    step_close
    stop_sampler
    write_metrics
    [[ -n "$RUN_STATUS_FILE" ]] || return 0
    local total=$(( $(date +%s) - RUN_START_EPOCH ))
    local pico=0
    [[ -s "${SAMPLES_FILE:-/dev/null}" ]] && pico=$(awk -F'\t' \
        'BEGIN{m=0} $2+0>m{m=$2+0} END{printf "%.0f", m/1024}' "$SAMPLES_FILE")
    {
        printf 'exit_code\t%s\n'        "$code"
        # Tres estados y no dos: "interrumpido" no es lo mismo que "error". Un
        # error significa que el pipeline falló y hay que arreglar algo; una
        # interrupción significa que no se sabe nada, y sobre todo que no se
        # puede dar por hecha.
        printf 'status\t%s\n'           "$(if [[ -n "$RUN_SIGNAL" ]]; then echo interrumpido
                                             elif [[ "$code" -eq 0 ]]; then echo success
                                             else echo error; fi)"
        [[ -n "$RUN_SIGNAL" ]] && printf 'senal\t%s\n' "$RUN_SIGNAL"
        printf 'finished_at\t%s\n'      "$(date '+%Y-%m-%d %H:%M:%S')"
        printf 'duration_seconds\t%s\n' "$total"
        printf 'peak_rss_mb\t%s\n'      "$pico"
    } > "$RUN_STATUS_FILE"
}
trap 'write_exit_status $?' EXIT

check_command() {
    if ! command -v "$1" >/dev/null 2>&1; then
        log "Error: '$1' no encontrado en PATH. Instala el paquete correspondiente o ajusta tu PATH."
        exit 1
    fi
}

validate_salmon() {
    if ! salmon --version >/dev/null 2>&1; then
        local version_output
        version_output=$(salmon --version 2>&1 || true)
        log "Error: el comando 'salmon' parece ser un paquete Python o una instalación incorrecta."
        log "Salida detectada: $version_output"
        log "Asegúrate de instalar el binario Salmon correcto (por ejemplo desde bioconda) y no el paquete Python 'salmon'."
        exit 1
    fi
}

salmon_version() {
    salmon --version 2>&1 | awk '/^[Ss]almon/ { for (i=1;i<=NF;i++) if ($i ~ /^[0-9]+\.[0-9]+\.[0-9]+$/) print $i; exit }'
}

salmon_index_type() {
    local version
    version=$(salmon_version)
    if [[ "$version" =~ ^1\.[0-9]+\.[0-9]+$ ]]; then
        echo "puff"
    else
        echo "quasi"
    fi
}

# El estado de salida se habilita ANTES de validar las herramientas: la causa
# más frecuente de fallo temprano es justamente que falte un binario en PATH, y
# esa ejecución también tiene que dejar constancia de por qué murio.
mkdir -p "$OUTPUT"
RUN_STATUS_FILE="${OUTPUT}/exit_status.tsv"
SAMPLES_FILE="${OUTPUT}/.rss_samples.tsv"
STEPS_FILE="${OUTPUT}/.steps.tsv"
METRICS_FILE="${OUTPUT}/metrics.tsv"
: > "$SAMPLES_FILE"; : > "$STEPS_FILE"
start_sampler

# Validate required tools early
step preparacion "Validando herramientas en PATH..."
log "PATH=$PATH"
check_command fastqc
check_command fastp
check_command multiqc
if [[ "$ALIGNMENT_TYPE" == "bowtie2" ]]; then
    check_command bowtie2
    check_command samtools
    check_command featureCounts
elif [[ "$ALIGNMENT_TYPE" == "subjunc" ]]; then
    # subjunc y subread-buildindex vienen en el mismo paquete que featureCounts,
    # así que esta ruta no añade ninguna dependencia nueva.
    check_command subjunc
    check_command subread-buildindex
    check_command samtools
    check_command featureCounts
elif [[ "$ALIGNMENT_TYPE" == "salmon" ]]; then
    check_command salmon
    validate_salmon
elif [[ "$ALIGNMENT_TYPE" == "kallisto" ]]; then
    check_command kallisto
fi

# Output subdirectories
QC="${OUTPUT}/01_quality"
TRIMMED="${OUTPUT}/02_trimmed_reads"
ALIGNMENTS="${OUTPUT}/03_alignments/${ALIGNMENT_TYPE}"
COUNTS="${OUTPUT}/04_counts"

# Index paths for different aligners
# Los índices se cachean FUERA del directorio de la ejecución, en una carpeta
# indexada por el md5 de la referencia. Antes vivian en ${OUTPUT}/índices, y como
# la app crea un directorio nuevo por ejecución, la lógica de "índice existente,
# se reutiliza" no se activaba nunca: se reconstruia el índice en cada corrida,
# que en genomas grandes es la parte más lenta del pipeline.
#
# El md5 de la referencia como clave garantiza que un cambio en el genoma o el
# transcriptoma inválida el índice automaticamente.
GENOME_MD5="$(md5_of "$GENOME_FILE")"
INDEX_CACHE="${INDEX_CACHE_DIR:-$(dirname "$OUTPUT")/.index_cache}/${GENOME_MD5}"
mkdir -p "$INDEX_CACHE" 2>/dev/null || INDEX_CACHE="${OUTPUT}/indices"
BOWTIE_INDEX="${INDEX_CACHE}/bowtie2/ref"
SALMON_INDEX="${INDEX_CACHE}/salmon/ref"
KALLISTO_INDEX="${INDEX_CACHE}/kallisto/ref.idx"
SUBJUNC_INDEX="${INDEX_CACHE}/subjunc/ref"

mkdir -p "$QC" "$TRIMMED" "$ALIGNMENTS" "$COUNTS"
mkdir -p "$(dirname "$BOWTIE_INDEX")" "$(dirname "$SALMON_INDEX")" "$(dirname "$KALLISTO_INDEX")"


FASTQ_FILES=( "$INPUT"/*.fastq.gz "$INPUT"/*.fastq )
R1_FILES=( "$INPUT"/*_1.fastq.gz "$INPUT"/*_R1.fastq.gz "$INPUT"/*_1.fastq "$INPUT"/*_R1.fastq )
SINGLE_FILES=()
for fq in "${FASTQ_FILES[@]}"; do
    [[ -e "$fq" ]] || continue
    if [[ "$fq" == *_2.fastq.gz || "$fq" == *_R2.fastq.gz || "$fq" == *_2.fastq || "$fq" == *_R2.fastq ]]; then
        continue
    fi
    SINGLE_FILES+=( "$fq" )
done

if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
    log "Error: no se encontraron FASTQ en $INPUT (*.fastq.gz o *.fastq)."
    exit 1
fi
if [[ "$READ_TYPE" == "pe" && ${#R1_FILES[@]} -eq 0 ]]; then
    log "Error: no se encontraron archivos R1 en $INPUT (*_1.fastq[.gz] o *_R1.fastq[.gz])."
    exit 1
fi
if [[ "$READ_TYPE" == "se" && ${#SINGLE_FILES[@]} -eq 0 ]]; then
    log "Error: no se encontraron FASTQ single-end en $INPUT."
    exit 1
fi

# ── Integridad de los FASTQ ────────────────────────────────────────────────
# Un .gz truncado (transferencia interrumpida, disco lleno) fallaba tarde y con
# un error criptico de fastp o del alineador, a veces después de horas de
# ejecución. Comprobarlo aquí cuesta segundos y el mensaje dice que fichero es.
log "Comprobando integridad de los FASTQ..."
corruptos=0
for fq in "${FASTQ_FILES[@]}"; do
    [[ -e "$fq" ]] || continue
    if [[ "$fq" == *.gz ]]; then
        if ! gzip -t "$fq" 2>/dev/null; then
            log "! FASTQ corrupto o truncado: $fq"
            corruptos=$((corruptos + 1))
        fi
    elif [[ ! -s "$fq" ]]; then
        log "! FASTQ vacio: $fq"
        corruptos=$((corruptos + 1))
    fi
done
if [[ $corruptos -gt 0 ]]; then
    log "Error: $corruptos fichero(s) FASTQ no superan la comprobación de integridad."
    exit 1
fi
log "+ ${#FASTQ_FILES[@]} FASTQ integros"

log "Entrada: $INPUT"
log "Salida: $OUTPUT"
log "Genoma/transcriptoma: $GENOME_FILE"
log "Anotación: $ANNOTATION_FILE"
log "Tipo de alineamiento: $ALIGNMENT_TYPE"
log "Tipo de lectura: $READ_TYPE"

# ── Atributo de conteo: comprobarlo contra la anotación de verdad ──────────
# FEATURE_ATTR viene por defecto como `locus_tag` porque el trabajo empezó con
# procariotas de NCBI, donde ese atributo es la norma. Las anotaciones de
# Ensembl para eucariotas NO lo traen: el GTF de S. cerevisiae tiene cero
# apariciones de locus_tag. featureCounts con `-g locus_tag` sobre ese fichero
# no cuenta nada, y el fallo salia como una matriz vacía después de haber
# alineado todas las muestras.
#
# La interfaz tampoco pasa --FEATURE_ATTR, así que desde la aplicación era
# imposible corregirlo. Aquí se comprueba si el atributo pedido existe y, si no,
# se cae a la primera alternativa que si esté, dejandolo dicho en el log y en
# run_params.tsv.
if [[ -f "$ANNOTATION_FILE" ]]; then
    # El tipo de feature va primero y NO se sustituye. Un tipo equivocado da
    # cero conteos sin que ninguna alternativa sea evidentemente la buena, así
    # que es mejor pararse aquí que entregar una matriz vacía. Comprobarlo antes
    # que el atributo tambien da el mensaje correcto: con un tipo que no existe,
    # ningún atributo aparece, y el error del atributo despistaria.
    if ! awk -F'\t' -v tipo="$FEATURE_TYPE" '
            /^#/ { next } $3 == tipo { encontrado = 1; exit }
            END { exit(encontrado ? 0 : 1) }' "$ANNOTATION_FILE"; then
        log "Error: la anotación no contiene ninguna linea de tipo '$FEATURE_TYPE'."
        log "  Tipos presentes: $(awk -F'\t' '!/^#/ {print $3}' "$ANNOTATION_FILE" | sort -u | head -12 | tr '\n' ' ')"
        exit 1
    fi

    attr_presente() {
        awk -F'\t' -v tipo="$2" -v attr="$3" '
            /^#/ { next }
            $3 == tipo && $9 ~ attr" [\"=]" { encontrado = 1; exit }
            END { exit(encontrado ? 0 : 1) }
        ' "$1"
    }
    if ! attr_presente "$ANNOTATION_FILE" "$FEATURE_TYPE" "$FEATURE_ATTR"; then
        FEATURE_ATTR_PEDIDO="$FEATURE_ATTR"
        FEATURE_ATTR=""
        for cand in locus_tag gene_id ID Name gene; do
            if attr_presente "$ANNOTATION_FILE" "$FEATURE_TYPE" "$cand"; then
                FEATURE_ATTR="$cand"; break
            fi
        done
        if [[ -z "$FEATURE_ATTR" ]]; then
            log "Error: en las lineas de tipo '$FEATURE_TYPE' la anotación no trae"
            log "  ninguno de los atributos conocidos (locus_tag, gene_id, ID, Name, gene)."
            exit 1
        fi
        log "Atributo de conteo: '$FEATURE_ATTR_PEDIDO' no existe en la anotación; se usa '$FEATURE_ATTR'."
    fi
    log "Conteo: -t $FEATURE_TYPE -g $FEATURE_ATTR"
fi

# Parámetros en formato legible por la app. Sin esto, una ejecución guardada no
# deja rastro de con que anotación se hizo, y la app no puede construir el mapa
# transcrito-gen para tximport ni identificar los genes de rRNA.
RUN_ID="$(date '+%Y%m%d_%H%M%S')_$$"
{
    printf 'run_id\t%s\n'          "$RUN_ID"
    printf 'input_dir\t%s\n'       "$INPUT"
    printf 'output_dir\t%s\n'      "$OUTPUT"
    printf 'genome_file\t%s\n'     "$GENOME_FILE"
    printf 'annotation_file\t%s\n' "$ANNOTATION_FILE"
    printf 'alignment_type\t%s\n'  "$ALIGNMENT_TYPE"
    printf 'tool\t%s\n'            "$ALIGNMENT_TYPE"
    printf 'read_type\t%s\n'       "$READ_TYPE"
    printf 'fragment_length\t%s\n' "$FRAGMENT_LENGTH"
    printf 'fragment_sd\t%s\n'     "$FRAGMENT_SD"
    printf 'inferential_reps\t%s\n' "$INFERENTIAL_REPS"
    printf 'feature_type\t%s\n'     "$FEATURE_TYPE"
    printf 'feature_attr\t%s\n'     "$FEATURE_ATTR"
    printf 'threads\t%s\n'          "$THREADS"
    printf 'n_fastq\t%s\n'          "${#FASTQ_FILES[@]}"
    # Quien, donde y con que versión del pipeline. Es el mínimo de un registro
    # de auditoria y lo exigen las buenas prácticas de laboratorio clínico para
    # los componentes informáticos.
    printf 'user\t%s\n'             "$(whoami 2>/dev/null || echo '—')"
    printf 'host\t%s\n'             "$(hostname 2>/dev/null || echo '—')"
    printf 'workflow_md5\t%s\n'     "$(md5_of "${BASH_SOURCE[0]}")"
    printf 'workflow_git_commit\t%s\n' \
        "$(git -C "$SCRIPT_DIR" rev-parse --short HEAD 2>/dev/null || echo '—')"
    printf 'started_at\t%s\n'      "$(date '+%Y-%m-%d %H:%M:%S')"
} > "${OUTPUT}/run_params.tsv"

# ── Versiones de las herramientas ──────────────────────────────────────────
# Sin esto es imposible reconstruir que produjo exactamente un resultado. Se ha
# demostrado que cambiar la versión de una herramienta del pipeline altera la
# lista de genes diferenciales (Wessels Perelo et al., NAR Genom Bioinform 2024),
# así que la versión es parte del resultado, no un detalle de instalación.
log "Registrando versiones de las herramientas..."
{
    printf 'tool\tversion\tpath\n'
    for t in fastqc fastp multiqc bowtie2 subjunc samtools featureCounts salmon kallisto; do
        if command -v "$t" >/dev/null 2>&1; then
            # Dos cuidados aquí, los dos aprendidos ejecutando.
            #
            # Primero, NO se canaliza hacía `head`: con `set -o pipefail`, una
            # herramienta que sigue escribiendo después de la primera línea
            # recibe SIGPIPE cuando head cierra, la tuberia devuelve 141 y
            # `set -e` aborta el script. Ocurría con bowtie2, que imprime nueve
            # líneas: la ejecución moria al registrar versiones, antes de
            # alinear nada.
            #
            # Segundo, no todas aceptan `--version`: featureCounts usa `-v` y
            # kallisto usa `version`. Y bowtie2 antepone dos líneas de aviso a
            # la versión real. Coger la primera línea a ciegas registraba un
            # aviso o un mensaje de error en el fichero de procedencia, que es
            # justo el fichero que debe poder creerse.
            case "$t" in
                featureCounts|subjunc) v_full="$("$t" -v 2>&1 || true)" ;;
                kallisto)      v_full="$("$t" version 2>&1 || true)" ;;
                *)             v_full="$("$t" --version 2>&1 || true)" ;;
            esac
            # Primera línea no vacia que no sea un aviso ni un error, en UNA
            # sola pasada de awk. Encadenar dos filtros donde el segundo sale
            # pronto reproduce el mismo SIGPIPE que se acaba de corregir.
            v="$(printf '%s\n' "$v_full" \
                 | awk 'tolower($0) !~ /^(\[warning\]|warning|error)/ && NF { print; exit }' \
                 || true)"
            v="$(printf '%s' "$v" | tr -d '\r' | sed 's/\t/ /g')"
            printf '%s\t%s\t%s\n' "$t" "${v:-—}" "$(command -v "$t")"
        else
            printf '%s\t(no instalado)\t—\n' "$t"
        fi
    done
    printf 'R\t%s\t%s\n' \
        "$(Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || echo '—')" \
        "$(command -v Rscript || echo '—')"
} > "${OUTPUT}/versions.tsv"

# ── Huella de los ficheros de entrada ──────────────────────────────────────
# Permite demostrar que dos ejecuciones partieron de los mismos datos, y
# detectar que una entrada cambio después de analizarla.
log "Calculando checksums de las entradas..."
{
    printf 'file\tsize_bytes\tmd5\n'
    for f in "$GENOME_FILE" "$ANNOTATION_FILE"; do
        [[ -f "$f" ]] || continue
        printf '%s\t%s\t%s\n' "$f" "$(file_size "$f")" "$(md5_of "$f")"
    done
    for fq in "${FASTQ_FILES[@]}"; do
        [[ -e "$fq" ]] || continue
        printf '%s\t%s\t%s\n' "$fq" "$(file_size "$fq")" "$(md5_of "$fq")"
    done
} > "${OUTPUT}/checksums.tsv"
log "+ versions.tsv y checksums.tsv escritos"
if [[ "$READ_TYPE" == "se" && "$ALIGNMENT_TYPE" == "kallisto" ]]; then
    log "Fragment length for kallisto single-end: mean=$FRAGMENT_LENGTH sd=$FRAGMENT_SD"
fi
log "FASTQ detectados: ${#FASTQ_FILES[@]}"
if [[ "$READ_TYPE" == "pe" ]]; then
    log "R1 detectados: ${#R1_FILES[@]}"
else
    log "FASTQ single-end detectados: ${#SINGLE_FILES[@]}"
fi

# Build genome index based on alignment type
step indice "Building ${ALIGNMENT_TYPE} index..."

if [[ "$ALIGNMENT_TYPE" == "bowtie2" ]]; then
    if [[ -s "${BOWTIE_INDEX}.1.bt2" && -s "${BOWTIE_INDEX}.2.bt2" && -s "${BOWTIE_INDEX}.3.bt2" && -s "${BOWTIE_INDEX}.4.bt2" && -s "${BOWTIE_INDEX}.rev.1.bt2" && -s "${BOWTIE_INDEX}.rev.2.bt2" ]]; then
        log "Índice Bowtie2 existente detectado; se reutiliza."
    else
        rm -f "${BOWTIE_INDEX}".*.bt2.tmp "${BOWTIE_INDEX}".*.bt2
        log "Construyendo índice Bowtie2. Este paso puede tardar varios minutos y puede no imprimir progreso continuo."
        run_cmd bowtie2-build "$GENOME_FILE" "$BOWTIE_INDEX"
    fi
    INDEX="$BOWTIE_INDEX"
elif [[ "$ALIGNMENT_TYPE" == "salmon" ]]; then
    if [[ -d "$SALMON_INDEX" && -s "${SALMON_INDEX}/versionInfo.json" ]]; then
        log "Índice Salmon existente detectado; se reutiliza."
    else
        rm -rf "$SALMON_INDEX"
        SALMON_INDEX_TYPE=$(salmon_index_type)
        log "Using Salmon index type: $SALMON_INDEX_TYPE"
        run_cmd salmon index -t "$GENOME_FILE" -i "$SALMON_INDEX" --type "$SALMON_INDEX_TYPE" -k 31
    fi
elif [[ "$ALIGNMENT_TYPE" == "kallisto" ]]; then
    if [[ -s "$KALLISTO_INDEX" ]]; then
        log "Índice kallisto existente detectado; se reutiliza."
    else
        rm -f "$KALLISTO_INDEX"
        run_cmd kallisto index -i "$KALLISTO_INDEX" "$GENOME_FILE"
    fi
elif [[ "$ALIGNMENT_TYPE" == "subjunc" ]]; then
    # subread-buildindex deja varios ficheros con el mismo prefijo; el .00.b.tab
    # es el último que escribe, así que su presencia indica un índice completo.
    if [[ -s "${SUBJUNC_INDEX}.00.b.tab" && -s "${SUBJUNC_INDEX}.reads" ]]; then
        log "Índice Subread existente detectado; se reutiliza."
    else
        mkdir -p "$(dirname "$SUBJUNC_INDEX")"
        rm -f "${SUBJUNC_INDEX}".*
        log "Construyendo índice Subread. Puede tardar varios minutos."
        run_cmd subread-buildindex -o "$SUBJUNC_INDEX" "$GENOME_FILE"
    fi
    INDEX="$SUBJUNC_INDEX"
fi

# Initial quality control
step qc_inicial "Control de calidad inicial con FastQC..."
run_cmd fastqc "${FASTQ_FILES[@]}" -t "$THREADS" -o "$QC"

if [[ "$READ_TYPE" == "pe" ]]; then
    SAMPLE_FILES=( "${R1_FILES[@]}" )
else
    SAMPLE_FILES=( "${SINGLE_FILES[@]}" )
fi

# Process each sample
# Es el paso que domina el tiempo total: recorte con fastp más alineamiento o
# cuantificación, una vez por muestra.
step procesado "Procesando ${#SAMPLE_FILES[@]} muestras (recorte y ${ALIGNMENT_TYPE})..."
for READ1 in "${SAMPLE_FILES[@]}"; do

    [ -e "$READ1" ] || continue

    SAMPLE=$(basename "$READ1")
    SAMPLE=${SAMPLE%%_1.fastq.gz}
    SAMPLE=${SAMPLE%%_R1.fastq.gz}
    SAMPLE=${SAMPLE%%.fastq.gz}
    SAMPLE=${SAMPLE%%_1.fastq}
    SAMPLE=${SAMPLE%%_R1.fastq}
    SAMPLE=${SAMPLE%%.fastq}

    if [[ "$READ_TYPE" == "pe" ]]; then
        if [[ "$READ1" == *_1.fastq.gz ]]; then
            READ2="${READ1%_1.fastq.gz}_2.fastq.gz"
        elif [[ "$READ1" == *_R1.fastq.gz ]]; then
            READ2="${READ1%_R1.fastq.gz}_R2.fastq.gz"
        elif [[ "$READ1" == *_1.fastq ]]; then
            READ2="${READ1%_1.fastq}_2.fastq"
        else
            READ2="${READ1%_R1.fastq}_R2.fastq"
        fi

        if [[ ! -s "$READ2" ]]; then
            log "Error: falta R2 para $SAMPLE: $READ2"
            exit 1
        fi
    fi

    log "Processing sample: $SAMPLE"

    # Adapter trimming and low-quality read filtering
    if [[ "$READ_TYPE" == "pe" ]]; then
        run_cmd fastp \
            -i "$READ1" \
            -I "$READ2" \
            -o "${TRIMMED}/${SAMPLE}_R1_trimmed.fastq.gz" \
            -O "${TRIMMED}/${SAMPLE}_R2_trimmed.fastq.gz" \
            --detect_adapter_for_pe \
            --thread "$THREADS" \
            --html "${TRIMMED}/${SAMPLE}_fastp.html" \
            --json "${TRIMMED}/${SAMPLE}_fastp.json"
    else
        run_cmd fastp \
            -i "$READ1" \
            -o "${TRIMMED}/${SAMPLE}_trimmed.fastq.gz" \
            --thread "$THREADS" \
            --html "${TRIMMED}/${SAMPLE}_fastp.html" \
            --json "${TRIMMED}/${SAMPLE}_fastp.json"
    fi

    # Perform alignment based on alignment type
    if [[ "$ALIGNMENT_TYPE" == "bowtie2" ]]; then
        log "  Aligning with Bowtie2..."
        if [[ "$READ_TYPE" == "pe" ]]; then
            log "+ bowtie2 -x $INDEX -1 ${TRIMMED}/${SAMPLE}_R1_trimmed.fastq.gz -2 ${TRIMMED}/${SAMPLE}_R2_trimmed.fastq.gz -p $THREADS | samtools sort -@ $THREADS -o ${ALIGNMENTS}/${SAMPLE}.bam"
            bowtie2 \
                -x "$INDEX" \
                -1 "${TRIMMED}/${SAMPLE}_R1_trimmed.fastq.gz" \
                -2 "${TRIMMED}/${SAMPLE}_R2_trimmed.fastq.gz" \
                -p "$THREADS" \
            | samtools sort \
                -@ "$THREADS" \
                -o "${ALIGNMENTS}/${SAMPLE}.bam"
        else
            log "+ bowtie2 -x $INDEX -U ${TRIMMED}/${SAMPLE}_trimmed.fastq.gz -p $THREADS | samtools sort -@ $THREADS -o ${ALIGNMENTS}/${SAMPLE}.bam"
            bowtie2 \
                -x "$INDEX" \
                -U "${TRIMMED}/${SAMPLE}_trimmed.fastq.gz" \
                -p "$THREADS" \
            | samtools sort \
                -@ "$THREADS" \
                -o "${ALIGNMENTS}/${SAMPLE}.bam"
        fi

        # Index BAM file for downstream analysis
        run_cmd samtools index "${ALIGNMENTS}/${SAMPLE}.bam"

    elif [[ "$ALIGNMENT_TYPE" == "subjunc" ]]; then
        # subjunc SÍ es splice-aware, a diferencia de bowtie2: detecta las
        # uniones exón-exón y por tanto vale para eucariotas con intrones, donde
        # bowtie2 perdería las lecturas que cruzan una unión y subestimaría los
        # conteos de forma sistemática.
        #
        # Escribe el BAM directamente en lugar de por tubería, porque no tiene
        # salida a stdout: se ordena después con samtools, igual que bowtie2.
        log "  Aligning with Subjunc (splice-aware)..."
        if [[ "$READ_TYPE" == "pe" ]]; then
            run_cmd subjunc \
                -i "$INDEX" \
                -r "${TRIMMED}/${SAMPLE}_R1_trimmed.fastq.gz" \
                -R "${TRIMMED}/${SAMPLE}_R2_trimmed.fastq.gz" \
                -o "${ALIGNMENTS}/${SAMPLE}.unsorted.bam" \
                -T "$THREADS"
        else
            run_cmd subjunc \
                -i "$INDEX" \
                -r "${TRIMMED}/${SAMPLE}_trimmed.fastq.gz" \
                -o "${ALIGNMENTS}/${SAMPLE}.unsorted.bam" \
                -T "$THREADS"
        fi
        run_cmd samtools sort -@ "$THREADS" \
            -o "${ALIGNMENTS}/${SAMPLE}.bam" "${ALIGNMENTS}/${SAMPLE}.unsorted.bam"
        rm -f "${ALIGNMENTS}/${SAMPLE}.unsorted.bam"
        run_cmd samtools index "${ALIGNMENTS}/${SAMPLE}.bam"

    elif [[ "$ALIGNMENT_TYPE" == "salmon" ]]; then
        log "  Quasi-aligning with Salmon..."
        if [[ "$READ_TYPE" == "pe" ]]; then
            run_cmd salmon quant \
                -i "$SALMON_INDEX" \
                -l A \
                -1 "${TRIMMED}/${SAMPLE}_R1_trimmed.fastq.gz" \
                -2 "${TRIMMED}/${SAMPLE}_R2_trimmed.fastq.gz" \
                -p "$THREADS" \
                --numGibbsSamples "$INFERENTIAL_REPS" \
                -o "${ALIGNMENTS}/${SAMPLE}"
        else
            run_cmd salmon quant \
                -i "$SALMON_INDEX" \
                -l A \
                -r "${TRIMMED}/${SAMPLE}_trimmed.fastq.gz" \
                -p "$THREADS" \
                --numGibbsSamples "$INFERENTIAL_REPS" \
                -o "${ALIGNMENTS}/${SAMPLE}"
        fi

    elif [[ "$ALIGNMENT_TYPE" == "kallisto" ]]; then
        log "  Quasi-aligning with kallisto..."
        if [[ "$READ_TYPE" == "pe" ]]; then
            run_cmd kallisto quant \
                -i "$KALLISTO_INDEX" \
                -o "${ALIGNMENTS}/${SAMPLE}" \
                -t "$THREADS" \
                -b "$INFERENTIAL_REPS" \
                "${TRIMMED}/${SAMPLE}_R1_trimmed.fastq.gz" \
                "${TRIMMED}/${SAMPLE}_R2_trimmed.fastq.gz"
        else
            run_cmd kallisto quant \
                -i "$KALLISTO_INDEX" \
                -o "${ALIGNMENTS}/${SAMPLE}" \
                -t "$THREADS" \
                -b "$INFERENTIAL_REPS" \
                --single \
                -l "$FRAGMENT_LENGTH" \
                -s "$FRAGMENT_SD" \
                "${TRIMMED}/${SAMPLE}_trimmed.fastq.gz"
        fi

    fi

done

# Quality control after trimming
TRIMMED_FASTQ=( "$TRIMMED"/*.fastq.gz )
if [[ ${#TRIMMED_FASTQ[@]} -eq 0 ]]; then
    log "Error: no se encontraron FASTQ recortados en $TRIMMED."
    exit 1
fi
step qc_final "Control de calidad post-trimming con FastQC..."
if [[ ${#TRIMMED_FASTQ[@]} -eq 0 ]]; then
    log "Error: el recorte no dejó ningún FASTQ en $TRIMMED."
    exit 1
fi
run_cmd fastqc "${TRIMMED_FASTQ[@]}" -t "$THREADS" -o "$QC"

# ── Orientación de la libreria ─────────────────────────────────────────────
# featureCounts contaba SIEMPRE como no orientada (-s 0). La mayoria de los
# protocolos actuales son stranded (dUTP, que corresponde a -s 2), y contar una
# libreria stranded como no orientada suma las lecturas antisentido: en genomas
# de alta densidad genica, como los procariotas, eso infla los conteos de genes
# solapantes y degrada la especificidad (Zhao et al., BMC Genomics 2015).
#
# En modo "auto" se infiere de los propios datos contando un subconjunto de
# lecturas con las tres orientaciones y quedandose con la que asigna más.
infer_strandedness() {
    local bam="$1" s0 s1 s2 best
    local tmp="${COUNTS}/.strand_check"
    local -a pe_flags=()
    # featureCounts ABORTA con "Paired-end reads were detected in single-end read
    # library" si recibe un BAM pareado sin -p. La inferencia lo llamaba sin esos
    # flags y con `|| true`, así que el error se tragaba, el fichero de resumen no
    # se escribía y las tres orientaciones salian 0: la función devolvía siempre
    # 0 (sin orientar) creyendo haberlo medido.
    [[ "$READ_TYPE" == "pe" ]] && pe_flags=(-p --countReadPairs)
    mkdir -p "$tmp"
    for s in 0 1 2; do
        # `${pe_flags[@]+"${pe_flags[@]}"}` y no `"${pe_flags[@]}"`: macOS trae
        # bash 3.2, donde expandir un array VACÍO bajo `set -u` aborta con
        # "unbound variable". Con datos single-end pe_flags está vacío, así que
        # la forma directa mata la ejecución justo aquí, después de haber
        # alineado todo. No es cosmético y no se puede "simplificar".
        featureCounts -T "$THREADS" ${pe_flags[@]+"${pe_flags[@]}"} -s "$s" \
            -t "$FEATURE_TYPE" -g "$FEATURE_ATTR" \
            -a "$ANNOTATION_FILE" -o "${tmp}/s${s}.txt" "$bam" >/dev/null 2>&1 || true
    done
    s0=$(awk 'NR>1 && $1=="Assigned" {print $2}' "${tmp}/s0.txt.summary" 2>/dev/null || echo 0)
    s1=$(awk 'NR>1 && $1=="Assigned" {print $2}' "${tmp}/s1.txt.summary" 2>/dev/null || echo 0)
    s2=$(awk 'NR>1 && $1=="Assigned" {print $2}' "${tmp}/s2.txt.summary" 2>/dev/null || echo 0)
    s0=${s0:-0}; s1=${s1:-0}; s2=${s2:-0}
    # El diagnóstico va a STDERR: esta función se invoca con $(...) y su stdout es
    # el valor de retorno. Escribiendolo en stdout, el mensaje acababa DENTRO de
    # la variable y `-s` recibía una línea de log entera.
    log "  asignadas por orientación -> sin orientar: $s0 | directa: $s1 | inversa: $s2" >&2
    # La comparación NO puede ser "la que más asigna". El modo sin orientar cuenta
    # las lecturas en ambos sentidos, así que su recuento es por construcción mayor
    # o igual que el de cualquiera de los dos modos orientados: ganaria siempre, y
    # la inferencia devolveria "sin orientar" para cualquier libreria, también para
    # una dUTP perfectamente orientada.
    #
    # El criterio correcto es el que usan las herramientas al uso (RSeQC): mirar
    # como se REPARTEN entre los dos sentidos las lecturas que si son asignables
    # por hebra. Si el reparto es muy asimétrico, la libreria está orientada en ese
    # sentido; si está repartido, no lo está.
    local total_orientado=$(( s1 + s2 ))
    best=0
    if [[ "$total_orientado" -lt 1000 ]]; then
        log "  Muy pocas lecturas asignables por hebra ($total_orientado) para decidir; se usa 0." >&2
    else
        local pct_inversa=$(( 100 * s2 / total_orientado ))
        log "  Reparto entre sentidos: inversa ${pct_inversa} %, directa $(( 100 - pct_inversa )) %." >&2
        if   [[ "$pct_inversa" -ge 80 ]]; then best=2
        elif [[ "$pct_inversa" -le 20 ]]; then best=1
        else best=0
        fi
    fi
    rm -rf "$tmp"
    # Guarda final: si nada se pudo medir, no se puede afirmar una orientación.
    if [[ "$s0" -eq 0 && "$s1" -eq 0 && "$s2" -eq 0 ]]; then
        log "  ADVERTENCIA: no se pudo medir la orientación; se usa 0 (sin orientar)." >&2
    fi
    printf '%s' "$best"
}

if [[ ( "$ALIGNMENT_TYPE" == "bowtie2" || "$ALIGNMENT_TYPE" == "subjunc" ) && "$STRANDEDNESS" == "auto" ]]; then
    FIRST_BAM=$( (shopt -s nullglob; b=( "$ALIGNMENTS"/*.bam ); echo "${b[0]:-}") )
    if [[ -n "$FIRST_BAM" && -f "$FIRST_BAM" ]]; then
        log "Infiriendo la orientación de la libreria a partir de $(basename "$FIRST_BAM")..."
        STRANDEDNESS=$(infer_strandedness "$FIRST_BAM")
        log "+ Orientación inferida: -s $STRANDEDNESS"
    else
        STRANDEDNESS=0
    fi
elif [[ "$STRANDEDNESS" == "auto" ]]; then
    # En pseudoalineamiento la detecta el propio cuantificador (salmon -l A).
    STRANDEDNESS=0
fi
printf 'strandedness\t%s\n' "$STRANDEDNESS" >> "${OUTPUT}/run_params.tsv"

# Generate count matrix based on alignment type
step conteo "Generando la matriz de conteos..."
if [[ "$ALIGNMENT_TYPE" == "bowtie2" || "$ALIGNMENT_TYPE" == "subjunc" ]]; then
    log "Generando la matriz de conteos a partir de los alineamientos de ${ALIGNMENT_TYPE}..."
    # Count reads per gene
    BAM_FILES=( "$ALIGNMENTS"/*.bam )
    if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
        log "Error: no se generaron BAM en $ALIGNMENTS."
        exit 1
    fi
    if [[ "$READ_TYPE" == "pe" ]]; then
        run_cmd featureCounts \
            -T "$THREADS" \
            -p \
            --countReadPairs \
            -s "$STRANDEDNESS" \
            -t "$FEATURE_TYPE" \
            -g "$FEATURE_ATTR" \
            -a "$ANNOTATION_FILE" \
            -o "${COUNTS}/gene_counts.txt" \
            "${BAM_FILES[@]}"
    else
        run_cmd featureCounts \
            -T "$THREADS" \
            -s "$STRANDEDNESS" \
            -t "$FEATURE_TYPE" \
            -g "$FEATURE_ATTR" \
            -a "$ANNOTATION_FILE" \
            -o "${COUNTS}/gene_counts.txt" \
            "${BAM_FILES[@]}"
    fi

    # Create count matrix in TSV format
    log "+ Crear matriz de conteos: ${COUNTS}/count_matrix.tsv"
    # La cabecera se limpia CAMPO A CAMPO, no con un sed sobre la línea entera.
    # El `sed '1s|.*/||g'` que había aquí es codicioso sobre toda la línea: con
    # rutas absolutas colapsaba `Geneid\t/ruta/A.bam\t/ruta/B.bam` en un único
    # campo `B`, dejando una cabecera de 1 columna frente a N+1 de datos. El
    # lector de la app no encontraba el patron `^Geneid\t`, read.delim abortaba y
    # devolvía NULL, y la interfaz anunciaba "Workflow finalizado correctamente"
    # sin matriz. Con awk cada campo se procesa por separado.
    awk 'BEGIN { FS=OFS="\t" } !/^#/ { print }' "${COUNTS}/gene_counts.txt" \
        | cut -f1,7- \
        | awk 'BEGIN { FS=OFS="\t" }
               NR==1 { for (i=2; i<=NF; i++) { sub(/.*\//, "", $i); sub(/\.bam$/, "", $i) } }
               { print }' \
        > "${COUNTS}/count_matrix.tsv"

elif [[ "$ALIGNMENT_TYPE" == "salmon" || "$ALIGNMENT_TYPE" == "kallisto" ]]; then
    # Resumen de transcritos a genes con tximport, reutilizando la misma función
    # que usa la aplicación (scripts/build_count_matrix.R). Antes se escribía un
    # fichero con SOLO la cabecera y un comentario diciendo que hacía falta un
    # script de R: quien se llevaba la carpeta de resultados obtenia una matriz
    # vacia presentada como el entregable principal.
    log "Resumiendo transcritos a genes con tximport (${ALIGNMENT_TYPE})..."
    matrix_script="${SCRIPT_DIR}/scripts/build_count_matrix.R"
    if [[ ! -f "$matrix_script" ]]; then
        log "! No se encuentra ${matrix_script}: no se genera la matriz por gen."
    elif Rscript "$matrix_script" "$OUTPUT" "$ALIGNMENT_TYPE" \
            "${COUNTS}/count_matrix.tsv" "$ANNOTATION_FILE"; then
        log "+ Matriz de conteos por gen: ${COUNTS}/count_matrix.tsv"
    else
        # Sin anotación utilizable no hay resumen a gen posible. Se avisa y NO se
        # deja un fichero a medias: es preferible que falte el artefacto a que
        # exista uno que parece válido y no lo es.
        rm -f "${COUNTS}/count_matrix.tsv"
        log "! No se ha podido construir la matriz por gen (revisa la anotación)."
        log "  Las cuantificaciones por muestra quedan en ${COUNTS}/ y la app"
        log "  puede resumirlas al cargar la ejecución."
    fi

    # Copia de las cuantificaciones por muestra, con el nombre de la muestra en
    # el fichero. Antes se copiaban con `cp "$ALIGNMENTS"/*/*.sf "${COUNTS}/"`:
    # todas se llaman quant.sf (o abundance.tsv), así que cada copia sobrescribia
    # la anterior y solo sobrevivia la última muestra, en silencio por el
    # `|| true`.
    if [[ "$ALIGNMENT_TYPE" == "salmon" ]]; then quant_name="quant.sf"
    else quant_name="abundance.tsv"; fi
    for dir in "$ALIGNMENTS"/*/; do
        [[ -d "$dir" ]] || continue
        sample=$(basename "$dir")
        if [[ -f "${dir}${quant_name}" ]]; then
            cp "${dir}${quant_name}" "${COUNTS}/${sample}_${quant_name}"
        fi
    done

fi

# Generate combined quality and alignment report
step multiqc "Generando informe MultiQC..."
run_cmd multiqc "$OUTPUT" -o "$OUTPUT"

log "Analysis completed successfully."
log "Alignment type: $ALIGNMENT_TYPE"
log "Results saved in: $OUTPUT"
log "Quality control: $QC"
log "Trimmed reads: $TRIMMED"
log "Alignments: $ALIGNMENTS"
log "Count matrices: $COUNTS"
