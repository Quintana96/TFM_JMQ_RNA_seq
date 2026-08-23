#!/bin/bash

# -E (errtrace) hace que el trap ERR se herede en funciones y subshells. Sin el,
# un fallo dentro de run_cmd o del pipe bowtie2|samtools terminaba el script
# (correcto) pero sin el mensaje que dice en que linea, que es justo lo que se
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
# Replicas inferenciales de la cuantificacion (salmon: --numGibbsSamples,
# kallisto: -b). Las necesita Swish para propagar la incertidumbre de asignacion
# de lecturas entre transcritos que comparten secuencia. 0 las desactiva.
INFERENTIAL_REPS=20
# Tipo de feature que cuenta featureCounts (-t) y atributo que agrupa (-g).
# El default se mantiene en gene/locus_tag para no cambiar el comportamiento
# existente. Las buenas practicas para procariotas (Genome Biology 2021) hacen el
# analisis diferencial a nivel de CDS; con FEATURE_TYPE=CDS se obtiene eso, a
# cambio de excluir los genes no codificantes (rRNA, tRNA) del recuento.
FEATURE_TYPE="gene"
FEATURE_ATTR="locus_tag"
# Orientacion de la libreria para featureCounts (-s): 0 sin orientar, 1 directa,
# 2 inversa (el caso de los protocolos dUTP, que son la mayoria hoy). "auto"
# la infiere de los propios datos.
STRANDEDNESS="auto"
# Hilos. Antes estaba cableado a 8, lo que sobresuscribe maquinas mas pequenas
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
if [[ ! "$ALIGNMENT_TYPE" =~ ^(bowtie2|salmon|kallisto)$ ]]; then
    echo "Error: Invalid ALIGNMENT_TYPE. Choose: bowtie2, salmon, or kallisto"
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

# Numero de hilos: el indicado, o los nucleos disponibles menos uno para no
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

trap 'log "ERROR: fallo en la linea ${BASH_LINENO[0]:-$LINENO}${FUNCNAME[0]:+ (funcion ${FUNCNAME[0]})}. Revisa el comando anterior."' ERR

# Estado de salida como FICHERO, no solo como texto en el log.
#
# La app deducia si una ejecucion habia terminado bien buscando la frase
# "Analysis completed successfully" en el log y clasificaba como error cualquier
# tail que contuviera "Error". Eso es fragil por partida doble: cambiar el texto
# del log rompe la deteccion, y el aviso de una herramienta que mencione "Error"
# marca como fallida una ejecucion correcta. Con esto queda un dato explicito.
RUN_STATUS_FILE=""
write_exit_status() {
    local code=$1
    [[ -n "$RUN_STATUS_FILE" ]] || return 0
    {
        printf 'exit_code\t%s\n'   "$code"
        printf 'status\t%s\n'      "$([[ "$code" -eq 0 ]] && echo success || echo error)"
        printf 'finished_at\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S')"
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
# mas frecuente de fallo temprano es justamente que falte un binario en PATH, y
# esa ejecucion tambien tiene que dejar constancia de por que murio.
mkdir -p "$OUTPUT"
RUN_STATUS_FILE="${OUTPUT}/exit_status.tsv"

# Validate required tools early
log "Validando herramientas en PATH..."
log "PATH=$PATH"
check_command fastqc
check_command fastp
check_command multiqc
if [[ "$ALIGNMENT_TYPE" == "bowtie2" ]]; then
    check_command bowtie2
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
# Los indices se cachean FUERA del directorio de la ejecucion, en una carpeta
# indexada por el md5 de la referencia. Antes vivian en ${OUTPUT}/indices, y como
# la app crea un directorio nuevo por ejecucion, la logica de "indice existente,
# se reutiliza" no se activaba nunca: se reconstruia el indice en cada corrida,
# que en genomas grandes es la parte mas lenta del pipeline.
#
# El md5 de la referencia como clave garantiza que un cambio en el genoma o el
# transcriptoma invalida el indice automaticamente.
GENOME_MD5="$(md5_of "$GENOME_FILE")"
INDEX_CACHE="${INDEX_CACHE_DIR:-$(dirname "$OUTPUT")/.index_cache}/${GENOME_MD5}"
mkdir -p "$INDEX_CACHE" 2>/dev/null || INDEX_CACHE="${OUTPUT}/indices"
BOWTIE_INDEX="${INDEX_CACHE}/bowtie2/ref"
SALMON_INDEX="${INDEX_CACHE}/salmon/ref"
KALLISTO_INDEX="${INDEX_CACHE}/kallisto/ref.idx"

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
# un error criptico de fastp o del alineador, a veces despues de horas de
# ejecucion. Comprobarlo aqui cuesta segundos y el mensaje dice que fichero es.
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
    log "Error: $corruptos fichero(s) FASTQ no superan la comprobacion de integridad."
    exit 1
fi
log "+ ${#FASTQ_FILES[@]} FASTQ integros"

log "Entrada: $INPUT"
log "Salida: $OUTPUT"
log "Genoma/transcriptoma: $GENOME_FILE"
log "Anotacion: $ANNOTATION_FILE"
log "Tipo de alineamiento: $ALIGNMENT_TYPE"
log "Tipo de lectura: $READ_TYPE"

# Parametros en formato legible por la app. Sin esto, una ejecucion guardada no
# deja rastro de con que anotacion se hizo, y la app no puede construir el mapa
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
    # Quien, donde y con que version del pipeline. Es el minimo de un registro
    # de auditoria y lo exigen las buenas practicas de laboratorio clinico para
    # los componentes informaticos.
    printf 'user\t%s\n'             "$(whoami 2>/dev/null || echo '—')"
    printf 'host\t%s\n'             "$(hostname 2>/dev/null || echo '—')"
    printf 'workflow_md5\t%s\n'     "$(md5_of "${BASH_SOURCE[0]}")"
    printf 'workflow_git_commit\t%s\n' \
        "$(git -C "$SCRIPT_DIR" rev-parse --short HEAD 2>/dev/null || echo '—')"
    printf 'started_at\t%s\n'      "$(date '+%Y-%m-%d %H:%M:%S')"
} > "${OUTPUT}/run_params.tsv"

# ── Versiones de las herramientas ──────────────────────────────────────────
# Sin esto es imposible reconstruir que produjo exactamente un resultado. Se ha
# demostrado que cambiar la version de una herramienta del pipeline altera la
# lista de genes diferenciales (Wessels Perelo et al., NAR Genom Bioinform 2024),
# asi que la version es parte del resultado, no un detalle de instalacion.
log "Registrando versiones de las herramientas..."
{
    printf 'tool\tversion\tpath\n'
    for t in fastqc fastp multiqc bowtie2 samtools featureCounts salmon kallisto; do
        if command -v "$t" >/dev/null 2>&1; then
            # Dos cuidados aqui, los dos aprendidos ejecutando.
            #
            # Primero, NO se canaliza hacia `head`: con `set -o pipefail`, una
            # herramienta que sigue escribiendo despues de la primera linea
            # recibe SIGPIPE cuando head cierra, la tuberia devuelve 141 y
            # `set -e` aborta el script. Ocurria con bowtie2, que imprime nueve
            # lineas: la ejecucion moria al registrar versiones, antes de
            # alinear nada.
            #
            # Segundo, no todas aceptan `--version`: featureCounts usa `-v` y
            # kallisto usa `version`. Y bowtie2 antepone dos lineas de aviso a
            # la version real. Coger la primera linea a ciegas registraba un
            # aviso o un mensaje de error en el fichero de procedencia, que es
            # justo el fichero que debe poder creerse.
            case "$t" in
                featureCounts) v_full="$("$t" -v 2>&1 || true)" ;;
                kallisto)      v_full="$("$t" version 2>&1 || true)" ;;
                *)             v_full="$("$t" --version 2>&1 || true)" ;;
            esac
            # Primera linea no vacia que no sea un aviso ni un error, en UNA
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
# detectar que una entrada cambio despues de analizarla.
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
log "Building ${ALIGNMENT_TYPE} index..."

if [[ "$ALIGNMENT_TYPE" == "bowtie2" ]]; then
    if [[ -s "${BOWTIE_INDEX}.1.bt2" && -s "${BOWTIE_INDEX}.2.bt2" && -s "${BOWTIE_INDEX}.3.bt2" && -s "${BOWTIE_INDEX}.4.bt2" && -s "${BOWTIE_INDEX}.rev.1.bt2" && -s "${BOWTIE_INDEX}.rev.2.bt2" ]]; then
        log "Indice Bowtie2 existente detectado; se reutiliza."
    else
        rm -f "${BOWTIE_INDEX}".*.bt2.tmp "${BOWTIE_INDEX}".*.bt2
        log "Construyendo indice Bowtie2. Este paso puede tardar varios minutos y puede no imprimir progreso continuo."
        run_cmd bowtie2-build "$GENOME_FILE" "$BOWTIE_INDEX"
    fi
    INDEX="$BOWTIE_INDEX"
elif [[ "$ALIGNMENT_TYPE" == "salmon" ]]; then
    if [[ -d "$SALMON_INDEX" && -s "${SALMON_INDEX}/versionInfo.json" ]]; then
        log "Indice Salmon existente detectado; se reutiliza."
    else
        rm -rf "$SALMON_INDEX"
        SALMON_INDEX_TYPE=$(salmon_index_type)
        log "Using Salmon index type: $SALMON_INDEX_TYPE"
        run_cmd salmon index -t "$GENOME_FILE" -i "$SALMON_INDEX" --type "$SALMON_INDEX_TYPE" -k 31
    fi
elif [[ "$ALIGNMENT_TYPE" == "kallisto" ]]; then
    if [[ -s "$KALLISTO_INDEX" ]]; then
        log "Indice kallisto existente detectado; se reutiliza."
    else
        rm -f "$KALLISTO_INDEX"
        run_cmd kallisto index -i "$KALLISTO_INDEX" "$GENOME_FILE"
    fi
fi

# Initial quality control
log "Control de calidad inicial con FastQC..."
run_cmd fastqc "${FASTQ_FILES[@]}" -t "$THREADS" -o "$QC"

if [[ "$READ_TYPE" == "pe" ]]; then
    SAMPLE_FILES=( "${R1_FILES[@]}" )
else
    SAMPLE_FILES=( "${SINGLE_FILES[@]}" )
fi

# Process each sample
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
log "Control de calidad post-trimming con FastQC..."
run_cmd fastqc "${TRIMMED_FASTQ[@]}" -t "$THREADS" -o "$QC"

# ── Orientacion de la libreria ─────────────────────────────────────────────
# featureCounts contaba SIEMPRE como no orientada (-s 0). La mayoria de los
# protocolos actuales son stranded (dUTP, que corresponde a -s 2), y contar una
# libreria stranded como no orientada suma las lecturas antisentido: en genomas
# de alta densidad genica, como los procariotas, eso infla los conteos de genes
# solapantes y degrada la especificidad (Zhao et al., BMC Genomics 2015).
#
# En modo "auto" se infiere de los propios datos contando un subconjunto de
# lecturas con las tres orientaciones y quedandose con la que asigna mas.
infer_strandedness() {
    local bam="$1" s0 s1 s2 best
    local tmp="${COUNTS}/.strand_check"
    local -a pe_flags=()
    # featureCounts ABORTA con "Paired-end reads were detected in single-end read
    # library" si recibe un BAM pareado sin -p. La inferencia lo llamaba sin esos
    # flags y con `|| true`, asi que el error se tragaba, el fichero de resumen no
    # se escribia y las tres orientaciones salian 0: la funcion devolvia siempre
    # 0 (sin orientar) creyendo haberlo medido.
    [[ "$READ_TYPE" == "pe" ]] && pe_flags=(-p --countReadPairs)
    mkdir -p "$tmp"
    for s in 0 1 2; do
        featureCounts -T "$THREADS" "${pe_flags[@]}" -s "$s" \
            -t "$FEATURE_TYPE" -g "$FEATURE_ATTR" \
            -a "$ANNOTATION_FILE" -o "${tmp}/s${s}.txt" "$bam" >/dev/null 2>&1 || true
    done
    s0=$(awk 'NR>1 && $1=="Assigned" {print $2}' "${tmp}/s0.txt.summary" 2>/dev/null || echo 0)
    s1=$(awk 'NR>1 && $1=="Assigned" {print $2}' "${tmp}/s1.txt.summary" 2>/dev/null || echo 0)
    s2=$(awk 'NR>1 && $1=="Assigned" {print $2}' "${tmp}/s2.txt.summary" 2>/dev/null || echo 0)
    s0=${s0:-0}; s1=${s1:-0}; s2=${s2:-0}
    # El diagnostico va a STDERR: esta funcion se invoca con $(...) y su stdout es
    # el valor de retorno. Escribiendolo en stdout, el mensaje acababa DENTRO de
    # la variable y `-s` recibia una linea de log entera.
    log "  asignadas por orientacion -> sin orientar: $s0 | directa: $s1 | inversa: $s2" >&2
    # La comparacion NO puede ser "la que mas asigna". El modo sin orientar cuenta
    # las lecturas en ambos sentidos, asi que su recuento es por construccion mayor
    # o igual que el de cualquiera de los dos modos orientados: ganaria siempre, y
    # la inferencia devolveria "sin orientar" para cualquier libreria, tambien para
    # una dUTP perfectamente orientada.
    #
    # El criterio correcto es el que usan las herramientas al uso (RSeQC): mirar
    # como se REPARTEN entre los dos sentidos las lecturas que si son asignables
    # por hebra. Si el reparto es muy asimetrico, la libreria esta orientada en ese
    # sentido; si esta repartido, no lo esta.
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
    # Guarda final: si nada se pudo medir, no se puede afirmar una orientacion.
    if [[ "$s0" -eq 0 && "$s1" -eq 0 && "$s2" -eq 0 ]]; then
        log "  ADVERTENCIA: no se pudo medir la orientacion; se usa 0 (sin orientar)." >&2
    fi
    printf '%s' "$best"
}

if [[ "$ALIGNMENT_TYPE" == "bowtie2" && "$STRANDEDNESS" == "auto" ]]; then
    FIRST_BAM=$( (shopt -s nullglob; b=( "$ALIGNMENTS"/*.bam ); echo "${b[0]:-}") )
    if [[ -n "$FIRST_BAM" && -f "$FIRST_BAM" ]]; then
        log "Infiriendo la orientacion de la libreria a partir de $(basename "$FIRST_BAM")..."
        STRANDEDNESS=$(infer_strandedness "$FIRST_BAM")
        log "+ Orientacion inferida: -s $STRANDEDNESS"
    else
        STRANDEDNESS=0
    fi
elif [[ "$STRANDEDNESS" == "auto" ]]; then
    # En pseudoalineamiento la detecta el propio cuantificador (salmon -l A).
    STRANDEDNESS=0
fi
printf 'strandedness\t%s\n' "$STRANDEDNESS" >> "${OUTPUT}/run_params.tsv"

# Generate count matrix based on alignment type
if [[ "$ALIGNMENT_TYPE" == "bowtie2" ]]; then
    log "Generating count matrix from Bowtie2 alignments..."
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
    # La cabecera se limpia CAMPO A CAMPO, no con un sed sobre la linea entera.
    # El `sed '1s|.*/||g'` que habia aqui es codicioso sobre toda la linea: con
    # rutas absolutas colapsaba `Geneid\t/ruta/A.bam\t/ruta/B.bam` en un unico
    # campo `B`, dejando una cabecera de 1 columna frente a N+1 de datos. El
    # lector de la app no encontraba el patron `^Geneid\t`, read.delim abortaba y
    # devolvia NULL, y la interfaz anunciaba "Workflow finalizado correctamente"
    # sin matriz. Con awk cada campo se procesa por separado.
    awk 'BEGIN { FS=OFS="\t" } !/^#/ { print }' "${COUNTS}/gene_counts.txt" \
        | cut -f1,7- \
        | awk 'BEGIN { FS=OFS="\t" }
               NR==1 { for (i=2; i<=NF; i++) { sub(/.*\//, "", $i); sub(/\.bam$/, "", $i) } }
               { print }' \
        > "${COUNTS}/count_matrix.tsv"

elif [[ "$ALIGNMENT_TYPE" == "salmon" || "$ALIGNMENT_TYPE" == "kallisto" ]]; then
    # Resumen de transcritos a genes con tximport, reutilizando la misma funcion
    # que usa la aplicacion (scripts/build_count_matrix.R). Antes se escribia un
    # fichero con SOLO la cabecera y un comentario diciendo que hacia falta un
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
        # Sin anotacion utilizable no hay resumen a gen posible. Se avisa y NO se
        # deja un fichero a medias: es preferible que falte el artefacto a que
        # exista uno que parece valido y no lo es.
        rm -f "${COUNTS}/count_matrix.tsv"
        log "! No se ha podido construir la matriz por gen (revisa la anotacion)."
        log "  Las cuantificaciones por muestra quedan en ${COUNTS}/ y la app"
        log "  puede resumirlas al cargar la ejecucion."
    fi

    # Copia de las cuantificaciones por muestra, con el nombre de la muestra en
    # el fichero. Antes se copiaban con `cp "$ALIGNMENTS"/*/*.sf "${COUNTS}/"`:
    # todas se llaman quant.sf (o abundance.tsv), asi que cada copia sobrescribia
    # la anterior y solo sobrevivia la ultima muestra, en silencio por el
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
log "Generando informe MultiQC..."
run_cmd multiqc "$OUTPUT" -o "$OUTPUT"

log "Analysis completed successfully."
log "Alignment type: $ALIGNMENT_TYPE"
log "Results saved in: $OUTPUT"
log "Quality control: $QC"
log "Trimmed reads: $TRIMMED"
log "Alignments: $ALIGNMENTS"
log "Count matrices: $COUNTS"
