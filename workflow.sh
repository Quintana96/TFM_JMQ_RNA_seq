#!/bin/bash

set -euo pipefail
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

# Usage: workflow.sh --INPUT <path> --OUTPUT <path> --GENOME_FILE <file> --ANNOTATION_FILE <file> [--ALIGNMENT_TYPE bowtie2|salmon|kallisto] [--READ_TYPE pe|se] [--FRAGMENT_LENGTH <mean>] [--FRAGMENT_SD <sd>]

# Default values
ALIGNMENT_TYPE="bowtie2"
READ_TYPE="pe"
FRAGMENT_LENGTH=200
FRAGMENT_SD=20

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
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
done

# Validate required arguments
if [[ -z "${INPUT:-}" || -z "${OUTPUT:-}" || -z "${GENOME_FILE:-}" || -z "${ANNOTATION_FILE:-}" ]]; then
    echo "Error: Missing required arguments."
    echo "Usage: pipeline_ecoli.sh --INPUT <path> --OUTPUT <path> --GENOME_FILE <file> --ANNOTATION_FILE <file> [--ALIGNMENT_TYPE bowtie2|salmon|kallisto] [--READ_TYPE pe|se] [--FRAGMENT_LENGTH <mean>] [--FRAGMENT_SD <sd>]"
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

# Number of threads to accelerate analysis
THREADS=8

log() {
    printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"
}

run_cmd() {
    log "+ $*"
    "$@"
}

trap 'log "ERROR: fallo en la linea ${LINENO}. Revisa el comando anterior."' ERR

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
BOWTIE_INDEX="${OUTPUT}/indices/bowtie2/ecoli"
SALMON_INDEX="${OUTPUT}/indices/salmon/ecoli"
KALLISTO_INDEX="${OUTPUT}/indices/kallisto/ecoli.idx"

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

log "Entrada: $INPUT"
log "Salida: $OUTPUT"
log "Genoma/transcriptoma: $GENOME_FILE"
log "Anotacion: $ANNOTATION_FILE"
log "Tipo de alineamiento: $ALIGNMENT_TYPE"
log "Tipo de lectura: $READ_TYPE"

# Parametros en formato legible por la app. Sin esto, una ejecucion guardada no
# deja rastro de con que anotacion se hizo, y la app no puede construir el mapa
# transcrito-gen para tximport ni identificar los genes de rRNA.
{
    printf 'input_dir\t%s\n'       "$INPUT"
    printf 'output_dir\t%s\n'      "$OUTPUT"
    printf 'genome_file\t%s\n'     "$GENOME_FILE"
    printf 'annotation_file\t%s\n' "$ANNOTATION_FILE"
    printf 'alignment_type\t%s\n'  "$ALIGNMENT_TYPE"
    printf 'read_type\t%s\n'       "$READ_TYPE"
    printf 'started_at\t%s\n'      "$(date '+%Y-%m-%d %H:%M:%S')"
} > "${OUTPUT}/run_params.tsv"
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
                -o "${ALIGNMENTS}/${SAMPLE}"
        else
            run_cmd salmon quant \
                -i "$SALMON_INDEX" \
                -l A \
                -r "${TRIMMED}/${SAMPLE}_trimmed.fastq.gz" \
                -p "$THREADS" \
                -o "${ALIGNMENTS}/${SAMPLE}"
        fi

    elif [[ "$ALIGNMENT_TYPE" == "kallisto" ]]; then
        log "  Quasi-aligning with kallisto..."
        if [[ "$READ_TYPE" == "pe" ]]; then
            run_cmd kallisto quant \
                -i "$KALLISTO_INDEX" \
                -o "${ALIGNMENTS}/${SAMPLE}" \
                -t "$THREADS" \
                "${TRIMMED}/${SAMPLE}_R1_trimmed.fastq.gz" \
                "${TRIMMED}/${SAMPLE}_R2_trimmed.fastq.gz"
        else
            run_cmd kallisto quant \
                -i "$KALLISTO_INDEX" \
                -o "${ALIGNMENTS}/${SAMPLE}" \
                -t "$THREADS" \
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
            -s 0 \
            -t gene \
            -g locus_tag \
            -a "$ANNOTATION_FILE" \
            -o "${COUNTS}/gene_counts.txt" \
            "${BAM_FILES[@]}"
    else
        run_cmd featureCounts \
            -T "$THREADS" \
            -s 0 \
            -t gene \
            -g locus_tag \
            -a "$ANNOTATION_FILE" \
            -o "${COUNTS}/gene_counts.txt" \
            "${BAM_FILES[@]}"
    fi

    # Create count matrix in TSV format
    log "+ Crear matriz de conteos: ${COUNTS}/count_matrix.tsv"
    awk 'BEGIN { FS=OFS="\t" } !/^#/ { print }' "${COUNTS}/gene_counts.txt" \
        | cut -f1,7- \
        | sed '1s|.bam||g; 1s|.*/||g' \
        > "${COUNTS}/count_matrix.tsv"

elif [[ "$ALIGNMENT_TYPE" == "salmon" ]]; then
    log "Generating count matrix from Salmon quantifications..."
    # Extract TPM and counts from Salmon outputs
    {
        # Create header with sample names
        echo -n "gene_id"
        for dir in "$ALIGNMENTS"/*; do
            sample=$(basename "$dir")
            echo -n -e "\t${sample}_TPM\t${sample}_Counts"
        done
        echo
        
        # Extract counts from quant.sf files
        # This assumes Salmon produced quantifications
        # For proper count matrix, you might need custom R script
    } > "${COUNTS}/count_matrix.tsv"
    
    # Copy individual Salmon outputs for reference
    cp -r "$ALIGNMENTS"/*/*.sf "${COUNTS}/" 2>/dev/null || true

elif [[ "$ALIGNMENT_TYPE" == "kallisto" ]]; then
    log "Generating count matrix from kallisto abundances..."
    # Extract counts from kallisto outputs
    {
        # Create header with sample names
        echo -n "gene_id"
        for dir in "$ALIGNMENTS"/*; do
            sample=$(basename "$dir")
            echo -n -e "\t${sample}_TPM\t${sample}_Counts"
        done
        echo
        
        # Extract abundance data from abundance.tsv files
        # For proper count matrix, you might need custom R script
    } > "${COUNTS}/count_matrix.tsv"
    
    # Copy individual kallisto outputs for reference
    cp -r "$ALIGNMENTS"/*/abundance.tsv "${COUNTS}/" 2>/dev/null || true

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
