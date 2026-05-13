

#!/bin/bash

set -euo pipefail

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

# Usage: workflow.sh --INPUT <path> --OUTPUT <path> --GENOME_FILE <file> --ANNOTATION_FILE <file> [--ALIGNMENT_TYPE bowtie2|salmon|kallisto]

# Default values
ALIGNMENT_TYPE="bowtie2"

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --INPUT) INPUT="$2"; shift 2 ;;
        --OUTPUT) OUTPUT="$2"; shift 2 ;;
        --GENOME_FILE) GENOME_FILE="$2"; shift 2 ;;
        --ANNOTATION_FILE) ANNOTATION_FILE="$2"; shift 2 ;;
        --ALIGNMENT_TYPE) ALIGNMENT_TYPE="$2"; shift 2 ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
done

# Validate required arguments
if [[ -z "${INPUT:-}" || -z "${OUTPUT:-}" || -z "${GENOME_FILE:-}" || -z "${ANNOTATION_FILE:-}" ]]; then
    echo "Error: Missing required arguments."
    echo "Usage: pipeline_ecoli.sh --INPUT <path> --OUTPUT <path> --GENOME_FILE <file> --ANNOTATION_FILE <file> [--ALIGNMENT_TYPE bowtie2|salmon|kallisto]"
    exit 1
fi

# Validate alignment type
if [[ ! "$ALIGNMENT_TYPE" =~ ^(bowtie2|salmon|kallisto)$ ]]; then
    echo "Error: Invalid ALIGNMENT_TYPE. Choose: bowtie2, salmon, or kallisto"
    exit 1
fi

# Number of threads to accelerate analysis
THREADS=8

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

# Build genome index based on alignment type
echo "Building ${ALIGNMENT_TYPE} index..."

if [[ "$ALIGNMENT_TYPE" == "bowtie2" ]]; then
    bowtie2-build "$GENOME_FILE" "$BOWTIE_INDEX"
    INDEX="$BOWTIE_INDEX"
elif [[ "$ALIGNMENT_TYPE" == "salmon" ]]; then
    salmon index -t "$GENOME_FILE" -i "$SALMON_INDEX" --type quasi -k 31
elif [[ "$ALIGNMENT_TYPE" == "kallisto" ]]; then
    kallisto index -i "$KALLISTO_INDEX" "$GENOME_FILE"
fi

# Initial quality control
fastqc "$INPUT"/*.fastq.gz -t "$THREADS" -o "$QC"

# Process each sample
for R1 in "$INPUT"/*_1.fastq.gz "$INPUT"/*_R1.fastq.gz; do

    [ -e "$R1" ] || continue

    SAMPLE=$(basename "$R1")
    SAMPLE=${SAMPLE%%_1.fastq.gz}
    SAMPLE=${SAMPLE%%_R1.fastq.gz}

    R2="${R1/_1.fastq.gz/_2.fastq.gz}"
    R2="${R2/_R1.fastq.gz/_R2.fastq.gz}"

    echo "Processing sample: $SAMPLE"

    # Adapter trimming and low-quality read filtering
    fastp \
        -i "$R1" \
        -I "$R2" \
        -o "${TRIMMED}/${SAMPLE}_R1_trimmed.fastq.gz" \
        -O "${TRIMMED}/${SAMPLE}_R2_trimmed.fastq.gz" \
        --detect_adapter_for_pe \
        --thread "$THREADS" \
        --html "${TRIMMED}/${SAMPLE}_fastp.html" \
        --json "${TRIMMED}/${SAMPLE}_fastp.json"

    # Perform alignment based on alignment type
    if [[ "$ALIGNMENT_TYPE" == "bowtie2" ]]; then
        echo "  Aligning with Bowtie2..."
        # Alignment against reference genome
        bowtie2 \
            -x "$INDEX" \
            -1 "${TRIMMED}/${SAMPLE}_R1_trimmed.fastq.gz" \
            -2 "${TRIMMED}/${SAMPLE}_R2_trimmed.fastq.gz" \
            -p "$THREADS" \
        | samtools sort \
            -@ "$THREADS" \
            -o "${ALIGNMENTS}/${SAMPLE}.bam"

        # Index BAM file for downstream analysis
        samtools index "${ALIGNMENTS}/${SAMPLE}.bam"

    elif [[ "$ALIGNMENT_TYPE" == "salmon" ]]; then
        echo "  Quasi-aligning with Salmon..."
        salmon quant \
            -i "$SALMON_INDEX" \
            -l A \
            -1 "${TRIMMED}/${SAMPLE}_R1_trimmed.fastq.gz" \
            -2 "${TRIMMED}/${SAMPLE}_R2_trimmed.fastq.gz" \
            -p "$THREADS" \
            -o "${ALIGNMENTS}/${SAMPLE}"

    elif [[ "$ALIGNMENT_TYPE" == "kallisto" ]]; then
        echo "  Quasi-aligning with kallisto..."
        kallisto quant \
            -i "$KALLISTO_INDEX" \
            -o "${ALIGNMENTS}/${SAMPLE}" \
            -t "$THREADS" \
            "${TRIMMED}/${SAMPLE}_R1_trimmed.fastq.gz" \
            "${TRIMMED}/${SAMPLE}_R2_trimmed.fastq.gz"

    fi

done

# Quality control after trimming
fastqc "$TRIMMED"/*.fastq.gz -t "$THREADS" -o "$QC"

# Generate count matrix based on alignment type
if [[ "$ALIGNMENT_TYPE" == "bowtie2" ]]; then
    echo "Generating count matrix from Bowtie2 alignments..."
    # Count reads per gene
    featureCounts \
        -T "$THREADS" \
        -p \
        --countReadPairs \
        -s 0 \
        -t gene \
        -g locus_tag \
        -a "$ANNOTATION_FILE" \
        -o "${COUNTS}/gene_counts.txt" \
        "$ALIGNMENTS"/*.bam

    # Create count matrix in TSV format
    cut -f1,7- "${COUNTS}/gene_counts.txt" \
        | sed '1s|.bam||g; 1s|.*/||g' \
        > "${COUNTS}/count_matrix.tsv"

elif [[ "$ALIGNMENT_TYPE" == "salmon" ]]; then
    echo "Generating count matrix from Salmon quantifications..."
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
    echo "Generating count matrix from kallisto abundances..."
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
multiqc "$OUTPUT" -o "$OUTPUT"

echo "Analysis completed successfully."
echo "Alignment type: $ALIGNMENT_TYPE"
echo "Results saved in: $OUTPUT"
echo "Quality control: $QC"
echo "Trimmed reads: $TRIMMED"
echo "Alignments: $ALIGNMENTS"
echo "Count matrices: $COUNTS"