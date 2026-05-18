#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

cat <<'EOF'
------------------------------------------------------------
TFM_JMQ_RNA_seq - comprobación de requerimientos
------------------------------------------------------------
EOF

# ---------------------------------------------
# 1) Comprobar programas CLI requeridos
# ---------------------------------------------
cli_programs=(
  Rscript
  bash
  bowtie2
  samtools
  fastqc
  fastp
  multiqc
  featureCounts
  salmon
  kallisto
)

missing_cli=()
for prog in "${cli_programs[@]}"; do
  if ! command -v "$prog" >/dev/null 2>&1; then
    missing_cli+=("$prog")
  fi
done

if [ ${#missing_cli[@]} -gt 0 ]; then
  printf "\nFaltan programas de sistema:\n"
  for prog in "${missing_cli[@]}"; do
    printf "  - %s\n" "$prog"
  done
  printf '\nInstala estos programas con tu gestor de paquetes preferido (conda/bioconda, brew, apt, etc.).\n'
else
  printf "\nTodos los programas CLI requeridos están disponibles.\n"
fi

# ---------------------------------------------
# 2) Comprobar e instalar paquetes R
# ---------------------------------------------
cat <<'EOF'

Comprobando paquetes R...
EOF

Rscript - <<'RS'
required_cran <- c(
  "shiny",
  "bslib",
  "shinyjs",
  "DT",
  "ggplot2",
  "plotly",
  "pheatmap",
  "ggrepel",
  "RColorBrewer",
  "processx"
)
required_bioc <- c(
  "DESeq2",
  "tximport",
  "clusterProfiler",
  "enrichplot"
)
installed <- rownames(installed.packages())
missing_cran <- setdiff(required_cran, installed)
missing_bioc <- setdiff(required_bioc, installed)

if (length(missing_cran) > 0) {
  message("Instalando paquetes CRAN: ", paste(missing_cran, collapse=", "))
  install.packages(missing_cran, repos="https://cloud.r-project.org")
} else {
  message("Todos los paquetes CRAN requeridos ya están instalados.")
}

if (length(missing_bioc) > 0) {
  if (!requireNamespace("BiocManager", quietly=TRUE)) {
    install.packages("BiocManager", repos="https://cloud.r-project.org")
  }
  message("Instalando paquetes Bioconductor: ", paste(missing_bioc, collapse=", "))
  BiocManager::install(missing_bioc, ask=FALSE, update=FALSE)
} else {
  message("Todos los paquetes Bioconductor requeridos ya están instalados.")
}

message("\nRevisión de paquetes R completada.")
RS

# ---------------------------------------------
# 3) Comprobar que workflow.sh existe
# ---------------------------------------------
if [ ! -f "$ROOT_DIR/workflow.sh" ]; then
  printf "\nERROR: no se encontró workflow.sh en %s\n" "$ROOT_DIR"
  exit 1
fi

printf "\nworkflow.sh encontrado.\n"
printf "\nListo. Puedes ejecutar la aplicación con: Rscript app.R\n"
