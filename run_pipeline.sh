#!/usr/bin/env bash
set -euo pipefail

# Default python executables (can override via options)
PYTHON=${PYTHON:-python}
PYTHON3=${PYTHON3:-python3}

usage() {
  cat <<EOF
Usage: $0 --prokka_input <path> --prokka_out <path> --uniref90 <path> --uniref50 <path> --index_dir <path> --index <name> --tax_infofile <path> --gmb_pkl <path> --gmb_tax_tsv <path> --gmb_fasta <path> --bt2_threads <num> [--genome_ext fna|fasta] [--threads_uniref <num>] [--threads_bowtie <num>] [--py <python_exec>] [--py3 <python3_exec>]

Required:
  --prokka_input     Path passed to 1-prokka_batch.py -i
  --prokka_out       Base output directory for prokka (will contain ffn_outputs, uniref-annotated-result, core-gene-select, bowtie2-result, filtered-core-gene, etc.)
  --uniref90         Path to uniref90 database file (dmnd) for script 2
  --uniref50         Path to uniref50 database file (dmnd) for script 2
  --index_dir        Path to bowtie2 index directory used by script 5
  --index            Index name used by script 5 (e.g. test-20251116_NEW)
  --tax_infofile     Taxonomy info file used by script 6 (and as new_tax_tsv for script 7)
  --gmb_pkl          Path to gmb pkl for script 7
  --gmb_tax_tsv      Path to gmb tax tsv for script 7
  --gmb_fasta        Path to gmb fasta for script 7
  --bt2_threads      Number for --bt2_threads in script 7

Optional:
  --genome_ext       Genome file extension in prokka_input (fna or fasta). Default: fna
  --threads_uniref   Threads for script 2 (default 60)
  --threads_bowtie   Threads for script 5 (default 30)
  --py               python executable to use (overrides PYTHON)
  --py3              python3 executable to use (overrides PYTHON3)

Example:
  $0 --prokka_input ../test-4-genome/ --prokka_out ../test-4-genome/prokka-result/ \\
     --uniref90 /data/database/humann3/uniref/uniref90/uniref90_201901b.dmnd \\
     --uniref50 /data/database/humann3/uniref/uniref50/uniref50_201901b_full.dmnd \\
     --index_dir /mnt/ME5012-Vol01/chenc/MetaphIan/7-seventh-all-maker-gene-db/20251116-filtered-gene/ \\
     --index test-20251116_NEW --tax_infofile ../test-4-genome/genus1_species1.tsv \\
     --gmb_pkl ../../GMB-db/test-20251116_NEW.pkl --gmb_tax_tsv ../../GMB-db/new-all-taxinfo.modified-20251111.tsv \\
     --gmb_fasta ../../GMB-db/test-20251116.fasta --bt2_threads 200
EOF
  exit 1
}

# parse args
if [ $# -eq 0 ]; then usage; fi

# defaults
THREADS_UNIREF=60
THREADS_BOWTIE=30
GENOME_EXT="fna"

while [ "$#" -gt 0 ]; do
  case "$1" in
    --prokka_input) PROKKA_INPUT="$2"; shift 2;;
    --prokka_out) PROKKA_OUT="$2"; shift 2;;
    --uniref90) UNIREF90="$2"; shift 2;;
    --uniref50) UNIREF50="$2"; shift 2;;
    --index_dir) INDEX_DIR="$2"; shift 2;;
    --index) INDEX_NAME="$2"; shift 2;;
    --tax_infofile) TAX_INFOFILE="$2"; shift 2;;
    --gmb_pkl) GMB_PKL="$2"; shift 2;;
    --gmb_tax_tsv) GMB_TAX_TSV="$2"; shift 2;;
    --gmb_fasta) GMB_FASTA="$2"; shift 2;;
    --bt2_threads) BT2_THREADS="$2"; shift 2;;
    --threads_uniref) THREADS_UNIREF="$2"; shift 2;;
    --threads_bowtie) THREADS_BOWTIE="$2"; shift 2;;
    --genome_ext) GENOME_EXT="$2"; shift 2;;
    --py) PYTHON="$2"; shift 2;;
    --py3) PYTHON3="$2"; shift 2;;
    --help|-h) usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

# check required
: "${PROKKA_INPUT:?missing --prokka_input}"
: "${PROKKA_OUT:?missing --prokka_out}"
: "${UNIREF90:?missing --uniref90}"
: "${UNIREF50:?missing --uniref50}"
: "${INDEX_DIR:?missing --index_dir}"
: "${INDEX_NAME:?missing --index}"
: "${TAX_INFOFILE:?missing --tax_infofile}"
: "${GMB_PKL:?missing --gmb_pkl}"
: "${GMB_TAX_TSV:?missing --gmb_tax_tsv}"
: "${GMB_FASTA:?missing --gmb_fasta}"
: "${BT2_THREADS:?missing --bt2_threads}"

# validate genome ext
if [ "${GENOME_EXT}" != "fna" ] && [ "${GENOME_EXT}" != "fasta" ]; then
  echo "Error: --genome_ext must be 'fna' or 'fasta'. Got: ${GENOME_EXT}"
  exit 1
fi

# normalize paths (remove trailing slash for PROKKA_OUT for consistent concatenation)
PROKKA_OUT=${PROKKA_OUT%/}
FFN_DIR="${PROKKA_OUT}/ffn_outputs"
UNIREF_OUT="${PROKKA_OUT}/uniref-annotated-result"
CORE_SELECT_OUT="${PROKKA_OUT}/core-gene-select"
BOWTIE_OUT="${PROKKA_OUT}/bowtie2-result"
FILTERED_CORE_OUT="${PROKKA_OUT}/filtered-core-gene"
MERGED_FOLDER="${BOWTIE_OUT}/${INDEX_NAME}_bowtie_result/"

mkdir -p "${PROKKA_OUT}" "${FFN_DIR}" "${UNIREF_OUT}" "${CORE_SELECT_OUT}" "${BOWTIE_OUT}" "${FILTERED_CORE_OUT}"

echo "===== Pipeline run started ====="
echo "Prokka input: ${PROKKA_INPUT}"
echo "Prokka base out: ${PROKKA_OUT}"
echo "Genome extension check: ${GENOME_EXT}"
echo "Uniref90: ${UNIREF90}"
echo "Uniref50: ${UNIREF50}"
echo "Bowtie2 index dir: ${INDEX_DIR}"
echo "Index name: ${INDEX_NAME}"
echo "Tax info file: ${TAX_INFOFILE}"
echo "GMB pkl: ${GMB_PKL}"
echo "GMB tax tsv: ${GMB_TAX_TSV}"
echo "GMB fasta: ${GMB_FASTA}"
echo "bt2_threads (for script7): ${BT2_THREADS}"
echo "threads_uniref (script2): ${THREADS_UNIREF}"
echo "threads_bowtie (script5): ${THREADS_BOWTIE}"
echo

# Check file formats in prokka_input directory
# Count chosen ext, and detect presence of the other ext to enforce only-one-format rule
OTHER_EXT="fasta"
if [ "${GENOME_EXT}" = "fasta" ]; then
  OTHER_EXT="fna"
fi

# Use find -maxdepth 1 to count files in top-level of PROKKA_INPUT. Adjust if recursive wanted.
COUNT_CHOSEN=$(find "${PROKKA_INPUT%/}" -maxdepth 1 -type f -iname "*.${GENOME_EXT}" | wc -l)
COUNT_OTHER=$(find "${PROKKA_INPUT%/}" -maxdepth 1 -type f -iname "*.${OTHER_EXT}" | wc -l)

echo "Found ${COUNT_CHOSEN} *.${GENOME_EXT} files in ${PROKKA_INPUT}"
if [ "${COUNT_OTHER}" -gt 0 ]; then
  echo "Error: Found ${COUNT_OTHER} *.${OTHER_EXT} files in ${PROKKA_INPUT} while --genome_ext is set to ${GENOME_EXT}."
  echo "Please ensure only one genome file format exists or set --genome_ext appropriately."
  exit 1
fi

# Decide path flow based on genome count
if [ "${COUNT_CHOSEN}" -gt 2 ]; then
  echo "Genome count > 2 (${COUNT_CHOSEN}) => run full pipeline (steps 1..7)."
  RUN_FULL=1
else
  echo "Genome count <= 2 (${COUNT_CHOSEN}) => run steps 1..3, then concatenate ffn -> simulate reads (alternative), then steps 5..7."
  RUN_FULL=0
fi

# 1) run 1-prokka_batch.py
echo "----- Step 1: Prokka batch -----"
CMD1="${PYTHON} 1-prokka_batch.py -i ${PROKKA_INPUT} -o ${PROKKA_OUT}"
echo "$CMD1"
eval "$CMD1"
echo "Step 1 finished"
echo

# 2) run 2-uniref_annotate_batch.py
echo "----- Step 2: Uniref annotate batch -----"
CMD2="${PYTHON} 2-uniref_annotate_batch.py -i ${FFN_DIR} -o ${UNIREF_OUT} --threads ${THREADS_UNIREF} --uniref90 ${UNIREF90} --uniref50 ${UNIREF50}"
echo "$CMD2"
eval "$CMD2"
echo "Step 2 finished"
echo

# 3) run 3-mmseqs_modify_perform-new-round-uniclustering.py
echo "----- Step 3: MMseqs modify / perform clustering -----"
CMD3="${PYTHON} 3-mmseqs_modify_perform-new-round-uniclustering.py -f ${UNIREF_OUT}"
echo "$CMD3"
eval "$CMD3"
echo "Step 3 finished"
echo

if [ "${RUN_FULL}" -eq 1 ]; then
  # 4) run 4-select-coreness_and_simulate.py
  echo "----- Step 4: Select coreness and simulate (original) -----"
  CMD4="${PYTHON} 4-select-coreness_and_simulate.py -i ${UNIREF_OUT} -o ${CORE_SELECT_OUT} --paired"
  echo "$CMD4"
  eval "$CMD4"
  echo "Step 4 finished"
  echo
else
  # When <=2 genomes: construct core_clusters.fasta by concatenating all .ffn files under UNIREF_OUT subfolders
  echo "----- Special small-sample flow: concatenate FFNs to core_clusters.fasta -----"
  # Find all .ffn files under UNIREF_OUT recursively
  CORE_FASTA="${CORE_SELECT_OUT}/core_clusters.fasta"
  rm -f "${CORE_FASTA}"
  mkdir -p "${CORE_SELECT_OUT}"
  # cat all *.ffn files from subdirectories
  # If none found, exit with error
  FFN_COUNT=$(find "${UNIREF_OUT}" -type f -iname "*.ffn" | wc -l)
  if [ "${FFN_COUNT}" -eq 0 ]; then
    echo "Error: No .ffn files found under ${UNIREF_OUT}. Cannot create core_clusters.fasta."
    exit 1
  fi
  echo "Found ${FFN_COUNT} .ffn files. Concatenating into ${CORE_FASTA} ..."
  # Use a deterministic ordering
  find "${UNIREF_OUT}" -type f -iname "*.ffn" | sort | while read -r f; do
    echo "cat ${f} >> ${CORE_FASTA}"
    cat "${f}" >> "${CORE_FASTA}"
  done
  echo "Concatenation finished. core_clusters.fasta saved to ${CORE_FASTA}"
  echo

  # Run simulate reads alternative script
  echo "----- Step 4-alt: Simulate reads (alternative) -----"
  CMD4A="${PYTHON} 4-simulated_reads_alternative.py -i ${CORE_SELECT_OUT}"
  echo "$CMD4A"
  eval "$CMD4A"
  echo "Step 4-alt finished"
  echo
fi

# 5) run 5-bowtie2_pipeline.py
echo "----- Step 5: Bowtie2 pipeline -----"
CMD5="${PYTHON} 5-bowtie2_pipeline.py --fq_dir ${CORE_SELECT_OUT} --index_dir ${INDEX_DIR} --index ${INDEX_NAME} --txt_out_dir ${BOWTIE_OUT} --sam_out_dir ${BOWTIE_OUT} --threads_bowtie ${THREADS_BOWTIE}"
echo "$CMD5"
eval "$CMD5"
echo "Step 5 finished"
echo

# 6) run 6-filter_core_gene_pipeline_dedup_rename-test.py
echo "----- Step 6: Filter core gene pipeline -----"
CMD6="${PYTHON3} 6-filter_core_gene_pipeline_dedup_rename-test.py --merged_folder ${MERGED_FOLDER} --core_fasta_dir ${CORE_SELECT_OUT} --tax_infofile ${TAX_INFOFILE} --outdir ${FILTERED_CORE_OUT}"
echo "$CMD6"
eval "$CMD6"
echo "Step 6 finished"
echo

# 7) run 7-add_markers.py
echo "----- Step 7: Add markers and build new db -----"
CMD7="${PYTHON} 7-add_markers.py --new_tax_tsv ${TAX_INFOFILE} --ffn_dir ${FILTERED_CORE_OUT} --gmb_pkl ${GMB_PKL} --gmb_tax_tsv ${GMB_TAX_TSV} --gmb_fasta ${GMB_FASTA} --out_dir ../../GMB-db/new-db-test-1 --bt2_threads ${BT2_THREADS}"
echo "$CMD7"
eval "$CMD7"
echo "Step 7 finished"
echo

echo "===== Pipeline completed successfully ====="
