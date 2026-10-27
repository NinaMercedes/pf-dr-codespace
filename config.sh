# =====================================================================
#  Workshop settings. Every step script reads this file.
#  You can change values here, or override any of them for one run:
#     SUBSAMPLE_PAIRS=2000000 bash steps/01_get_reads.sh
# =====================================================================

# --- Sample -----------------------------------------------------------
SAMPLE="${SAMPLE:-ERR018972}"
READ1_URL="${READ1_URL:-https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR018/ERR018972/ERR018972_1.fastq.gz}"
READ2_URL="${READ2_URL:-https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR018/ERR018972/ERR018972_2.fastq.gz}"

# Set to a number (e.g. 2000000) to keep only that many read pairs.
# Empty = use all reads. Subsampling makes every step faster but lowers
# coverage, which matters for mixed-infection and copy-number calls.
SUBSAMPLE_PAIRS="${SUBSAMPLE_PAIRS:-}"

# --- Reference (P. falciparum 3D7 from PlasmoDB) ----------------------
# Leave empty to download the current PlasmoDB release automatically.
REF_FASTA_URL="${REF_FASTA_URL:-}"
REF_GFF_URL="${REF_GFF_URL:-}"
PLASMODB_FALLBACK_RELEASE="${PLASMODB_FALLBACK_RELEASE:-68}"
# Or point at local files (used by the instructor self-test)
REF_FASTA_LOCAL="${REF_FASTA_LOCAL:-}"
REF_GFF_LOCAL="${REF_GFF_LOCAL:-}"

# --- Host removal -----------------------------------------------------
# Name of a standard hostile index, or a path to a custom Bowtie2 index
HOSTILE_INDEX="${HOSTILE_INDEX:-human-t2t-hla}"

# --- Variant calling --------------------------------------------------
# "dr"     = call only in the drug-resistance genes (+ padding) - fast
# "genome" = call across the whole genome - slower, closer to real life
CALL_REGIONS="${CALL_REGIONS:-dr}"
INTERVAL_PADDING="${INTERVAL_PADDING:-1000}"
# Pf blood-stage parasites are haploid, but patient samples are often
# mixed infections. Calling with ploidy 2 lets GATK report the second clone
# as a "heterozygous" genotype (this is what MalariaGEN does).
PLOIDY="${PLOIDY:-2}"

# --- Reporting --------------------------------------------------------
MIN_DP="${MIN_DP:-5}"            # min depth to make a genotype call at a codon
MIXED_MIN_AF="${MIXED_MIN_AF:-0.10}"  # minor allele fraction to report as mixed

# --- Compute ----------------------------------------------------------
THREADS="${THREADS:-$(nproc)}"
JAVA_MEM="${JAVA_MEM:-6g}"

# --- Paths (no need to change) ----------------------------------------
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF_DIR="${REF_DIR:-$REPO_DIR/ref}"
DATA_DIR="${DATA_DIR:-$REPO_DIR/data}"
RES_DIR="${RES_DIR:-$REPO_DIR/results}"
REF="$REF_DIR/Pf3D7.fasta"
GFF="$REF_DIR/Pf3D7.gff"
DR_GENES_TSV="${DR_GENES_TSV:-$REPO_DIR/resources/dr_genes.tsv}"
DR_MARKERS_TSV="${DR_MARKERS_TSV:-$REPO_DIR/resources/dr_markers.tsv}"
DR_BED="$REF_DIR/dr_genes.bed"
DR_CDS_JSON="$REF_DIR/dr_genes_cds.json"
