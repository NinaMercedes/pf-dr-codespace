# *Plasmodium falciparum* drug-resistance from WGS — a step-by-step workshop

A small, self-contained teaching pipeline that takes **raw Illumina reads from one
*P. falciparum* sample** and walks — one command at a time — through the whole
path to a **drug-resistance genotype**:

```
FastQC → trim (fastp) → remove human reads (hostile) → map (bwa mem -M)
      → mark duplicates → coverage (mosdepth) → call variants (GATK)
      → hard-filter → interpret resistance markers → plot mdr1 coverage → MultiQC
```

It is built to run **in GitHub Codespaces with zero local installation**. Every
tool is pre-installed when the Codespace is created. **You type each command
yourself**, with its real parameters, so you understand every stage — the point is
to learn the steps, not to hide them behind a workflow manager.

The example sample is **ERR018972** (a *P. falciparum* run on the ENA), mapped to
the **3D7 reference genome** from PlasmoDB.

---

## 1. Start the Codespace

1. On the GitHub repo page: **Code ▸ Codespaces ▸ Create codespace on main**.
2. Wait a few minutes. On first creation it installs micromamba and builds the
   `pfdr` conda environment with every tool (FastQC, fastp, MultiQC, hostile,
   bwa, samtools, bcftools, mosdepth, GATK4, seqtk). You’ll see it working in the
   terminal; the message `Setup complete` means it’s done.
3. **Open a new terminal** (``Terminal ▸ New Terminal``). You should see the
   prompt prefixed with `(pfdr)` — the environment is active automatically.

> Check you’re ready: `bwa 2>&1 | head -3` and `gatk --version` should both work.

No Codespaces? See [Running locally](#running-locally) at the bottom.

---

## 2. Run the pipeline, one command at a time

Work through the steps **in order**, reading the output of each before moving on.
First set two variables you’ll reuse in every command (threads: a Codespace has
~4 cores):

```bash
S=ERR018972     # sample name
T=4             # threads
```

> **Steps 00 and 01 stay as scripts** — they only download the reference + reads
> and build indexes (boilerplate, nothing to learn by typing out). **Every other
> step you run as the real commands below.** If you ever want to check your typing,
> the identical commands also live in `steps/NN_*.sh`, and `bash run_all.sh` runs
> the whole thing end to end.

---

### Step 00 — Reference genome *(script)*

```bash
bash steps/00_reference.sh
```

Downloads the 3D7 genome from PlasmoDB and builds the companion files every later
tool needs: a **`.fai`** (`samtools faidx`), a **`.dict`** (`gatk
CreateSequenceDictionary`), the **BWA index** (`bwa index`), and — by parsing the
PlasmoDB **GFF** — a BED of the drug-resistance genes (`ref/dr_genes.bed`) plus a
CDS model (`ref/dr_genes_cds.json`) used later to translate variants into codons.

*P. falciparum* is extremely **AT-rich (~81% AT)** — remember that when you see the
GC plots; it’s normal, not a problem with your data.

### Step 01 — Download the reads *(script)*

```bash
bash steps/01_get_reads.sh
```

Pulls the paired FASTQs from the ENA into `data/` (`${S}_1.fastq.gz`,
`${S}_2.fastq.gz`).

**To keep a classroom session fast**, subsample first — e.g. 2 million read pairs:

```bash
SUBSAMPLE_PAIRS=2000000 bash steps/01_get_reads.sh
```

The later commands below use the raw files. **If you subsampled**, swap
`data/${S}_1.fastq.gz` → `data/${S}.sub_1.fastq.gz` (and `_2`) in step 02/03.
(Subsampling lowers coverage, which weakens copy-number and mixed-infection calls.)

---

### Step 02 — FastQC on the raw reads

Look at the data *before* changing it: per-base quality, adapter content,
duplication, and the GC distribution.

```bash
mkdir -p results/01_fastqc_raw

fastqc -t $T --outdir results/01_fastqc_raw \
    data/${S}_1.fastq.gz data/${S}_2.fastq.gz
```

**Look at:** per-base sequence quality, adapter content, and the GC peak. A second
GC peak around ~41 % is a hint of **human contamination** (removed in step 04).
Open `results/01_fastqc_raw/*_fastqc.html` from the file explorer.

### Step 03 — Trim adapters & low-quality bases (fastp), then re-QC

`fastp` removes adapters and low-quality read ends and drops reads that become too
short; `--detect_adapter_for_pe` auto-detects adapters for paired-end data. Then
run FastQC again to compare **before vs after**.

```bash
mkdir -p results/02_trim

fastp \
    --in1 data/${S}_1.fastq.gz --in2 data/${S}_2.fastq.gz \
    --out1 results/02_trim/${S}.trim_1.fastq.gz \
    --out2 results/02_trim/${S}.trim_2.fastq.gz \
    --detect_adapter_for_pe \
    --cut_front --cut_tail --cut_mean_quality 20 \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --thread $T \
    --json results/02_trim/fastp.json \
    --html results/02_trim/fastp.html

fastqc -t $T --outdir results/02_trim \
    results/02_trim/${S}.trim_1.fastq.gz \
    results/02_trim/${S}.trim_2.fastq.gz
```

**Look at:** `results/02_trim/fastp.html` — the one-file summary of what was
trimmed and the before/after quality.

### Step 04 — Remove human (host) reads with hostile

Clinical malaria samples come from **patient blood**, so the data mixes parasite
and human reads. `hostile` aligns reads to a human reference (Bowtie2 index
`human-t2t-hla`, **downloaded & cached on first use, ~3.5 GB**) and keeps only the
reads that *don’t* map.

```bash
mkdir -p results/03_dehost

hostile clean \
    --fastq1 results/02_trim/${S}.trim_1.fastq.gz \
    --fastq2 results/02_trim/${S}.trim_2.fastq.gz \
    --aligner bowtie2 \
    --index human-t2t-hla \
    --threads $T \
    --out-dir results/03_dehost \
    --force > results/03_dehost/hostile_log.json

# hostile names outputs *.clean_1.fastq.gz / *.clean_2.fastq.gz —
# give them tidy names the next step expects:
cp results/03_dehost/*clean_1.fastq.gz results/03_dehost/${S}.clean_1.fastq.gz
cp results/03_dehost/*clean_2.fastq.gz results/03_dehost/${S}.clean_2.fastq.gz
```

**Look at:** `results/03_dehost/hostile_log.json` — `reads_removed_proportion` is
the fraction flagged as human.

### Step 05 — Map to the reference with `bwa mem -M`

Two flags matter:
- **`-M`** marks shorter split hits as *secondary* — the conventional, safe choice
  for Picard/GATK compatibility (the flag this workshop specifies).
- **`-R`** adds a **read group**. GATK **refuses to run** without `@RG`, because it
  groups reads by sample (`SM`) and library (`LB`).

Pipe straight into `samtools sort` (no giant intermediate SAM), then index.

```bash
mkdir -p results/04_map

bwa mem -M -t $T \
    -R "@RG\tID:${S}\tSM:${S}\tLB:${S}_lib1\tPL:ILLUMINA" \
    ref/Pf3D7.fasta \
    results/03_dehost/${S}.clean_1.fastq.gz \
    results/03_dehost/${S}.clean_2.fastq.gz \
  | samtools sort -@ $T -o results/04_map/${S}.sorted.bam -

samtools index results/04_map/${S}.sorted.bam
samtools flagstat results/04_map/${S}.sorted.bam
```

**Look at:** the `flagstat` output — what fraction of reads mapped, and how many
are properly paired?

### Step 06 — Mark duplicates

PCR/optical duplicates are the *same* fragment sequenced many times; they must not
each count as independent evidence. `MarkDuplicates` **flags** them (doesn’t delete
them) so the caller ignores them.

```bash
mkdir -p results/05_markdup

gatk MarkDuplicates \
    -I results/04_map/${S}.sorted.bam \
    -O results/05_markdup/${S}.markdup.bam \
    -M results/05_markdup/${S}.dup_metrics.txt \
    --CREATE_INDEX true
```

**Look at:** `PERCENT_DUPLICATION` in `${S}.dup_metrics.txt`.

> **BQSR?** A human-style pipeline would run `BaseRecalibrator` + `ApplyBQSR` here.
> It needs a trusted known-variants set, which isn’t standard for *P. falciparum*,
> so this workshop skips it — see [Scaling to a cohort](#scaling-to-a-cohort).

### Step 07 — Coverage with mosdepth

Two runs: **genome-wide in 500 bp windows** (the track plotted across mdr1 in step
11), and **restricted to the DR genes** (a per-gene mean-depth table). The mosdepth
output prefix can include a path, so no `cd` is needed.

```bash
mkdir -p results/06_coverage

# (a) genome-wide, 500 bp windows
mosdepth --no-per-base --fast-mode -t $T --by 500 \
    results/06_coverage/${S}.genome \
    results/05_markdup/${S}.markdup.bam

# (b) over the drug-resistance genes
mosdepth --no-per-base -t $T --by ref/dr_genes.bed \
    results/06_coverage/${S}.dr \
    results/05_markdup/${S}.markdup.bam

# peek at per-gene mean depth
zcat results/06_coverage/${S}.dr.regions.bed.gz
```

**Look at:** the per-gene depths — mdr1 noticeably higher than the others is the
amplification signal you’ll plot in step 11.

### Step 08 — Call variants (GATK, the GVCF way)

Two moves, exactly as in the GATK germline short-variant best practice:
1. **HaplotypeCaller** in **GVCF mode** (`-ERC GVCF`) → a record for *every*
   position (variant or reference block).
2. **GenotypeGVCFs** → the final genotyped VCF.

Why two moves for one sample? Because **this is what scales to a cohort** (see
below). `-L ref/dr_genes.bed` restricts calling to the DR genes to keep it fast.
**Ploidy 2:** blood-stage *P. falciparum* is haploid, but patient samples are often
**mixed infections** — calling with `--sample-ploidy 2` lets a second clone appear
as a “heterozygous” genotype (the MalariaGEN convention).

```bash
mkdir -p results/07_variants

gatk HaplotypeCaller \
    -R ref/Pf3D7.fasta \
    -I results/05_markdup/${S}.markdup.bam \
    -O results/07_variants/${S}.g.vcf.gz \
    -ERC GVCF \
    --sample-ploidy 2 \
    -L ref/dr_genes.bed --interval-padding 1000

gatk GenotypeGVCFs \
    -R ref/Pf3D7.fasta \
    -V results/07_variants/${S}.g.vcf.gz \
    -O results/07_variants/${S}.vcf.gz

# how many raw variants?
zcat results/07_variants/${S}.vcf.gz | grep -vc '^#'
```

**Whole genome instead?** Drop the `-L … --interval-padding …` line (slower).

### Step 09 — Hard-filter the variants

“Hard filtering” flags variants on fixed thresholds — the right approach when you
**don’t have enough data to train** VQSR (the single-sample case). SNPs and indels
need **different** thresholds, so: split → filter each → merge. Filtered records are
**labelled** in the FILTER column, not deleted.

```bash
mkdir -p results/08_filter

# 1) split into SNPs and indels
gatk SelectVariants -R ref/Pf3D7.fasta -V results/07_variants/${S}.vcf.gz \
    --select-type-to-include SNP   -O results/08_filter/${S}.snps.raw.vcf.gz
gatk SelectVariants -R ref/Pf3D7.fasta -V results/07_variants/${S}.vcf.gz \
    --select-type-to-include INDEL -O results/08_filter/${S}.indels.raw.vcf.gz

# 2a) SNP hard filters (GATK-recommended starting thresholds)
gatk VariantFiltration -R ref/Pf3D7.fasta \
    -V results/08_filter/${S}.snps.raw.vcf.gz \
    -O results/08_filter/${S}.snps.filt.vcf.gz \
    --filter-expression "QD < 2.0"              --filter-name "QD2" \
    --filter-expression "FS > 60.0"             --filter-name "FS60" \
    --filter-expression "MQ < 40.0"             --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum" \
    --filter-expression "SOR > 3.0"             --filter-name "SOR3"

# 2b) indel hard filters (looser strand/position thresholds)
gatk VariantFiltration -R ref/Pf3D7.fasta \
    -V results/08_filter/${S}.indels.raw.vcf.gz \
    -O results/08_filter/${S}.indels.filt.vcf.gz \
    --filter-expression "QD < 2.0"               --filter-name "QD2" \
    --filter-expression "FS > 200.0"             --filter-name "FS200" \
    --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum" \
    --filter-expression "SOR > 10.0"             --filter-name "SOR10"

# 3) merge back together
gatk MergeVcfs \
    -I results/08_filter/${S}.snps.filt.vcf.gz \
    -I results/08_filter/${S}.indels.filt.vcf.gz \
    -O results/08_filter/${S}.filtered.vcf.gz

# PASS vs flagged breakdown
zcat results/08_filter/${S}.filtered.vcf.gz | grep -v '^#' | cut -f7 | sort | uniq -c
```

**Look at:** the FILTER breakdown. See the [Annotation glossary](#annotation-glossary)
for what QD / FS / MQ / SOR / ReadPosRankSum mean and why indels get looser cutoffs.

### Step 10 — Interpret the resistance markers

This step uses the repo’s helper (`scripts/report_dr.py`): for each known marker it
finds the codon’s genomic position from the CDS model (**handling + and − strand
genes**), **self-checks** that the reference codon really translates to the
documented wild-type amino acid (catches a reference/GFF release mismatch), applies
the PASS variants, translates, and classifies the codon as **wild-type / resistant
marker / other mutation**, with allele fractions so **mixed infections** show up.

```bash
mkdir -p results/09_report

python scripts/report_dr.py \
    --fasta ref/Pf3D7.fasta \
    --cds ref/dr_genes_cds.json \
    --markers resources/dr_markers.tsv \
    --vcf results/08_filter/${S}.filtered.vcf.gz \
    --coverage results/06_coverage/${S}.dr.regions.bed.gz \
    --sample $S \
    --outdir results/09_report \
    --min-dp 5 \
    --mixed-min-af 0.10
```

**Outputs:** `results/09_report/${S}.dr_variants.tsv` (every marker),
`${S}.dr_calls.tsv` (non-wild-type only), and a JSON. The marker list is the
editable TSV **`resources/dr_markers.tsv`** (covers `crt`, `mdr1`, `dhfr`, `dhps`,
`kelch13`, `cytb`). It is a **teaching set, not a clinical-grade panel**.

### Step 11 — Plot mdr1 coverage & estimate copy number

`mdr1` resistance has **two** mechanisms: point mutations (step 10) **and gene
amplification**. Copy number ≈ *median depth over mdr1 / genome-wide median depth*.

```bash
mkdir -p results/10_plots

python scripts/plot_coverage.py \
    --windows results/06_coverage/${S}.genome.regions.bed.gz \
    --summary results/06_coverage/${S}.genome.mosdepth.summary.txt \
    --dr-regions results/06_coverage/${S}.dr.regions.bed.gz \
    --cds ref/dr_genes_cds.json \
    --sample $S \
    --gene mdr1 --gene-id PF3D7_0523000 \
    --out results/10_plots/${S}.mdr1_coverage.png
```

**Look at:** `results/10_plots/${S}.mdr1_coverage.png`. A value clearly above 1 =
amplification, linked to **reduced lumefantrine/mefloquine susceptibility**. The
same logic applies to **plasmepsin 2/3** amplification and **piperaquine**.

### Step 12 — Aggregate everything with MultiQC

Stitches FastQC, fastp, samtools, mosdepth and GATK outputs into one interactive
HTML report.

```bash
mkdir -p results/11_multiqc

multiqc -f -o results/11_multiqc -n ${S}_multiqc.html results
```

**Look at:** `results/11_multiqc/${S}_multiqc.html` — the whole run on one page.

---

All settings (sample, URLs, threads, ploidy, regions to call, filtering depth) also
live in **`config.sh`** with comments, which is what the `steps/*.sh` scripts and
`run_all.sh` read.

---

## 4. The drug-resistance biology in one table

| Gene | PlasmoDB ID | Mechanism | Drugs affected |
|------|-------------|-----------|----------------|
| **crt** | PF3D7_0709000 | SNPs (K76T …; plus PPQ-associated SNPs) | chloroquine, amodiaquine, piperaquine |
| **mdr1** | PF3D7_0523000 | SNPs **and copy number** | lumefantrine, mefloquine, amodiaquine, chloroquine |
| **dhfr** | PF3D7_0417200 | SNPs (N51I, C59R, S108N, I164L) | pyrimethamine (SP) |
| **dhps** | PF3D7_0810800 | SNPs (A437G, K540E, A581G) | sulfadoxine (SP) |
| **kelch13** | PF3D7_1343700 | propeller-domain SNPs (e.g. C580Y, R561H) | artemisinin (partial resistance) |
| **cytb** | PF3D7_MIT02300 | SNP at codon 268 | atovaquone |
| **plasmepsin 2/3** | PF3D7_1408000/…8100 | copy number | piperaquine |

Full marker-by-marker detail is in `resources/dr_markers.tsv`.

---

## Scaling to a cohort

This workshop does one sample. For a real study you’d make **one GVCF per sample**
(step 08 as-is), then genotype them **together** so every site is evaluated in
every sample:

```bash
# 1) per sample (parallel): steps 00-08 -> sample1.g.vcf.gz, sample2.g.vcf.gz, ...

# 2) consolidate the GVCFs
gatk GenomicsDBImport \
     -V sample1.g.vcf.gz -V sample2.g.vcf.gz ... \
     --genomicsdb-workspace-path cohort_db \
     -L ref/dr_genes.bed

# 3) joint genotyping
gatk GenotypeGVCFs -R ref/Pf3D7.fasta -V gendb://cohort_db -O cohort.vcf.gz
```

### BQSR and VQSR — why this workshop skips them

- **BQSR** (Base Quality Score Recalibration, between steps 06 and 08) corrects
  systematic errors in the sequencer’s quality scores. It needs a trusted
  **known-variants** VCF to mask real variation. For *P. falciparum* there isn’t a
  single standard known-sites set, and for a teaching run the benefit is small, so
  we skip it. For a production cohort you can run
  `gatk BaseRecalibrator`/`ApplyBQSR` using a high-confidence site list (e.g. one
  derived from MalariaGEN’s Pf releases).
- **VQSR** (Variant Quality Score Recalibration) replaces the step-09 hard filters
  with a model **trained** on truth/training resource sets. It needs many samples
  and good truth sets to be stable — so it’s a **cohort** tool, not a single-sample
  one. With enough samples and a training set, swap step 09 for
  `VariantRecalibrator` + `ApplyVQSR`. Until then, hard-filtering is the correct
  choice.

---

## Annotation glossary

Thresholds used in step 09 (GATK-recommended starting points):

| Annotation | Meaning | SNP filter | Indel filter |
|-----------|---------|-----------|--------------|
| **QD** | Quality normalised by depth | `< 2.0` | `< 2.0` |
| **FS** | Fisher strand bias (Phred) | `> 60` | `> 200` |
| **SOR** | Strand odds ratio (another strand-bias measure) | `> 3.0` | `> 10.0` |
| **MQ** | Root-mean-square mapping quality | `< 40` | — |
| **MQRankSum** | Mapping-quality difference ref vs alt reads | `< −12.5` | — |
| **ReadPosRankSum** | Are alt alleles biased to read ends? | `< −8.0` | `< −20.0` |

Indels get looser strand/position thresholds because they are intrinsically harder
to align, which inflates strand and position metrics even for real indels.

---

## Running locally

You don’t need Codespaces. With [micromamba](https://mamba.readthedocs.io) (or
conda/mamba):

```bash
micromamba create -y -f environment.yml
micromamba activate pfdr
# then the same commands as above, starting with:
bash steps/00_reference.sh
```

Needs internet (PlasmoDB + ENA + the hostile human index), ~4 CPUs and ~16 GB RAM
for comfort, and a few GB of disk.

---

## Repository layout

```
.devcontainer/      Codespaces setup (installs the pfdr environment)
environment.yml     all tools + versions
config.sh           every setting, with comments
steps/              00–12, one script per stage (same commands as the README)
scripts/            Python helpers (region builder, DR report, coverage plot)
resources/          dr_genes.tsv, dr_markers.tsv  (editable reference tables)
tests/              synthetic end-to-end self-test with known answers
run_all.sh          run every step in order
```

### Instructor self-test
`bash tests/selftest.sh` builds a tiny synthetic genome with **planted** mutations
(crt K76T on the + strand, dhfr S108N on the − strand) and a duplicated mdr1
locus, runs the real step scripts, and asserts the report and copy-number estimate
are correct. Runs in well under a minute and needs no downloads — a good smoke test
before a workshop.

---

## Going further

- **MalariaGEN Pf7** — the reference resource for *P. falciparum* genomic
  epidemiology and curated resistance classifications.
- **Dedicated DR callers** for curated, clinical-grade output (GRC/amplicon
  workflows). Use those for reporting; use this workshop to understand what they do
  under the hood.

> This pipeline is for **training and research**, not clinical diagnosis.
