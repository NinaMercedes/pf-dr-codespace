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
tool is pre-installed when the Codespace is created. You run the steps yourself,
read what each one prints, and look at the outputs — the point is to understand
each stage, not to hide it behind a workflow manager.

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

## 2. Run the pipeline, one step at a time

Run these **in order**, reading the output of each before moving on. Every script
prints the exact command it runs (after `$`) so you can see — and copy — the real
tool invocation.

| Step | Command | What it does |
|------|---------|--------------|
| 00 | `bash steps/00_reference.sh` | Download 3D7 reference; build `.fai`, `.dict`, BWA index, and the DR-gene region files |
| 01 | `bash steps/01_get_reads.sh` | Download the paired FASTQ reads from the ENA |
| 02 | `bash steps/02_fastqc_raw.sh` | **FastQC** on the raw reads |
| 03 | `bash steps/03_trim.sh` | **fastp** adapter/quality trimming, then FastQC again |
| 04 | `bash steps/04_remove_host.sh` | **hostile** removes human reads |
| 05 | `bash steps/05_map.sh` | **bwa mem -M** → sorted BAM (+ flagstat) |
| 06 | `bash steps/06_markdup.sh` | **GATK MarkDuplicates** |
| 07 | `bash steps/07_coverage.sh` | **mosdepth** genome-wide + per-DR-gene coverage |
| 08 | `bash steps/08_call_variants.sh` | **GATK HaplotypeCaller** (GVCF) → **GenotypeGVCFs** |
| 09 | `bash steps/09_filter.sh` | **GATK** hard-filtering (SNPs + indels) |
| 10 | `bash steps/10_report.sh` | Translate variants → **resistance marker report** |
| 11 | `bash steps/11_plot_mdr1.sh` | Coverage plot across **mdr1** + copy-number estimate |
| 12 | `bash steps/12_multiqc.sh` | **MultiQC** — one HTML summary of everything |

Outputs land in `results/`, numbered to match the steps. Open the `.html` reports
and the `.png` plot straight from the VS Code file explorer (right-click an HTML
file ▸ *Open with Live Preview*; click a PNG to view it).

**In a hurry / want a full pass?** `bash run_all.sh` runs everything in order —
but do the steps by hand the first time.

### Make it faster

The first run downloads the full read set and (in step 04) a ~3.5 GB human index.
To keep a classroom session snappy, subsample the reads before step 02:

```bash
# edit config.sh and set SUBSAMPLE_PAIRS=2000000, OR for one run:
SUBSAMPLE_PAIRS=2000000 bash steps/01_get_reads.sh
```

Everything downstream picks up the subsampled files automatically. Note that
subsampling lowers coverage, which weakens copy-number and mixed-infection calls.

All knobs live in **`config.sh`** (sample, URLs, threads, ploidy, which regions to
call, filtering depth, etc.) with comments explaining each.

---

## 3. What each step is really doing (and what to look for)

### 00 — Reference genome
Mapping and GATK need three companion files beside the FASTA: a **`.fai`** index
(`samtools faidx`) for random access, a **`.dict`** (`gatk CreateSequenceDictionary`)
listing contigs, and the **BWA index** (`bwa index`). We also parse the PlasmoDB
**GFF** to get the coding coordinates of the drug-resistance genes and write a BED
(used to focus coverage and calling) and a CDS model (used to translate variants
into codon changes).

*P. falciparum* is extremely **AT-rich (~81% AT)** — keep that in mind when you see
the GC plots later; it’s normal, not a problem with your data.

### 02 / 03 — QC and trimming
**FastQC** first, so you look at the data before changing it: per-base quality,
adapter content, duplication, and the GC distribution. **fastp** then removes
adapters and low-quality read ends and drops reads that become too short; we run
FastQC again so you can compare *before vs after*. A second GC peak around ~41% in
the raw data is a hint of **human contamination**, which the next step removes.

### 04 — Host (human) read removal
Clinical malaria samples are sequenced from **patient blood**, so the data is a
mix of parasite and human reads. **hostile** aligns reads to a human reference and
keeps only the reads that *don’t* map. Removing human reads speeds up mapping,
avoids spurious alignments, and is good practice for data sharing (no human
genetic data left behind). The JSON log reports how many reads were removed.

### 05 — Mapping with `bwa mem -M`
`bwa mem` aligns the cleaned reads to 3D7. The flags that matter:
- **`-M`** marks shorter split hits as *secondary* — the conventional, safe choice
  for compatibility with Picard/GATK (this is the flag the workshop specifies).
- **`-R`** adds a **read group** (`@RG … SM:<sample>`). GATK **refuses to run**
  without it, because it groups reads by sample and library.

We pipe straight into `samtools sort` and index. Read the `flagstat` output: what
fraction of reads mapped, and how many are properly paired?

### 06 — Mark duplicates
PCR/optical duplicates are the *same* original fragment sequenced multiple times;
they must not each count as independent evidence. `MarkDuplicates` **flags** them
(doesn’t delete them) so the caller ignores them. Check the duplicate rate in the
metrics file.

### 07 — Coverage with mosdepth
Two runs: **genome-wide in 500 bp windows** (the track we plot across mdr1), and
**restricted to the DR genes** (a quick per-gene mean-depth table). Depth is the
basis of the copy-number estimate in step 11.

### 08 — Variant calling (GATK, the GVCF way)
Two moves, exactly as in the GATK germline short-variant best practice:
1. **HaplotypeCaller** in **GVCF mode** (`-ERC GVCF`) → a record for *every*
   position (variant or reference block).
2. **GenotypeGVCFs** → the final genotyped VCF.

Why two moves for one sample? Because **this is what scales to a cohort**: one
GVCF per sample, combined and genotyped together (see
[Scaling to a cohort](#scaling-to-a-cohort)). By default we call only inside the
DR genes (fast); set `CALL_REGIONS=genome` in `config.sh` for the whole genome.

**Ploidy.** Blood-stage *P. falciparum* is **haploid**, but patient samples are
often **mixed infections** (more than one clone). We call with `--sample-ploidy 2`
so a second clone shows up as a “heterozygous” genotype — the convention
MalariaGEN uses. The report flags such sites as *mixed*.

### 09 — Hard-filtering
“Hard filtering” flags variants using fixed thresholds on VCF annotations. It’s
the right approach when you **don’t have enough data to train** the machine-learning
filters (VQSR / newer tools) — which is exactly the single-sample case. We split
into **SNPs** and **indels** (they get different thresholds), apply the
GATK-recommended starting thresholds, and merge back. Filtered records are
**labelled** in the FILTER column, not deleted, so you can see *why* each was
dropped. The DR report keeps only `PASS`.

> See [Annotation glossary](#annotation-glossary) for what QD / FS / MQ / SOR /
> ReadPosRankSum mean and why the thresholds differ for indels.

### 10 — Resistance interpretation
`scripts/report_dr.py` takes each known marker (gene, codon, wild-type and
resistant amino acids), finds the codon’s genomic position from the CDS model
(**handling + and − strand genes**), **self-checks** that the reference codon
really translates to the documented wild-type amino acid (guards against a
reference/GFF release mismatch), applies the PASS variants, translates, and
classifies the codon as **wild-type / resistant marker / other mutation**, with
allele depths so **mixed infections are visible**. Outputs:
`results/09_report/<sample>.dr_variants.tsv` (every marker),
`…dr_calls.tsv` (non-wild-type only), and a JSON.

The marker list is in **`resources/dr_markers.tsv`** — a readable TSV you can
extend. It covers `crt`, `mdr1`, `dhfr`, `dhps`, `kelch13` and `cytb`. It is a
teaching set, **not** a clinical-grade panel; for curated calls use a dedicated
tool (see [Going further](#going-further)).

### 11 — mdr1 coverage & copy number
`mdr1` resistance has **two** mechanisms: point mutations (handled in step 10)
**and gene amplification** (extra copies). We estimate:

```
copy number ≈ median depth over mdr1 / genome-wide median depth
```

The plot shows the depth profile across mdr1 (with the 1× and 2× baselines) and a
per-DR-gene depth panel. A value clearly above 1 indicates amplification, linked
to **reduced lumefantrine and mefloquine susceptibility**. The same logic applies
to **plasmepsin 2/3** amplification and **piperaquine** resistance.

### 12 — MultiQC
Stitches FastQC, fastp, samtools, mosdepth and GATK outputs into one interactive
HTML report — handy for a final eyeball and, with many samples, for comparing them.

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
     -L resources/… (intervals)

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
bash steps/00_reference.sh
# ... and so on
```

Needs internet (PlasmoDB + ENA + the hostile human index), ~4 CPUs and ~16 GB RAM
for comfort, and a few GB of disk.

---

## Repository layout

```
.devcontainer/      Codespaces setup (installs the pfdr environment)
environment.yml     all tools + versions
config.sh           every setting, with comments
steps/              00–12, one script per stage (read these!)
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
- **Dedicated DR callers** for curated, clinical-grade output: *Pf-HaplotypeCaller
  / GRC* workflows, **miniature** amplicon callers, and tools like
  **ampseq**/**[resmarker]** pipelines. Use those for reporting; use this workshop
  to understand what they do under the hood.

> This pipeline is for **training and research**, not clinical diagnosis.
