## 1. Overview

This repository contains a complete RNA-seq preprocessing and quantification pipeline, implemented using Snakemake for workflow automation and reproducibility.

It runs the full workflow:

- FastQC on raw FASTQ
- Adapter trimming (Cutadapt – paired-end)
- FastQC post-trimming
- STAR genome index generation
- STAR paired-end alignment
- BAM processing (fixmate, sort, dedup) — samtools + Picard
- Filtering + featureCounts
- DESeq2 normalization 

**All steps are organized into Snakemake rules with explicit dependencies, automatic folder creation, file tracking, and full reproducibility.**

You can download the full pipeline by cloning this repository:

```code
git clone https://github.com/robirritteri/uca_m2bi_snakemake.git
```

## 2. Pipeline Structure
```text
m2bi_snakemake/
│
├── workflow/
│   ├── snakefile.smk
│   
├── config/
│   └── config.yaml
│
├── scripts/
│   ├── all_counts.R
│   └── Analyse.R
│
├── results/
│   ├── FastQC/
│   ├── Trimmed/
│   ├── FastQC_trim/
│   ├── STAR/
│   │   ├── index/
│   │   ├── alignments/
│   │   └── clean/
│   └── counts/
│
└── rnaseq.yaml
```
## 3. Conda Environment

Create the required environment:
```code
conda env create -f rnaseq.yaml
conda activate rnaseq
```

This environment includes:

**Cutadapt, samtools, STAR, picard, subread (featureCounts), R + DESeq2, tidyverse.**

It is fully compatible with all scripts.

## 4. Configuration File

All paths, references, and sample names are defined in:

```code
config/config.yaml
```

Example:
```code
paths:
  raw_dir: "/home/users/shared/data/stategra/RNAseq/raw"

  fastqc_raw_dir: "results/FastQC"
  trim_dir: "results/Trimmed"
  fastqc_trim_dir: "results/FastQC_trim"

  star_index_dir: "results/STAR/index"
  star_align_dir: "results/STAR/alignments"
  clean_bam_dir: "results/STAR/clean"
  counts_dir: "results/STAR/counts"

reference:
  genome_fasta: "/home/users/shared/databanks/Mus_musculus_GRCm39/fasta/all.fasta"
  gtf: "/home/users/shared/databanks/Mus_musculus_GRCm39/flat/genomic.gtf"

samples:
  - control_0h_B1_R1
  - control_0h_B1_R2
```

This file can be edited to customize input folders, genome references, or sample lists.

## 5. Running the Pipeline

To preview all planned steps:

```code
snakemake --snakefile workflow/snakefile.smk --dry-run
```

To execute the workflow:

```code
snakemake --snakefile workflow/snakefile.smk --cores 8 --printshellcmds
```

Snakemake will:
- create missing directories
- run tasks in the correct order
- skip outputs already up-to-date
- automatically call R scripts after quantification

## 6. Pipeline Steps
### 01 — FastQC on raw reads

Runs FastQC on all raw FASTQ files using an sbatch array.
Also generates a MultiQC report.

### 02 — Adapter trimming (paired-end, Cutadapt)

Trims adapters and low-quality bases with Cutadapt using the adapters you provided.

Output:
results/RNAseq/trimmed/paired/*.trim.fastq.gz

### 03 — FastQC after trimming

FastQC + MultiQC on trimmed FASTQ.

### 04 — STAR genome index generation

Uses the Mus musculus GRCm39 FASTA + GTF:
```code
/home/users/shared/databanks/Mus_musculus_GRCm39/fasta/all.fasta
/home/users/shared/databanks/Mus_musculus_GRCm39/flat/genomic.gtf
```
Index is created in:
```code
$HOME
```

### 05 — STAR alignment (paired-end)

Aligns trimmed paired-end reads.

Produces:
```text
$HOME/results/RNAseq/alignment/paired/<sample>/
   ├── <sample>_Aligned.sortedByCoord.out.bam
   ├── <sample>_Log.final.out
   └── <sample>_ReadsPerGene.out.tab
```

### 06 — Picard + Samtools cleaning

Operations for paired-end BAMs:
```text
- fixmate
- sort
- filter quality
- remove duplicates
- write metrics
```
Output:
```code
results/RNAseq/picard/bamTraite/paired/dedup_<sample>.bam
```

### 07 — featureCounts quantification

Counts reads with:
```
-s 2 (reverse stranded)
-p (paired-end)
-a (GTF="/home/users/shared/databanks/Mus_musculus_GRCm39/flat/genomic.gtf")
```

Output:

results/RNAseq/samtools/<sample>/<sample>.featureCounts.txt

### 08 — R analysis

- imports merged featureCounts table
- computes normalization factors
- outputs VST and normalized counts

Outputs:
```code
$HOME/results/RNAseq/counts/deseq2_normalized_counts.csv
$HOME/results/RNAseq/counts/deseq2_vst.csv
$HOME/results/RNAseq/counts/dds.rds
```

## 7. Pre-run checklist

- Check that your environment is loaded
```code
module load conda
conda activate rnaseq
```

- Check that you are inside your /workflow/ directory

- Check that input FASTQ exist
```code
ls /home/users/shared/data/stategra/RNAseq/raw
```

- Preview a single rule:
```code
snakemake --snakefile workflow/snakefile.smk fastqc_raw -n
```

- Run a specific output:
```code
snakemake results/Trimmed/control_0h_B1_R1_1.trim.fastq.gz
```

## 8. License

MIT License — feel free to reuse and adapt.
