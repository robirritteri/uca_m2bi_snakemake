###############################################
# RNA-seq pipeline for all samples
###############################################

# Load configuration
configfile: "/home/users/student14/m2bi_snakemake/config/config.yaml"

# ----- PATHS -----
PATHS           = config["paths"]
RAW_DIR         = PATHS["raw_dir"]
FASTQC_RAW_DIR  = PATHS["fastqc_raw_dir"]
TRIM_DIR        = PATHS["trim_dir"]
FASTQC_TRIM_DIR = PATHS["fastqc_trim_dir"]
STAR_INDEX_DIR  = PATHS["star_index_dir"]
STAR_ALIGN_DIR  = PATHS["star_align_dir"]
CLEAN_BAM_DIR   = PATHS["clean_bam_dir"]
COUNTS_DIR      = PATHS["counts_dir"]

# ----- REFERENCE -----
GENOME_FASTA = config["reference"]["genome_fasta"]
GTF_ANNOT    = config["reference"]["gtf"]

# ----- SAMPLES -----
SAMPLES = config["samples"]


###############################################
# RULE
###############################################

rule all:
    input:
        # Raw FastQC
        expand(f"{FASTQC_RAW_DIR}/{{sample}}_1_fastqc.html", sample=SAMPLES),
        expand(f"{FASTQC_RAW_DIR}/{{sample}}_2_fastqc.html", sample=SAMPLES),

        # Trimmed FASTQ
        expand(f"{TRIM_DIR}/{{sample}}_1.trim.fastq.gz", sample=SAMPLES),
        expand(f"{TRIM_DIR}/{{sample}}_2.trim.fastq.gz", sample=SAMPLES),

        # Trimmed FastQC
        expand(f"{FASTQC_TRIM_DIR}/{{sample}}_1.trim_fastqc.html", sample=SAMPLES),
        expand(f"{FASTQC_TRIM_DIR}/{{sample}}_2.trim_fastqc.html", sample=SAMPLES),

        # STAR alignment outputs
        expand(f"{STAR_ALIGN_DIR}/{{sample}}_Aligned.sortedByCoord.out.bam", sample=SAMPLES),

        # Clean BAM
        expand(f"{CLEAN_BAM_DIR}/{{sample}}_ready.bam", sample=SAMPLES),

        # FeatureCounts
        expand(f"{COUNTS_DIR}/{{sample}}_counts.txt", sample=SAMPLES)

		# Final merged + R analysis result
        "results/R/analysis/DE_results_mydata.tsv"


###############################################
# 1) FastQC raw
###############################################

rule fastqc_raw:
    input:
        lambda wc: f"{RAW_DIR}/{wc.sample}_{wc.pair}.fastq.gz"
    output:
        f"{FASTQC_RAW_DIR}/{{sample}}_{{pair}}_fastqc.zip",
        f"{FASTQC_RAW_DIR}/{{sample}}_{{pair}}_fastqc.html"
    shell:
        """
        module load conda
        conda activate rnaseq

        mkdir -p {FASTQC_RAW_DIR}
        fastqc --outdir {FASTQC_RAW_DIR} {input.r}
        """


rule fastqc_raw_r2:
    input:
        r = lambda wc: f"{RAW_DIR}/{wc.sample}_2.fastq.gz"
    output:
        zip  = f"{FASTQC_RAW_DIR}/{{sample}}_2_fastqc.zip",
        html = f"{FASTQC_RAW_DIR}/{{sample}}_2_fastqc.html"
    shell:
        """
        module load conda
        conda activate rnaseq

        mkdir -p {FASTQC_RAW_DIR}
        fastqc --outdir {FASTQC_RAW_DIR} {input.r}
        """


###############################################
# 2) Cutadapt trimming
###############################################

rule cutadapt:
    input:
        r1 = lambda wc: f"{RAW_DIR}/{wc.sample}_1.fastq.gz",
        r2 = lambda wc: f"{RAW_DIR}/{wc.sample}_2.fastq.gz"
    output:
        r1 = f"{TRIM_DIR}/{{sample}}_1.trim.fastq.gz",
        r2 = f"{TRIM_DIR}/{{sample}}_2.trim.fastq.gz"
    shell:
        """
        module load conda
        conda activate rnaseq

        mkdir -p {TRIM_DIR}

        cutadapt \
            -q 20,20 \
            -m 30 \
            -a AGATGTGTATAAGAGACAG \
            -A AGATGTGTATAAGAGACAG \
            -o {output.r1} \
            -p {output.r2} \
            {input.r1} {input.r2}
        """


###############################################
# 3) FastQC trimmed
###############################################

rule fastqc_trim:
    input:
        lambda wc: f"{TRIM_DIR}/{wc.sample}_{wc.pair}.trim.fastq.gz"
    output:
        f"{FASTQC_TRIM_DIR}/{{sample}}_{{pair}}.trim_fastqc.zip",
        f"{FASTQC_TRIM_DIR}/{{sample}}_{{pair}}.trim_fastqc.html"
    shell:
        """
        module load conda
        conda activate rnaseq

        mkdir -p {FASTQC_TRIM_DIR}
        fastqc --outdir {FASTQC_TRIM_DIR} {input}
        """


###############################################
# 4) STAR INDEX
###############################################

rule star_index:
    input:
        fasta = GENOME_FASTA,
        gtf   = GTF_ANNOT
    output:
        f"{STAR_INDEX_DIR}/genomeParameters.txt"
    threads: 8
    shell:
        """
        module load conda
        conda activate rnaseq

        mkdir -p {STAR_INDEX_DIR}

        STAR \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {STAR_INDEX_DIR} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang 99
        """


###############################################
# 5) STAR ALIGNMENT
###############################################

rule star_align:
    input:
        r1    = f"{TRIM_DIR}/{{sample}}_1.trim.fastq.gz",
        r2    = f"{TRIM_DIR}/{{sample}}_2.trim.fastq.gz",
        index = f"{STAR_INDEX_DIR}/genomeParameters.txt"
    output:
        bam    = f"{STAR_ALIGN_DIR}/{{sample}}_Aligned.sortedByCoord.out.bam",
        counts = f"{STAR_ALIGN_DIR}/{{sample}}_ReadsPerGene.out.tab"
    threads: 8
    shell:
        """
        module load conda
        conda activate rnaseq

        mkdir -p {STAR_ALIGN_DIR}

        STAR \
            --runThreadN {threads} \
            --runMode alignReads \
            --genomeDir {STAR_INDEX_DIR} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {STAR_ALIGN_DIR}/{{wildcards.sample}}_ \
            --sjdbOverhang 99 \
            --quantMode TranscriptomeSAM GeneCounts
        """


###############################################
# 6) BAM CLEANING
###############################################

rule clean_bam:
    input:
        bam = f"{STAR_ALIGN_DIR}/{{sample}}_Aligned.sortedByCoord.out.bam"
    output:
        ready   = f"{CLEAN_BAM_DIR}/{{sample}}_ready.bam",
        metrics = f"{CLEAN_BAM_DIR}/{{sample}}_dedup.metrics.txt"
    threads: 4
    shell:
        """
        module load conda
        conda activate rnaseq

        mkdir -p {CLEAN_BAM_DIR}

        fixmate="{CLEAN_BAM_DIR}/{{wildcards.sample}}_fixmate.bam"
        sorted="{CLEAN_BAM_DIR}/{{wildcards.sample}}_sorted.bam"

        samtools fixmate -m {input.bam} "$fixmate"
        samtools view -b -q 28 "$fixmate" | samtools sort -@ {threads} -o "$sorted"

        picard MarkDuplicates \
            I="$sorted" \
            O="{output.ready}" \
            M="{output.metrics}" \
            REMOVE_DUPLICATES=true \
            VALIDATION_STRINGENCY=LENIENT

        rm "$sorted" "$fixmate"
        """

###############################################
# 7) FEATURECOUNTS
###############################################

rule featurecounts:
    input:
        bam = f"{CLEAN_BAM_DIR}/{{sample}}_ready.bam"
    output:
        counts = f"{COUNTS_DIR}/{{sample}}_counts.txt"
    params:
        gtf = GTF_ANNOT
    threads: 4
    shell:
        """
        module load conda
        conda activate rnaseq

        mkdir -p {COUNTS_DIR}

        featureCounts \
            -T {threads} \
            -a {params.gtf} \
            -o {output.counts} \
            {input.bam}
        """

###############################################
# 8) Merge featureCounts outputs
###############################################

R_RESULTS_DIR = "results/R"

rule merge_counts:
    input:
        expand(f"{COUNTS_DIR}/{{sample}}_counts.txt", sample=SAMPLES)
    output:
        f"{R_RESULTS_DIR}/all_featureCounts_counts.txt"
    params:
        counts_dir = COUNTS_DIR
    shell:
        """
        module load conda
        conda activate rnaseq

        mkdir -p {R_RESULTS_DIR}

        Rscript /home/users/student14/m2bi_snakemake/scripts/all_counts.R \
            {params.counts_dir} \
            {output}
        """


###############################################
# 9) R Analysis
###############################################

rule r_analysis:
    input:
        "results/R/all_featureCounts_counts.txt"
    output:
        "results/R/analysis/DE_results_mydata.tsv"
    params:
        rdir = "results/R"
    shell:
        """
        module load conda
        conda activate rnaseq

        mkdir -p results/R/analysis

        Rscript config/Analyse.R {params.rdir}
        """
