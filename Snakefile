import os
import pandas as pd

configfile: "config/config.yaml"

SAMPLES_DF = pd.read_csv(config["samples"], sep="\t")
if not {"sample", "r1", "r2"}.issubset(SAMPLES_DF.columns):
    raise ValueError("samples TSV 必须包含列: sample, r1, r2")

SAMPLES = SAMPLES_DF["sample"].tolist()
SAMPLE_TO_R1 = dict(zip(SAMPLES_DF["sample"], SAMPLES_DF["r1"]))
SAMPLE_TO_R2 = dict(zip(SAMPLES_DF["sample"], SAMPLES_DF["r2"]))
CONTRASTS = config["contrasts"]
CONTRAST_FILES = [c.replace(":", "_") for c in CONTRASTS]

STAR_INDEX = config["star_index"]
SALMON_INDEX = config["salmon_index"]
GTF = config["gtf"]
FASTA = config.get("fasta", "")
TX2GENE = config["tx2gene"]
METADATA = config["metadata"]
STRANDNESS = config.get("strandness", "A")
THREADS = config.get("threads", 8)

rule all:
    input:
        "results/multiqc_qc/multiqc_report.html",
        "results/multiqc_final/multiqc_report.html",
        expand("results/deseq2/{contrast}_DEG.tsv", contrast=CONTRAST_FILES),
        "results/deseq2/PCA_plot.pdf",
        "results/deseq2/normalized_counts.tsv"

rule fastqc_raw:
    input:
        r1=lambda wc: SAMPLE_TO_R1[wc.sample],
        r2=lambda wc: SAMPLE_TO_R2[wc.sample]
    output:
        html_r1="results/raw_fastqc/{sample}_R1_fastqc.html",
        zip_r1="results/raw_fastqc/{sample}_R1_fastqc.zip",
        html_r2="results/raw_fastqc/{sample}_R2_fastqc.html",
        zip_r2="results/raw_fastqc/{sample}_R2_fastqc.zip"
    log:
        "logs/fastqc_raw_{sample}.log"
    threads: 2
    conda:
        "envs/qc.yaml"
    shell:
        """
        fastqc -t {threads} -o results/raw_fastqc {input.r1} {input.r2} &> {log}
        """

rule fastp:
    input:
        r1=lambda wc: SAMPLE_TO_R1[wc.sample],
        r2=lambda wc: SAMPLE_TO_R2[wc.sample]
    output:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz",
        html="results/trimmed/{sample}.fastp.html",
        json="results/trimmed/{sample}.fastp.json"
    params:
        extra=config.get("fastp_extra", "")
    log:
        "logs/fastp_{sample}.log"
    threads:
        THREADS
    conda:
        "envs/qc.yaml"
    shell:
        """
        fastp \
          -i {input.r1} -I {input.r2} \
          -o {output.r1} -O {output.r2} \
          --thread {threads} \
          --html {output.html} \
          --json {output.json} \
          {params.extra} &> {log}
        """

rule fastqc_trimmed:
    input:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz"
    output:
        html_r1="results/trimmed_fastqc/{sample}_R1_fastqc.html",
        zip_r1="results/trimmed_fastqc/{sample}_R1_fastqc.zip",
        html_r2="results/trimmed_fastqc/{sample}_R2_fastqc.html",
        zip_r2="results/trimmed_fastqc/{sample}_R2_fastqc.zip"
    log:
        "logs/fastqc_trimmed_{sample}.log"
    threads: 2
    conda:
        "envs/qc.yaml"
    shell:
        """
        fastqc -t {threads} -o results/trimmed_fastqc {input.r1} {input.r2} &> {log}
        """

rule multiqc_qc:
    input:
        expand("results/raw_fastqc/{sample}_R1_fastqc.zip", sample=SAMPLES),
        expand("results/raw_fastqc/{sample}_R2_fastqc.zip", sample=SAMPLES),
        expand("results/trimmed_fastqc/{sample}_R1_fastqc.zip", sample=SAMPLES),
        expand("results/trimmed_fastqc/{sample}_R2_fastqc.zip", sample=SAMPLES),
        expand("results/trimmed/{sample}.fastp.json", sample=SAMPLES)
    output:
        "results/multiqc_qc/multiqc_report.html"
    log:
        "logs/multiqc_qc.log"
    conda:
        "envs/qc.yaml"
    shell:
        """
        multiqc results/raw_fastqc results/trimmed_fastqc results/trimmed \
          -o results/multiqc_qc -f &> {log}
        """

rule star_align:
    input:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz"
    output:
        bam="results/star/{sample}.Aligned.sortedByCoord.out.bam",
        log_final="results/star/{sample}.Log.final.out",
        sj="results/star/{sample}.SJ.out.tab"
    params:
        prefix="results/star/{sample}.",
        extra=config.get("star_extra", "")
    threads:
        THREADS
    log:
        "logs/star_{sample}.log"
    conda:
        "envs/alignment.yaml"
    shell:
        """
        STAR \
          --runThreadN {threads} \
          --genomeDir {STAR_INDEX} \
          --readFilesIn {input.r1} {input.r2} \
          --readFilesCommand zcat \
          --outSAMtype BAM SortedByCoordinate \
          --outFileNamePrefix {params.prefix} \
          {params.extra} &> {log}
        """

rule samtools_index:
    input:
        bam="results/star/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        bam="results/samtools/{sample}.bam",
        bai="results/samtools/{sample}.bam.bai"
    params:
        outbam="results/samtools/{sample}.bam"
    threads: 2
    log:
        "logs/samtools_{sample}.log"
    conda:
        "envs/alignment.yaml"
    shell:
        """
        cp {input.bam} {output.bam}
        samtools index -@ {threads} {output.bam} {output.bai} &> {log}
        """

rule salmon_quant:
    input:
        r1="results/trimmed/{sample}_R1.trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2.trimmed.fastq.gz"
    output:
        quant="results/salmon/{sample}/quant.sf"
    params:
        outdir="results/salmon/{sample}",
        libtype=config.get("salmon_libtype", "A"),
        extra=config.get("salmon_extra", "")
    threads:
        THREADS
    log:
        "logs/salmon_{sample}.log"
    conda:
        "envs/salmon.yaml"
    shell:
        """
        salmon quant \
          -i {SALMON_INDEX} \
          -l {params.libtype} \
          -1 {input.r1} -2 {input.r2} \
          -p {threads} \
          -o {params.outdir} \
          {params.extra} &> {log}
        """

rule rseqc:
    input:
        bam="results/samtools/{sample}.bam",
        bai="results/samtools/{sample}.bam.bai"
    output:
        infer="results/rseqc/{sample}.infer_experiment.txt",
        bamstat="results/rseqc/{sample}.bam_stat.txt",
        genebody="results/rseqc/{sample}.geneBodyCoverage.txt"
    params:
        bed=config["rseqc_bed"],
        prefix="results/rseqc/{sample}"
    log:
        "logs/rseqc_{sample}.log"
    conda:
        "envs/qc.yaml"
    shell:
        """
        infer_experiment.py -i {input.bam} -r {params.bed} > {output.infer} 2>> {log}
        bam_stat.py -i {input.bam} > {output.bamstat} 2>> {log}
        geneBody_coverage.py -r {params.bed} -i {input.bam} -o {params.prefix} >> {log} 2>&1
        cp {params.prefix}.geneBodyCoverage.txt {output.genebody}
        """

rule qualimap:
    input:
        bam="results/samtools/{sample}.bam",
        bai="results/samtools/{sample}.bam.bai"
    output:
        html="results/qualimap/{sample}/qualimapReport.html",
        stats="results/qualimap/{sample}/rnaseq_qc_results.txt"
    params:
        outdir="results/qualimap/{sample}",
        stranded=STRANDNESS
    threads: 4
    log:
        "logs/qualimap_{sample}.log"
    conda:
        "envs/qc.yaml"
    shell:
        """
        qualimap rnaseq \
          -bam {input.bam} \
          -gtf {GTF} \
          -outdir {params.outdir} \
          -p {params.stranded} \
          --java-mem-size=8G &> {log}
        """

rule multiqc_final:
    input:
        expand("results/star/{sample}.Log.final.out", sample=SAMPLES),
        expand("results/salmon/{sample}/quant.sf", sample=SAMPLES),
        expand("results/rseqc/{sample}.bam_stat.txt", sample=SAMPLES),
        expand("results/qualimap/{sample}/rnaseq_qc_results.txt", sample=SAMPLES)
    output:
        "results/multiqc_final/multiqc_report.html"
    log:
        "logs/multiqc_final.log"
    conda:
        "envs/qc.yaml"
    shell:
        """
        multiqc results/star results/salmon results/rseqc results/qualimap \
          -o results/multiqc_final -f &> {log}
        """

rule deseq2:
    input:
        quant_files=expand("results/salmon/{sample}/quant.sf", sample=SAMPLES),
        tx2gene=TX2GENE,
        metadata=METADATA
    output:
        pca="results/deseq2/PCA_plot.pdf",
        norm="results/deseq2/normalized_counts.tsv",
        deg=expand("results/deseq2/{contrast}_DEG.tsv", contrast=CONTRAST_FILES)
    params:
        outdir="results/deseq2",
        contrasts=",".join(CONTRASTS),
        reference_level=config.get("reference_level", "")
    log:
        "logs/deseq2.log"
    conda:
        "envs/deseq2.yaml"
    shell:
        """
        Rscript scripts/run_deseq2.R \
          --salmon_dir results/salmon \
          --metadata {input.metadata} \
          --tx2gene {input.tx2gene} \
          --design '{config[design]}' \
          --condition_col {config[condition_col]} \
          --contrasts {params.contrasts} \
          --outdir {params.outdir} \
          --pca {output.pca} \
          --norm {output.norm} \
          --reference_level '{params.reference_level}' &> {log}
        """
