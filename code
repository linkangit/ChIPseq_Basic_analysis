# Generic ChIP-seq Analysis Pipeline
# Simple, clean, and easy to use

configfile: "config.yaml"

# Define your samples here or in config.yaml
SAMPLES = config.get("samples", {
    "input_sample": "input",      # Your input/control sample
    "chip_sample1": "H3K4me3",    # Your ChIP samples
    "chip_sample2": "H3K27ac"
})

# All sample names
ALL_SAMPLES = list(SAMPLES.keys())

rule all:
    input:
        # Quality control
        expand("results/qc/{sample}_fastqc.html", sample=ALL_SAMPLES),
        expand("results/qc/{sample}.flagstat", sample=ALL_SAMPLES),
        
        # Alignments
        expand("results/bam/{sample}.bam", sample=ALL_SAMPLES),
        expand("results/bam/{sample}.bam.bai", sample=ALL_SAMPLES),
        
        # Peak calling (ChIP samples vs input)
        expand("results/peaks/{chip}_peaks.narrowPeak", 
               chip=[s for s in ALL_SAMPLES if s != "input_sample"])

# Step 1: Quality control of raw reads
rule fastqc:
    input:
        "data/{sample}.fastq.gz"
    output:
        html="results/qc/{sample}_fastqc.html",
        zip="results/qc/{sample}_fastqc.zip"
    threads: 2
    shell:
        """
        fastqc -t {threads} -o results/qc {input}
        """

# Step 2: Align reads to genome
rule align:
    input:
        "data/{sample}.fastq.gz"
    output:
        temp("results/temp/{sample}.unsorted.bam")
    threads: 8
    params:
        genome=config.get("genome_index", "genome/bowtie2_index")
    shell:
        """
        bowtie2 -x {params.genome} -U {input} -p {threads} \
        | samtools view -Sb - > {output}
        """

# Step 3: Sort BAM files
rule sort_bam:
    input:
        "results/temp/{sample}.unsorted.bam"
    output:
        "results/bam/{sample}.bam"
    threads: 4
    shell:
        """
        samtools sort -@ {threads} -o {output} {input}
        """

# Step 4: Index BAM files
rule index_bam:
    input:
        "results/bam/{sample}.bam"
    output:
        "results/bam/{sample}.bam.bai"
    shell:
        """
        samtools index {input}
        """

# Step 5: Get alignment statistics
rule flagstat:
    input:
        "results/bam/{sample}.bam"
    output:
        "results/qc/{sample}.flagstat"
    shell:
        """
        samtools flagstat {input} > {output}
        """

# Step 6: Call peaks (ChIP vs Input)
rule call_peaks:
    input:
        chip="results/bam/{chip}.bam",
        control="results/bam/input_sample.bam"
    output:
        peaks="results/peaks/{chip}_peaks.narrowPeak",
        summits="results/peaks/{chip}_summits.bed"
    params:
        name="{chip}",
        genome_size=config.get("genome_size", "hs")
    shell:
        """
        macs2 callpeak \
            -t {input.chip} \
            -c {input.control} \
            -f BAM \
            -g {params.genome_size} \
            -n {params.name} \
            --outdir results/peaks \
            -q 0.01
        """
