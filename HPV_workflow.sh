#!/bin/bash

# Step-by-step setup for HPV16 Nanopore Analysis Pipeline

# 1. Create project structure
mkdir -p hpv_project/{{raw,trimmed,reference,results/{alignments,consensus,variants,igv}}}
cd hpv_project

# 2. Place files:
for file in {01..24}; do 
cd barcode$file/ && unpigz * && cat *.fastq > barcode$file.fastq && pigz * && rm FA* && cd ..;
done
# - Place raw reads (FASTQ) as barcode01.fastq ... barcodeN.fastq into raw/
# - Place the HPV16 reference FASTA into reference/hpv16_reference.fasta

# 3. Save the Snakefile
cat << 'EOF' > Snakefile
import glob
SAMPLES = [s.split("/")[-1].replace(".fastq", "") for s in glob.glob("raw/*.fastq")]

rule all:
    input:
        expand("results/consensus/{sample}_consensus.fasta", sample=SAMPLES),
        expand("results/variants/{sample}.vcf", sample=SAMPLES),
        expand("results/alignments/{sample}.sorted.bam.bai", sample=SAMPLES)

rule trim_primers:
    input:
        fq="raw/{sample}.fastq"
    output:
        fq_trimmed="trimmed/{sample}.trimmed.fastq"
    shell:
        """
        cutadapt -g TATAGTTCCAGGGTCTCCAC -a TTTACAAGCACACATACAAGCA \
                -o {output.fq_trimmed} {input.fq} > trimmed/{wildcards.sample}.cutadapt.log
        """

rule align_reads:
    input:
        fq="trimmed/{sample}.trimmed.fastq",
        ref="reference/hpv16_reference.fasta"
    output:
        bam="results/alignments/{sample}.sorted.bam",
        bai="results/alignments/{sample}.sorted.bam.bai"
    threads: 4
    shell:
        """
        minimap2 -ax map-ont {input.ref} {input.fq} | \
        samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

rule medaka_consensus:
    input:
        fq="trimmed/{sample}.trimmed.fastq",
        ref="reference/hpv16_reference.fasta"
    output:
        "results/consensus/{sample}_consensus.fasta"
    params:
        outdir="results/consensus/{sample}"
    threads: 4
    shell:
        """
        medaka_consensus -i {input.fq} -d {input.ref} -o {params.outdir} -t {threads} -m r941_min_high_g360
        cp {params.outdir}/consensus.fasta {output}
        """

rule call_variants:
    input:
        ref="reference/hpv16_reference.fasta",
        bam="results/alignments/{sample}.sorted.bam"
    output:
        vcf="results/variants/{sample}.vcf"
    shell:
        """
        bcftools mpileup -f {input.ref} {input.bam} | \
        bcftools call -mv -Ov -o {output.vcf}
        """

rule create_igv_index:
    input:
        bam="results/alignments/{sample}.sorted.bam"
    output:
        bai="results/alignments/{sample}.sorted.bam.bai"
    shell:
        "samtools index {input.bam}"

rule prepare_igv_ready:
    input:
        bam="results/alignments/{sample}.sorted.bam",
        bai="results/alignments/{sample}.sorted.bam.bai"
    output:
        tdf="results/igv/{sample}.tdf"
    params:
        ref="reference/hpv16_reference.fasta"
    shell:
        "igvtools count -z 5 -w 25 -e 0 {input.bam} {output.tdf} {params.ref}"

EOF

# 4. Save environment YAML
cat << 'EOF' > environment.yml
name: hpv_env
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - snakemake
  - minimap2
  - samtools
  - bcftools
  - medaka
  - porechop
  - seqtk
  - bedtools
  - cutadapt
  - python=3.9
  - igvtools
  - pandas
EOF

# 5. Create and activate environment
mamba env create -f environment.yml -n hpv_env || conda env create -f environment.yml -n hpv_env
conda activate hpv_env

# 6. Run pipeline
# snakemake --cores 4 --use-condasnakemake --cores 4 --use-conda --conda-prefix ~/.snakemake/conda


# You're done. Results will be in hpv_project/results/
