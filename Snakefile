import pandas as pd
import numpy as np
import os
import re
import itertools
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("6.15.0")

##### Editable variables #####

configfile: "bin/config.yaml"

bamCoverage_binsize = 10


ref_fa = config['ref']['sequence']
ref_fai = config['ref']['fai']
gtf = config['ref']['annotation']
bwa_index = config['ref']['index']

mito_chrom=config['ref']['mito_chr']
fai_parsed = pd.read_table(ref_fai, names=['chr','len','offset','bases_per_line','bytes_per_line'])
if mito_chrom not in fai_parsed['chr'].values:
    raise Exception('{mito} not found in reference fai file.'.format(mito=mito_chrom))

chroms_no_mito = ' '.join(fai_parsed[fai_parsed['chr'] != mito_chrom]['chr'].values)


##### load config and sample sheets #####
samplesheet="bin/samples.tsv"
units = pd.read_table(samplesheet, dtype={"sample" : str})
units['se_or_pe'] = ["SE" if x else "PE" for x in units['fq2'].isnull()]

samples = units[["sample","se_or_pe"]].drop_duplicates()
if not samples['sample'].is_unique:
    raise Exception('A sample has more than one combination of and se_or_pe.')


# if SE data, error out
if (units['fq2'].isnull().any()):
    raise Exception('SE not supported.')

snakemake_dir = os.getcwd() + "/"

# make a tmp directory for analyses
tmp_dir = os.path.join(snakemake_dir, "tmp")
if not os.path.exists(tmp_dir):
    os.mkdir(tmp_dir)

shared_snakemake_dir=config['shared_snakemake_repo'] + "/rules/"

include: os.path.join(shared_snakemake_dir, "post_alignment/flagstat")
include: os.path.join(shared_snakemake_dir, "post_alignment/CollectInsertSizeMetrics")
include: os.path.join(shared_snakemake_dir, "post_alignment/CollectAlignmentSummaryMetrics")

rule all:
    input:
        "analysis/multiqc/multiqc_report.html",

def get_orig_fastq(wildcards):
    if wildcards.read == "R1":
            fastq = expand("raw_data/{fq}", fq = units[units["sample"] == wildcards.sample]["fq1"].values)
    elif wildcards.read == "R2":
            fastq = expand("raw_data/{fq}", fq = units[units["sample"] == wildcards.sample]["fq2"].values)
    return fastq 

rule rename_fastqs:
    """
    Rename fastqs by biologically meaningful name. Concatenate different runs of same library.
    """
    input:
        get_orig_fastq
    output:
        "analysis/renamed_data/{sample}_{read}.fastq.gz"
    log:
        stdout="logs/rename_fastqs/{sample}_{read}.o",
        stderr="logs/rename_fastqs/{sample}_{read}.e",
    benchmark:
        "benchmarks/rename_fastqs/{sample}_{read}.txt"
    params:
    threads: 1
    resources:
        mem_gb=8
    envmodules:
    shell:
        """
        if [ `printf '{input}' | wc -w` -gt 1 ]
        then
            cat {input} > {output}
        else
            ln -sr {input} {output}
        fi

        """

rule fastqc:
    """
    Run fastqc on raw_data/ files.
    """
    input:
        "analysis/renamed_data/{fq_pref}.fastq.gz"
    output:
        html="analysis/fastqc/{fq_pref}_fastqc.html",
        zip="analysis/fastqc/{fq_pref}_fastqc.zip"
    params:
        outdir="analysis/fastqc/"
    log:
        stdout="logs/fastqc/{fq_pref}.o",
        stderr="logs/fastqc/{fq_pref}.e"
    benchmark:
        "benchmarks/fastqc/{fq_pref}.txt"
    envmodules:
        config['modules']['fastqc']
    threads: 1
    resources:
        mem_gb = 32
    shell:
        """
        fastqc --outdir {params.outdir} {input}
        """

rule fastq_screen:
    """
    Run fastq_screen to detect any contamination from other species or excessive rRNA.
    """
    input:
        "analysis/renamed_data/{fq_pref}.fastq.gz"
    output:
        html = "analysis/fastq_screen/{fq_pref}_screen.html",
        txt = "analysis/fastq_screen/{fq_pref}_screen.txt",
    params:
    log:
        stdout="logs/fastq_screen/{fq_pref}.o",
        stderr="logs/fastq_screen/{fq_pref}.e"
    benchmark:
        "benchmarks/fastq_screen/{fq_pref}.txt"
    envmodules:
        config['modules']['fastq_screen']
    threads: 8
    resources:
        mem_gb = 32
    shell:
        """
        fastq_screen --threads {threads} --outdir analysis/fastq_screen/ {input}
        """

rule trim_galore_PE:
    """
    Run trim_galore on paired-end reads.
    """
    input:
        expand("analysis/renamed_data/{{sample}}_R{read}.fastq.gz", read=["1","2"])
    output:
        temp(expand("analysis/trim_galore/{{sample}}_R{ext}", ext=["1_val_1.fq.gz","2_val_2.fq.gz"])),
        expand("analysis/trim_galore/{{sample}}_R1{ext}", ext=[".fastq.gz_trimming_report.txt","_val_1_fastqc.html"]),
        expand("analysis/trim_galore/{{sample}}_R2{ext}", ext=[".fastq.gz_trimming_report.txt","_val_2_fastqc.html"]),
    params:
    log:
        stdout="logs/trim_galore/{sample}.o",
        stderr="logs/trim_galore/{sample}.e"
    benchmark:
        "benchmarks/trim_galore/{sample}.txt"
    envmodules:
        config['modules']['trim_galore']
    threads: 4
    resources:
        mem_gb = 64
    shell:
        """
        trim_galore --paired {input} --output_dir analysis/trim_galore/ --cores {threads} --fastqc
        """

rule trim_galore_SE:
    """
    Run trim_galore on single-end reads.
    """
    input:
        "analysis/renamed_data/{sample}_R1.fastq.gz"
    output:
        temp("analysis/trim_galore/{sample}_R1_trimmed.fq.gz"),
        "analysis/trim_galore/{sample}_R1.fastq.gz_trimming_report.txt",
        expand("analysis/trim_galore/{{sample}}_R1_trimmed_fastqc{ext}", ext=['.html','.zip']),
    params:
    log:
        stdout="logs/trim_galore/{sample}.o",
        stderr="logs/trim_galore/{sample}.e"
    benchmark:
        "benchmarks/trim_galore/{sample}.txt"
    envmodules:
        config['modules']['trim_galore']
    threads: 4
    resources:
        mem_gb = 64
    shell:
        """
        trim_galore {input} --output_dir analysis/trim_galore/ --cores {threads} --fastqc
        """

rule bwamem:
    input:
        lambda wildcards: expand("analysis/trim_galore/{sample}_R{read}_val_{read}.fq.gz", read=["1","2"], sample=wildcards.sample) if samples[samples['sample']==wildcards.sample]['se_or_pe'].values=="PE" else "analysis/trim_galore/{sample}_R1_trimmed.fq.gz"
    output:
        outbam=temp("analysis/bwamem/{sample}.bam"),
        outbai="analysis/bwamem/{sample}.bam.bai",
        idxstat="analysis/bwamem/{sample}.bam.idxstat",
        samblaster_err="analysis/bwamem/{sample}.samblaster.e",
    log:
        stdout="logs/bwamem/{sample}.o",
        stderr="logs/bwamem/{sample}.e",
    benchmark:
        "benchmarks/bwamem/{sample}.txt"
    params:
        bwa_RG='"@RG\tID:{sample}\tLB:{sample}\tPL:ILLUMINA\tSM:{sample}"',
        bwa_idx=bwa_index,
        samblaster_params=lambda wildcards: "" if samples[samples['sample']==wildcards.sample]['se_or_pe'].values=="PE" else "--ignoreUnmated"
    threads: 16
    envmodules:
        config['modules']['bwa'],
        config['modules']['samblaster'],
        config['modules']['samtools'],
    resources:
        mem_gb=180
    shell:
        """
        bwa mem \
        -t {threads} \
        -R {params.bwa_RG} \
        {params.bwa_idx} \
        {input} | \
        samblaster {params.samblaster_params} 2>{output.samblaster_err} | \
        samtools sort \
        -m 6G \
        -@ {threads} \
        -O "BAM" \
        -o {output.outbam} \
        -

        echo "END bwamem"
        echo "END bwamem" 1>&2

        samtools index -@ {threads} {output.outbam}

        echo "END indexing"
        echo "END indexing" 1>&2
        
        samtools idxstats {output.outbam} > {output.idxstat}

        echo "END idxstats"
        echo "END idxstats" 1>&2
 
        """

rule qualimap:
    """
    Run qualimap on unfiltered alignments.
    """
    input:
        "analysis/bwamem/{sample}.bam",
    output:
        touch("analysis/qualimap/{sample}/done")
    log:
        stdout="logs/qualimap/{sample}.o",
        stderr="logs/qualimap/{sample}.e"
    benchmark:
        "benchmarks/qualimap/{sample}.txt"
    envmodules:
        config['modules']['qualimap']
    params:
    resources:
        mem_gb=100
    threads: 8
    shell:
        """
        qualimap bamqc -bam {input} --java-mem-size={resources.mem_gb}G --paint-chromosome-limits -outdir analysis/qualimap/{wildcards.sample} -nt {threads}

        """

rule multiqc:
    input:
        expand("analysis/fastqc/{sample.sample}_R1_fastqc.html", sample=samples.itertuples()),
        expand("analysis/fastqc/{sample.sample}_R2_fastqc.html", sample=samples[samples['se_or_pe']=="PE"].itertuples()),
        expand("analysis/trim_galore/{sample.sample}_R{read}_val_{read}_fastqc.html", sample=samples[samples['se_or_pe']=="PE"].itertuples(), read=["1","2"]),
        expand("analysis/trim_galore/{sample.sample}_R1_trimmed_fastqc.html", sample=samples[samples['se_or_pe']=="SE"].itertuples()),
        expand("analysis/bwamem/flagstat/{sample.sample}.flagstat", sample=samples.itertuples()),
        expand("analysis/bwamem/CollectInsertSizeMetrics/{sample.sample}.insert_size_metrics.txt", sample=samples[samples['se_or_pe']=="PE"].itertuples()),
        expand("analysis/bwamem/{sample.sample}.bam.idxstat", sample=samples.itertuples()),
        expand("analysis/bwamem/CollectAlignmentSummaryMetrics/{sample.sample}.aln_metrics.txt", sample=samples.itertuples()),
        expand("analysis/bwamem/{sample.sample}.samblaster.e", sample=samples.itertuples()),
        expand("analysis/fastq_screen/{sample.sample}_R1_screen.html", sample=samples.itertuples()),
        expand("analysis/fastq_screen/{sample.sample}_R2_screen.html", sample=samples[samples['se_or_pe']=="PE"].itertuples()),
        expand("analysis/qualimap/{sample.sample}/done", sample=samples.itertuples()),
    output:
        "analysis/multiqc/multiqc_report.html",
    log:
        stdout="logs/multiqc/multiqc.o",
        stderr="logs/multiqc/multiqc.e"
    benchmark:
        "benchmarks/multiqc/multiqc.txt"
    params:
        workdir="analysis/multiqc/",
        dirs=["analysis/trim_galore/",
        "analysis/fastqc/",
        "analysis/qualimap/",
        "analysis/preseq_complexity/",
        "analysis/bwamem/flagstat/",
        "analysis/bwamem/",
        "analysis/fastq_screen/",
        "analysis/bwamem/CollectAlignmentSummaryMetrics/"],
        PE_dirs=["analysis/filt_bams/CollectInsertSizeMetrics/"] if not all(x == "SE" for x in samples['se_or_pe'].values)  else [],
        outfile="multiqc_report"
    envmodules:
        config['modules']['multiqc']
    threads: 4
    resources:
        mem_gb=100
    shell:
        """
        multiqc \
        --force \
        --config bin/multiqc_config.yaml \
        --outdir {params.workdir} \
        --filename {params.outfile} \
        {params.dirs} {params.PE_dirs}

        """

