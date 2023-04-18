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

contigs_file = "bin/grouped_contigs.tsv"

# read in file with col1 as contig group name and col2 as a commma separated list of contigs. This was written for use with GATK for splitting the variant calling steps by chromosome/contig. We group the unplaced contigs together since those tend to be small.
contig_groups = pd.read_table(contigs_file)

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
        #expand("analysis/final/06_filter_vcf/all.merged.filt.PASS.{var_type}.vcf.gz.vt_peek.txt", var_type=["SNP", "INDEL"]),
        "analysis/final/07_snpEff/all.merged.filt.PASS.snpEff.vcf.gz"

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
        expand("analysis/trim_galore/{{sample}}_R{ext}", ext=["1_val_1.fq.gz","2_val_2.fq.gz"]),
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
        outbam="analysis/bwamem/{sample}.bam",
        outbai="analysis/bwamem/{sample}.bam.bai",
        idxstat="analysis/bwamem/{sample}.bam.idxstat",
        samblaster_err="analysis/bwamem/{sample}.samblaster.e",
    log:
        stdout="logs/bwamem/{sample}.o",
        stderr="logs/bwamem/{sample}.e",
    benchmark:
        "benchmarks/bwamem/{sample}.txt"
    params:
        bwa_RG='"@RG\\tID:{sample}\\tLB:{sample}\\tPL:ILLUMINA\\tSM:{sample}"',
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

rule haplotypecaller:
    input:
        bam=lambda wildcards: "analysis/bwamem/{sample}.bam" if wildcards.round == "bqsr" else "analysis/final/00_BQSR/{sample}.bqsr.bam"
    output:
        "analysis/{round}/01_haplotypecaller/{sample}.{contig_group}.g.vcf.gz"
    log:
        stdout="logs/{round}/01_haplotypecaller/{sample}.{contig_group}.o",
        stderr="logs/{round}/01_haplotypecaller/{sample}.{contig_group}.e"
    benchmark:
        "benchmarks/{round}/01_haplotypecaller/{sample}.{contig_group}.txt"
    params:
        dbsnp=lambda wildcards: f'--dbsnp {config["ref"]["known_snps"]}' if config["ref"]["known_snps"] != "" else "",
        ref_fasta=config["ref"]["sequence"],
        contigs = lambda wildcards: "-L " + contig_groups[contig_groups.name == wildcards.contig_group]['contigs'].values[0].replace(",", " -L "),
        
    envmodules:
        config["modules"]["gatk"]
    threads: 4
    resources: 
        mem_gb = 80
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        HaplotypeCaller \
        -R {params.ref_fasta} \
        -I {input.bam} \
        -O {output} \
        -ERC GVCF \
        --native-pair-hmm-threads {threads} \
        {params.dbsnp} \
        {params.contigs}
        """

rule combinevar:
    input:
        lambda wildcards: expand("analysis/{{round}}/01_haplotypecaller/{sample}.{contig_group}.g.vcf.gz", sample=samples['sample'].unique(), contig_group=wildcards.contig_group)

    output:
        touch=touch("analysis/{round}/02_combinevar/{contig_group}.done"),
        genomicsdb=directory("analysis/{round}/02_combinevar/{contig_group}.genomicsdb"),
    log:
        stdout="logs/{round}/02_combinevar/all.{contig_group}.o",
        stderr="logs/{round}/02_combinevar/all.{contig_group}.e"
    benchmark:
        "benchmarks/{round}/02_combinevar/{contig_group}.txt"
    params:
        sample_gvcfs = lambda wildcards, input: list(map("-V {}".format, input)),
        contigs = lambda wildcards: "-L " + contig_groups[contig_groups.name == wildcards.contig_group]['contigs'].values[0].replace(",", " -L "),
    envmodules:
        config["modules"]["gatk"]
    threads: 4
    resources: 
        mem_gb = 80
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        GenomicsDBImport \
        {params.sample_gvcfs} \
        --genomicsdb-workspace-path {output.genomicsdb} \
        {params.contigs}
        """

rule jointgeno:
    input:
        "analysis/{round}/02_combinevar/{contig_group}.done"
    output:
        vcf="analysis/{round}/03_jointgeno/all.{contig_group}.vcf.gz",
    log:
        stdout="logs/{round}/03_jointgeno/all.{contig_group}.o",
        stderr="logs/{round}/03_jointgeno/all.{contig_group}.e"
    benchmark:
        "benchmarks/{round}/03_jointgeno/all.{contig_group}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        genomicsdb="analysis/{round}/02_combinevar/{contig_group}.genomicsdb"
    envmodules:
        config["modules"]["gatk"]
    threads: 4
    resources: 
        mem_gb = 80
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        GenotypeGVCFs \
        -R {params.ref_fasta} \
        -V gendb://{params.genomicsdb} \
        -O {output.vcf}
        """

rule sortVCF:
    """
    Sort the output VCFs from joint genotyping. Merging errors out sometimes if we do not do this step.
    """
    input:
        vcf="analysis/{round}/03_jointgeno/all.{contig_group}.vcf.gz",
    output:
        sorted_vcf="analysis/{round}/04_sortvcf/all.{contig_group}.sort.vcf.gz"
    log:
        stdout="logs/{round}/04_sortvcf/all.{contig_group}.o",
        stderr="logs/{round}/04_sortvcf/all.{contig_group}.e"
    benchmark:
        "benchmarks/{round}/04_sortvcf/all.{contig_group}.txt"
    params:
        dictionary=config['ref']['dict'],
    envmodules:
        config["modules"]["gatk"]
    threads: 4
    resources: 
        mem_gb = 80
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        SortVcf \
        -I {input.vcf} \
        -O {output.sorted_vcf} \
        -SD {params.dictionary} 
        """

rule merge_vcf:
    """
    Merge the contig group VCFs into one unified VCF.
    """
    input:
        expand("analysis/{{round}}/04_sortvcf/all.{contig_grp}.sort.vcf.gz", contig_grp=contig_groups.name)
    output:
        raw="analysis/{round}/05_merge_vcf/all.merged.vcf.gz",
        raw_tbi="analysis/{round}/05_merge_vcf/all.merged.vcf.gz.tbi",
        vt_peek_raw="analysis/{round}/05_merge_vcf/all.merged.vcf.gz.vt_peek.txt",
    log:
        stdout="logs/{round}/05_merge_vcf/out.o",
        stderr="logs/{round}/05_merge_vcf/err.e"
    benchmark:
        "benchmarks/{round}/05_merge_vcf/benchmark.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        dictionary=config['ref']['dict'],
        in_vcfs = lambda wildcards, input: ' '.join(['--INPUT ' + vcf for vcf in input]) 
    envmodules:
        config["modules"]["gatk"],
        config["modules"]["vt"],
    threads: 4
    resources: 
        mem_gb = 80
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        MergeVcfs \
        {params.in_vcfs} \
        --SEQUENCE_DICTIONARY {params.dictionary} \
        --OUTPUT {output.raw} 
        
        vt peek -r {params.ref_fasta} {output.raw} 2> {output.vt_peek_raw} 1>>{log.stdout}
        """

def get_filt_params (wildcards):
    # parameters adapted from https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/
    if (wildcards.var_type == "SNP"):
        return """--cluster-window-size 10 --cluster-size 3 \
        --filter-name 'QD' --filter 'QD < 2.0' \
        --filter-name 'FS' --filter 'FS > 60.0' \
        --filter-name 'MQ_filter' --filter 'MQ < 40.0' \
        --filter-name 'SOR_filter' --filter 'SOR > 4.0' \
        --filter-name 'MQRankSum_filter' --filter 'MQRankSum < -12.5' \
        --filter-name 'ReadPosRankSum_filter' --filter 'ReadPosRankSum < -8.0'"""
    elif (wildcards.var_type == "INDEL"):
        return """--filter-name 'QD' --filter 'QD < 2.0' \
        --filter-name 'FS' --filter 'FS > 200.0' \
        --filter-name 'SOR_filter' --filter 'SOR > 10.0'"""
    else:
        raise Exception("var_type wildcard must be either SNP or INDEL.")

rule filter_vcf:
    """
    Do quality filters. Use different paramters depending on SNPs verus indels ('SNP' or 'INDEL').
    """
    input:
        "analysis/{round}/05_merge_vcf/all.merged.vcf.gz"
    output:
        raw="analysis/{round}/06_filter_vcf/all.merged.{var_type}.vcf.gz",
        filt="analysis/{round}/06_filter_vcf/all.merged.filt.{var_type}.vcf.gz",
        pass_only="analysis/{round}/06_filter_vcf/all.merged.filt.PASS.{var_type}.vcf.gz",
        vt_peek_pass="analysis/{round}/06_filter_vcf/all.merged.filt.PASS.{var_type}.vcf.gz.vt_peek.txt"
    log:
        stdout="logs/{round}/06_filter_vcf/{var_type}.o",
        stderr="logs/{round}/06_filter_vcf/{var_type}.e"
    benchmark:
        "benchmarks/{round}/06_filter_vcf/{var_type}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        filt_params=get_filt_params
    envmodules:
        config["modules"]["gatk"],
        config["modules"]["vt"]
    threads: 4
    resources: 
        mem_gb = 80
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        SelectVariants \
        -R {params.ref_fasta} \
        -V {input} \
        --select-type-to-include {wildcards.var_type} \
        -O {output.raw}

        echo "SelectVariants 1 done." >> {log.stdout}
        echo "SelectVariants 1 done." >> {log.stderr}
        
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        VariantFiltration \
        --R {params.ref_fasta} \
        --V {output.raw} \
        {params.filt_params} \
        -O {output.filt} 
        
        echo "VariantFiltration done." >> {log.stdout}
        echo "VariantFiltration done." >> {log.stderr}

        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        SelectVariants \
        -R {params.ref_fasta} \
        -V {output.filt} \
        --exclude-filtered \
        -O {output.pass_only} 
        
        echo "SelectVariants 2 done." >> {log.stdout}
        echo "SelectVariants 2 done." >> {log.stderr}

        vt peek -r {params.ref_fasta} {output.pass_only} 2> {output.vt_peek_pass} 1>>{log.stdout}
        """

rule BQSR:
    """
    Base quality score recalibrator
    """
    input:
        bam="analysis/bwamem/{sample}.bam",
        known_snps_vcf=config["ref"]["known_snps"] if config["ref"]["known_snps"] else "analysis/bqsr/06_filter_vcf/all.merged.filt.PASS.SNP.vcf.gz",
        known_indels_vcf=config["ref"]["known_indels"] if config["ref"]["known_indels"] else "analysis/bqsr/06_filter_vcf/all.merged.filt.PASS.INDEL.vcf.gz",
    output:
        recal_table="analysis/final/00_BQSR/{sample}.bqsr.table",
        recal_bam="analysis/final/00_BQSR/{sample}.bqsr.bam"
    log:
        stdout="logs/final/00_BQSR/{sample}.o",
        stderr="logs/final/00_BQSR/{sample}.e"
    benchmark:
        "benchmarks/final/00_BQSR/{sample}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
    envmodules:
        config["modules"]["gatk"],
    threads: 4
    resources: 
        mem_gb = 80
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        BaseRecalibrator \
        -R {params.ref_fasta} \
        -I {input.bam} \
        --known-sites {input.known_snps_vcf} \
        --known-sites {input.known_indels_vcf} \
        -O {output.recal_table} 
 
        echo "BaseRecalibrator done." >> {log.stdout}
        echo "BaseRecalibrator done." >> {log.stderr}
        
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        ApplyBQSR \
        -R {params.ref_fasta} \
        -I {input.bam} \
        -bqsr {output.recal_table} \
        -O {output.recal_bam}

        echo "ApplyBQSR done." >> {log.stdout}
        echo "ApplyBQSR done." >> {log.stderr}

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

rule snpEff:
    """
    Concatenate SNPs and indels and run snpEff.
    """
    input:
        vcfs=expand("analysis/final/06_filter_vcf/all.merged.filt.PASS.{var_type}.vcf.gz", var_type=["SNP", "INDEL"]),
    output:
        html="analysis/final/07_snpEff/snpEff.html",
        genes="analysis/final/07_snpEff/snpEff.genes.txt",
        merged_vcf="analysis/final/07_snpEff/all.merged.filt.PASS.vcf.gz",
        annot_vcf="analysis/final/07_snpEff/all.merged.filt.PASS.snpEff.vcf.gz"
    log:
        stdout="logs/final/07_snpEff/out.o",
        stderr="logs/final/07_snpEff/err.e"
    benchmark:
        "benchmarks/final/07_snpEff/bench.txt"
    envmodules:
        config['modules']['snpEff'],
        config['modules']['bcftools']
    params:
        snpEff_dataDir=f'-dataDir {config["ref"]["snpEff_dataDir"]}' if config['ref']['snpEff_dataDir'] else "",
        snpEff_genome_id=config['ref']['snpEff_genome_ID']
    resources:
        mem_gb=100
    threads: 8
    shell:
        """
        # merge
        bcftools concat --threads {threads} --allow-overlaps -O z -o {output.merged_vcf} {input.vcfs}
        bcftools index --tbi --threads {threads} {output.merged_vcf}
        
        # annotate
        java -Xms8g -Xmx80g -Djava.io.tmpdir=./tmp -jar $SNPEFF/snpEff.jar {params.snpEff_dataDir} -v -stats {output.html} {params.snpEff_genome_id} {output.merged_vcf} | bcftools view --threads {threads} -O z -o {output.annot_vcf} -
        bcftools index --tbi --threads {threads} {output.annot_vcf}

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
        "analysis/bwamem/flagstat/",
        "analysis/bwamem/",
        "analysis/fastq_screen/",
        "analysis/bwamem/CollectAlignmentSummaryMetrics/"],
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
        --outdir {params.workdir} \
        --filename {params.outfile} \
        {params.dirs} 

        """

