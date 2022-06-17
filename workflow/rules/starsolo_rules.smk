# use STARsolo to map the raw reads
import snakemake.utils
from snakemake.io import glob_wildcards, expand
import numpy as np


# Find the input fastq files with the cDNA reads
# setup of standard 10X runs is assumed with scDNA in R2
# fastq files are assumed to begin with the sample name, contain "R2", and end with "fastq.gz"
def find_cdna_read_fastqs(wildcards):
    sample_name = wildcards.sample
    mids = glob_wildcards(config["inputOutput"]["input_fastqs"] + sample_name + "{middle}_R2_{end}.fastq.gz").middle
    ends = glob_wildcards(config["inputOutput"]["input_fastqs"] + sample_name + "{middle}_R2_{end}.fastq.gz").end
    fastqs = expand(config["inputOutput"]["input_fastqs"] + '{sample_name}{middle}_R2_{end}.fastq.gz',
            sample_name = sample_name,
            middle = mids,
            end = ends)
    fast_np = np.array(fastqs)
    fastqs_uniq = np.unique(fast_np)
    fastqs_comma = ','.join(fastqs_uniq)
    return fastqs_comma

# Find the input fastq files with the barcode reads
# setup of standard 10X runs is assumed with barcode information in R1
# fastq files are assumed to begin with the sample name, contain "R1", and end with "fastq.gz"
def find_barcode_read_fastqs(wildcards):
    sample_name = wildcards.sample
    mids = glob_wildcards(config["inputOutput"]["input_fastqs"] + sample_name + "{middle}_R1_{end}.fastq.gz").middle
    ends = glob_wildcards(config["inputOutput"]["input_fastqs"] + sample_name + "{middle}_R1_{end}.fastq.gz").end
    fastqs = expand(config["inputOutput"]["input_fastqs"] + '{sample_name}{middle}_R1_{end}.fastq.gz',
            sample_name = sample_name,
            middle = mids,
            end = ends)
    fast_np = np.array(fastqs)
    fastqs_uniq = np.unique(fast_np)
    fastqs_comma = ','.join(fastqs_uniq)
    return fastqs_comma

# rule to run STARsolo for mapping the reads
# About formatting of input files see input functions above
rule starsolo:
    input:
        indir = config["inputOutput"]["input_fastqs"]
    output:
        features_file = 'results/starsolo/{sample}.features.tsv',
        matrix_file = 'results/starsolo/{sample}.matrix.mtx',
        barcodes_file = 'results/starsolo/{sample}.barcodes.tsv'
    params:
        cdna_reads = find_cdna_read_fastqs,
        barcode_reads = find_barcode_read_fastqs,
        sample = '{sample}',
        starsolo_dir = 'results/starsolo/',
        outdir = 'results/starsolo/{sample}',
        genomeDir = config['resources']['genome_index_starsolo'],
        soloCBwhitelist = config['resources']['soloCBwhitelist'],
        soloUMIlen = config['tools']['starsolo']['soloUMIlen'],
        variousParams = config['tools']['starsolo']['variousParams'],
    resources:
        mem_mb = config['tools']['starsolo']['mem'],
        time_min = config['tools']['starsolo']['time']
    threads:
        config['tools']['starsolo']['threads']
    log:
        "logs/starsolo/{sample}.log"
    benchmark:
        'logs/starsolo/{sample}.benchmark'
    shell:
        'mkdir -p {params.outdir} ; '
        '(cd {params.outdir} ; '
        '{config[tools][starsolo][call]} '
        '--genomeDir {params.genomeDir} '
        '--readFilesIn {params.cdna_reads} {params.barcode_reads} '
        '--soloCBwhitelist {params.soloCBwhitelist} '
        '--soloUMIlen {params.soloUMIlen} '
        '--runThreadN {config[tools][starsolo][threads]} '
        '{params.variousParams} ) ; '
        'pwd ; '
        'ln -sr {params.outdir}/Solo.out/Gene/filtered/features.tsv {output.features_file} ; '
        'ln -sr {params.outdir}/Solo.out/Gene/filtered/matrix.mtx {output.matrix_file} ; '
        'ln -sr {params.outdir}/Solo.out/Gene/filtered/barcodes.tsv {output.barcodes_file} '

