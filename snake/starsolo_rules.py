# use STARsolo to map the raw reads
import snakemake.utils
from snakemake.io import glob_wildcards, expand
import numpy as np


# With STARSOLO_IN find the input fastq files with the cDNA reads
# the setup of standard 10X runs is assumed with scDNA in R2
# the fastq files are assumed to begin with the sample name, contain "R2", and end with "fastq.gz"
def find_cdna_read_fastqs(wildcards):
    sample_name = wildcards.sample
    directory = STARSOLO_IN
    mids = glob_wildcards(STARSOLO_IN + sample_name + "{middle}_R2_{end}.fastq.gz").middle
    ends = glob_wildcards(STARSOLO_IN + sample_name + "{middle}_R2_{end}.fastq.gz").end
    fastqs = expand(STARSOLO_IN + '{sample_name}{middle}_R2_{end}.fastq.gz',
            sample_name = sample_name,
            middle = mids,
            end = ends)
    fast_np = np.array(fastqs)
    fastqs_uniq = np.unique(fast_np)
    fastqs_comma = ','.join(fastqs_uniq)
    return fastqs_comma

# With STARSOLO_IN find the input fastq files with the barcode reads
# the setup of standard 10X runs is assumed with barcode information in R1
# the fastq files are assumed to begin with the sample name, contain "R1", and end with "fastq.gz"
def find_barcode_read_fastqs(wildcards):
    sample_name = wildcards.sample
    directory = STARSOLO_IN
    mids = glob_wildcards(STARSOLO_IN + sample_name + "{middle}_R1_{end}.fastq.gz").middle
    ends = glob_wildcards(STARSOLO_IN + sample_name + "{middle}_R1_{end}.fastq.gz").end
    fastqs = expand(STARSOLO_IN + '{sample_name}{middle}_R1_{end}.fastq.gz',
            sample_name = sample_name,
            middle = mids,
            end = ends)
    fast_np = np.array(fastqs)
    fastqs_uniq = np.unique(fast_np)
    fastqs_comma = ','.join(fastqs_uniq)
    return fastqs_comma

if not 'STARSOLO_IN' in globals():
    STARSOLO_IN = INPUTDIR
if not 'STARSOLO_OUT' in globals():
    STARSOLO_OUT = OUTDIR + 'starsolo_mapping/'

# rule to run STARsolo for mapping the reads
# Regarding formatting of input files see input functions above
rule starsolo:
    input:
        indir = STARSOLO_IN
    output:
        features_file = STARSOLO_OUT + '{sample}.features.tsv',
        matrix_file = STARSOLO_OUT +'{sample}.matrix.mtx',
        barcodes_file = STARSOLO_OUT + '{sample}.barcodes.tsv'
    params:
        cdna_reads = find_cdna_read_fastqs,
        barcode_reads = find_barcode_read_fastqs,
        sample = '{sample}',
        starsolo_dir = STARSOLO_OUT,
        outdir = STARSOLO_OUT + '{sample}',
        genomeDir = config['resources']['genome_index_starsolo'],
        soloCBwhitelist = config['resources']['soloCBwhitelist'],
        soloUMIlen = config['tools']['starsolo']['soloUMIlen'],
        variousParams = config['tools']['starsolo']['variousParams'],
        lsfoutfile = STARSOLO_OUT + '{sample}.starsolo.lsfout.log',
        lsferrfile = STARSOLO_OUT + '{sample}.starsolo.lsferr.log',
        scratch = config['tools']['starsolo']['scratch'],
        mem = config['tools']['starsolo']['mem'],
        time = config['tools']['starsolo']['time'],
    threads:
        config['tools']['starsolo']['threads']
    benchmark:
        STARSOLO_OUT + '{sample}.starsolo.benchmark'
    shell:
        'mkdir -p {params.outdir} ; '
        'cd {params.outdir} ; '
        '{config[tools][starsolo][call]} '
        '--genomeDir {params.genomeDir} '
        '--readFilesIn {params.cdna_reads} {params.barcode_reads} '
        '--soloCBwhitelist {params.soloCBwhitelist} '
        '--soloUMIlen {params.soloUMIlen} '
        '--runThreadN {config[tools][starsolo][threads]} '
        '{params.variousParams} ; '
        'ln -s {params.outdir}/Solo.out/Gene/filtered/features.tsv {output.features_file} ; '
        'ln -s {params.outdir}/Solo.out/Gene/filtered/matrix.mtx {output.matrix_file} ; '
        'ln -s {params.outdir}/Solo.out/Gene/filtered/barcodes.tsv {output.barcodes_file} '


