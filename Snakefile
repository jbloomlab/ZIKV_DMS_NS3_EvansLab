"""``snakemake`` pipeline that runs analysis."""


import os

import pandas as pd


configfile: 'config.yml'

wildcard_constraints:
    tile="tile_\d+"

rule all:
    input:
        expand("results/{tile}",
               tile=config['tiles']),
        expand("results/summary/dms_{tile}_analysis.md",
               tile=config['tiles']),


rule jupnb_to_md:
    """Convert Jupyter notebook to Markdown format."""
    input: notebook="results/notebooks/{notebook}.ipynb"
    output: markdown="results/summary/{notebook}.md"
    params: outdir=lambda wildcards, output: os.path.dirname(output.markdown)
    conda: 'environment.yml'
    shell: 
        """
        jupyter nbconvert \
            --output-dir {params.outdir} \
            --to markdown \
            {input.notebook}
        """

rule dms_tile_analysis:
    """Analyze DMS data for a tile."""
    input:
        amplicon="data/{tile}_amplicon.fasta",
        alignspecs="data/{tile}_subamplicon_alignspecs.txt",
        samplelist="data/{tile}_samplelist.csv",
    output:
        resultsdir=directory("results/{tile}"),
    params:
        errpre=lambda wc: config['tiles'][wc.tile]['errpre'],
        site_number_offset=lambda wc: config['tiles'][wc.tile]['site_number_offset']
    threads: config['max_cpus']
    conda: 'environment.yml'
    log: notebook='results/notebooks/dms_{tile}_analysis.ipynb'
    notebook: 'dms_tile_analysis.py.ipynb'
    