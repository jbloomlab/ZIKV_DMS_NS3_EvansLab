"""
Pipeline that runs the ZIKV NS2B/NS3 analysis for each tile.
Authors: Caroline Kikawa, Will Hannon, and David Bascik
"""

#### ----------------------- Imports ----------------------- ####

import os
import pandas as pd

#### -------------------- Configuration -------------------- ####

configfile: 'config.yml'

#### ----------------------- Targets ----------------------- ####

wildcard_constraints:
    tile="tile_\d+"

rule all:
    input:
        "results/summary/all_tiles_effects_and_preferences.csv",
        "results/summary/all_tiles_effects_and_preferences_with_stops.csv",
        expand("results/summary/dms_{tile}_analysis.md",
               tile=config['tiles']),
        expand("results/summary/dms_{tile}_analysis.html",
            tile=config['tiles']),

#### ------------------------ Rules ------------------------ ####

rule clean:
    shell:
        """
        rm -rf logs/
        rm -rf tmp/
        rm -f slurm*.out
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


rule combine_all_tiles:
    input:  expand("results/{tile}", tile=config['tiles'])
    output: without_stops_csv = "results/summary/all_tiles_effects_and_preferences.csv",
            with_stops_csv = "results/summary/all_tiles_effects_and_preferences_with_stops.csv",
    params: tiles = config['tiles']
    conda: 'environment.yml'
    log: notebook='results/notebooks/combine_all_tiles.ipynb'
    notebook: 'combine_all_tiles.ipynb'


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
    
rule jupnb_to_html:
    """Convert Jupyter notebook to HTML format."""
    input: notebook="results/notebooks/{notebook}.ipynb"
    output: html="results/summary/{notebook}.html"
    params: outdir=lambda wildcards, output: os.path.dirname(output.html)
    conda: 'environment.yml'
    shell: 
        """
        jupyter nbconvert \
            --output-dir {params.outdir} \
            --to html \
            {input.notebook}
        """
