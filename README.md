# Deep mutational scanning of ZIKV NS2B/NS3 protease

Experiments by **Blake Richardson** and **Matt Evans**
Analysis by **Caroline Kikawa** and **Will Hannon**

## Analysis

A deep mutational scanning library was created in 3 discrete non-overlapping 'tiles' across the ZIKV NS2B/NS3 protease. The resulting pool of virus particles expressing these variants was selected by passaging on cells. The counts of each variant in the plasmid library were compared to their counts after passaging using a pipeline adapted from [a `Snakemake` pipeline]() written by Jesse Bloom and David Bacsik. The pipeline takes deep sequencing reads as input and processes them using [`dms_tools2`](https://jbloomlab.github.io/dms_tools2/). The **amino acid preferences** or `prefs` are calculated based on count of variants pre- and post- selection by passaging on cells. We then calculate **mutational effects** or `muteffects` for each variant by taking the log ratio of the variant mutation versus wild-type.

## Results

All results are located in [results/](results), although large files are not tracked in the GitHub repo. Results files are sub-divided by tile and analysis (e.g., [results/tile_1/prefs](results/tile_1/prefs), [results/tile_1/muteffects](results/tile_1/muteffects), etc).

See [results/summary/](results/summary/) for markdown summaries of the analysis for each tile (e.g., [results/summary/dms_tile_1_analysis.md](results/summary/dms_tile_1_analysis.md), etc).

## Running analysis

First activate the *ZIKV_DMS_NS3_EvansLab* [conda](https://docs.conda.io/projects/conda/en/latest/index.html) environment for the analysis.
If you have not already created this environment, build it from [environment.yml](environment.yml) with:

```bash
conda env create -f environment.yml
```
Then activate the environment with:

```bash
conda activate ZIKV_DMS_NS3_EvansLab
```

The analysis is run by the [`snakemake`](https://snakemake.readthedocs.io/) pipeline in [Snakefile](Snakefile). Essentially, this pipeline runs the Jupyter notebook [dms_tile_analysis.ipynb](dms_tile_analysis.ipynb) for each deep mutational scanning tile, with the tile information specified in [config.yml](config.yml).

To run the pipeline using 36 jobs, use the command:

```bash
snakemake --keep-incomplete -j 36 
```
To run on the Hutch cluster using `slurm`, do:

```bash
sbatch -c 36 run_analysis.bash
```

## Input data

The input data are in [data/](data):

 - `data/tile_*_amplicon.fasta`: amplicons for each tile of the barcoded-subamplicon sequencing.

 - `data/tile_*_subamplicon_alignspecs.txt`: the alignment specs for the [barcoded subamplicon sequencing](https://jbloomlab.github.io/dms_tools2/bcsubamp.html) for each amplicon.

 - `data/tile_*_samplelist.csv`: all the samples that we sequenced and the locations of the associated deep-sequencing data for each amplicon.

