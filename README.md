# Deep mutational scanning of ZIKV NS3 protein
Experiments by Blake Richardson and Matt Evans.
Analysis by Caroline Kikawa, adapting a pipeline by Jesse Bloom and David Bacsik.

## Analysis pipeline overview
This `Snakemake` takes input deep sequencing reads and processes them using the `dms_tools2` package written by the Bloom lab. For this project, a deep mutational scanning library was created in discrete 'tiles' across the NS2B/NS3 protein from Zika virus. A pool of virus particles expressing these variants was selected by passaging on cells. The pre-selection counts of each variant (e.g., their frequency in the plasmid library) were compared to their counts in the cell-passaging selected libraries. The *amino acid preferences* or `prefs` are the result of this comparison. We then calculate *mutational effects* or `muteffects` for each variant by taking the log ratio of the variant mutation versus wild-type.

## Results
For a summary of the results, see [results/summary/](results/summary/), which has Markdown summaries for the analysis of each tile (e.g., [results/summary/dms_tile_1_analysis.md](results/summary/dms_tile_1_analysis.md), etc).

Other results are placed in [./results/](results), although not all files are tracked in the GitHub repo. Again, these files are sub-divided by tile and analysis (e.g., [results/tile_1/prefs](results/tile_1/prefs), [results/tile_1/muteffects](results/tile_1/muteffects), etc).

## Running analysis
First activate the *ZIKV_DMS_NS5_EvansLab* [conda](https://docs.conda.io/projects/conda/en/latest/index.html) environment for the analysis.
If you have not already created this environment, build it from [environment.yml](ZIKV_DMS_NS3_EvansLab) with:

    conda env create -f environment.yml

Then activate the environment with:

    conda activate ZIKV_DMS_NS3_EvansLab

The analysis is run by the [snakemake](https://snakemake.readthedocs.io/) pipeline in [Snakefile](Snakefile).
Essentially, this pipeline runs the Jupyter notebook [dms_tile_analysis.ipynb](dms_tile_analysis.ipynb) for each deep mutational scanning tile, with the tile information specified in [config.yml](config.yml).
To run the pipeline using 36 jobs, use the command:

    snakemake -j 36

Add the `--keep-incomplete` flag if you don't want to delete results on an error.
To run on the Hutch cluster using `slurm`, do:

    sbatch -c 36 run_Snakemake.bash


## Input data
The input data are in [./data/](data):

 - `./data/tile_*_amplicon.fasta`: amplicons for each tile of the barcoded-subamplicon sequencing.

 - `./data/tile_*_subamplicon_alignspecs.txt`: the alignment specs for the [barcoded subamplicon sequencing](https://jbloomlab.github.io/dms_tools2/bcsubamp.html) for each amplicon.

 - `./data/tile_*_samplelist.csv`: all the samples that we sequenced and the locations of the associated deep-sequencing data for each amplicon.


