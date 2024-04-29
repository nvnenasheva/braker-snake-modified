# braker-snake

Simple snakemake workflows for handling BRAKER on large data sets to prepare training data for DeepFinder. 

Ultimately it should do this:

    1. Download the available assemblies for taxa from NCBI
    2. Prioritize in case of species duplications (first choice: annotated refseq, second choice: max N50)
    3. Separate into annotated and un-annotated genomes
    4. Download the respective data sets from NCBI datasets (either genome only, or genome, annotation, proteins)
    5. Download OrthoDB partitions
    6. Check availability of RNA-Seq data for all downloaded genomes
    7. If less than N libraries, full download, alignment, sorting (N to be determined later) <- currently implemented up to here
    8. otherwise run VARUS <- this is a problem because of the environment modification that we need to perform for VARUS, not sure how to do this in a container, yet
    9. Run BRAKER3 on the un-annotated genomes with RNA-Seq <- this will go into a separate Snakefile because it should launch one job per species, not per taxon, runtime issue otherwise
    10. Run BRAKER2 on the un-annotated genomes without RNA-Seq
    11. Run BUSCO on all the protein data sets and compile a summary

(For Clara's project, we do not need other steps, but the pipeline could serve as template for further expansion in the future.)

## Developing

First rule of development: **always commit and push your changes! Do not keep local unsynchronized changes for longer than 2 hours.**

If you want to develop on this pipeline, please at the moment, communicate with Katharina, Stepan, Clara. Many rules are still completely blank, please let the others know what rule file you are working on to avoid git conflicts.

If you want to fix a bug in an implemented set of rules, no need to communicate. Create a branch, fix the bug, make a pull request. If it has been tested properly, you can merge it yourself. Otherwise wait for Katharina to test the fix prior merging.

Once we have a running end-to-end pipeline, you can create branch for your own development. Please communicate with others to avoid conflicts down the line.

## Installation

Go do a directory where you have space for many GBs of data. Clone the repository:

```git clone https://github.com/KatharinaHoff/braker-snake.git```

Python dependencies (this workflow cannot execute with a singularity container of snakemake due to problems with SLURM communication):

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc # to activate miniconda
conda install -n base -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake
mamba activate snakemake
pip install snakemake-executor-plugin-slurms
```

Bash dependencies (are usually available on a cluster): singularity, curl, zcat, unzip, rm, echo, mkdir, cat, sed, ...

## Configuration

Contents of the config.ini file are important, adjust before running! Do not set the orthodb_path into your home directory, output is huge!

## Input data

The input data is expected to be in the following format:

```
<taxon> <odb_partition>
```

You find a small example in input.csv.  You may want to modify input.csv

Beware: currently, specifying a taxon that includes another taxon in the input is dangerous! Don't do this!

## Running

Execute this only in a place where you have space for many GBs of data, output is written to a subdirectory data and will be very, very big!

Due to a weird binding issue with current snakemake/SLURM (see Issues), you have to create a config file and fill it with the bindings directory for your singularity job. I hope this will be fixed, eventually.

```
mkdir -p ~/profile/apptainer/
touch ~/profile/apptainer/config.v8+.yaml
```

Add the following content to the file (adapt to your own working directory):

```
use-singularity: True
singularity-args: "\"--bind /home/hoffk83/git/braker-snake\""
```

Run the pipeline:

```
mamba activate snakemake
module load singularity
cd braker-snake
snakemake -s Snakefile_dataprep --executor slurm --default-resources slurm_account=none slurm_partition=batch --jobs=100 --use-apptainer
```

## Current DAG with example data

(can always be generated with `snakemake -s Snakefile_dataprep --dag | dot -Tpng > dag.png`, potentially not on BRAIN because dot may not be installed)

![DAG](dag.png)