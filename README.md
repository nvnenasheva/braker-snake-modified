# braker-snake

Simple snakemake workflows for handling BRAKER on large data sets to prepare training data for DeepFinder. 

Ultimately it should do this:

    1. Download the available assemblies for taxa from NCBI
    2. Prioritize in case of species duplications (first choice: annotated refseq, second choice: max N50)
    3. Separate into annotated and un-annotated genomes
    4. Download the respective data sets from NCBI datasets (either genome only, or genome, annotation, proteins)
    5. Download OrthoDB partitions
    6. Check availability of RNA-Seq data for all downloaded genomes <- currently implemented up to here
    7. If less than N libraries, full download, alignment, sorting (N to be determined later)
    8. otherwise run VARUS
    9. Run BRAKER3 on the un-annotated genomes with RNA-Seq
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

Python dependencies (on BRAIN, you may first have to install conda, then install pip with conda, then install pandas with pip):

```
pip install pandas
```

Bash dependencies (are usually available on a cluster): singularity, curl, zcat, unzip, rm, echo, mkdir, ...

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

```
# only on BRAIN cluster:
module load singularity
cd braker-snake
singularity build snakemake.sif docker://snakemake/snakemake:latest
# once the image was built, only execute the following line:
singularity exec -B $PWD:$PWD snakemake.sif snakemake --cores 2 --use-singularity
```

