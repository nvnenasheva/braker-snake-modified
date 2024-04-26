# braker-snake

Simple snakemake workflows for handling BRAKER on large data sets. Ultimately it should do this:

    1. Download the available assemblies for taxa from NCBI
    2. Prioritize in case of species duplications (first choice: annotated refseq, second choice: max N50)
    3. Separate into annotated and un-annotated genomes
    4. Download the repspective data sets from NCBI datasets
    5. Download OrthoDB partitions <- currently implemented up to here
    6. Check availability of RNA-Seq data for all downloaded genomes
    7. If less than 4 libraries, full download, otherwise run VARUS
    8. Run BRAKER3 on the un-annotated genomes with RNA-Seq
    9. Run BRAKER2 on the un-annotated genomes without RNA-Seq
    10. Run BUSCO on all the protein data sets and compile a summary

(For Clara's project, we do not need other steps, but the pipeline could serve as template for further expansion in the future.)

## Installation

```git clone https://github.com/KatharinaHoff/braker-snake.git```

Python dependencies:

```
pip install pandas
pip install snakemake
```

Bash dependencies (are usually available on a cluster):

   * singularity

   * curl

   * zcat

   * unzip

## Configuration

Contents of the config.ini file are important, adjust before running! Do not set the orthodb_path into your home directory, output is huge!

## Input data

The input data is expected to be in the following format:

```
<taxon> <odb_partition>
```

You find a small example in input.csv.  You may want to modify input.csv

## Running

Call with (currently only dry runs, but I tested that the singularity container is pulled)

```
cd braker-snake
snakemake --cores 1 --use-singularity
```

Beware: currently, specifying a taxon that includes another taxon in the input is dangerous! Don't do this!


Bash dependencies:

   * Singuarlity

   * curl

   * zcat

   * unzip

