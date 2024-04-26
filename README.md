# braker-snake
Simple snakemake workflows for handling BRAKER on large data sets

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

