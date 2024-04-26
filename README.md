# braker-snake
Simple snakemake workflows for handling BRAKER on large data sets

Contents of the config.ini file are important, adjust before running! Do not set the orthodb_path into your home directory, output is huge!

Call with (currently only dry runs, but I tested that the singularity container is pulled)

```
snakemake --cores 1 --use-singularity
```

Beware: currently, specifying a taxon that includes another taxon in the input is dangerous! Don't do this!


Bash dependencies:

   * Singuarlity

   * curl

   * zcat

   * unzip

