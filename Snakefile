"""
braker-snake

Katharina J. Hoff, Stepan Saenko, Clara Pitzschel

University of Greifswald

A joint effort to automated bulk genome annotation
"""

__author__ = "Katharina J. Hoff, Stepan Saenko, Clara Pitzschel"

import configparser
import pandas as pd
from pathlib import Path
import json
import re

# Load and parse the config file
config = configparser.ConfigParser()
config.read('config.ini')
input_csv = config['INPUT']['input_csv']

# Read the input CSV file to get taxon and odb partition names
# Assuming the CSV file is tab-separated
data = pd.read_csv(input_csv, header=None, sep=' ', names=['taxa', 'odb_partition'])

# Create separate lists for taxa and unique odb partitions
taxa_list = data['taxa'].tolist()
unique_odb_partitions = data['odb_partition'].unique().tolist()

# Include other rule files (assuming they define their own targets without using wildcards inappropriately)
include: "rules/genome_download.smk"
include: "rules/odb_download.smk"
include: "rules/rnaseq_download.smk"

# Main rule to process each taxon
rule all:
    input:
        expand("data/{taxon}_download.done", taxon=taxa_list),
        expand(config['BRAKER']['orthodb_path'] + "/{odb_partition}.fa", odb_partition=unique_odb_partitions),
        expand("data/{taxon}_rnaseq_info.done", taxon=taxa_list),
        expand("data/{taxon}_hisat2_index.done", taxon=taxa_list)
        # This is the place where you have to expand when you are waiting for more targets! For example,
        # we will implement RNA-Seq download, and this is where you have to add the RNA-Seq targets.

