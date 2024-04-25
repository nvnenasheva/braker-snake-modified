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

# Load and parse the config file
config = configparser.ConfigParser()
config.read('config.ini')
print(config.sections())
input_csv = config['INPUT']['input_csv']
print(input_csv)

# Read the input CSV file to get clade & odb partition names
taxa_list = pd.read_csv(input_csv, header=None, names=['taxa'])
print(taxa_list)

# Include other rule files (assuming they define their own targets without using wildcards inappropriately)
include: "rules/genome_download.smk"

# Main rule to process each taxon
rule all:
    input:
        expand("data/{taxon}_processed.json", taxon=taxa_list['taxa'])

