"""
braker-snake

Katharina J. Hoff, Stepan Saenko, Clara Pitzschel

University of Greifswald

A joint effort to automated bulk genome annotation
"""

__author__ = "Katharina J. Hoff, Stepan Saenko, Clara Pitzschel"

import configparser
import pandas as pd

# Load and parse the config file
config = configparser.ConfigParser()
config.read('config.ini')
input_csv = config['INPUT']['input_csv']

# Read the input CSV file to get species
species_list = pd.read_csv(input_csv, header=None, names=['species'])
species_dirs = [s.replace(' ', '_') for s in species_list['species']]

# Generate a complete list of expected output files
expected_outputs = [f"{species_dir}/genome/genome.fa" for species_dir in species_dirs]
print(expected_outputs)

# Include other rule files (assuming they define their own targets without using wildcards inappropriately)
include: "rules/genome_download.smk"

# Main rule to ensure all directories and initial files are created
rule all:
    input: expected_outputs
