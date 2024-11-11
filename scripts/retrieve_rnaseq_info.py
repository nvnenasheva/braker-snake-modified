#!/usr/bin/env python3

import argparse
import logging
from Bio import Entrez
import pandas as pd
# generate a seed for consistent random sampling
import random
random.seed(42)

__author__ = "Katharina J. Hoff"

"""
    Search the NCBI SRA (Sequence Read Archive) for RNA libraries of a specific species that use the Illumina platform and paired-end layout. 
    
    This function constructs a query to search for transcriptomic data from the specified species that was sequenced on Illumina equipment with paired-end reads. It retrieves the count of available datasets and the accession numbers of each dataset.
    
    Parameters:
    species_name (str): The scientific name of the species to query. This name must be recognized by NCBI.
    email (str): Email address used to identify the user making requests to NCBI's servers. This is a requirement by NCBI to use their services.
    
    Returns:
    tuple:
        int: The count of RNA libraries found matching the query.
        list of str: A list of accession numbers for the RNA libraries found, if any.
    
    Usage:
    count, accessions = search_sra("Thalassiosira pseudonana", "your.email@example.com")

    Note:
    The function fetches detailed run information, which is then processed to extract only the accession numbers, omitting unnecessary details. This includes handling large textual data from NCBI's output and parsing it accordingly.
    
    The email provided should be a valid email, as NCBI might use it to contact in case of excessive usage or other issues.
    """
def search_sra(species_name, email):
    Entrez.email = email  # Set the email for NCBI access

    # Define the search query, we are on purpose not searching for field species name because a lot of diatom data comes from mixed cultures with bacteria
    if '_sp_' in species_name:
        print(species_name)
        query = treat_sp(species_name)
    else:
        query = f'"{species_name}" AND "transcriptomic"[Source] AND "Illumina"[Platform] AND "PAIRED"[Layout] AND PolyA'

    # Use Entrez.esearch to search SRA
    handle = Entrez.esearch(db="sra", term=query)
    record = Entrez.read(handle)
    handle.close()

    # Retrieve detailed information about each run
    ids = record['IdList']
    run_info = []
    if ids:
        handle = Entrez.efetch(db="sra", id=','.join(ids), rettype="runinfo", retmode="text")
        details = handle.read()
        handle.close()
        # the details are a byte object, convert to string:
        details = details.decode('utf-8')
        # parse the details to filter for the accession numbers of all libraries:
        # this data structure contains info on one library per line, the lines themselves are comma-separated, we need
        # only the first field of each line, we can skip the header
        details = details.split('\n')
        details = [line.split(',')[0] for line in details[1:]]
        run_info.append(details)
    
    if record['Count'] == '0':
        run_info = []
    else:
        run_info = run_info[0]
        # remove empty list entries:
        run_info = [x for x in run_info if x]
    return record['Count'], run_info

def treat_sp(string):
    if "_sp_" in string:
        # Split the string at '_sp_' and take the part after it
        strain_extension = string.split("_sp_")[-1].strip()
        logging.debug(f"Extracted strain extension '{strain_extension}' from '{string}'")
        string = strain_extension
    return string

def main(args):
    # Read species from args.species_table with pandas to get species names
    species_table = pd.read_csv(args.species_tables, sep='\t')
    logging.info(f"Read species table from '{args.species_tables}'")
    species_list = species_table['species'].tolist()
    logging.debug(f"Original species list: {species_list}")
    species_list = treat_sp(species_list)
    print(species_list)
    logging.debug(f"Modified species list: {species_list}")

    # Search SRA for each species
    all_data = {}
    for species in species_list:
        nRecords, accessions_list = search_sra(species, args.email)
        # randomly sort the accessions_list contents
        random.shuffle(accessions_list)
        ## if we want to return to VARUS, we need to delete the next two lines
        # store only the first n accessions in accessions_list
        if len(accessions_list) > int(args.n_threshold):
            accessions_list = accessions_list[:int(args.n_threshold)]
            logging.debug(f"Trimmed accessions list for '{species}' to {args.n_threshold} entries")
        all_data[species] = {'nRecords': nRecords, 'accessions': accessions_list}


    # Print species with less than n_threshold records, if they have more than 0
    # including the accession numbers for direct fastq-dump processing
    try:
        with open(args.fastqdump_out_list, 'w') as fastq_list:
            for species in species_list:
                nRecords = all_data[species]['nRecords']
                accessions_list = all_data[species]['accessions']
                if int(nRecords) != 0:
                    fastq_list.write(f"{species}\t")
                    for acc in accessions_list:
                        if acc != accessions_list[-1]:
                            fastq_list.write(f"{acc},")
                        else:
                            fastq_list.write(f"{acc}\n")
                else:
                    logging.warning(f"No RNA-Seq data for species '{species}' in SRA")
                    #print("Info: There is no RNA-Seq data for species " + species + " in SRA")

    except IOError:
        print("Could not write to file " + args.fastqdump_out_list)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Search NCBI SRA for RNA libraries of specified species.')
    parser.add_argument('-e', '--email', required=True, help='Email address for NCBI access')
    parser.add_argument('-t', '--species_tables', required=True, help='File path to tab separated file that contains species information, format specific to braker-snake')
    parser.add_argument('-f', '--fastqdump_out_list', required=True, help='List of species ')
    parser.add_argument('-n', '--n_threshold', required=False, default=6, help='Threshold for number of records to be considered for further processing with VARUS (default: 6), those with fewer records will be output into a separate file for direct fastq-dump processing')

    args = parser.parse_args()
    main(args)