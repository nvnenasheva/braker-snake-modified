#!/usr/bin/env python3

import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser(description='Build annotation command')
parser.add_argument('-c', '--csv', type=str, required=True, help='csv file with all info about all species')
parser.add_argument('-s', '--species', type=str, required=True, help='Species name in data structure')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file name with annotation command')
parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use')
parser.add_argument('-w', '--workingdir', type=str, default='.', help='Working directory')
parser.add_argument('-l', '--logfile', type=str, default='braker.log', help='Log file name')
args = parser.parse_args()

# read csv file
df = pd.read_csv(args.csv)

# find the row where column "species" has the value of args.species
row = df.loc[df['species'] == args.species]
print(row)
cmd = ""
# if annotation_file is empty
if pd.isna(row['annotation_file'].values[0]):
    "I am in the building phase"
    cmd += 'braker.pl --genome=' + row['genome_file'].values[0] + ' '
    cmd += '--prot_seq=' + row['odb_file'].values[0]
    # if legacy_prot_file is not empty
    if not(pd.isna(row['legacy_prot_file'].values[0])):
        cmd += ',' + row['legacy_prot_file'].values[0]
    # if rnaseq_file is not empty
    if not(pd.isna(row['rnaseq_file'].values[0])):
        cmd += ' --bam=' + row['rnaseq_file'].values[0]
    cmd += ' --busco_lineage=' + row['busco_lineage'].values[0]
    cmd += ' --workingdir=' + args.workingdir
    cmd += ' --threads=' + str(args.threads)
    cmd += ' &> ' + args.logfile

# write command to output file
try:
    with open(args.output, 'w') as f:
        if len(cmd) > 0:
            f.write(cmd)
            f.write('\n')
        else:
            f.write('echo "No annotation command needed for species ' + str(args.species)+ '" &> ' + args.logfile)

except IOError:
    print("Error: can't write to file" + args.output)
    exit(1)

