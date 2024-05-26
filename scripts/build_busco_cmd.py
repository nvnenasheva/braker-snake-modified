#!/usr/bin/env python3

import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser(description='Build busco command')
parser.add_argument('-c', '--csv', type=str, required=True, help='csv file with all info about all species')
parser.add_argument('-s', '--species', type=str, required=True, help='Species name in data structure')
parser.add_argument('-o', '--output', type=str, required=True, help='Output file name with annotation command')
parser.add_argument('-b', '--braker', type=str, required=True, help='Output Directory of Braker Run')
parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use')
parser.add_argument('-w', '--workingdir', type=str, default='.', help='Working directory')
parser.add_argument('-l', '--logfile', type=str, default='braker.log', help='Log file name')
args = parser.parse_args()

# read csv file
df = pd.read_csv(args.csv)

# find the row where column "species" has the value of args.species
row = df.loc[df['species'] == args.species]
cmd = ""
cmd += 'busco --in=' + str(args.braker)
cmd += ' --mode=genome'
cmd += ' --lineage_dataset=' + row['busco_lineage'].values[0]
cmd += ' --cpu=' + str(args.threads)
cmd += ' --out=busco_' + str(args.species)
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
