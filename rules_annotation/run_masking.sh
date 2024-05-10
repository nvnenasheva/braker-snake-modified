#!/bin/bash

#get Output folder
output_folder = $OUTPUT_FOLDER

wd=${{PWD}}
mkdir $output_folder
cp data/species/{params.spid}/genome/genome.fa $output_folder/genome.fa
cd $output_folder
log = $output_folder/{params_spid}_masking.log
touch $log
echo "BuildDatabase -name {params.spid} -dir data/species/{params.spid}/genome" &>> $log
BuildDatabase -name {params.spid} -pa {params.threads} genome.fa
echo "RepeatModeler -database {params.spid} -threads {params.threads} -LTRStruct" &>> $log
RepeatModeler -database {params.spid} -pa {params.threads} -LTRStruct
echo "RepeatMasker -threads {params.threads} -xsmall -lib {params.spid}-families.fa -dir data/species/{params.spid}/genome/" &>> $log
RepeatMasker -pa {params.threads} -xsmall -lib {params.spid}-families.fa genome.fa
cp genome.fa.masked data/species/{params.spid}/genome/genome.fa.masked
cd $wd
rm -rf $output_folder
touch {output}
