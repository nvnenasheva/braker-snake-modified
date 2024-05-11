#!/bin/bash

OPTSTRING=":s:t:f:o:"
spid=""
threads=""
output_folder=""
output=""

while getopts ${OPTSTRING} opt; do
	case ${opt} in
		s) spid="${OPTARG}";;
		t) threads="${OPTARG}";;
		f) output_folder="${OPTARG}";;
		o) output="${OPTARG}";;
		?) echo "Invalid option"; exit 1;;
	esac
done 

wd=${{PWD}}
mkdir ${output_folder}
cp data/species/${spid}/genome/genome.fa ${output_folder}/genome.fa
cd ${output_folder}
log = ${output_folder}/${spid}_masking.log
touch $log
echo "BuildDatabase -name {params.spid} -dir data/species/{params.spid}/genome" &>> $log
BuildDatabase -name ${spid} -pa ${threads} genome.fa
echo "RepeatModeler -database {params.spid} -threads {params.threads} -LTRStruct" &>> $log
RepeatModeler -database ${spid} -pa ${threads} -LTRStruct
echo "RepeatMasker -threads {params.threads} -xsmall -lib {params.spid}-families.fa -dir data/species/{params.spid}/genome/" &>> $log
RepeatMasker -pa ${threads} -xsmall -lib ${spid}-families.fa genome.fa
cp genome.fa.masked $wd/data/species/${spid}/genome/genome.fa.masked
cd $wd
rm -rf ${output_folder}
touch ${output}
