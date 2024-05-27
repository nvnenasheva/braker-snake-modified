#!/bin/bash

OPTSTRING=":s:t:c:o:"
spid=""
threads=""
csv=""
output=""

while getopts ${OPTSTRING} opt; do
	case ${opt} in
		s) spid="${OPTARG}";;
		t) threads="${OPTARG}";;
		f) csv="${OPTARG}";;
		o) output="${OPTARG}";;
		?) echo "Invalid option"; exit 1;;
	esac
done 

conda activate busco_env
cmd_file=data/checkpoints_annotate/${spid}_busco.cmd
log_file=data/checkpoints_annotate/${spid}_busco.log
wdir=data/species/${spid}/busco
python3 scripts/build_busco_cmd.py -c ${csv} -s ${spid} -t ${threads} -o $cmd_file -b data/species/${spid}/braker/braker.aa -l $log_file -w $wdir
bash $cmd_file 
touch ${output}
