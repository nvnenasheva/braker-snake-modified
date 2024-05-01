# add a rule here that performs a BUSCO assessment on the braker.aa / annotation protein files
# look at how this was solved for building the braker commands in braker.smk / build_annotation_cmd.py
# BUSCO is available in the response container, you can directly activate the conda environment
# in the rule's shell command.