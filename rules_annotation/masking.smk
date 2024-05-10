# We had orignally planned to use RED for masking but that has no data for our
# clade at all, and we do not know how to generate the required data.
# It is really bad that we cannot change directories :-(
# This rule requires a cleanup rule, later, to remove all the
# weird species-specific files that are generated in the current directory.
# These include folders with a prefix RM_ and files with a prefix {spid}.
# We also want to move, not delete, the species specific repeat library files.
rule mask_repeats:
    output:
        output = "data/checkpoints_annotate/{spid}_repeats.done"
    params:
        spid="{spid}",
        threads = config['SLURM_ARGS']['cpus_per_task']
    singularity:
        "docker://dfam/tetools:latest"
    threads: int(config['SLURM_ARGS']['cpus_per_task'])
    resources:
        mem_mb=int(config['SLURM_ARGS']['mem_of_node']),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"
	wd=${{PWD}}
	export OUTPUT_FOLDER={{config['TARGET']['output_folder']}}
        log=data/checkpoints_annotate/{params.spid}_repeats.log
        touch $log
        echo "$wd/run_masking.sh" &>> $log
	$wd/run_masking.sh
        touch {output}
        """
