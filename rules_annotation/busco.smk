# add a rule here that performs a BUSCO assessment on the braker.aa / annotation protein files
# look at how this was solved for building the braker commands in braker.smk / build_annotation_cmd.py
# BUSCO is available in the response container, you can directly activate the conda environment
# in the rule's shell command.
rule run_busco:
    input:
        input = "data/checkpoints_annotate/{spid}_braker.done"
    output:
        done = "data/checkpoints_annotate/{spid}_busco.done"
    params:
        spid="{spid}",
        threads = config['SLURM_ARGS']['cpus_per_task'],
        csv = config['INPUT']['species_csv']
    singularity:
        "docker://katharinahoff/response:devel"
    threads: int(config['SLURM_ARGS']['cpus_per_task'])
    resources:
        mem_mb=int(config['SLURM_ARGS']['mem_of_node']),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}";
        cmd_file=data/checkpoints_annotate/{params.spid}_busco.cmd
        log_file=data/checkpoints_annotate/{params.spid}_busco.log
        wdir=data/species/{params.spid}/busco
        python3 scripts/build_busco_cmd.py -c {params.csv} -s {params.spid} -t {params.threads} -o $cmd_file -b data/species/{params.spid}/braker -l $log_file -w $wdir
        bash $cmd_file 
        touch {output.done}
        """
