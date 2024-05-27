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
    threads: int(config['SLURM_ARGS']['cpus_per_task'])
    resources:
        mem_mb=int(config['SLURM_ARGS']['mem_of_node']),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}";
        wd=${{PWD}}
        log=data/checkpoints_annotate/{params.spid}_runningbusco.log
        touch $log
        echo "$wd/scripts/run_busco.sh -s {params.spid} -t {params.threads} -c {params.csv} -o {output}" &>> $log
        $wd/scripts/run_busco.sh -s {params.spid} -t {params.threads} -c {params.csv} -o {output}
        """
