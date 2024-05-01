rule run_braker:
    output:
        done = "data/checkpoints_annotate/{id}_braker.done"
    params:
        id="{id}",
        threads = config['SLURM_ARGS']['cpus_per_task'],
        csv = config['INPUT']['species_csv']
    singularity:
        "docker://teambraker/braker3:latest"
    threads: int(config['SLURM_ARGS']['cpus_per_task'])
    resources:
        mem_mb=int(config['SLURM_ARGS']['mem_of_node']),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}";
        cmd_file=data/checkpoints_annotate/{params.id}_braker.cmd
        log_file=data/checkpoints_annotate/{params.id}_braker.log
        wdir=data/species/{params.id}/braker
        python3 scripts/build_annotation_cmd.py -c {params.csv} -s {params.id} -t {params.threads} -o $cmd_file -l $log_file -w $wdir
        bash $cmd_file
        touch {output.done}
        """