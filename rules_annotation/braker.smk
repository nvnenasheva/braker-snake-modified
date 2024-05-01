rule build_repeats:
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
        mkdir /dev/shm/dfam
        #cp data/species/{params.spid}/genome/genome.fa /dev/shm/dfam/genome.fa
        #cd /dev/shm/dfam
        #BuildDatabase -name {params.spid} -pa {params.threads} genome.fa
        #RepeatModeler -database {params.id} -pa {params.threads} -LTRStruct
        #RepeatMasker -pa {params.threads} -xsmall -lib {params.spid}-families.fa genome.fa
        #cp genome.fa.masked data/species/{params.spid}/genome/genome.fa.masked
        #cd $wd
        rm -rf /dev/shm/dfam
        touch {output}
        """

# This rule determines the BRAKER run mode from the input
# file species_csv, automatically. It runs either BRAKER3,
# or BRAKER2 (if present including legacy proteins from 
# a previously annotated assembly of that species).
# It is not running BRAKER1, we hope that this will not be
# necessary / i.e. will resort to BRAKER2 if BRAKER1 fails.
# Failure is not caught/resumed.
rule run_braker:
    input:
        input = "data/checkpoints_annotate/{spid}_repeats.done"
    output:
        done = "data/checkpoints_annotate/{spid}_braker.done"
    params:
        spid="{spid}",
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
        cmd_file=data/checkpoints_annotate/{params.spid}_braker.cmd
        log_file=data/checkpoints_annotate/{params.spid}_braker.log
        wdir=data/species/{params.spid}/braker
        python3 scripts/build_annotation_cmd.py -c {params.csv} -s {params.spid} -t {params.threads} -o $cmd_file -l $log_file -w $wdir
        #bash $cmd_file
        touch {output.done}
        """