# We had orignally planned to use RED for masking but that has no data for our
# clade at all, and we do not know how to generate the required data.
# Therefore resorting to RepeatModeler2/RepeatMasker
# Implementation is fragile because if anything goes wrong, we do not see the log
# and the directory on /dev/shm is possibly not cleaned up.
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
        mkdir /dev/shm/dfam
        #cp data/species/{params.spid}/genome/genome.fa /dev/shm/dfam/genome.fa
        #cd /dev/shm/dfam
        #BuildDatabase -name {params.spid} -pa {params.threads} genome.fa
        #RepeatModeler -database {params.spid} -pa {params.threads} -LTRStruct
        #RepeatMasker -pa {params.threads} -xsmall -lib {params.spid}-families.fa genome.fa
        #cp genome.fa.masked data/species/{params.spid}/genome/genome.fa.masked
        #cd $wd
        rm -rf /dev/shm/dfam
        touch {output}
        """