rule run_braker3:
    output:
        done = "data/checkpoints_annotate/braker3.done"
    singularity:
        "docker://teambraker/braker3:latest"
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        wd=braker3;
        braker.pl --genome=/opt/BRAKER/example/genome.fa \
            --prot_seq=/opt/BRAKER/example/proteins.fa --bam=/opt/BRAKER/example/RNAseq.bam \
            --workingdir=$wd \
	        --threads 8 --gm_max_intergenic 10000 --skipOptimize \
            --busco_lineage eukaryota_odb10 &> braker.log
        touch {output.done}
        """
