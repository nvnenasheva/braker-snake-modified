rule run_braker1:
    output:
        done = "data/checkpoints_annotate/braker3.done"
    singularity:
        "docker://teambraker/braker3:latest"
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        wd=braker2;
        braker.pl --genome=/opt/BRAKER/example/genome.fa \
            --workingdir=$wd \
	        --threads 8 --gm_max_intergenic 10000 --skipOptimize \
            --busco_lineage eukaryota_odb10 &> braker.log
        touch {output.done}
        """