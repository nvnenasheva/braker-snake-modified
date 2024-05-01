rule run_braker2:
    input:
        csv = config['INPUT']['species_csv']
    output:
        done = "data/checkpoints_annotate/{id}_braker2.done"
    params:
        id="{id}"  
    singularity:
        "docker://teambraker/braker3:latest"
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        echo {params.id}


        #wd=braker2;
        #braker.pl --genome=/opt/BRAKER/example/genome.fa \
        #    --prot_seq=/opt/BRAKER/example/proteins.fa \
        #    --workingdir=$wd \
	    #    --threads 8 --gm_max_intergenic 10000 --skipOptimize \
        #    --busco_lineage eukaryota_odb10 &> braker.log
        touch {output.done}
        """