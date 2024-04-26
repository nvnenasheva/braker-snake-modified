rule retrieve_rnaseq_info_from_sra:
    input:
        download_script = "scripts/retrieve_rnaseq_info.py",
        unannotated_species = "data/{taxon}_blank.tbl"
    output:
        fastqdump_lst = "data/{taxon}_rnaseq_for_fastqdump.lst",
        varus_list = "data/{taxon}_rnaseq_for_varus.lst",
        done = "data/{taxon}_rnaseq_info.done"
    params:
        taxon = lambda wildcards: wildcards.taxon,
        email = config['ESEARCH']['email']
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://teambraker/braker3:latest"
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        python3 {input.download_script} -e {params.email} -t {input.unannotated_species} -l {output.varus_list} -f {output.fastqdump_lst}; \
        touch {output.done}
        """

rule download_fastq:
    input:
        fastqdump_lst = "data/{taxon}_rnaseq_for_fastqdump.lst"
    output:
        done = "data/{taxon}_fastqdump.done"
    singularity:
        "docker://teambraker/braker3:latest"
    shell:
        """
        # Ensure the PWD is bound correctly if using singularity
        export APPTAINER_BIND="${{PWD}}:${{PWD}}";

        # Read the input file and process it
        while IFS=$'\\t' read -r species sra_ids; do
            species_fixed=$(echo "$species" | sed 's/ /_/g')  # Replace space with underscore
            mkdir -p data/$species_fixed/fastq  # Create a directory for the species

            # Convert comma-separated string to array
            IFS=',' read -ra ids <<< "$sra_ids"

            # Process each SRA ID locally using the array
            for id in "${{ids[@]}}"; do
                fastq-dump --split-files --outdir data/$species_fixed/fastq $id
                gzip data/$species_fixed/fastq/${{id}}_1.fastq
                gzip data/$species_fixed/fastq/${{id}}_2.fastq
            done
        done < {input.fastqdump_lst} &> data/{wildcards.taxon}_fastqdump.log
        
        touch {output.done}
        """

