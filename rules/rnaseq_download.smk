# Rule: retrieve_rnaseq_info_from_sra
# Purpose:
#   This rule automates the retrieval of RNA sequencing (RNA-seq) information for specified taxa
#   from the NCBI Sequence Read Archive (SRA). It utilizes a Python script to query NCBI for RNA-seq
#   data related to unannotated species listed in a provided table. The rule processes each taxon
#   listed in the 'data/{taxon}_blank.tbl', extracting necessary information for downstream RNA-seq
#   analysis pipelines.
#
# Inputs:
#   - download_script: A Python script ('scripts/retrieve_rnaseq_info.py') that handles the querying
#     of NCBI databases and prepares output lists of SRA accession numbers and other relevant data.
#   - unannotated_species: A table listing species without annotated genomes; used to determine which
#     species' RNA-seq data to retrieve.
#
# Outputs:
#   - fastqdump_lst: A list file ('data/{taxon}_rnaseq_for_fastqdump.lst') containing accession numbers
#     and other parameters formatted for the `fastq-dump` utility, which can be used to download sequence data.
#   - varus_list: A list file ('data/{taxon}_rnaseq_for_varus.lst') used for VARUS, a pipeline that automatically
#     optimizes parameters for RNA-seq data retrieval and assembly.
#   - done: A simple checkpoint file ('data/{taxon}_rnaseq_info.done') indicating successful completion of data retrieval.
#
# Parameters:
#   - taxon: A wildcard parameter defining the specific taxon being processed. Extracted from the filename pattern.
#   - email: An email address registered with NCBI, required for making API requests to NCBI databases.
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


# Rule: download_fastq
# Warning: 
#   This rule takes extremely long! It is not run in SLURM parallel because i/o may be a bottleneck.
#   It must regrettably expected that this rule fails, repeatedly. Do not worry, just restart the
#   workflow. The rule is implemented in a way that it will pick up where it left off. 
# To Do:
#   The failures are due to the fact that the NCBI servers are not always available.
#   There's a tweet by Heng Li that may help to eventually debug this:
#   https://twitter.com/lh3lh3/status/1779876367200387172
#   This suggests to replace fastq-dump by fasterq-dump. This is not yet implemented in the rule.
#   It suggest to first download a prefetch file, and then execute the dump.
# Purpose:
#   This rule is designed to process a list of SRA accession numbers and use `fastq-dump` to download
#   paired-end fastq files for each accession. It handles the preparation of directory structures
#   and the organization of downloaded files within those directories. This is crucial for downstream
#   analysis such as RNA-seq data processing or genomic assemblies.
#
# Inputs:
#   - fastqdump_lst: A file ('data/{taxon}_rnaseq_for_fastqdump.lst') that lists species and their
#     corresponding SRA accession numbers. Each line contains a species name and a comma-separated list
#     of SRA IDs.
#
# Outputs:
#   - done: A file ('data/{taxon}_fastqdump.done') that signals the successful completion of the downloads
#     and processing of all SRA data listed in the input file.
#
# Operations:
#   - The rule first sets up the environment to ensure correct directory binding when using Singularity.
#   - It reads each line from the input file, processes the species name to replace spaces with underscores
#     (for file naming consistency), and creates a directory for storing the downloaded fastq files.
#   - For each SRA ID, `fastq-dump` is executed to retrieve paired-end reads, and the resultant files are
#     immediately compressed using gzip to save space.
#   - A log file is generated to capture the output and errors of the download process, helping in troubleshooting
#     and ensuring transparency of the operation.
#   - Upon successful execution of all commands, a 'done' file is created as a checkpoint indicating
#     the completion of the task.
rule download_fastq:
    input:
        fastqdump_lst = "data/{taxon}_rnaseq_for_fastqdump.lst",
        genome_done = "data/{taxon}_download.done" # This is a dummy file to ensure that the download rule is executed after the retrieval rule
    output:
        done = "data/{taxon}_fastqdump.done"
    singularity:
        "docker://teambraker/braker3:latest"
    shell:
        """
        echo "I am done with the SIF file, now starting, will take very long..."
        export APPTAINER_BIND="${{PWD}}:${{PWD}}";
        logfile=$PWD/data/{wildcards.taxon}_fastqdump.log
        # Read the input file and process it
        while IFS=$'\\t' read -r species sra_ids; do
            echo "$PWD" &>> $logfile
            species_fixed=$(echo "$species" | sed 's/ /_/g')  # Replace space with underscore
            echo "cd data/species/$species_fixed" &>> $logfile
            cd data/species/$species_fixed
            # if the directory fastq does not exist, yet
            if [ ! -d "fastq" ]; then
                echo "mkdir fastq" &>> $logfile
                mkdir fastq  # Create a directory for the species
            fi
            echo "cd ../../../" &>> $logfile
            cd ../../../ 
            echo "$PWD" &>> $logfile
            # Convert comma-separated string to array
            IFS=',' read -ra ids <<< "$sra_ids"

            # Process each SRA ID locally using the array
            for id in "${{ids[@]}}"; do
                # this is such a long and expensive process that we do not want it to execute if fastq.gz files already exist
                if [ ! -f "data/species/$species_fixed/fastq/${{id}}_1.fastq.gz" ] && [ ! -f "data/species/$species_fixed/fastq/${{id}}_2.fastq.gz" ]; then
                    echo "fastq-dump --split-files --outdir data/species/$species_fixed/fastq $id" &>> $logfile
                    fastq-dump --split-files --outdir data/species/$species_fixed/fastq $id &>> $logfile
                    echo "gzip data/species/$species_fixed/fastq/${{id}}_1.fastq" &>> $logfile
                    gzip data/species/$species_fixed/fastq/${{id}}_1.fastq &>> $logfile
                    gzip data/species/$species_fixed/fastq/${{id}}_2.fastq
                fi
            done
        done < {input.fastqdump_lst} &>> $logfile
        
        echo "touch {output.done}" &>> $logfile
        touch {output.done}
        """

