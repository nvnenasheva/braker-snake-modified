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

rule write_hisat_index_script:
    input:
        fastqdump_lst = "data/{taxon}_rnaseq_for_fastqdump.lst",
        download_done = "data/{taxon}_fastqdump.done",
        genome_done = "data/{taxon}_download.done"
    output:
        done = "data/{taxon}_hisat2_index_scripts.done"
    params:
        taxon = lambda wildcards: wildcards.taxon,
        threads = config['SLURM_ARGS']['cpus_per_task']
    wildcard_constraints:
        taxon="[^_]+"
    run:
        try:
            with open(input.fastqdump_lst, "r") as f:
                lines = f.readlines()
        except IOError:
            raise FileNotFoundError("Could not read the input file " + input.fastqdump_lst)
        # construct commands, one index per species will be needed
        cmds = []
        for line in lines:
            species, sra_ids = line.strip().split("\t")
            species_fixed = species.replace(" ", "_")
            # if species folder for scripts does not exist, create it
            script_dir = f"data/species/{species_fixed}/scripts_hisat2"
            if not os.path.exists(script_dir):
                os.makedirs(script_dir, exist_ok=True)
            # open script file, this does not need a random name since it is exactly one per species
            script_file = f"{script_dir}/hisat2_index.sh"
            try:
                with open(script_file, "w") as f:
                    # construct a genome index
                    f.write(f"hisat2-build -p {params.threads} data/species/{species_fixed}/genome.fa data/species/{species_fixed}/genome/genome.fa.idx\n")
            except IOError:
                raise FileNotFoundError("Could not write the script file " + script_file)
        # create the file that signals that the scripts are written
        try:
            with open(output.done, "w") as f:
                f.write("done\n")
        except IOError:
            raise FileNotFoundError("Could not write the output file " + output.done)


rule run_hisat2_index:
    input:
        fastqdump_lst = "data/{taxon}_rnaseq_for_fastqdump.lst",
        download_done = "data/{taxon}_hisat2_index_scripts.done",
        genome_done = "data/{taxon}_download.done"
    output:
        done = "data/{taxon}_hisat2_index.done"
    params:
        taxon=lambda wildcards: wildcards.taxon,
        singularity_image="docker://teambraker/braker3:latest",
        threads = config['SLURM_ARGS']['cpus_per_task']
    threads: 48
    resources:
        mem_mb=int(config['SLURM_ARGS']['mem_of_node']),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        log="data/{wildcards.taxon}_hisat2_index.log"
        echo "" > $log
        readarray -t species_list < <(cat {input.fastqdump_lst} | sed 's/ /_/' | cut -d$'\t' -f1)
        for species in "${{species_list[@]}}"; do
            if [ ! -d "data/species/$species/hisat2" ]; then
                mkdir -p "data/species/$species/hisat2"
            fi
            echo "hisat2-build -p {params.threads} data/species/$species/genome/genome.fa data/species/$species/genome/genome.fa.idx" >> $log
            which hisat2-build &>> $log
            hisat2-build -p {params.threads} data/species/$species/genome/genome.fa data/species/$species/genome/genome.fa.idx &>> $log
        done
        touch {output.done}
        """

'''
rule run_hisat2_index:
    input:
        fastqdump_lst = "data/{taxon}_rnaseq_for_fastqdump.lst",
        download_done = "data/{taxon}_hisat2_index_scripts.done",
        genome_done = "data/{taxon}_download.done"
    output:
        done = "data/{taxon}_hisat2_index.done"
    params:
        taxon=lambda wildcards: wildcards.taxon,
        singularity_image="docker://teambraker/braker3:latest",
        partition = config['SLURM_ARGS']['partition'],
        threads = config['SLURM_ARGS']['cpus_per_task']
    wildcard_constraints:
        taxon="[^_]+"
    shell:
        """
        while IFS=$'\t' read -r species sra_ids; do
            species_fixed=$(echo "$species" | sed 's/ /_/g')  # Replace space with underscore
            echo $species_fixed
            slurm_script="data/species/$species_fixed/scripts_hisat2/run_hisat2_index.sh"
            # Create the SLURM script with proper quoting
            echo '#!/bin/bash' > "$slurm_script"                             # Use single quotes here
            echo 'module load singularity' >> "$slurm_script"
            # load singularity bindings
            echo "export APPTAINER_BIND=\"${{PWD}}:${{PWD}}\"" >> "$slurm_script" # Properly quote variables
            # execute the actual index script
            echo "singularity exec docker://teambraker/braker3:latest bash data/species/$species_fixed/scripts_hisat2/hisat2_index.sh" >> "$slurm_script"
            # submit this job with sbatch
            sbatch --partition=snowball \
                --cpus-per-task=72 \
                --job-name="${species_fixed}_hisat2_index" \
                --output="data/species/$species_fixed/hisat2_index.out" \
                "$slurm_script"
        done < data/Mediophyceae_rnaseq_for_fastqdump.lst
        
        # Create an output file to signal that the job is done
        touch {output.done}
        """



rule write_hisat2_scripts:
    input:
        fastqdump_lst = "data/{taxon}_rnaseq_for_fastqdump.lst",
        download_done = "data/{taxon}_fastqdump.done",
        genome_done = "data/{taxon}_download.done"
    output:
        done = "data/{taxon}_hisat2_commands.lst"
    params:
        taxon = lambda wildcards: wildcards.taxon,
        threads = config['SLURM_ARGS']['cpus_per_task']
    wildcard_constraints:
        taxon="[^_]+"
    run:
        try:
            with open(input.fastqdump_lst, "r") as f:
                lines = f.readlines()
        except IOError:
            raise FileNotFoundError("Could not read the input file " + input.fastqdump_lst)
        # construct commands
        cmds = []
        for line in lines:
            species, sra_ids = line.strip().split("\t")
            species_fixed = species.replace(" ", "_")
            # construct a genome index
            cmd = f"hisat2-build -p {params.threads} data/species/{species_fixed}/genome.fa data/species/{species_fixed}/genome; "
            # make hisat2 folder
            cmd += f"mkdir -p data/species/{species_fixed}/hisat2; "
            for sra_id in sra_ids.split(","):
                # perform hisat2 alignment
                cmd += f"hisat2 -x data/species/{species_fixed}/genome -1 data/species/{species_fixed}/fastq/{sra_id}_1.fastq.gz -2 data/species/{species_fixed}/fastq/{sra_id}_2.fastq.gz -S data/species/{species_fixed}/hisat2/{sra_id}.sam"
                # convert sam 2 bam:
                cmd += f"samtools view --threads {params.threads} -bS data/species/{species_fixed}/hisat2/{sra_id}.sam > data/species/{species_fixed}/hisat2/{sra_id}.bam; "
                # sort bam file:
                cmd += f"samtools sort --threads {params.threads} data/species/{species_fixed}/hisat2/{sra_id}.bam -o data/species/{species_fixed}/hisat2/{sra_id}.sorted.bam; "
                # index the sorted bam file
                cmd += f"samtools index data/species/{species_fixed}/hisat2/{sra_id}.sorted.bam; "
                # remove the sam file
                cmd += f"rm data/species/{species_fixed}/hisat2/{sra_id}.sam; "
                # remove the unsorted bam file
                cmd += f"rm data/species/{species_fixed}/hisat2/{sra_id}.bam; "
            cmd += "\n"
            cmds.append(cmd)
        # create a folder for the scripts in data/species/{species_fixed}/hisat2_scripts
        script_dir = f"data/scripts_{params.taxon}"
        os.makedirs(script_dir, exist_ok=True)
        # write the commands script files
        for i, cmd in enumerate(cmds):
            # 


rule run_hisat2:
    input:
        cmds= "data/{taxon}_hisat2_commands.lst"
    output:
        done= "data/{taxon}_hisat2.done"
    params:
        taxon=lambda wildcards: wildcards.taxon,
        singularity_image="docker://teambraker/braker3:latest",
        partition = config['SLURM_ARGS']['partition'],
        threads = config['SLURM_ARGS']['cpus_per_task'],
        script_dir="data/scripts_{taxon}"  # Directory for scripts
    shell:
        """
        mkdir -p {params.script_dir}
        # Split commands into separate scripts based on library ID
        awk '/^hisat2/{{printf "\\n"}} /SRR/{{print > "{params.script_dir}/job_"substr($1, match($1, /SRR[0-9]+/), 10)".sh"}}' {input.cmds}

        # Add headers to each script and prepare for execution
        for script in {params.script_dir}/*.sh; do
            mv $script {script}.tmp  # Temporarily rename for safe appending
            echo '#!/bin/bash' > $script
            echo 'module load singularity' >> $script
            echo 'export APPTAINER_BIND="${{PWD}}:${{PWD}}"' >> $script
            cat {script}.tmp >> $script
            rm {script}.tmp  # Remove temporary file
            chmod +x $script
        done

        # Submit the array job
        sbatch --partition={params.partition} \
               --cpus-per-task={params.threads} \
               --job-name={params.taxon}_hisat2 \
               --array=1-$(ls {params.script_dir}/*.sh | wc -l) \
               --output={params.taxon}_hisat2_%A_%a.out \
               --wrap="singularity exec {params.singularity_image} bash {params.script_dir}/job_%a.sh"

        # Mark as done
        touch {output.done}
        """

'''

