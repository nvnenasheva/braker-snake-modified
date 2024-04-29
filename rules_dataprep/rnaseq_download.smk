# Rule: retrieve_rnaseq_info_from_sra
# Purpose:
#   This rule automates the retrieval of RNA sequencing (RNA-seq) information for specified taxa
#   from the NCBI Sequence Read Archive (SRA). It utilizes a Python script to query NCBI for RNA-seq
#   data related to unannotated species listed in a provided table. The rule processes each taxon
#   listed in the 'data/checkpoints_dataprep/{taxon}_blank.tbl', extracting necessary information for downstream RNA-seq
#   analysis pipelines.
#
# Inputs:
#   - download_script: A Python script ('scripts/retrieve_rnaseq_info.py') that handles the querying
#     of NCBI databases and prepares output lists of SRA accession numbers and other relevant data.
#   - unannotated_species: A table listing species without annotated genomes; used to determine which
#     species' RNA-seq data to retrieve.
#
# Outputs:
#   - fastqdump_lst: A list file ('data/checkpoints_dataprep/{taxon}_rnaseq_for_fastqdump.lst') containing accession numbers
#     and other parameters formatted for the `fastq-dump` utility, which can be used to download sequence data.
#   - varus_list: A list file ('data/checkpoints_dataprep/{taxon}_rnaseq_for_varus.lst') used for VARUS, a pipeline that automatically
#     optimizes parameters for RNA-seq data retrieval and assembly.
#   - done: A simple checkpoint file ('data/checkpoints_dataprep/{taxon}_rnaseq_info.done') indicating successful completion of data retrieval.
#
# Parameters:
#   - taxon: A wildcard parameter defining the specific taxon being processed. Extracted from the filename pattern.
#   - email: An email address registered with NCBI, required for making API requests to NCBI databases.
#
# ToDo:
#   - make the number of libraries to download random & configurable in terms of number by config.ini
#   - make sure the libraries are mRNA and do not match microbiome or metagenome
rule retrieve_rnaseq_info_from_sra:
    input:
        download_script = "scripts/retrieve_rnaseq_info.py",
        unannotated_species = "data/checkpoints_dataprep/{taxon}_blank.tbl"
    output:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_rnaseq_for_fastqdump.lst",
        varus_list = "data/checkpoints_dataprep/{taxon}_rnaseq_for_varus.lst",
        done = "data/checkpoints_dataprep/{taxon}_rnaseq_info.done"
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
#   This rule takes extremely long! It runs on SLURM but cannot take full advantage of a larger
#   number of threads since i/o is the bottleneck. We are using /dev/shm to speed up as much
#   possible.
#   It must regrettably expected that this rule fails, repeatedly. Do not worry, just restart the
#   workflow. The rule is implemented in a way that it will pick up where it left off. 
#   This rule is designed to process a list of SRA accession numbers and use `fasterq-dump` to download
#   paired-end fastq files for each accession. It handles the preparation of directory structures
#   and the organization of downloaded files within those directories. This is crucial for downstream
#   analysis such as RNA-seq data processing or genomic assemblies.
#
# Inputs:
#   - fastqdump_lst: A file ('data/checkpoints_dataprep/{taxon}_rnaseq_for_fastqdump.lst') that lists species and their
#     corresponding SRA accession numbers. Each line contains a species name and a comma-separated list
#     of SRA IDs.
#
# Outputs:
#   - done: A file ('data/checkpoints_dataprep/{taxon}_fastqdump.done') that signals the successful completion of the downloads
#     and processing of all SRA data listed in the input file.
#
# Operations:
#   - The rule first sets up the environment to ensure correct directory binding when using Singularity.
#   - It reads each line from the input file, processes the species name to replace spaces with underscores
#     (for file naming consistency), and creates a directory for storing the downloaded fastq files.
#   - For each SRA ID, `fasterq-dump` is executed to retrieve paired-end reads, and the resultant files are
#     immediately compressed using gzip to save space.
#   - A log file is generated to capture the output and errors of the download process, helping in troubleshooting
#     and ensuring transparency of the operation.
#   - Upon successful execution of all commands, a 'done' file is created as a checkpoint indicating
#     the completion of the task.
rule download_fastq:
    input:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_rnaseq_for_fastqdump.lst",
        genome_done = "data/checkpoints_dataprep/{taxon}_download.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_fastqdump.done"
    params:
        threads = int(int(config['SLURM_ARGS']['cpus_per_task'])/4),
        tmpdir = config['SLURM_ARGS']['tmp_dir']
    singularity:
        "docker://katharinahoff/varus-notebook:v0.0.5"
    threads: 12
    resources:
        mem_mb=int(int(config['SLURM_ARGS']['mem_of_node'])/int(config['SLURM_ARGS']['cpus_per_task'])*12),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        echo "I am done with the SIF file, now starting, will take very long..."
        export APPTAINER_BIND="${{PWD}}:${{PWD}}";
        logfile=$PWD/data/checkpoints_dataprep/{wildcards.taxon}_fastqdump.log
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
                    echo "prefetch -O data/species/$species_fixed/sra $id" &>> $logfile
                    prefetch -O data/species/$species_fixed/sra $id &>> $logfile
                    echo "fasterq-dump data/species/$species_fixed/sra/$id/$id.sra --split-files -O data/species/$species_fixed/fastq -t {params.tmpdir} -e {params.threads}" &>> $logfile
                    fasterq-dump data/species/$species_fixed/sra/$id/$id.sra --split-files -O data/species/$species_fixed/fastq -t {params.tmpdir} -e {params.threads} &>> $logfile
                    rm data/species/$species_fixed/sra/$id/$id.sra
                    echo "gzip data/species/$species_fixed/fastq/${{id}}_1.fastq" &>> $logfile
                    gzip data/species/$species_fixed/fastq/${{id}}_1.fastq &>> $logfile
                    gzip data/species/$species_fixed/fastq/${{id}}_2.fastq
                fi
            done
            rm -rf data/species/$species_fixed/sra
        done < {input.fastqdump_lst} &>> $logfile
        
        echo "touch {output.done}" &>> $logfile
        touch {output.done}
        """

# Rule: run_hisat2
# Description: This rule is responsible for building a HISAT2 index for RNA-seq data analysis, aligning RNA-seq reads to the genome, and processing the alignment outputs.
# It performs several operations:
# 1. Reads a list of species and their corresponding SRA IDs from a specified file.
# 2. For each species, it checks if the necessary directories exist, and if not, it creates them.
# 3. It then proceeds to build the HISAT2 index for the genome of each species.
# 4. For each SRA ID associated with a species, the rule runs HISAT2 to align the RNA-seq reads to the genome, converts the output SAM to BAM format, sorts and indexes the BAM files.
# 5. It merges all BAM files for a given species into a single sorted BAM file, indexes it, and cleans up individual BAM files to save space.
# 6. Finally, it touches a file to indicate completion of the process.
# This rule utilizes significant computational resources (as specified in SLURM_ARGS and passed to the rule) and runs in a Singularity container to ensure environment consistency.
# WARNING: this rule executes once for each taxon. Depending on genome size and the number of SRA IDs, this rule may take a long time to complete.
#          this may exceed the runtime limit on our cluster!
#          If we figure out that it exceeds indeed the runtime limit, we may have to further take this apart (e.g. one rule building the index,
#          one executing hisat, one converting to bam, one sorting, one merging).
rule run_hisat2:
    input:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_rnaseq_for_fastqdump.lst",
        download_done = "data/checkpoints_dataprep/{taxon}_fastqdump.done",
        genome_done = "data/checkpoints_dataprep/{taxon}_download.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_hisat2.done"
    params:
        taxon=lambda wildcards: wildcards.taxon,
        threads = config['SLURM_ARGS']['cpus_per_task']
    singularity:
        "docker://teambraker/braker3:latest"
    threads: int(config['SLURM_ARGS']['cpus_per_task'])
    resources:
        mem_mb=int(config['SLURM_ARGS']['mem_of_node']),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        log="data/checkpoints_dataprep/{params.taxon}_hisat2.log"
        echo "" > $log
        readarray -t lines < <(cat {input.fastqdump_lst})
        for line in "${{lines[@]}}"; do
            # Replace the first space with an underscore in the species name part of the line
            modified_line=$(echo "$line" | sed 's/\\([^\\t]*\\) /\\1_/')
            species=$(echo "$modified_line" | cut -f1)
            if [ ! -d "data/species/$species/hisat2" ]; then
                mkdir -p "data/species/$species/hisat2"
            fi
            echo "hisat2-build -p {params.threads} data/species/$species/genome/genome.fa data/species/$species/genome/genome" >> $log
            which hisat2-build &>> $log
            hisat2-build -p {params.threads} data/species/$species/genome/genome.fa data/species/$species/genome/genome &>> $log
            sra_ids=$(echo "$modified_line" | cut -f2)
            IFS=',' read -r -a sra_array <<< "$sra_ids"
            for sra_id in "${{sra_array[@]}}"; do
                if [ ! -f "data/species/$species/hisat2/${{sra_id}}.sorted.bam" ]; then
                    echo "hisat2 -p {params.threads} -x data/species/$species/genome/genome -1 data/species/$species/fastq/${{sra_id}}_1.fastq.gz -2 data/species/$species/fastq/${{sra_id}}_2.fastq.gz -S data/species/$species/hisat2/${{sra_id}}.sam" >> $log
                    hisat2 -p {params.threads} -x data/species/$species/genome/genome -1 data/species/$species/fastq/${{sra_id}}_1.fastq.gz -2 data/species/$species/fastq/${{sra_id}}_2.fastq.gz -S data/species/$species/hisat2/${{sra_id}}.sam &>> $log
                    echo "samtools view --threads {params.threads} -bS data/species/$species/hisat2/${{sra_id}}.sam > data/species/$species/hisat2/${{sra_id}}.bam" &>> $log
                    samtools view --threads {params.threads} -bS data/species/$species/hisat2/${{sra_id}}.sam > data/species/$species/hisat2/${{sra_id}}.bam &>> $log
                    echo "samtools sort --threads {params.threads} data/species/$species/hisat2/${{sra_id}}.bam -o data/species/$species/hisat2/${{sra_id}}.sorted.bam" &>> $log
                    samtools sort --threads {params.threads} data/species/$species/hisat2/${{sra_id}}.bam -o data/species/$species/hisat2/${{sra_id}}.sorted.bam &>> $log
                    echo "samtools index data/species/$species/hisat2/${{sra_id}}.sorted.bam" &>> $log
                    samtools index data/species/$species/hisat2/${{sra_id}}.sorted.bam &>> $log
                    echo "rm data/species/$species/hisat2/${{sra_id}}.sam" &>> $log
                    rm data/species/$species/hisat2/${{sra_id}}.sam &>> $log
                    echo "rm data/species/$species/hisat2/${{sra_id}}.bam" &>> $log
                    rm data/species/$species/hisat2/${{sra_id}}.bam &>> $log
                fi
            done
            # merge all bam files of the same species with samtools merge 
            if [ ! -f "data/species/$species/hisat2/${{species}}.sorted.bam" ]; then
                echo "samtools merge --threads {params.threads} data/species/$species/hisat2/${{species}}.sorted.bam data/species/$species/hisat2/*.sorted.bam" &>> $log
                samtools merge --threads {params.threads} data/species/$species/hisat2/${{species}}.sorted.bam data/species/$species/hisat2/*.sorted.bam &>> $log
                echo "samtools index data/species/$species/hisat2/${{species}}.sorted.bam &>> $log"
                samtools index data/species/$species/hisat2/${{species}}.sorted.bam &>> $log
            fi
            # delete the individual bam files
            for sra_id in "${{sra_array[@]}}"; do
                echo "rm data/species/$species/hisat2/${{sra_id}}.sorted.bam" &>> $log
                rm data/species/$species/hisat2/${{sra_id}}.sorted.bam &>> $log
            done
        done
        touch {output.done}
        """


rule run_varus:
    input:
        varus_lst = "data/checkpoints_dataprep/{taxon}_rnaseq_for_varus.lst",
        genome_done = "data/checkpoints_dataprep/{taxon}_download.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_varus.done"
    params:
        threads = int(config['SLURM_ARGS']['cpus_per_task']),
        taxon=lambda wildcards: wildcards.taxon
    singularity:
        "docker://katharinahoff/varus-notebook:v0.0.5"
    threads: 12
    resources:
        mem_mb=int(int(config['SLURM_ARGS']['mem_of_node'])/int(config['SLURM_ARGS']['cpus_per_task'])*12),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        log="data/checkpoints_dataprep/{params.taxon}_varus.log"
        echo "" > $log
        readarray -t lines < <(cat {input.varus_lst})
        for line in "${{lines[@]}}"; do
            echo "Original line: $line" >> $log
            modified_line=$(echo "$line" | sed 's/\\([^\\t]*\\) /\\1_/')
            species=$(echo "$modified_line" | cut -f1)
            echo "Modified species line: $species" >> $log
            pwd >> $log
            if [ ! -d "data/species/$species/varus" ]; then
                mkdir -p "data/species/$species/varus"
            fi
            IFS='_' read -r genus spec <<< "$species"
            echo "runVARUS.pl --aligner=HISAT --runThreadN={params.threads} --speciesGenome=data/species/$species/genome/genome.fa --readFromTable=0 --createindex=1 --verbosity=5 --latinGenus=$genus --latinSpecies=$spec --varusParameters=VARUSparameters.txt --genomeDir ./ --outFileNamePrefix data/species/$species/varus/ --logfile=data/species/$species/varus/varus.log 2> data/species/$species/varus/varus.err" &>> $log
            runVARUS.pl --aligner=HISAT --runThreadN={params.threads} \
                --speciesGenome=data/species/$species/genome/genome.fa \
                --readFromTable=0 --createindex=1 --verbosity=5 \
                --latinGenus=$genus --latinSpecies=$spec \
                --varusParameters=VARUSparameters.txt \
                --outFileDir data/species/$species/varus/ \
                --logfile=data/species/$species/varus/varus.log 2> data/species/$species/varus/varus.err
        done
        touch {output.done}
        """


rule cleanup_data:
    input:
        hisat_done = "data/checkpoints_dataprep/{taxon}_hisat2.done",
        varus_done = "data/checkpoints_dataprep/{taxon}_varus.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_cleanup.done"
    shell:
        """
        echo "Hello world"
        touch {output.done}
        """