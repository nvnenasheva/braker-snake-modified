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
#     and other parameters formatted for the `prefetch` and `fasterq-dump` utility, which can be used to download sequence data.
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
        unannotated_species = "data/checkpoints_dataprep/{taxon}_A03_blank.tbl"
    output:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B01_rnaseq_for_fastqdump.lst",
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
        python3 {input.download_script} -e {params.email} -t {input.unannotated_species} -f {output.fastqdump_lst};
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
#
# This job uses at most 12 threads but occupies the entire RAM of a node because it stores data in RAM
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
#
# This job uses at most 12 threads but occupies the entire RAM of a node because it stores data in RAM
rule run_download_fastq:
    input:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B01_rnaseq_for_fastqdump.lst",
        genome_done = "data/checkpoints_dataprep/{taxon}_A09_shorten_genomic_headers.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_B03_fastqdump.done"
    params:
        threads = 10,
        taxon = lambda wildcards: wildcards.taxon
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://katharinahoff/varus-notebook:v0.0.5"
    threads: 12
    resources:
        mem_mb=int(config['SLURM_ARGS']['mem_of_node']),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        echo "I am done with the SIF file, now starting, will take very long..."
        logfile=$PWD/data/checkpoints_dataprep/{params.taxon}_B03_fastqdump.log
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
                    echo "fasterq-dump data/species/$species_fixed/sra/$id/$id.sra --split-files -O data/species/$species_fixed/fastq -e {params.threads}" &>> $logfile
                    fasterq-dump data/species/$species_fixed/sra/$id/$id.sra --split-files -O data/species/$species_fixed/fastq -e {params.threads} &>> $logfile
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


rule run_hisat2_index:
    input:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B01_rnaseq_for_fastqdump.lst",
        download_done = "data/checkpoints_dataprep/{taxon}_B03_fastqdump.done",
        genome_done = "data/checkpoints_dataprep/{taxon}_A09_shorten_genomic_headers.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_B04_hisat2_index.done"
    params:
        taxon=lambda wildcards: wildcards.taxon,
        threads = config['SLURM_ARGS']['cpus_per_task']
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://teambraker/braker3:latest"
    threads: int(config['SLURM_ARGS']['cpus_per_task'])
    resources:
        mem_mb=int(config['SLURM_ARGS']['mem_of_node']),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        log="data/checkpoints_dataprep/{params.taxon}_B04_hisat2_index.log"
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
        done
        touch {output.done}
        """

rule run_hisat2:
    input:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B01_rnaseq_for_fastqdump.lst",
        genome_done = "data/checkpoints_dataprep/{taxon}_B04_hisat2_index.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_B05_hisat2.done"
    params:
        taxon=lambda wildcards: wildcards.taxon,
        threads = config['SLURM_ARGS']['cpus_per_task']
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://teambraker/braker3:latest"
    threads: int(config['SLURM_ARGS']['cpus_per_task'])
    resources:
        mem_mb=int(config['SLURM_ARGS']['mem_of_node']),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}";
        log=data/checkpoints_dataprep/{params.taxon}_B05_hisat2.log
        # if file $log does not exist yet
        if [ ! -f $log ]; then
            touch $log
        fi
        readarray -t lines < <(cat {input.fastqdump_lst})
        for line in "${{lines[@]}}"; do
            # Replace the first space with an underscore in the species name part of the line
            modified_line=$(echo "$line" | sed 's/\\([^\\t]*\\) /\\1_/')
            species=$(echo "$modified_line" | cut -f1)
            sra_ids=$(echo "$modified_line" | cut -f2)
            IFS=',' read -r -a sra_array <<< "$sra_ids"
            for sra_id in "${{sra_array[@]}}"; do
                if [ ! -f "data/species/$species/hisat2/${{sra_id}}.sam" ]; then
                    echo "hisat2 -p {params.threads} -x data/species/$species/genome/genome -1 data/species/$species/fastq/${{sra_id}}_1.fastq.gz -2 data/species/$species/fastq/${{sra_id}}_2.fastq.gz -S data/species/$species/hisat2/${{sra_id}}.sam" >> $log
                    hisat2 -p {params.threads} -x data/species/$species/genome/genome -1 data/species/$species/fastq/${{sra_id}}_1.fastq.gz -2 data/species/$species/fastq/${{sra_id}}_2.fastq.gz -S data/species/$species/hisat2/${{sra_id}}.sam &>> $log
                else
                    echo "data/species/$species/hisat2/${{sra_id}}.sam already exists" &>> $log
                fi
            done
        done
        touch {output.done}
        """

# remove bad libraries
rule remove_bad_libraries:
    input:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B01_rnaseq_for_fastqdump.lst",
        hisat2_log = "data/checkpoints_dataprep/{taxon}_B05_hisat2.log"
    output:
        done = "data/checkpoints_dataprep/{taxon}_B06_remove_bad_libraries.done",
        new_fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B06_rnaseq_for_fastqdump.lst",
    params:
        taxon=lambda wildcards: wildcards.taxon,
        mapping_threshold = float(config['FASTERQDUMP']['minimial_aligned_reads_percent']),
    wildcard_constraints:
        taxon="[^_]+"
    shell:
        '''
        log=data/checkpoints_dataprep/{params.taxon}_B06_remove_bad_libraries.log
        cmd="python3 scripts/remove_bad_sra.py -t {params.mapping_threshold} -l {input.hisat2_log} -o {output.new_fastqdump_lst} -p data/species -i {input.fastqdump_lst}"
        echo $cmd &> $log
        $cmd &>> $log
        touch {output.done}
        '''

rule run_sam_to_bam:
    input:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B06_rnaseq_for_fastqdump.lst",
        remove_bad_done = "data/checkpoints_dataprep/{taxon}_B06_remove_bad_libraries.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_B06_sam2bam.done"
    params:
        taxon=lambda wildcards: wildcards.taxon,
        threads = config['SLURM_ARGS']['cpus_per_task']
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://teambraker/braker3:latest"
    threads: int(config['SLURM_ARGS']['cpus_per_task'])
    resources:
        mem_mb=int(config['SLURM_ARGS']['mem_of_node']),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        log="data/checkpoints_dataprep/{params.taxon}_B06_sam2bam.log"
        echo "" > $log
        readarray -t lines < <(cat {input.fastqdump_lst})
        for line in "${{lines[@]}}"; do
            # Replace the first space with an underscore in the species name part of the line
            modified_line=$(echo "$line" | sed 's/\\([^\\t]*\\) /\\1_/')
            species=$(echo "$modified_line" | cut -f1)
            sra_ids=$(echo "$modified_line" | cut -f2)
            IFS=',' read -r -a sra_array <<< "$sra_ids"
            for sra_id in "${{sra_array[@]}}"; do
                if [ ! -f "data/species/$species/hisat2/${{sra_id}}.bam" ] && [ -f data/species/$species/hisat2/${{sra_id}}.sam ]; then
                    echo "samtools view --threads {params.threads} -bS data/species/$species/hisat2/${{sra_id}}.sam > data/species/$species/hisat2/${{sra_id}}.bam" &>> $log
                    samtools view --threads {params.threads} -bS data/species/$species/hisat2/${{sra_id}}.sam > data/species/$species/hisat2/${{sra_id}}.bam 2>> $log
                else
                    echo "data/species/$species/hisat2/${{sra_id}}.bam already exists" &>> $log
                fi
            done
        done
        touch {output.done}
        """

rule run_samtools_sort_single:
    input:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B06_rnaseq_for_fastqdump.lst",
        genome_done = "data/checkpoints_dataprep/{taxon}_B06_sam2bam.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_B07_samtools_sort_single.done"
    params:
        taxon=lambda wildcards: wildcards.taxon,
        threads = config['SLURM_ARGS']['cpus_per_task']
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://teambraker/braker3:latest"
    threads: int(config['SLURM_ARGS']['cpus_per_task'])
    resources:
        mem_mb=int(config['SLURM_ARGS']['mem_of_node']),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        log="data/checkpoints_dataprep/{params.taxon}_B07_samtools_sort_single.log"
        echo "" > $log
        readarray -t lines < <(cat {input.fastqdump_lst})
        for line in "${{lines[@]}}"; do
            # Replace the first space with an underscore in the species name part of the line
            modified_line=$(echo "$line" | sed 's/\\([^\\t]*\\) /\\1_/')
            species=$(echo "$modified_line" | cut -f1)
            sra_ids=$(echo "$modified_line" | cut -f2)
            IFS=',' read -r -a sra_array <<< "$sra_ids"
            for sra_id in "${{sra_array[@]}}"; do
                if [ ! -f "data/species/$species/hisat2/${{sra_id}}.sorted.bam" ] && [ -f data/species/$species/hisat2/${{sra_id}}.bam ]; then
                    echo "samtools sort --threads {params.threads} data/species/$species/hisat2/${{sra_id}}.bam -o data/species/$species/hisat2/${{sra_id}}.sorted.bam" &>> $log
                    samtools sort --threads {params.threads} data/species/$species/hisat2/${{sra_id}}.bam -o data/species/$species/hisat2/${{sra_id}}.sorted.bam &>> $log
                else
                    echo "data/species/$species/hisat2/${{sra_id}}.sorted.bam already exists" &>> $log
                fi
            done
        done
        touch {output.done}
        """

rule run_samtools_index_single:
    input:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B06_rnaseq_for_fastqdump.lst",
        genome_done = "data/checkpoints_dataprep/{taxon}_B07_samtools_sort_single.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_B08_samtools_index_single.done"
    params:
        taxon=lambda wildcards: wildcards.taxon,
        threads = config['SLURM_ARGS']['cpus_per_task']
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://teambraker/braker3:latest"
    threads: int(config['SLURM_ARGS']['cpus_per_task'])
    resources:
        mem_mb=int(config['SLURM_ARGS']['mem_of_node']),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        log="data/checkpoints_dataprep/{params.taxon}_B08_samtools_index_single.log"
        echo "" > $log
        readarray -t lines < <(cat {input.fastqdump_lst})
        for line in "${{lines[@]}}"; do
            # Replace the first space with an underscore in the species name part of the line
            modified_line=$(echo "$line" | sed 's/\\([^\\t]*\\) /\\1_/')
            species=$(echo "$modified_line" | cut -f1)
            sra_ids=$(echo "$modified_line" | cut -f2)
            IFS=',' read -r -a sra_array <<< "$sra_ids"
            for sra_id in "${{sra_array[@]}}"; do
                echo "samtools index data/species/$species/hisat2/${{sra_id}}.sorted.bam" &>> $log
                samtools index data/species/$species/hisat2/${{sra_id}}.sorted.bam &>> $log
            done
        done
        touch {output.done}
        """

rule cleanup_sam_bam_unsorted_files:
    input:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B06_rnaseq_for_fastqdump.lst",
        genome_done = "data/checkpoints_dataprep/{taxon}_B08_samtools_index_single.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_B09_cleanup_sam_bam_unsorted.done"
    params:
        taxon=lambda wildcards: wildcards.taxon,
        threads = config['SLURM_ARGS']['cpus_per_task']
    wildcard_constraints:
        taxon="[^_]+"
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        log="data/checkpoints_dataprep/{params.taxon}_B09_cleanup_sam_bam_unsorted.log"
        echo "" > $log
        readarray -t lines < <(cat {input.fastqdump_lst})
        for line in "${{lines[@]}}"; do
            # Replace the first space with an underscore in the species name part of the line
            modified_line=$(echo "$line" | sed 's/\\([^\\t]*\\) /\\1_/')
            species=$(echo "$modified_line" | cut -f1)
            sra_ids=$(echo "$modified_line" | cut -f2)
            IFS=',' read -r -a sra_array <<< "$sra_ids"
            for sra_id in "${{sra_array[@]}}"; do
                if [ -f data/species/$species/hisat2/${{sra_id}}.sam ] && [ data/species/$species/hisat2/${{sra_id}}.bam ]; then
                    echo "rm data/species/$species/hisat2/${{sra_id}}.sam" &>> $log
                    rm data/species/$species/hisat2/${{sra_id}}.sam &>> $log
                    echo "rm data/species/$species/hisat2/${{sra_id}}.bam" &>> $log
                    rm data/species/$species/hisat2/${{sra_id}}.bam &>> $log
                else
                    echo "data/species/$species/hisat2/${{sra_id}}.sam or data/species/$species/hisat2/${{sra_id}}.bam does not exist" &>> $log
                fi
            done
        done
        touch {output.done}
        """


rule run_merge_bam:
    input:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B06_rnaseq_for_fastqdump.lst",
        genome_done = "data/checkpoints_dataprep/{taxon}_B09_cleanup_sam_bam_unsorted.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_B10_merge_bam.done"
    params:
        taxon=lambda wildcards: wildcards.taxon,
        threads = config['SLURM_ARGS']['cpus_per_task']
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://teambraker/braker3:latest"
    threads: int(config['SLURM_ARGS']['cpus_per_task'])
    resources:
        mem_mb=int(config['SLURM_ARGS']['mem_of_node']),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        log="data/checkpoints_dataprep/{params.taxon}_B10_samtools_merge.log"
        echo "" > $log
        readarray -t lines < <(cat {input.fastqdump_lst})
        for line in "${{lines[@]}}"; do
            # Replace the first space with an underscore in the species name part of the line
            modified_line=$(echo "$line" | sed 's/\\([^\\t]*\\) /\\1_/')
            species=$(echo "$modified_line" | cut -f1)
            # count number of files matching the pattern data/species/$species/hisat2/*.sorted.bam, excluding file data/species/$species/hisat2/${{species}}.sorted.bam
            num_files=$(ls -1 data/species/$species/hisat2/*.sorted.bam | grep -v data/species/$species/hisat2/${{species}}.sorted.bam | wc -l)
            # merge all bam files of the same species with samtools merge 
            if [ $num_files -gt 1 ] && [ ! -f "data/species/$species/hisat2/${{species}}.bam" ]; then
                echo "samtools merge --threads {params.threads} data/species/$species/hisat2/${{species}}.sorted.bam data/species/$species/hisat2/*.sorted.bam" &>> $log
                samtools merge --threads {params.threads} data/species/$species/hisat2/${{species}}.bam data/species/$species/hisat2/*.sorted.bam &>> $log
            elif [ $numfiles -eq 1 ] && [ ! -f "data/species/$species/hisat2/${{species}}.bam" ] ; then
                cp data/species/$species/hisat2/*.sorted.bam data/species/$species/hisat2/${{species}}.bam
            else
                echo "data/species/$species/hisat2/${{species}}.bam already exists" &>> $log
            fi
        done
        touch {output.done}
        """

rule cleanup_sorted_bam_files:
    input:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B06_rnaseq_for_fastqdump.lst",
        genome_done = "data/checkpoints_dataprep/{taxon}_B10_merge_bam.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_B11_cleanup_sorted_bam.done"
    params:
        taxon=lambda wildcards: wildcards.taxon,
        threads = config['SLURM_ARGS']['cpus_per_task']
    wildcard_constraints:
        taxon="[^_]+"
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        log="data/checkpoints_dataprep/{params.taxon}_B11_cleanup_single_bam.log"
        echo "" > $log
        readarray -t lines < <(cat {input.fastqdump_lst})
        for line in "${{lines[@]}}"; do
            # Replace the first space with an underscore in the species name part of the line
            modified_line=$(echo "$line" | sed 's/\\([^\\t]*\\) /\\1_/')
            species=$(echo "$modified_line" | cut -f1)
            sra_ids=$(echo "$modified_line" | cut -f2)
            IFS=',' read -r -a sra_array <<< "$sra_ids"
            # delete the individual bam files
            for sra_id in "${{sra_array[@]}}"; do
                if [ -f "data/species/$species/hisat2/${{sra_id}}.sorted.bam" ]; then
                    echo "rm data/species/$species/hisat2/${{sra_id}}.sorted.bam" &>> $log
                    rm data/species/$species/hisat2/${{sra_id}}.sorted.bam &>> $log
                else
                    echo "data/species/$species/hisat2/${{sra_id}}.sorted.bam does not exist" &>> $log
                fi
            done
        done
        touch {output.done}
        """

rule run_sort_merged_bam:
    input:
        fastqdump_lst = "data/checkpoints_dataprep/{taxon}_B06_rnaseq_for_fastqdump.lst",
        genome_done = "data/checkpoints_dataprep/{taxon}_B10_merge_bam.done"
    output:
        done = "data/checkpoints_dataprep/{taxon}_B12_sort_merged_bam.done"
    params:
        taxon=lambda wildcards: wildcards.taxon,
        threads = config['SLURM_ARGS']['cpus_per_task']
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://teambraker/braker3:latest"
    threads: int(config['SLURM_ARGS']['cpus_per_task'])
    resources:
        mem_mb=int(config['SLURM_ARGS']['mem_of_node']),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}"; \
        log="data/checkpoints_dataprep/{params.taxon}_B12_samtools_sort_merged.log"
        echo "" > $log
        readarray -t lines < <(cat {input.fastqdump_lst})
        for line in "${{lines[@]}}"; do
            # Replace the first space with an underscore in the species name part of the line
            modified_line=$(echo "$line" | sed 's/\\([^\\t]*\\) /\\1_/')
            species=$(echo "$modified_line" | cut -f1)
            # merge all bam files of the same species with samtools merge 
            if [ ! -f "data/species/$species/hisat2/${{species}}.sorted.bam" ] && [ -f  data/species/$species/hisat2/${{species}}.bam ]; then
                echo "samtools sort --threads {params.threads} data/species/$species/hisat2/${{species}}.bam -o data/species/$species/hisat2/${{species}}.sorted.bam &>> $log" &>> $log
                samtools sort --threads {params.threads} data/species/$species/hisat2/${{species}}.bam -o data/species/$species/hisat2/rnaseq.s.bam &>> $log
            else
                echo "data/species/$species/hisat2/${{species}}.sorted.bam already exists" &>> $log
            fi
        done
        touch {output.done}
        """

'''
# A perfectly fine VARUS rule that does not work simply because VARUS tried to chdir
# chdir is not allowed in snakemake.
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
            echo "runVARUS.pl --aligner=HISAT --runThreadN={params.threads} --speciesGenome=../genome/genome.fa --readFromTable=0 --createindex=1 --verbosity=5 --latinGenus=$genus --latinSpecies=$spec --varusParameters=VARUSparameters.txt --outFileDir data/species/$species/varus 2> data/species/$species/varus/varus.err" &>> $log
            ( runVARUS.pl --aligner=HISAT --runThreadN={params.threads} \
                --speciesGenome=../genome/genome.fa \
                --readFromTable=0 --createindex=1 --verbosity=5 \
                --latinGenus=$genus --latinSpecies=$spec \
                --varusParameters=VARUSparameters.txt \
                --outFileDir data/species/$species/varus \
                2> data/species/$species/varus/varus.err )
        done
        touch {output.done}
        """
'''
