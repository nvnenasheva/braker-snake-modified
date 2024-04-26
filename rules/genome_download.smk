# Rule: download_assembly_info
# Purpose:
# This rule automates the process of downloading genome assembly info/ for a specified taxon
# from the NCBI database. It utilizes the NCBI datasets command-line tool within a Singularity
# container to ensure environment consistency and reproducibility. The downloaded data is
# initially compressed into a ZIP file, which is then extracted, and the relevant JSONL file
# containing assembly data is moved to a specified directory.
#
# Inputs:
#   - None directly specified, but the rule uses a parameter 'taxon' derived from wildcards,
#     typically provided by upstream rules or workflow orchestration.
#
# Outputs:
#   - raw_json: A JSONL file containing detailed assembly information for the specified taxon.
#     This file is used in downstream analysis and must be stored in the 'data' directory.
#
# Singularity:
#   - Uses a Docker container 'katharinahoff/varus-notebook:v0.0.1' converted to Singularity
#     image format, ensuring that the NCBI datasets tool and all its dependencies are correctly
#     configured and isolated from the host environment.
#
# Steps Executed:
# 1. Set up a bind mount to ensure current working directory ('${PWD}') is accessible inside the Singularity container.
# 2. Execute the 'datasets download genome' command specifying the taxon, source database (GenBank), and download format (dehydrated).
# 3. Unzip the downloaded file to a specified sub-directory.
# 4. Move the assembly data report (JSONL format) to the designated output location.
# 5. Clean up all intermediate files and directories to maintain a clean working environment.
rule download_assembly_info:
    output:
        raw_json = "data/{taxon}.json"
    params:
        taxon = lambda wildcards: wildcards.taxon
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://katharinahoff/varus-notebook:v0.0.1"
    shell:
        """
        export SINGULARITY_BIND="${{PWD}}:${{PWD}}"; \
        datasets download genome taxon "{params.taxon}" --assembly-source genbank --dehydrated --filename {params.taxon}_ncbi.zip; \
        unzip -o {params.taxon}_ncbi.zip -d {params.taxon}_ncbi_dataset; \
        mv {params.taxon}_ncbi_dataset/ncbi_dataset/data/assembly_data_report.jsonl {output.raw_json}; \
        rm -rf {params.taxon}_ncbi_dataset {params.taxon}_ncbi.zip; \
        mkdir -p data/species
        """


# Rule: assembly_json_to_tbl
# Purpose:
# This rule processes JSON files containing genome assembly metadata for different taxa.
# Each JSON file may contain one or multiple JSON objects representing different assembly reports.
# The rule reads the JSON file, parses multiple JSON objects if present, and then extracts
# relevant fields to produce a tab-separated values (TSV) file. The output file includes columns for
# accession, species name, assembly status, protein coding gene counts, contig N50 size, and the refSeq category.
# This rule handles JSON reading and parsing errors and will raise exceptions if files cannot be read or written.
#
# Input:
#   json_file - Path to the JSON file for a specific taxon, located in the 'data' directory.
# Output:
#   processed_tbl - Path to the output TSV file which will also be located in the 'data' directory.
# Processing Steps:
# 1. Read the entire JSON file content.
# 2. Split the content into individual JSON strings based on a specified delimiter.
# 3. Convert each JSON string back into a JSON object (dictionary) and handle any JSON parsing errors.
# 4. Extract the required information from each JSON object and write it to the output TSV file.
# 5. Handle and report any file input/output errors during reading and writing operations.
rule assembly_json_to_tbl:
    input:
        json_file = "data/{taxon}.json"
    output:
        processed_tbl = "data/{taxon}.tbl"
    params:
        taxon = lambda wildcards: wildcards.taxon
    wildcard_constraints:
        taxon="[^_]+"
    run:
        taxon = wildcards.taxon
        print(taxon)
        json_file_path = f"data/{taxon}.json"
        tbl_file_path = f"data/{taxon}.tbl"

        try:
            print(json_file_path)
            # Read the entire file into a single string
            with open(json_file_path, 'r') as file:
                file_content = file.read()
        except IOError:
            raise Exception(f"File not found: {json_file_path}")
        # Split the content by a delimiter that indicates the end of a JSON object, assuming each JSON object ends with '}' on a new line
        json_objects = file_content.split('}\n')
        # Filter out empty strings if any and add closing brace which was removed by split
        json_objects = [obj + '}' for obj in json_objects if obj.strip()]
        data = []
        for json_str in json_objects:
            try:
                # Load each JSON object as a dictionary
                data.append(json.loads(json_str))
            except json.JSONDecodeError as e:
                print(f"Failed to decode JSON: {e}")
            continue
        try:
            with open(tbl_file_path, "w") as f:
                f.write("accession\tspecies\tstatus\tproteins\tcontigN50\trefseqCategory\n")
                # If data is a list of entries, iterate through it
                # If it's a single dict, adjust the loop or wrap `data` in a list: [data]
                for entry in data:
                    accession = entry.get('accession', '')
                    species = entry['organism'].get('organismName', '')
                    status = entry['assemblyInfo'].get('assemblyStatus', '')
                    proteins = entry.get('annotationInfo', {}).get('stats', {}).get('geneCounts', {}).get('proteinCoding', 'N/A')
                    contigN50 = entry['assemblyStats'].get('contigN50', '')
                    refSeqCategory = entry['assemblyInfo'].get('refseqCategory', 'N/A')
                    f.write(f"{accession}\t{species}\t{status}\t{proteins}\t{contigN50}\t{refSeqCategory}\n")
        except IOError:
            raise Exception(f"Error writing to file: {tbl_file_path}")


# Rule: classify_species
# Purpose:
# This rule is designed to classify species based on their protein counts from genome assembly data.
# It separates the species into two categories: those that are well-annotated (with protein counts
# exceeding 1000) and those that are either poorly annotated or have missing protein data ('N/A').
#
# Inputs:
#   - tbl_file: A tab-separated file ('.tbl') for each taxon, containing genome assembly data.
#     This file includes various metrics among which are protein counts.
#
# Outputs:
#   - already_annotated_tbl: Outputs a file listing species that are considered well-annotated.
#   - not_annotated_tbl: Outputs a file listing species that are either not well-annotated or
#     whose protein count data is missing.
#
# Parameters:
#   - taxon: A wildcard parameter that dynamically accepts the name of the taxon from the workflow.
#
# Processing Steps:
# 1. Read the input TBL file into a pandas DataFrame. Non-numeric 'proteins' values are converted to NaN
#    to facilitate numerical operations and comparisons.
# 2. Filter the DataFrame to create two subsets:
#    - 'annotated_data': Contains entries with protein counts greater than 1000.
#      This subset is further processed to keep only the highest quality genome per species,
#      prioritized by 'refseqCategory' and the highest protein count.
#    - 'blank_data': Contains entries with protein counts less than or equal to 1000 or with missing protein data.
#      This subset is sorted by 'contigN50' to prioritize genomes based on genome assembly quality, keeping the
#      highest 'contigN50' value per species.
# 3. Both subsets are written to separate output files in a tab-separated format.
rule classify_species:
    input:
        tbl_file = "data/{taxon}.tbl"
    output:
        already_annotated_tbl = "data/{taxon}_annotated.tbl",
        not_annotated_tbl = "data/{taxon}_blank.tbl"
    params:
        taxon = lambda wildcards: wildcards.taxon
    wildcard_constraints:
        taxon="[^_]+"
    run:
        taxon = wildcards.taxon
        tbl_file_path = f"data/{taxon}.tbl"
        annotated_tbl_path = f"data/{taxon}_annotated.tbl"
        blank_tbl_path = f"data/{taxon}_blank.tbl"
        
        # Read the data, handling non-numeric proteins values
        try:
            data = pd.read_csv(tbl_file_path, sep="\t")
            data['proteins'] = pd.to_numeric(data['proteins'], errors='coerce')  # Convert 'N/A' to NaN
        except IOError:
            raise Exception(f"Error reading file: {tbl_file_path}")

        # Annotated data filter
        annotated_data = data[data['proteins'] > 1000]
        annotated_data = annotated_data.sort_values(by=['refseqCategory', 'proteins'], ascending=[False, False]).drop_duplicates(subset='species', keep='first')

        # Blank data filter
        blank_data = data[(data['proteins'] <= 1000) | (data['proteins'].isna())]
        blank_data = blank_data.sort_values(by=['contigN50'], ascending=False).drop_duplicates(subset='species', keep='first')

        # Output data to files
        try:
            annotated_data.to_csv(annotated_tbl_path, sep="\t", index=False)
            blank_data.to_csv(blank_tbl_path, sep="\t", index=False)
        except IOError:
            raise Exception(f"Error writing to file: {annotated_tbl_path} or {blank_tbl_path}")


# This Snakemake rule 'prepare_download_assemblies_from_ncbi' is designed to automate the preparation of shell scripts
# for downloading genomic assemblies from the NCBI database. The rule uses specific table files containing species and 
# accession numbers to generate commands for downloading and organizing genomic data into structured directories. The
# output is a shell script tailored for each taxon that handles directory creation, data retrieval, and data extraction.
rule prepare_download_assemblies_from_ncbi:
    input:
        annotated_tbl_file = "data/{taxon}_annotated.tbl",
        blank_tbl_file = "data/{taxon}_blank.tbl"
    output:
        download_script = "data/{taxon}_download.sh"
    params:
        taxon = lambda wildcards: wildcards.taxon
    wildcard_constraints:
        taxon="[^_]+"
    run:
        # Read the annotated table to get accessions and species names
        print(input.blank_tbl_file)
        df_blank = pd.read_csv(input.blank_tbl_file, sep="\t", usecols=['accession', 'species'])
        df_anno = pd.read_csv(input.annotated_tbl_file, sep="\t", usecols=['accession', 'species'])
        # Prepare a file to capture assembly information
        try:
            with open(output.download_script, 'w') as outfile:
                command = ""
                for index, row in df_blank.iterrows():
                    # take the row["species"] field and replace all spaces by _, otherwise remove all special characters
                    species = row["species"].replace(" ", "_")
                    species = re.sub(r'[/?,.*&\\;]+', '', species)
                    command += f"cd data/species; "
                    command += f"mkdir {species}; cd {species};"
                    command += f"datasets download genome accession {row['accession']} --filename {row['accession']}_assembly.zip; "
                    command += f"unzip -o {row['accession']}_assembly.zip; "
                    command += f"mkdir genome; "
                    command += f"mv ncbi_dataset/data/{row['accession']}/*.fna genome/genome.fa; "
                    command += f"rm -rf ncbi_dataset {row['accession']}_assembly.zip; cd ../../..;\n"
                for index, row in df_anno.iterrows():
                    # take the row["species"] field and replace all spaces by _, otherwise remove all special characters
                    species = row["species"].replace(" ", "_")
                    species = re.sub(r'[/?,.*&\\;]+', '', species)
                    command += f"cd data/species; "
                    command += f"mkdir {species}; cd {species};"
                    command += f"datasets download genome accession {row['accession']} --filename {row['accession']}_assembly.zip --include genome,gff3,protein,gtf; "
                    command += f"unzip -o {row['accession']}_assembly.zip; "
                    command += f"mkdir genome; "
                    command += f"mv ncbi_dataset/data/{row['accession']}/*.fna genome/genome.fa; "
                    command += f"mkdir annot;"
                    command += f"mv ncbi_dataset/data/{row['accession']}/*.gff annot/annot.gff3;"
                    command += f"mv ncbi_dataset/data/{row['accession']}/*.gtf annot/annot.gtf;"
                    command += f"mkdir prot;"
                    command += f"mv ncbi_dataset/data/{row['accession']}/*.faa prot/protein.faa;"
                    command += f"rm -rf ncbi_dataset {row['accession']}_assembly.zip; cd ../../..;\n"
                outfile.write(command)
        except IOError:
            raise Exception(f"Error writing to file: {output.download_script}")


# This rule executes the download commands specified in the download script generated by the rule 'prepare_download_assemblies_from_ncbi'. 
# The rule is not parallelized as the data input/output operations are likely to be the bottleneck.
rule run_download_commands:
    input:
        download_script = "data/{taxon}_download.sh"
    output:
        done = touch("data/{taxon}_download.done")
    params:
        taxon = lambda wildcards: wildcards.taxon
    wildcard_constraints:
        taxon="[^_]+"
    singularity:
        "docker://katharinahoff/varus-notebook:v0.0.1"
    shell:
        """
        export SINGULARITY_BIND="${{PWD}}:${{PWD}}"; \
        bash {input.download_script}; \
        touch {output.done}
        """

