rule download_assembly_info:
    output:
        raw_json = "data/{taxon}.json"
    params:
        taxon = lambda wildcards: wildcards.taxon
    singularity:
        "docker://katharinahoff/varus-notebook:v0.0.1"
    shell:
        """
        export SINGULARITY_BIND="${{PWD}}:${{PWD}}"; \
        datasets download genome taxon "{params.taxon}" --assembly-source genbank --dehydrated --filename {params.taxon}_ncbi.zip; \
        unzip -o {params.taxon}_ncbi.zip -d {params.taxon}_ncbi_dataset; \
        mv {params.taxon}_ncbi_dataset/ncbi_dataset/data/assembly_data_report.jsonl {output.raw_json}; \
        rm -rf {params.taxon}_ncbi_dataset {params.taxon}_ncbi.zip
        """

rule assembly_json_to_tbl:
    input:
        json_file = "data/{taxon}.json"
    output:
        processed_tbl = "data/{taxon}.tbl"
    params:
        taxon = lambda wildcards: wildcards.taxon
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

                    
