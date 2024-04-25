rule download_and_process_genomes:
    output:
        processed_json = "data/{taxon}_processed.json"
    params:
        taxon = lambda wildcards: wildcards.taxon
    singularity:
        "docker://katharinahoff/varus-notebook:devel"
    shell:
        """
        export SINGULARITY_BIND="${{PWD}}:${{PWD}}"; \
        ./datasets download genome taxon "{params.taxon}" --assembly-source genbank --dehydrated; \
        unzip -o ncbi_dataset.zip -d ncbi_dataset; \
        python3 -c "import json; import pandas as pd; \
            df = pd.read_json('ncbi_dataset/data/assembly_data_report.jsonl', lines=True); \
            df = df[df['assemblyStatus'] == 'current']; \
            records = df.apply(lambda x: {{'accession': x['accession'], 'organismName': x['organism']['organismName']}}, axis=1).to_dict(orient='records'); \
            with open('{output.processed_json}', 'w') as f: json.dump(records, f)"
        """