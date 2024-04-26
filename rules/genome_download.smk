rule download_assembly_info:
    output:
        processed_json = "data/{taxon}.json"
    params:
        taxon = lambda wildcards: wildcards.taxon
    singularity:
        "docker://katharinahoff/varus-notebook:v0.0.1"
    shell:
        """
        export SINGULARITY_BIND="${{PWD}}:${{PWD}}"; \
        datasets download genome taxon "{params.taxon}" --assembly-source genbank --dehydrated --filename {params.taxon}_ncbi.zip; \
        unzip -o {params.taxon}_ncbi.zip -d {params.taxon}_ncbi_dataset; \
        mv {params.taxon}_ncbi_dataset/ncbi_dataset/data/assembly_data_report.jsonl {output.processed_json}; \
        rm -rf {params.taxon}_ncbi_dataset {params.taxon}_ncbi.zip
        """
