# Define rules for processing that are expected to generate the files listed in expected_outputs
# This is where your processing happens, ensuring the expected outputs are created

rule process_species:
    output: # Ensure this is not defined with wildcards here
        touch("{species_dir}/genome/genome.fa")
    shell:
        """
        mkdir -p {output[0].rsplit('/', 1)[0]}
        touch {output[0]}
        """