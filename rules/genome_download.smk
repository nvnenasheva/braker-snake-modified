# Define rules for processing that are expected to generate the files listed in expected_outputs
# This is where your processing happens, ensuring the expected outputs are created

rule process_species:
    output:
        genome_fa = "{species_dir}/genome/genome.fa"
    shell:
        """
        mkdir -p {Path(output.genome_fa).parent}
        touch {output.genome_fa}
        """