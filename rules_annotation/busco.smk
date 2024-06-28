# rule that performs a BUSCO assessment on the braker.aa / annotation protein files and the genome.fa file
# works similar to the braker rule
rule run_busco:
    input:
        input = "data/checkpoints_annotate/{spid}_braker.done"
    output:
        done = "data/checkpoints_annotate/{spid}_busco.done"
    params:
        spid="{spid}",
        threads = config['SLURM_ARGS']['cpus_per_task'],
        csv = config['INPUT']['species_csv']
    threads: int(config['SLURM_ARGS']['cpus_per_task'])
    resources:
        mem_mb=int(config['SLURM_ARGS']['mem_of_node']),
        runtime=int(config['SLURM_ARGS']['max_runtime'])
    shell:
        """
        export APPTAINER_BIND="${{PWD}}:${{PWD}}";
        wd=${{PWD}}
        log=data/checkpoints_annotate/{params.spid}_runningbusco.log
        touch $log
        echo "$wd/scripts/run_busco.sh -s {params.spid} -t {params.threads} -c {params.csv} -o {output}" &>> $log
        $wd/scripts/run_busco.sh -s {params.spid} -t {params.threads} -c {params.csv} -o {output}
        """

# rule that extracts the busco scores 
rule extract_busco_scores:
    input:
        input = "data/checkpoints_annotate/{spid}_busco.done"
    output:
        table = "data/checkpoints_annotate/{spid}_busco.csv"
    params:
        spid = "{spid}"
    run:
        #get busco scores
        pattern = r'C:\d+(\.\d+)?%\[S:\d+(\.\d+)?%,D:\d+(\.\d+)?%\],F:\d+(\.\d+)?%,M:\d+(\.\d+)?%,n:\d+'
        scores = []

        subdirectories = ['genomefile', 'annotfile', 'brakerfile']
        for subdir in subdirectories:
            file_path = f"data/species/{params.spid}/busco/{subdir}/short_summary.specific.stramenopiles_odb10.{subdir}.txt"
            if os.path.exists(file_path):  # Check if file exists
                with open(file_path, 'r') as file:
                    for line in file:
                        match = re.search(pattern,line)
                        if match:
                            scores.append(match.group(0))
                            break
            else:
                scores.append('NA')

        #get a dataframe with the scores
        data = {'species': [params.spid], 'genome_score': scores[0], 'annot_score': scores[1], 'braker_score': scores[2]}
        df = pd.DataFrame(data)

        # Save the DataFrame to a CSV file
        df.to_csv(output.table, index=False) 
