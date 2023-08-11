rule sourmash_gather:
    input: 
        os.path.join(OUTPUT, "1-sketch", "{sample}.dna.sig.zip"),
        databases = lambda w: search_databases[f"k{w.ksize}"] 
    output:
        gather_csv=os.path.join(OUTPUT, '2-gather', '{sample}.k{ksize}.gather.csv'),
        gather_txt=os.path.join(OUTPUT, '2-gather', '{sample}.k{ksize}.gather.txt')
    params:
        threshold_bp = sourmash_params.get('threshold_bp', '50000')
    threads: 
        THREADS
    resources:
        mem_mb=BigJobMem,
        time=1320,
    log: 
        os.path.join(LOGS, "gather", "{sample}.k{ksize}.gather.log")
    benchmark: 
        os.path.join(BENCHMARKS, "gather", "{sample}.k{ksize}.gather.benchmark")
    conda: 
        os.path.join("..", "envs", "sourmash.yml")
    shell:
        # touch output to let workflow continue in cases where 0 results are found
        """
        echo "DB(s): {input.databases}"
        echo "DB(s): {input.databases}" > {log}

        sourmash gather {input.query} {input.databases} --dna --ksize {wildcards.ksize} \
                 --threshold-bp {params.threshold_bp} \
                 -o {output.gather_csv} > {output.gather_txt} 2>> {log}
        
        touch {output.gather_txt}
        touch {output.gather_csv}
        """


rule aggr_gather:
    input: 
        expand(os.path.join(OUTPUT, '2-gather', '{sample}.k{ksize}.gather.csv'), sample = SAMPLES),
        expand(os.path.join(OUTPUT, '2-gather', '{sample}.k{ksize}.gather.txt'), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "sourmash_gather.flag")
    threads: 
        1
    resources:
        mem_mb=SmallJobMem,
        time=3,
    shell:
        """
        touch {output[0]}
        """

