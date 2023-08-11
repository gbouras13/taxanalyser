rule sourmash_sketch_dna:
    input: 
        get_input_lr_fastqs
    output:
        os.path.join(OUTPUT, "1-sketch", "{sample}.dna.sig.zip")
    threads: 
        1
    resources:
        mem_mb=BigJobMem,
        time=240,
    log: 
        os.path.join(LOGS, "sketch", "{sample}.sketch_dna.log")
    benchmark: 
        os.path.join(BENCHMARKS, "sketch", "{sample}.sketch_dna.benchmark")
    conda: 
        os.path.join("..", "envs", "sourmash.yml")
    shell:
        """
        sourmash sketch dna {input} -p k=21,k=31,k=51,dna,scaled=1000,abund \
                                    --name {wildcards.sample} -o {output} 2> {log}
        """


rule aggr_sketch:
    input: 
        expand(os.path.join(OUTPUT, "1-sketch", "{sample}.dna.sig.zip"), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "sourmash_sketch.flag")
    threads: 
        1
    resources:
        mem_mb=SmallJobMem,
        time=3,
    shell:
        """
        touch {output[0]}
        """

