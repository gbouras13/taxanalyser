rule tax_annotate:
    input:
        gather = os.path.join(OUTPUT, '2-gather', '{sample}.k{ksize}.gather.csv'),
        lineages = sourmash_params['database_lineage_files'],
    output:
        os.path.join(OUTPUT, '2-gather', '{sample}.k{ksize}.gather.with-lineages.csv'),
    threads: 
        1
    resources:
        mem_mb=BigJobMem,
        time=240,
    params:
        outd= lambda w: os.path.join(out_dir, f'2-gather'),
    conda: 
        os.path.join("..", "envs", "sourmash.yml")
    log: 
        os.path.join(LOGS, "tax_annotate", "{sample}.k{ksize}.tax_annotate.log")
    benchmark: 
        os.path.join(BENCHMARKS, "tax_annotate", "{sample}.k{ksize}.tax_annotate.benchmark")
    shell:
        """
        mkdir -p {params.outd}
        sourmash tax annotate -g {input.gather} -t {input.lineages} -o {params.outd} 2> {log}
        """

rule aggr_tax_annotate:
    input: 
        expand(os.path.join(OUTPUT, '2-gather', '{sample}.k{ksize}.gather.with-lineages.csv'), sample = SAMPLES, ksize = KSIZE)
    output:
        os.path.join(FLAGS, "sourmash_tax_annotate.flag")
    threads: 
        1
    resources:
        mem_mb=SmallJobMem,
        time=3,
    shell:
        """
        touch {output[0]}
        """