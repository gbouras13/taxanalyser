
rule sourmash_tax_metagenome:
    input:
        gather = os.path.join(OUTPUT, '2-gather', '{sample}.k{ksize}.gather.csv'),
        lineages = sourmash_params['database_lineage_files'],
    output:
        os.path.join(OUTPUT, '3-taxprofile', '{sample}.k{ksize}.gather.genbank.krona.tsv'),
        os.path.join(OUTPUT, '3-taxprofile', '{sample}.k{ksize}.gather.genbank.summarized.csv'),
        os.path.join(OUTPUT, '3-taxprofile', '{sample}.k{ksize}.gather.genbank.kreport.txt'),
    threads: 
        1
    resources:
        mem_mb=BigJobMem,
        time=1000,
    params:
        outd= lambda w: os.path.join(OUTPUT, f'3-taxprofile'),
        out_base= lambda w: f'{w.sample}.k{w.ksize}.gather.genbank',
    conda: 
        os.path.join("..", "envs", "sourmash.yml")
    log: 
        os.path.join(LOGS, "tax_metagenome", "{sample}.k{ksize}.tax_metagenome.log")
    benchmark: 
        os.path.join(BENCHMARKS, "tax_metagenome", "{sample}.k{ksize}.tax_metagenome.benchmark")
    shell:
        """
        mkdir -p {params.outd}
        sourmash tax metagenome -g {input.gather} -t {input.lineages} -o {params.out_base} \
                                --output-dir {params.outd} --output-format krona csv_summary kreport \
                                --rank species 2> {log}
        """

rule aggr_tax_metagenome:
    input: 
        expand(os.path.join(OUTPUT, '3-taxprofile', '{sample}.k{ksize}.gather.genbank.krona.tsv'), sample = SAMPLES, ksize = KSIZE)
    output:
        os.path.join(FLAGS, "sourmash_tax_metagenome.flag")
    threads: 
        1
    resources:
        mem_mb=SmallJobMem,
        time=3,
    shell:
        """
        touch {output[0]}
        """