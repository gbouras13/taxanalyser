
rule run_sylph:
    input:
        reads = expand(os.path.join(dir.out.qc, "{sample}_filtlong.fastq.gz"), sample = SAMPLES)
    output:
        results=os.path.join(dir.out.slyph, 'results.tsv')
    threads: 
        config.resources.med.cpu
    params:
        gtdb_200 = config.sylph.gtdb_200
    resources:
        mem_mb=config.resources.med.mem,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.med.time,
    benchmark:
        os.path.join(dir.out.bench, "sylph", "sylph.benchmark")
    log:
        os.path.join(dir.out.stderr, "sylph", "sylph.log"),
    conda: 
        os.path.join(dir.env, "sylph.yaml")
    shell:
        """
        sylph query {input.reads} {params.gtdb_200} -t {threads}  > results.tsv
        """
