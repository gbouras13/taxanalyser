
"""

Taken from Rob's [atavide_lite](https://github.com/linsalrob/atavide_lite/blob/main/slurm/mmseqs_easy_taxonomy.slurm)

"""


rule run_mmseqs_easy_tax_uniref50:
    input: 
        fasta = os.path.join(dir.out.fastas, "{sample}.fasta")
    output:
        outdir=os.path.join(dir.out.mmseqs2, '{sample}'),
        outtouch=os.path.join(dir.out.mmseqs2, 'flags', '{sample}.done')
    params:
        db = os.path.join(config.mmseqs2.uniref_50_dir, 'UniRef50'),
        tmpdir = os.path.join(config.tmpdir, '{sample}_mmseqs2')
    threads: config.resources.big.cpu
    resources:
        mem_mb=config.resources.big.mem,
        mem=str(config.resources.big.mem) + "MB",
        time=config.resources.big.time,
    benchmark:
        os.path.join(dir.out.bench, "mmseqs2", "{sample}.taxonomy.benchmark")
    log:
        os.path.join(dir.out.stderr, "mmseqs2", "{sample}.taxonomy.log"),
    conda: 
        os.path.join(dir.env, "mmseqs2.yaml")
    shell:
        # touch output to let workflow continue in cases where 0 results are found
        """
        mmseqs easy-taxonomy {input.fasta} {params.db} {output.outdir} {params.tmpdir} -s 7 --threads {threads} --orf-filter 0 2>> {log}
        touch {output.outtouch}
        """


rule aggr_mmseqs2_easy_tax:
    input:
        expand(os.path.join(dir.out.mmseqs2, 'flags', '{sample}.done'), sample = SAMPLES)
    output:
        os.path.join(dir.out.flags, "aggr_mmseqs2_tax.flag")
    threads: config.resources.sml.cpu
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time
    shell:
        """
        touch {output[0]}
        """