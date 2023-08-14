
"""

Taken from Rob's [atavide_lite](https://github.com/linsalrob/atavide_lite/blob/main/slurm/mmseqs_easy_taxonomy.slurm)

"""


rule run_mmseqs_easy_tax:
    input: 
        fasta = os.path.join(FASTA, "{sample}.fasta")
    output:
        outdir=os.path.join(MMSEQS2, '{sample}')
    params:
        uniref50 = UNIREF50_DB,
        tmpdir = TMPDIR
    threads: 
        THREADS
    resources:
        mem_mb=BigJobMem,
        time=1320,
    log: 
        os.path.join(LOGS, "mmseqs2", "{sample}.taxonomy.log")
    benchmark: 
        os.path.join(BENCHMARKS, "mmseqs2", "{sample}.taxonomy.benchmark")
    conda: 
        os.path.join("..", "envs", "mmseqs2.yml")
    shell:
        # touch output to let workflow continue in cases where 0 results are found
        """
        mmseqs easy-taxonomy {input.fasta} {params.uniref50} {output.outdir} {params.tmpdir} --start-sens 1 --sens-steps 3 -s 7 --threads {threads}
        """


