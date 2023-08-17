
"""

Taken from Rob's [atavide_lite](https://github.com/linsalrob/atavide_lite/blob/main/slurm/mmseqs_easy_taxonomy.slurm)

"""


rule run_mmseqs_easy_tax:
    input: 
        fasta = os.path.join(FASTA, "{sample}.fasta")
    output:
        outdir=os.path.join(MMSEQS2, '{sample}'),
        outtouch=os.path.join(MMSEQS2, 'flags', '{sample}.done')
    params:
        gtdb = os.path.join(GTDB_DIR, 'GTDB'),
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
        mmseqs easy-taxonomy {input.fasta} {params.gtdb} {output.outdir} {params.tmpdir} --start-sens 1 --sens-steps 3 -s 7 --threads {threads}
        touch {output.outtouch}
        """


rule aggr_mmseqs2_easy_tax:
    input:
        expand(os.path.join(MMSEQS2, 'flags', '{sample}.done'), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "aggr_mmseqs2.flag")
    threads:
        1
    resources:
        mem_mb=SmallJobMem,
        time=SmallTime
    shell:
        """
        touch {output[0]}
        """