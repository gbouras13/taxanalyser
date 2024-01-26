rule fastq_to_fasta:
    input:
        get_input_lr_fastqs
    output:
        os.path.join(FASTA, "{sample}.fasta")
    threads:
        1
    log: 
        os.path.join(LOGS, "fastq_to_fastq", "{sample}.log")
    benchmark: 
        os.path.join(BENCHMARKS, "fastq_to_fastq", "{sample}.benchmark")
    conda: 
        os.path.join("..", "envs", "biopython.yml")
    resources:
        mem_mb=SmallJobMem,
        time=SmallTime
    script:
        '../scripts/fastq_to_fasta.py'


rule aggr_fastq_to_fasta:
    input:
        expand(os.path.join(FASTA, "{sample}.fasta"), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "aggr_fastq_to_fasta.flag")
    threads:
        1
    resources:
        mem_mb=SmallJobMem,
        time=SmallTime
    shell:
        """
        touch {output[0]}
        """