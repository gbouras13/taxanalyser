rule fastq_to_fasta:
    input:
        fastq=os.path.join(dir.out.qc, "{sample}_filtlong.fastq.gz")
    output:
        fasta=os.path.join(dir.out.fastas, "{sample}.fasta")
    threads: config.resources.sml.cpu
    benchmark:
        os.path.join(dir.out.bench, "fastq_to_fastq", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "fastq_to_fastq", "{sample}.log"),
    conda: 
        os.path.join(dir.env, "scripts.yaml")
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    script:
        '../scripts/fastq_to_fasta.py'


rule aggr_fastq_to_fasta:
    input:
        expand(os.path.join(dir.out.fastas, "{sample}.fasta"), sample = SAMPLES)
    output:
        flag=os.path.join(dir.out.flags, "aggr_fastq_to_fasta.flag"),
    threads: config.resources.sml.cpu
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    shell:
        """
        touch {output[0]}
        """