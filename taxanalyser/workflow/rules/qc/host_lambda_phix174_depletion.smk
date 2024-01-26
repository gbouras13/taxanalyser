"""
first 2 functions taken from and modified from trimnami  https://github.com/beardymcjohnface/Trimnami/blob/main/trimnami/workflow/rules/hostRemoval.smk

didn't just use trimnami as I don't want snaketool-ception
"""



rule index_host_genome:
    """Pre-index the host genome plus phage lambda plus phage  for mapping with minimap2"""
    input:
        t2t_fasta=config.databases.t2t_fasta
    output:
        index=os.path.join(dir.out.contaminant_index, "host.index"),
    resources:
        mem_mb=config.resources.med.mem,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.med.cpu
    conda:
        os.path.join(dir.env, "contaminants.yaml")
    benchmark:
        os.path.join(dir.out.bench, "index_host_genome.txt")
    log:
        os.path.join(dir.out.stderr, "index_host_genome.log"),
    shell:
        """
        minimap2 -t {threads} -I 8G -d {output.index} {input.t2t_fasta} &> {log}
        """

rule host_removal_mapping_long:
    """Map reads to host and return unmapped reads"""
    input:
        index=os.path.join(dir.out.contaminant_index, "host.index"),
        fastq=get_input_lr_fastqs,
    output:
        r1=os.path.join(
            dir.out.contaminant_removal, "{sample}", "{sample}.host_rm.fastq.gz"
        ),
        s=os.path.join(
            dir.out.contaminant_removal, "{sample}", "{sample}_s.host_rm.fastq.gz"
        ),
        o=os.path.join(
            dir.out.contaminant_removal, "{sample}", "{sample}_o.host_rm.fastq.gz"
        ),
        minimap2_version=os.path.join(dir.out.versions, "{sample}", "minimap2.version"),
        samtools_version=os.path.join(dir.out.versions, "{sample}", "samtools.version"),
    params:
        compression=config.qc.compression,
        minimap_mode=config.qc.minimapModel,
        flagFilt=config.qc.hostRemoveFlagstat,
    benchmark:
        os.path.join(dir.out.bench, "host_removal_mapping_long", "{sample}.txt")
    log:
        mm=os.path.join(dir.out.stderr, "host_removal_mapping_long", "{sample}.minimap.log"),
        sv=os.path.join(
            dir.out.stderr, "host_removal_mapping_long", "{sample}.samtoolsView.log"
        ),
        fq=os.path.join(
            dir.out.stderr, "host_removal_mapping_long", "{sample}.samtoolsFastq.log"
        ),
    resources:
        mem_mb=config.resources.med.mem,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.med.time,
    threads: 
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "contaminants.yaml")
    shell:
        (
            "minimap2 "
            "-ax {params.minimap_mode} "
            "-t {threads} "
            "--secondary=no "
            "{input.index} "
            "{input.fastq} "
            "2> {log.mm} "
            "| samtools view "
            "-h {params.flagFilt} "
            "2> {log.sv} "
            "| samtools fastq "
            "-n -O -c 1 "
            "-o {output.r1} "
            "-0 {output.o} "
            "-s {output.s} "
            "2> {log.fq}; "
            "cat {output.o} {output.s} >> {output.r1}; "
            "minimap2 --version > {output.minimap2_version} "
            "samtools --version  > {output.samtools_version} "
        )


rule  lambda_removal_mapping_long:
    """Map reads to phage lambda and return unmapped reads"""
    input:
        lambda_phage=os.path.join(dir.contaminant_genomes, "lambda.fasta"),
        fastq=os.path.join(
            dir.out.contaminant_removal, "{sample}", "{sample}.host_rm.fastq.gz"
        ),
    output:
        r1=os.path.join(
            dir.out.contaminant_removal, "{sample}", "{sample}.host_lambda_rm.fastq.gz"
        ),
        s=os.path.join(
            dir.out.contaminant_removal, "{sample}", "{sample}_s.host_lambda_rm.fastq.gz"
        ),
        o=os.path.join(
            dir.out.contaminant_removal, "{sample}", "{sample}_o.host_lambda_rm.fastq.gz"
        ),
        minimap2_version=os.path.join(dir.out.versions, "{sample}", "minimap2.version"),
        samtools_version=os.path.join(dir.out.versions, "{sample}", "samtools.version"),
    params:
        compression=config.qc.compression,
        minimap_mode=config.qc.minimapModel,
        flagFilt=config.qc.hostRemoveFlagstat,
    benchmark:
        os.path.join(dir.out.bench, "lambda_removal_mapping_long", "{sample}.txt")
    log:
        mm=os.path.join(dir.out.stderr, "lambda_removal_mapping_long", "{sample}.minimap.log"),
        sv=os.path.join(
            dir.out.stderr, "lambda_removal_mapping_long", "{sample}.samtoolsView.log"
        ),
        fq=os.path.join(
            dir.out.stderr, "lambda_removal_mapping_long", "{sample}.samtoolsFastq.log"
        ),
    resources:
        mem_mb=config.resources.med.mem,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.med.time,
    threads: 
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "contaminants.yaml")
    shell:
        (
            "minimap2 "
            "-ax {params.minimap_mode} "
            "-t {threads} "
            "--secondary=no "
            "{input.phage_lambda} "
            "{input.fastq} "
            "2> {log.mm} "
            "| samtools view "
            "-h {params.flagFilt} "
            "2> {log.sv} "
            "| samtools fastq "
            "-n -O -c 1 "
            "-o {output.r1} "
            "-0 {output.o} "
            "-s {output.s} "
            "2> {log.fq}; "
            "cat {output.o} {output.s} >> {output.r1}; "
            "minimap2 --version > {output.minimap2_version} "
            "samtools --version  > {output.samtools_version} "
        )

        # lambda_phage=os.path.join(dir.contaminant_genomes, "lambda.fasta")
        # phix174_phage=os.path.join(dir.contaminant_genomes, "phix174.fasta")



rule filtlong:
    """
    runs filtlong to filter quality and length
    """
    input:
        fastq=os.path.join(
            dir.out.contaminant_removal, "{sample}", "{sample}.host_lambda_rm.fastq.gz"
        ),
    output:
        fastq=os.path.join(dir.out.qc, "{sample}_filt.fastq.gz"),
        version=os.path.join(dir.out.versions, "{sample}", "filtlong.version"),
    conda:
        os.path.join(dir.env, "filtlong.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.sml.cpu
    params:
        qual=config.config.min_quality,
        length=config.config.min_length,
    benchmark:
        os.path.join(dir.out.bench, "filtlong", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "filtlong", "{sample}.log"),
    shell:
        """
        filtlong --min_mean_q {params.qual} --min_length {params.length} {input.fastq} | pigz > {output.fastq} 2> {log}
        filtlong --version > {output.version}
        """





# rule fastp:
#     """
#     runs fastp on the paired end short reads
#     """
#     input:
#         r1=get_input_r1,
#         r2=get_input_r2,
#     output:
#         r1=os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"),
#         r2=os.path.join(dir.out.fastp, "{sample}_2.fastq.gz"),
#         html=os.path.join(dir.out.fastp, "{sample}.html"),
#         json=os.path.join(dir.out.fastp, "{sample}.json"),
#         version=os.path.join(dir.out.versions, "{sample}", "fastp.version"),
#     conda:
#         os.path.join(dir.env, "fastp.yaml")
#     resources:
#         mem_mb=config.resources.med.mem,
#         mem=str(config.resources.med.mem) + "MB",
#         time=config.resources.med.time,
#     threads: config.resources.sml.cpu
#     benchmark:
#         os.path.join(dir.out.bench, "fastp", "{sample}.txt")
#     log:
#         os.path.join(dir.out.stderr, "fastp", "{sample}.log"),
#     shell:
#         """
#         fastp --in1 {input.r1} --in2 {input.r2} --out1 {output.r1} --out2 {output.r2} --html {output.html} --json {output.json} --thread {threads} 2> {log}
#         fastp --version 2> {output.version}
#         """


rule aggr_long_qc:
    """
    aggregates over all samples for long read qc
    """
    input:
        expand(os.path.join(dir.out.qc, "{sample}_filt.fastq.gz"), sample=SAMPLES),
    output:
        flag=os.path.join(dir.out.flags, "aggr_long_qc.flag"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        touch {output.flag}
        """


# rule aggr_short_qc:
#     """
#     aggregates over all samples
#     """
#     input:
#         expand(os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"), sample=SAMPLES),
#         expand(os.path.join(dir.out.fastp, "{sample}_2.fastq.gz"), sample=SAMPLES),
#     output:
#         flag=os.path.join(dir.out.flags, "aggr_short_qc.flag"),
#     resources:
#         mem_mb=config.resources.sml.mem,
#         mem=str(config.resources.sml.mem) + "MB",
#         time=config.resources.sml.time,
#     threads: config.resources.sml.cpu
#     shell:
#         """
#         touch {output.flag}
#         """
