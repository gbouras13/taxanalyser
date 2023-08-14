
rule trimnami:
    input:
        dir = INPUT
    output:
        dir = TRIMNAMI,
        out_trim_flag = os.path.join(FLAGS,"trimnami.done")
    params:
        host_genome = os.path.join(HOST_DIR, "human-t2t-hla.fa")
    resources:
        mem_mb=BigJobMem,
        time=BigTime
    threads:
        1
    shell:
        """
        trimnami run \
            --reads {input.dir} \
            --host {params.host_genome} \
            --minimap map-ont \
            --output {output.dir} \
            nanopore
        touch {output.out_trim_flag}
        """


