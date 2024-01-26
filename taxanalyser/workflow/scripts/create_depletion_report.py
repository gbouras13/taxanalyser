#!/usr/bin/env python3

import os
import subprocess as sp
import gzip

def count_fastq_reads(fastq_file):
    try:
        with gzip.open(fastq_file, 'rb') as f:
            line_count = sum(1 for line in f)
    except gzip.BadGzipFile:
        # If it's not a gzipped file, try opening it as a regular text file
        with open(fastq_file, 'r') as f:
            line_count = sum(1 for line in f)
    return line_count // 4  # Each read consists of 4 lines in a Fastq file



def create_depletion_report(sample, pre_depletion_fastq, post_lambda_fastq, post_chm13_fastq, summary_file_path):


    count_pre_depletion = count_fastq_reads(pre_depletion_fastq)
    count_post_lambda_depletion = count_fastq_reads(post_lambda_fastq)
    count_post_chm13_depletion = count_fastq_reads(post_chm13_fastq)


    summary = f"""Summary stats for {sample}:
        Significant reads (not host or lambda): {count_post_chm13_depletion}
        Total reads inc host (no lambda): {count_post_lambda_depletion}
        Proportion of reads mapped to chm13: {'N/A' if count_post_lambda_depletion == 0 else f'{round((count_post_lambda_depletion - count_post_chm13_depletion)/count_post_lambda_depletion*100, 2)}%'}
        Number of reads mapped to chm13 (and depleted): {count_post_lambda_depletion - count_post_chm13_depletion}
        Number of reads mapped to lambda phage (and depleted): {count_pre_depletion - count_post_lambda_depletion}
    """

    with open(summary_file_path, "w") as summary_file:
        summary_file.write(summary)



# to actually run the script
create_depletion_report(snakemake.wildcards.sample, snakemake.input.pre_depletion_fastq, snakemake.input.post_lambda_fastq, snakemake.input.post_chm13_fastq, snakemak.output.summary_file_path)
