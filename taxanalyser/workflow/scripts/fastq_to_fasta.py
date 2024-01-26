#!/usr/bin/env python3

import gzip
from Bio import SeqIO


def open_fastq_file(input_filename):
    if input_filename.endswith('.gz'):
        return gzip.open(input_filename, 'rt')
    else:
        return open(input_filename, 'r')

def fastq_to_fasta(input_filename, output_filename):
    try:
        with open_fastq_file(input_filename) as input_file, open(output_filename, 'w') as output_file:
            for record in SeqIO.parse(input_file, "fastq"):
                SeqIO.write(record, output_file, "fasta")
    except Exception as e:
        print(f"A parsing error occurred: {e}")

# to actually run the script
fastq_to_fasta(snakemake.input.fastq, snakemake.output.fasta)
