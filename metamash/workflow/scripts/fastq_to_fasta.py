#!/usr/bin/env python3

import gzip

def fastq_to_fasta(input_filename, output_filename):
    if input_filename.endswith('.gz'):
        open_func = gzip.open
    else:
        open_func = open

    with open_func(input_filename, 'rt') as input_file, open(output_filename, 'w') as output_file:
        line_number = 0
        sequence = ""

        for line in input_file:
            line = line.strip()

            if line_number % 4 == 1:
                sequence = line
                output_file.write(f'>{sequence}\n')
            elif line_number % 4 == 0:
                line_length = len(sequence)

                for i in range(0, line_length, 60):
                    output_file.write(sequence[i:i + 60] + '\n')

            line_number += 1

# to actually run the script
fastq_to_fasta(snakemake.input[0], snakemake.output[0])
