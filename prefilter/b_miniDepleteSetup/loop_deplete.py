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

your_data_directory = 'a_readsConcat'
run_name_list = [os.path.splitext(os.path.basename(file))[0].replace('.fastq', '').replace('.gz', '') for file in os.listdir(your_data_directory) if file.endswith('.fastq.gz')]

for run_name in run_name_list:
    print(f"""
    Running loop_deplete.py on {run_name}:
    """)

    chm13_fasta_filepath = "/Users/issyburdon/Library/CloudStorage/Box-Box/meta_pipe_box/x_contaminants/chm13_ref/GCF_009914755.1_T2T-CHM13v2.0_genomic.fasta"
    lambda_fasta_filepath = "/Users/issyburdon/Library/CloudStorage/Box-Box/meta_pipe_box/x_contaminants/phage_lambda.fasta"
    fastq_filepath = f"a_readsConcat/{run_name}.fastq.gz"
    output_filepath = f"b_depletedReads/{run_name}"
    metrics_filepath = f"b_readMetrics/"

    # Run minimap2 and samtools commands to deplete lambda reads
    command_map_lambda = f"minimap2 -t 4 -ax map-ont {lambda_fasta_filepath} {fastq_filepath} | samtools view -f 4 -b -o {output_filepath}_lambda.bam && samtools fastq -f 4 -0 {output_filepath}_lambda.fastq {output_filepath}_lambda.bam"

    print(f"""
            Depleting phage lambda reads... 
            """)
        
    try:
        lambda_process = sp.run(command_map_lambda, shell=True, check=True, stderr=sp.PIPE, text=True)
        lambda_minimap2_report = lambda_process.stderr
        print(f"""
            Phage lambda reads depleted. 
            Moving on to depletion of chm13 reads...
            """)
    except sp.CalledProcessError:
        print("An error occurred during processing. Please check your inputs and dependencies.")

    command_map_chm13 = f"minimap2 -t 4 -ax map-ont {chm13_fasta_filepath} {output_filepath}_lambda.fastq | samtools view -f 4 -b -o {output_filepath}_chm13.bam && samtools fastq -f 4 -0 {output_filepath}_chm13.fastq {output_filepath}_chm13.bam"
        
    try:
        chm13_process = sp.run(command_map_chm13, shell=True, check=True, stderr=sp.PIPE, text=True)
        chm13_minimap2_report = chm13_process.stderr
        print(f"""
            chm13 reads depleted, processing complete. 
            Writing summary...""")
    except sp.CalledProcessError:
        print("An error occurred during processing. Please check your inputs and dependencies.")

    # write a summary report
    count_pre_depletion = count_fastq_reads(fastq_filepath)
    count_post_lambda_depletion = count_fastq_reads(f"{output_filepath}_lambda.fastq")
    count_post_chm13_depletion = count_fastq_reads(f"{output_filepath}_chm13.fastq")

    summary = f"""Summary stats for {run_name}:
        Significant reads (not host or lambda): {count_post_chm13_depletion}
        Total reads inc host (no lambda): {count_post_lambda_depletion}
        Proportion of reads mapped to chm13: {'N/A' if count_post_lambda_depletion == 0 else f'{round((count_post_lambda_depletion - count_post_chm13_depletion)/count_post_lambda_depletion*100, 2)}%'}
        Number of reads mapped to chm13 (and depleted): {count_post_lambda_depletion - count_post_chm13_depletion}
        Number of reads mapped to lambda phage (and depleted): {count_pre_depletion - count_post_lambda_depletion}

        Paths submitted to program:
            chm13 fasta reference genome path = {chm13_fasta_filepath}
            lambda phage fasta reference genome path ={lambda_fasta_filepath}
            Fastq reads file path = {fastq_filepath}
            Output file path = {output_filepath}

    Lambda depletion minimap2 generated report:
    {lambda_minimap2_report}

    chm13 depletion minimap2 generated report:
    {chm13_minimap2_report}
    """

    print(summary)

    # last_slash_index = output_filepath.rfind("/")

    # if last_slash_index != -1:
    #     # Remove the last portion of the string
    #     output_filepath_cd = output_filepath[:last_slash_index]

    summary_file_path = f"{metrics_filepath}/{run_name}_summary.txt"
    with open(summary_file_path, "w") as summary_file:
        summary_file.write(summary)

    print(f"""  Summary txt file written.
        Saved to: {metrics_filepath}/{run_name}_summary.txt""")

    # Append summary to summary_compiled.txt
    summary_path = f"{metrics_filepath}/{run_name}_summary.txt"
    if os.path.exists(summary_path):
        with open("b_readMetrics/summary_compiled.txt", "a") as summary_compiled:
            with open(summary_path, "r") as summary:
                first_7_lines = [next(summary) for _ in range(7)]
                summary_compiled.writelines(first_7_lines)
    else:
        print(f"""
              Summary file not found for {run_name}""")

    print(f"""
          
    **** Depletion of {run_name} complete! ****

    """)