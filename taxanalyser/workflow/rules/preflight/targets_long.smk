"""
All target output files are declared here
"""

"""
long 
"""

TargetFilesLong = [
    os.path.join(dir.out.flags, "aggr_long_qc.flag"),
    os.path.join(dir.out.flags, "aggr_fastq_to_fasta.flag"),
    #os.path.join(dir.out.flags, "aggr_mmseqs2_tax.flag"),
    os.path.join(dir.out.slyph, 'results.tsv') #sylph
]


# TargetFilesSourmash = [
# os.path.join(FLAGS, "sourmash_sketch.flag"),
# os.path.join(FLAGS, "sourmash_gather.flag"),
# os.path.join(FLAGS, "sourmash_tax_annotate.flag"),
# os.path.join(FLAGS, "sourmash_tax_metagenome.flag")


# ]

# TargetFilesMMseqs2 = [
# os.path.join(FLAGS, "aggr_fastq_to_fasta.flag"),
# os.path.join(FLAGS, "aggr_mmseqs2.flag")

# ]

# #os.path.join(FLAGS, "aggr_fastq_to_fasta.txt")