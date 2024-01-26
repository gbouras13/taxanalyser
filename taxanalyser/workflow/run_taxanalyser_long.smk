import os
import glob
import attrmap as ap
import attrmap.utils as au
from pathlib import Path


# Concatenate Snakemake's own log file with the master log file
# log defined below
def copy_log_file():
    files = glob.glob(os.path.join(".snakemake", "log", "*.snakemake.log"))
    if not files:
        return None
    current_log = max(files, key=os.path.getmtime)
    shell("cat " + current_log + " >> " + LOG)

onsuccess:
    copy_log_file()

onerror:
    copy_log_file()



# config file
configfile: os.path.join(workflow.basedir, "../", "config", "config.yaml")

config = ap.AttrMap(config)

"""
start of pipeline
"""

include: os.path.join("rules", "preflight", "directories.smk")
# functions
include: os.path.join("rules", "preflight", "functions.smk")
# samples
include: os.path.join("rules", "preflight", "samples.smk")





### from config files
#  input as csv
INPUT = config.args.input
OUTPUT = config.args.output
LOG = os.path.join(OUTPUT, "taxanalyser.log")
THREADS = config.args.threads

# Parse the samples and read files
dictReads = parseSamples(INPUT, long_flag=True)  # long flag True
SAMPLES = list(dictReads.keys())

# mmseqs2 dirs
LAMBDA = os.path.join(dir.contaminant_genomes, "lambda.fasta")

GTDB_DIR = config.databases.mmseqs2.gtdb_dir
UNIREF50_DIR = config.databases.mmseqs2.uniref_50
TMPDIR = config.tmpdir



#####
# db search
#####



sourmash_params = config.sourmash
search_databases = config.sourmash.search_databases
KSIZE = config.sourmash.ksize


if not isinstance(KSIZE, list):
    KSIZE=[KSIZE]
for k in KSIZE:
    if k == 21 or k == 31 or k == 51:
        k_str = f"k{k}"
        if k_str not in search_databases.keys():
            raise ValueError(f"Database not specified for search ksize {k_str}. Please specify databases in `config.yaml` file.")
            sys.exit(-1)
    else:
        sys.exit("You have specified a k size for sourmash not 21, 31 or 51. Please check your config file.")


##############################
# Import rules and functions
##############################



# import dir

include: os.path.join("rules", "qc", "host_lambda_phix174_depletion.smk")
include: os.path.join("rules", "mmseqs2", "fastq_to_fasta.smk")
include: os.path.join("rules", "mmseqs2", "mmseqs2_easy_tax.smk")

#include: "rules/sourmash_sketch.smk"
#include: "rules/sourmash_gather.smk"
#include: "rules/sourmash_tax_annotate.smk"
#include: "rules/sourmash_tax_metagenome.smk"

#include: "rules/fastq_to_fasta.smk"
#include: "rules/mmseqs2_easy_tax.smk"


# import targets
include: os.path.join("rules", "preflight", "targets_long.smk")



rule all:
    input:
        TargetFilesLong


