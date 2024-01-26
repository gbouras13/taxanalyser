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
# targets
include: os.path.join("rules", "preflight", "targets_long.smk")




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

print(config)

sourmash_params = config.sourmash
search_databases = config.sourmash.search_databases
KSIZE = sourmash_params.get("ksize", [31, 51])

if not isinstance(KSIZE, list):
    KSIZE=[KSIZE]
for k in KSIZE:
    k_str = f"k{k}"
    if k_str not in search_databases.keys():
        raise ValueError(f"Database not specified for search ksize {k_str}. Please specify databases in `config.yaml` file.")
        sys.exit(-1)



##############################
# Import rules and functions
##############################



# import dir
include: "rules/directories.smk"

include: "rules/qc/remove_contaminants.smk"

#include: "rules/sourmash_sketch.smk"
#include: "rules/sourmash_gather.smk"
#include: "rules/sourmash_tax_annotate.smk"
#include: "rules/sourmash_tax_metagenome.smk"

#include: "rules/fastq_to_fasta.smk"
#include: "rules/mmseqs2_easy_tax.smk"


# import targets
include: "rules/targets.smk"



rule all:
    input:
        TargetFilesMMseqs2




# # Target file
# outTouch = os.path.join(config['output'], config['input'])


# # Mark target rules
# target_rules = []
# def targetRule(fn):
#     assert fn.__name__.startswith('__')
#     target_rules.append(fn.__name__[2:])
#     return fn


# @targetRule
# rule all:
#     input:
#         outTouch


# @targetRule
# rule print_targets:
#     run:
#         print("\nTop level rules are: \n", file=sys.stderr)
#         print("* " + "\n* ".join(target_rules) + "\n\n", file=sys.stderr)


