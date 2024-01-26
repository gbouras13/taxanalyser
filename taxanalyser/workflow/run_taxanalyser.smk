import glob


configfile: os.path.join(workflow.basedir, '../', 'config', 'config.yaml')


# Concatenate Snakemake's own log file with the master log file
def copy_log_file():
    files = glob.glob(os.path.join(".snakemake", "log", "*.snakemake.log"))
    if not files:
        return None
    current_log = max(files, key=os.path.getmtime)
    shell("cat " + current_log + " >> " + config['log'])

onsuccess:
    copy_log_file()

onerror:
    copy_log_file()



"""
start of pipeline
"""

# snakemake params 
BigJobMem = config["BigJobMem"]
BigJobCpu = config["BigJobCpu"]
SmallJobMem = config["SmallJobMem"]
SmallJobCpu = config["SmallJobCpu"]
SmallTime = config["SmallTime"]
BigTime = config["BigTime"]
MediumTime = config["MediumTime"]

INPUT = config['input']
OUTPUT = config['output']
THREADS = config['threads']

from metasnek import fastq_finder

# https://gist.github.com/beardymcjohnface/bb161ba04ae1042299f48a4849e917c8

# this parses the samples into a dictionary
sample_dict = fastq_finder.parse_samples_to_dictionary(INPUT)
SAMPLES = list(sample_dict.keys())

# get long reads
def get_input_lr_fastqs(wildcards):
    return sample_dict[wildcards.sample]['R1']


#####
# db search
#####

sourmash_params = config['sourmash']
search_databases = sourmash_params['search_databases'] # must be dictionary
KSIZE = sourmash_params.get("ksize", [31, 51])

if not isinstance(KSIZE, list):
    KSIZE=[KSIZE]
for k in KSIZE:
    k_str = f"k{k}"
    if k_str not in search_databases.keys():
        raise ValueError(f"Database not specified for search ksize {k_str}. Please specify databases in `config.yaml` file.")
        sys.exit(-1)

# mmseqs2 dirs
GTDB_DIR = config['gtdb_dir']
UNIREF50_DIR = config['uniref_50']
TMPDIR = config['tmp']

##############################
# Import rules and functions
##############################



# import dir
include: "rules/directories.smk"

#include: "rules/sourmash_sketch.smk"
#include: "rules/sourmash_gather.smk"
#include: "rules/sourmash_tax_annotate.smk"
#include: "rules/sourmash_tax_metagenome.smk"

include: "rules/fastq_to_fasta.smk"
include: "rules/mmseqs2_easy_tax.smk"


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


