"""
Ensures consistent variable names and file locations for the pipeline.

A lot taken and modified from hecatomb
"""

import attrmap as ap

dir = ap.AttrMap()

### DATABASE LOCATION


### OUTPUT LOCATION
try:
    assert (ap.utils.to_dict(config.args)["output"]) is not None
    dir.out.base = config.args.output
except (KeyError, AssertionError):
    dir.out.base = "taxanalyser_out"


### WORKFLOW DIRs
dir.env = os.path.join(workflow.basedir, "envs")
dir.rules = os.path.join(workflow.basedir, "rules")
dir.scripts = os.path.join(workflow.basedir, "scripts")


### OUTPUT DIRs

dir.out.flags = os.path.join(dir.out.base, "flags")
dir.out.versions = os.path.join(dir.out.base, "versions")
dir.out.processing = os.path.join(dir.out.base, "processing")
dir.out.qc = os.path.join(dir.out.base, "qc")
dir.out.qc_report = os.path.join(dir.out.base, "qc_report")


# logs and benchmarks
dir.out.bench = os.path.join(dir.out.base, "benchmarks")
dir.out.stderr = os.path.join(dir.out.base, "stderr")


# contaminants
dir.contaminant_genomes = os.path.join(workflow.basedir, "../", "contaminant_genomes")
dir.out.contaminant_index = os.path.join(dir.out.processing, "contaminant_index")
dir.out.contaminant_removal = os.path.join(dir.out.processing, "contaminant_removal")

# mmseqs2

dir.out.fastas = os.path.join(dir.out.processing, "fastas")
dir.out.mmseqs2 = os.path.join(dir.out.processing, "mmseqs2")

# sylph 

dir.out.slyph = os.path.join(dir.out.processing, "slyph")