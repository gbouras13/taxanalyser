# metamash
Snakemake and [Snaketool](https://github.com/beardymcjohnface/Snaketool) pipeline to taxonomically profile ONT long read metagenomics.

Snakemake and Snaketool pipeline to taxonomically profile ONT long read metagenomics with sourmash

Much code and ideas have been taken and adapted from the [pb-metagenomics-tools](https://github.com/PacificBiosciences/pb-metagenomics-tools) repository. 

Additionally this pipeline has been substantially informed by [this benchmarking study](https://doi.org/10.1186/s12859-022-05103-0).

Assumes that your reads are in a directory, `fastq.gz` suffixed.


Installation
=========

```
git clone "https://github.com/gbouras13/metamash.git"
cd metamash/
pip install -e .
metamash --help
metamash install --help
metamash run --help
```


Usage
=========

Steps:

1. Download the CHM13 `human-t2t-hla.fa` host genome (from [hostile](https://github.com/bede/hostile)) with the following command.

```
metamash install --database host_genome_db
```

2. Run [trimnami](https://github.com/beardymcjohnface/Trimnami) specifying the directory of FASTQ reads as `--reads` ensuring the mapping is for ONT

```
trimnami run --reads metamash/test_data --host host_genome_db/human-t2t-hla.fa --minimap map-ont --output trimnami_output nanopore
```

You will find the host depleted and trimmed reads in `<output>/nanopore` so here e.g. `trimnami_output/nanopore`.


3. Download the mmseqs2 v 13.45111 uniref50 database - this takes ages in Australia~

e.g.

```
# will write to 'uniref50' directory
mmseqs databases UniRef50 uniref50 tmp --threads 3  
```


4. Run [metamash]()