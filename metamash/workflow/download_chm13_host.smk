"""
Snakefile for downloading  
"""
import os

# load default config
configfile: os.path.join(workflow.basedir, '../', 'config', 'config.yaml')

# config

HostDir = config['database']

if not os.path.exists(os.path.join(HostDir)):
    os.makedirs(os.path.join(HostDir))


rule all:
    input:
        os.path.join(HostDir,"human-t2t-hla.fa")

rule get_db:
    """ 
    This can definitely be improved but hack code for now
    """
    params:
        host_db = HostDir
    conda:
        os.path.join( 'envs','gzip.yml')
    output:
        fasta = os.path.join(HostDir,"human-t2t-hla.fa"),
#        fin_flag = os.path.join(HostDir,"download.finished")
    shell:
        """
        cd {params.host_db}
        wget "https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/human-t2t-hla.fa.gz" -O human-t2t-hla.fa.gz
        gunzip human-t2t-hla.fa.gz
        """