https://github.com/wdecoster/nanoQC?tab=readme-ov-file

# install nanoQC (conda channel did not work)
pip install nanoQC

# python loop to get the gzipped input files 
python c_nanoQCSetup/filesfornanoQC.py

# usage
nanoQC -h

# to run on an individual file
nanoQC -o c_QCMetrics/barcode01 c_nanoQCInput/barcode01.fastq.gz

# to run on all files
python c_nanoQCSetup/loop_nanoQC.py

# install chopper and make env
conda create -n chopper_env
conda install -n chopper_env python=3.10
conda activate chopper_env
conda install -c bioconda chopper

# to trim first 80 and last 50 bases from the reads from one idividual file
conda activate chopper_env
gunzip -c c_nanoQCInput/barcode02.fastq.gz | chopper --headcrop 80 --tailcrop 50 | gzip > c_QCReadsChopped/QCbarcode02.fastq.gz
conda deactivate

# to trim first 80 and last 50 bases from the reads from multiple files
chmod +x c_nanoQCSetup/loop_chopper.sh #first time only
conda activate chopper_env
./c_nanoQCSetup/loop_chopper.sh
conda deactivate
