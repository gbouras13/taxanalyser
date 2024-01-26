# Create the output directory if it doesn't exist
mkdir -p c_QCReadsChopped

# Loop through all files ending with .fastq.gz in c_nanoQCInput
for file in c_nanoQCInput/*.fastq.gz; do
    # Extract the filename without the directory and extension
    filename=$(basename "$file" .gz)

    # Gunzip, chopper, and gzip each file
    gunzip -c "$file" | chopper --headcrop 80 --tailcrop 50 | gzip > "c_QCReadsChopped/QCd_$filename.gz"
done
