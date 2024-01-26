import os
import shutil
import gzip

input_directory = 'b_depletedReads'
output_directory = 'c_nanoQCInput'

# Ensure the output directory exists
os.makedirs(output_directory, exist_ok=True)

# List all files ending with '_chm13.fastq' in the input directory
input_files = [f for f in os.listdir(input_directory) if f.endswith('_chm13.fastq')]

# Process each file
for input_file in input_files:
    input_path = os.path.join(input_directory, input_file)

    # Remove '_chm13' from the filename
    output_file = os.path.splitext(input_file)[0].replace('_chm13', '') + '.fastq.gz'
    output_path = os.path.join(output_directory, output_file)

    # Print the file being converted
    print(f"Converting {input_file} to {output_file}")

    # Compress the file and move it to the output directory
    with open(input_path, 'rb') as f_in, gzip.open(output_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

print("Conversion completed.")
