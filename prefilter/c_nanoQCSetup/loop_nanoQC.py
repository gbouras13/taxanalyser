import os
import subprocess

input_directory = "c_nanoQCInput"
output_directory = "c_QCMetrics"

# Ensure the output directory exists
os.makedirs(output_directory, exist_ok=True)

# Loop through all files in the input directory
for input_file in os.listdir(input_directory):
    if input_file.endswith('.fastq.gz'):
        # Extract barcode from the file name
        barcode = os.path.splitext(os.path.basename(input_file))[0]

        # Run nanoQC command
        nanoQC_command = f'nanoQC -o "{os.path.join(output_directory, barcode)}" "{os.path.join(input_directory, input_file)}"'
        subprocess.run(nanoQC_command, shell=True)

print("Processing completed.")
