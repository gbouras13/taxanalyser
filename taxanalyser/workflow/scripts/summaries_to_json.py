#!/usr/bin/env python3
import re
import json
import glob
import os
import sys


def create_json(summaries_dir,  json_file_path):

    # Specify a pattern to match files (e.g., all .txt files)
    pattern = "*.txt"

    # Get a list of files that match the pattern in the directory
    summary_list = glob.glob(os.path.join(summaries_dir, pattern))

    # write all the summary dfs to a list
    summaries = []

    result_dict = {}

    for summary_file in summaries:
        with open(summary_file, 'r') as file:
            content = file.read()

            base_name = os.path.basename(summary_file)
            stripped_name = os.path.splitext(base_name)[0]  # remove extension
            if stripped_name.endswith('_depletion_report'):
                sample = stripped_name[:-len('_depletion_report')]
            else:
                sys.exit("Sample not found.")

  
            barcode_sections = re.split(r'\n\n', content.strip())


            # Process each barcode section
            for section in barcode_sections:
                lines = section.split('\n')

                # Extract summary stats for each barcode
                stats = {}
                for line in lines[1:]:
                    key, value = line.split(':')
                    stats[key.strip()] = value.strip()

                # Add the stats to the result dictionary
                result_dict[sample] = stats

    print(result_dict_sorted)
    
    # Reorganize barcodes by alphanumeric order
    result_dict_sorted = {k: result_dict[k] for k in sorted(result_dict)}

    # Convert the dictionary to JSON format
    json_data = json.dumps(result_dict_sorted, indent=4)

    # Write the JSON data to a new file
    with open(json_file_path, 'w') as json_file:
        json_file.write(json_data)

# to actually run the script
create_json(snakemake.params.summaries_dir, snakemake.output.json_file_path)
