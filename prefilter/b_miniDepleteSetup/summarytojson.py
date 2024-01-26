# to convert summary_compiled.txt to summary_metrics.json
import re
import json

# Read the input text file
with open('b_readMetrics/summary_compiled.txt', 'r') as file:
    content = file.read()

# Split the content into individual barcode sections
barcode_sections = re.split(r'\n\n', content.strip())

# Initialize an empty dictionary to store the results
result_dict = {}

# Process each barcode section
for section in barcode_sections:
    lines = section.split('\n')
    barcode = lines[0].split()[-1]

    # Extract summary stats for each barcode
    stats = {}
    for line in lines[1:]:
        key, value = line.split(':')
        stats[key.strip()] = value.strip()

    # Add the stats to the result dictionary
    result_dict[barcode] = stats

# Reorganize barcodes by alphanumeric order
result_dict_sorted = {k: result_dict[k] for k in sorted(result_dict)}

# Convert the dictionary to JSON format
json_data = json.dumps(result_dict_sorted, indent=4)

# Write the JSON data to a new file
with open('b_readMetrics/summary_metrics.json', 'w') as json_file:
    json_file.write(json_data)

print("Conversion completed. JSON file saved as 'summary_metrics.json' in 'b_readMetrics' folder.")
