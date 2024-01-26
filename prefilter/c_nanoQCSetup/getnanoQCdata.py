from bs4 import BeautifulSoup
import re
import json

html_file_path = 'c_QCMetrics/barcode01.fastq/nanoQC.html'
output_file_path = 'c_QCMetrics/barcode01.json'

# Read HTML content from the file
with open(html_file_path, 'r', encoding='utf-8') as html_file:
    html_code = html_file.read()

# Parse HTML code using BeautifulSoup
soup = BeautifulSoup(html_code, 'html.parser')

# Find the script tag with the specified ID
script_tag = soup.find('script', id='2293')

# Extract JSON content from the script tag
json_content = re.search(r'\{.*\}', script_tag.text).group()

# Load JSON data into a Python dictionary
data_dict = json.loads(json_content)

# Write the extracted JSON data to a file
with open(output_file_path, 'w', encoding='utf-8') as output_file:
    json.dump(data_dict, output_file, indent=2)

print(f"JSON data has been exported to {output_file_path}")
