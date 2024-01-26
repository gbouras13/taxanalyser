from scipy.stats import entropy
import json
import base64
import numpy as np
import pandas as pd
from bs4 import BeautifulSoup
import json

with open('nanoQC.html', 'r') as f:
    html_data = f.read()

# Parse the HTML content
soup = BeautifulSoup(html_data, 'html.parser')

# Find the <script> tag with the specified id
script_tag = soup.find('script', id='2293')

# Extract the content inside the <script> tag
if script_tag:
    json_content = script_tag.string.strip()
    
    # Load the JSON content
    try:
        json_data = json.loads(json_content)
        print(json_data)
    except json.JSONDecodeError as e:
        print("Error decoding JSON:", e)
else:
    print("Script tag with id '2293' not found.")

name = list(json_data)[0]


def find_indexes_with_data(json_obj, path=[]):
    index_list = []
    if isinstance(json_obj, dict):
        for key, value in json_obj.items():
            new_path = path + [key]
            if key == 'data' and 'attributes' in new_path:
                index_list.append(new_path[3])
            index_list.extend(find_indexes_with_data(value, new_path))
    elif isinstance(json_obj, list):
        for i, item in enumerate(json_obj):
            new_path = path + [i]
            index_list.extend(find_indexes_with_data(item, new_path))
    return index_list

# Assuming json_data is the loaded JSON data
index_list = find_indexes_with_data(json_data)

graph1_indices = []
graph2_indices = []
interim_index_list = []

for index in index_list:
    data = json_data[name]['roots']['references'][index]['attributes']['data']

    if 'y' in data:
        if isinstance(data['y'], dict) and '__ndarray__' in data['y']:
            graph1_indices.append(index)
        else:
            # Check for key error for 'left'
            try:
                if json_data[name]['roots']['references'][index]['attributes']['data']['left']:
                    # If 'left' key is present, do not add to graph 1 or 2
                    pass
            except KeyError:
                interim_index_list.append(index)

    # Now go through the interim_index_list
    if index in interim_index_list:
        # If key error occurs, check if 'y' has values greater than 2
        if any(i > 2 for i in data.get('y', [])):
            # If 'y' has values greater than 2, append to graph2_indices
            pass
        else:
            # If neither criteria met, append to graph2_indices
            graph2_indices.append(index)

# for graph one (start index)
y_1 = []
y_2 = []
y_3 = []
y_4 = []

for index in graph1_indices:
    json_string = json.dumps(json_data[name]['roots']['references'][index]['attributes']['data'], indent=2)
    data = json.loads(json_string)
    y_values = np.frombuffer(base64.b64decode(data["y"]["__ndarray__"]), dtype=data["y"]["dtype"])
    y_values = y_values.reshape(data["y"]["shape"])
    if index == graph1_indices[0]:
        y_1 = y_values
    elif index == graph1_indices[1]:
        y_2 = y_values
    elif index == graph1_indices[2]:
        y_3 = y_values
    elif index == graph1_indices[3]:
        y_4 = y_values

df_g1 = pd.DataFrame({'y_1': y_1, 'y_2': y_2, 'y_3': y_3, 'y_4': y_4})

def calculate_similarity(col1, col2, last_row):
    return abs(last_row[col1].values[0] - last_row[col2].values[0])

last_row_g1 = df_g1.tail(1)
similarities = [(col, calculate_similarity('y_1', col, last_row_g1)) for col in ['y_2', 'y_3', 'y_4']]
min_similarity_col = min(similarities, key=lambda x: x[1])[0]
df_g1.rename(columns={'y_1': 'g1_A1', min_similarity_col: 'g1_A2'}, inplace=True)

for col in df_g1.columns:
    if col not in ['g1_A1', 'g1_A2']:
        df_g1.rename(columns={col: 'g1_B1'}, inplace=True)
        break

for col in df_g1.columns:
    if col not in ['g1_A1', 'g1_A2', 'g1_B1']:
        df_g1.rename(columns={col: 'g1_B2'}, inplace=True)
        break

df_g1['diff_g1A'] = df_g1.apply(lambda row: abs(row['g1_A1'] - row['g1_A2']), axis=1)
df_g1['diff_g1B'] = df_g1.apply(lambda row: abs(row['g1_B1'] - row['g1_B2']), axis=1)

start_index = df_g1[(df_g1['diff_g1A'] < 0.01) & (df_g1['diff_g1B'] < 0.01)].index[0]

# for the graph 2 (tail index)
y_1 = []
y_2 = []
y_3 = []
y_4 = []

for index in graph2_indices:
    json_string = json.dumps(json_data[name]['roots']['references'][index]['attributes']['data'], indent=2)
    data = json.loads(json_string)
    y_values = data['y']
    if index == graph2_indices[0]:
        y_1 = y_values
    elif index == graph2_indices[1]:
        y_2 = y_values
    elif index == graph2_indices[2]:
        y_3 = y_values
    elif index == graph2_indices[3]:
        y_4 = y_values

df_g2 = pd.DataFrame({'y_1': y_1, 'y_2': y_2, 'y_3': y_3, 'y_4': y_4})

def calculate_similarity(col1, col2, last_row):
    return abs(last_row[col1].values[0] - last_row[col2].values[0])

last_row_g2 = df_g2.tail(1)
similarities = [(col, calculate_similarity('y_1', col, last_row_g2)) for col in ['y_2', 'y_3', 'y_4']]
min_similarity_col = min(similarities, key=lambda x: x[1])[0]
df_g2.rename(columns={'y_1': 'g2_A1', min_similarity_col: 'g2_A2'}, inplace=True)

for col in df_g2.columns:
    if col not in ['g2_A1', 'g2_A2']:
        df_g2.rename(columns={col: 'g2_B1'}, inplace=True)
        break

for col in df_g2.columns:
    if col not in ['g2_A1', 'g2_A2', 'g2_B1']:
        df_g2.rename(columns={col: 'g2_B2'}, inplace=True)
        break

df_g2['diff_A'] = df_g2.apply(lambda row: abs(row['g2_A1'] - row['g2_A2']), axis=1)
df_g2['diff_B'] = df_g2.apply(lambda row: abs(row['g2_B1'] - row['g2_B2']), axis=1)

end_index = df_g2[(df_g2['diff_A'] < 0.01) & (df_g2['diff_B'] < 0.01)].index[0]