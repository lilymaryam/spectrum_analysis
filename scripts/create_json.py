import json
import argparse
import pandas as pd

def parse_metadata_file(metadata_file):
    # Read the metadata file into a pandas DataFrame
    df = pd.read_csv(metadata_file, sep='\t')
    
    # Extract the min and max values for each spectrum column (columns 3 onwards)
    spectrum_columns = df.columns[2:]  # From 3rd column to the end
    min_max_dict = {col: {"min": df[col].min(), "max": df[col].max()} for col in spectrum_columns}
    
    return min_max_dict

def generate_config(min_max_dict):
    # Create the base structure of the config JSON
    config = { "colorRamps": {} }
    
    # Fill in the colorRamps with scales for each spectrum column
    for header_id, min_max in min_max_dict.items():
        min_val = min_max["min"]
        max_val = min_max["max"]
        config["colorRamps"]["meta_"+header_id] = {
            "scale": [
                [min_val, "#FF0000"],  # Red for min
                [max_val, "#0000FF"]   # Blue for max
            ]
        }
    
    return config

def write_json_config(config, output_file):
    # Write the generated config to a JSON file
    with open(output_file, 'w') as f:
        json.dump(config, f, indent=4)

def main():
    parser = argparse.ArgumentParser(description="Process metadata file and generate JSON config.")
    parser.add_argument("--metadata_file", help="Path to the input metadata file")
    parser.add_argument("--output_file", help="Path to save the output JSON config")
    
    args = parser.parse_args()
    
    # Step 1: Parse the metadata file and compute min/max for each spectrum column
    min_max_dict = parse_metadata_file(args.metadata_file)
    
    # Step 2: Generate the config JSON structure
    config = generate_config(min_max_dict)
    
    # Step 3: Write the JSON config to the output file
    write_json_config(config, args.output_file)

if __name__ == "__main__":
    main()