#!/usr/bin/env python3
"""
Combine all *_single_split_report.txt files from the reports directory.

This script:
1. Iterates through all *_single_split_report.txt files in reports/
2. Renames the first column to the filename (without the _single_split_report.txt suffix)
3. Combines all files into a single matrix with headers listed only once
4. Writes the combined matrix to reports/combined_single_splits.txt
"""

import pandas as pd
import glob
import os
from pathlib import Path

def combine_single_splits(reports_dir, output_file):
    """
    Combine all single_split report files into one matrix.
    
    Parameters:
    -----------
    reports_dir : str
        Path to the reports directory
    output_file : str
        Path to the output file
    """
    
    # Find all *_single_split_report.txt files
    pattern = os.path.join(reports_dir, "*_single_split_report.txt")
    files = sorted(glob.glob(pattern))
    
    if not files:
        print(f"No single_split files found in {reports_dir}")
        return
    
    print(f"Found {len(files)} single_split files")
    
    combined_df = None
    
    for i, file_path in enumerate(files):
        # Get the filename without path and suffix
        filename = Path(file_path).stem  # removes .txt
        sample_name = filename.replace("_single_split_report", "")
        
        print(f"Processing {i+1}/{len(files)}: {sample_name}")
        
        # Read the file
        df = pd.read_csv(file_path, sep="\t")
        
        # Replace the first column with the sample name for all rows
        first_col = df.columns[0]
        df[first_col] = sample_name
        
        # Rename the first column to "Sample"
        df = df.rename(columns={first_col: "Sample"})
        
        # Combine with previous dataframes by stacking vertically
        if combined_df is None:
            combined_df = df
        else:
            combined_df = pd.concat([combined_df, df], axis=0, ignore_index=False)
    
    # Write the combined matrix to output file
    print(f"\nWriting combined matrix to {output_file}")
    combined_df.to_csv(output_file, sep="\t", index=False)
    print(f"Done! Combined matrix has shape {combined_df.shape}")
    print(f"Columns: {list(combined_df.columns[:10])}... (showing first 10)")

if __name__ == "__main__":
    # Set paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    reports_dir = os.path.join(project_root, "reports")
    output_file = os.path.join(reports_dir, "combined_single_splits.txt")
    
    combine_single_splits(reports_dir, output_file)
