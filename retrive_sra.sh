#!/bin/bash

# Input CSV file
input_csv="/Users/siyuan/Library/CloudStorage/OneDrive-LSUHealthShreveport/pan_cell_line/ccle.csv"

# Output CSV file for metadata
output_csv="info.csv"

# Extract the "StrippedCellLineName" column from the CSV file
keywords=($(tail -n +2 "$input_csv" | cut -d ',' -f 4))

# Iterate through each keyword
for keyword in "${keywords[@]}"; do
    # Trim leading and trailing spaces from the keyword
    keyword=$(echo "$keyword" | xargs)

    # Check if the keyword is not empty
    if [ -n "$keyword" ]; then
        echo "Retrieving metadata for keyword: $keyword"
        
        # Use the provided code to retrieve metadata
        if esearch -db sra -query "$keyword" | efetch > metadata.csv; then
            # Append the keyword to each line in the metadata and save to the output CSV
            sed "s/$/,\"$keyword\"/" metadata.csv >> "$output_csv"
            echo "Metadata retrieved successfully for keyword: $keyword"
        else
            echo "Error retrieving metadata for keyword: $keyword. Skipping."
        fi
    fi
done
