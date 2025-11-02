#!/bin/bash

# StaphScope SCCmec Batch Processor
# Usage: ./run_sccmec_batch.sh *.fna

echo "=================================================="
echo "          StaphScope SCCmec Batch Processor"
echo "=================================================="

# Set relative paths based on current directory structure
SCRIPT_DIR="script_dir"
DATABASE_DIR="database" 
DB_CHOICE="reference"

# Check if input files are provided
if [ $# -eq 0 ]; then
    echo "Error: No input files provided"
    echo "Usage: ./run_sccmec_batch.sh *.fna"
    echo "Or:    ./run_sccmec_batch.sh *.fasta"
    exit 1
fi

# Check if required directories exist
if [ ! -d "$SCRIPT_DIR" ]; then
    echo "Error: script_dir not found in current directory"
    exit 1
fi

if [ ! -d "$DATABASE_DIR" ]; then
    echo "Error: database directory not found in current directory"
    exit 1
fi

# Count input files
file_count=$#
echo "Found $file_count input file(s) to process"
echo

# Process each file
counter=0
for file in "$@"; do
    counter=$((counter + 1))
    
    # Get filename without extension and path
    base_name=$(basename "$file")
    sample_name="${base_name%.*}"
    
    echo "--------------------------------------------------"
    echo "[$counter/$file_count] Processing: $base_name"
    echo "Sample name: $sample_name"
    echo "--------------------------------------------------"
    
    # Check if file exists
    if [ ! -f "$file" ]; then
        echo "Error: File $file not found, skipping..."
        continue
    fi
    
    # Create output directory name
    output_dir="s_${sample_name}"
    
    # Run SCCmecFinder
    python3 SCCmecFinder_v4.py \
        -iDb "$file" \
        -iKm "$file" \
        -o "sccmec_detailed_results.txt" \
        -d "$output_dir" \
        -db_dir "$DATABASE_DIR" \
        -sc_dir "$SCRIPT_DIR" \
        -db_choice "$DB_CHOICE"
    
    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "✓ Successfully processed $base_name"
        echo "  Output directory: $output_dir"
    else
        echo "✗ Failed to process $base_name"
    fi
    
    echo
done

echo "=================================================="
echo "           Batch processing completed!"
echo "            Processed $counter files"
echo "=================================================="
