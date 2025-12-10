```markdown
# Usage Guide

## Basic Usage

```bash
# Basic usage with single sample
staphscope -i sample.fna -o results

# Process multiple samples
staphscope -i "*.fasta" -o results

# Use wildcards for multiple files
staphscope -i "*.fna" -o results

# Use multiple threads for faster processing
staphscope -i "*.fna" -o results --threads 8

