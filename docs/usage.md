```markdown
# Usage Guide

## Basic Usage

```bash
# Basic usage with single sample
staphscope -i sample.fna -o results

# Process multiple samples
staphscope -i sample1.fna sample2.fna sample3.fna -o results

# Use wildcards for multiple files
staphscope -i *.fna -o results

# Use multiple threads for faster processing
staphscope -i *.fna -o results --threads 8

# Check environment and dependencies
staphscope --check
