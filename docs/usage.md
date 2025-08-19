#For questions and support, please open an issue on GitHub or contact brownbeckley94@gmail.com


## 7. docs/usage.md
```markdown
# Usage Guide

# Use multiple threads
staphscope -i *.fna -o results --threads 8

# Check environment
staphscope --check

# Update databases
staphscope --update
## Basic Usage

```bash
staphscope -i sample1.fna sample2.fna -o results


## 8. scripts/install_dependencies.sh
```bash
#!/bin/bash
# install_dependencies.sh

echo "Installing Staphscope Typing Tool dependencies..."

# Check if conda is available
if command -v conda &> /dev/null; then
    echo "Installing with conda..."
    conda create -n staphscope -c bioconda -c conda-forge mlst spatyper blast python=3.8 biopython pandas
    conda activate staphscope
    pip install staphscope-typing-tool
else
    echo "Conda not found. Attempting to install with pip..."
    pip install biopython pandas
    
    # Check for required tools
    for cmd in mlst spatyper blastn makeblastdb; do
        if ! command -v $cmd &> /dev/null; then
            echo "Warning: $cmd not found. Please install it manually."
        fi
    done
    
    pip install staphscope-typing-tool
fi

echo "Installation complete. Run 'staphscope --help' for usage information."
