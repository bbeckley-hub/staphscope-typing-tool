## Citation

If you use Staphscope Typing Tool in your research, please cite it as:

Brown, B. (2025). Staphscope Typing Tool: Unified MLST + spa + SCCmec typing for Staphylococcus aureus. GitHub repository. https://github.com/bbeckley-hub/staphscope-typing-tool

# Staphscope Typing Tool

A unified MLST + spa + SCCmec typing tool for *Staphylococcus aureus*.

## Features

- Multi-locus sequence typing (MLST)
- spa typing
- SCCmec typing using CGE SCCmecFinder
- Parallel processing for high-throughput analysis
- Comprehensive reporting in multiple formats



## Installation

### Option 1: Conda (Recommended)
```bash

conda create -n staphscope -c bioconda -c conda-forge -c bbeckley-hub blast staphscope
conda activate staphscope
conda install -c bbeckley-hub staphscope

**### Option 2: From Source**
```bash

git clone https://github.com/bbeckley-hub/staphscope-typing-tool.git
cd staphscope-typing-tool
pip install -e .

**### Option 3 :**

pip install staphscope
**
### Option 4 :**

sudo apt-get staphscope

# Usage
Staphscope is a unified typing tool for Staphylococcus aureus, supporting MLST, spa typing, and SCCmec typing. It also supports parallel processing for faster analysis of multiple genomes.
#Command-Line Interface (CLI)
staphscope -i <input_files> -o <output_dir> [options]
staphscope -i genomes/*.fna -o results --threads 4

# Install system dependencies (Ubuntu/Debian)
sudo apt-get install blastn makeblastdb any2fasta

# Or install BLAST via conda
conda install -c bioconda blast

#System Dependencies

The following system tools must be installed separately:

    BLAST+: blastn and makeblastdb

    any2fasta: For sequence format conversion

    Perl 5+ with JSON module
