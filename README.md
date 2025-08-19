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

### Using Conda (Recommended)
```bash
conda install -c bioconda staphscope-typing-tool

pip install staphscope-typing-tool

sudo apt install staphscope-typing-tool
# Basic Commands
staphscope -i /path/to/genomes/*.fasta -o results/
# Enable Multithreading
staphscope -i /path/to/genomes/*.fasta -o results/ --threads 4
# Check Environment
staphscope --check
# Update Databases
staphscope --update
# Display Help
staphscope --help
# Output
Organized output directory with MLST, spa, and SCCmec typing results.

Supports multiple output formats (CSV, JSON).

Includes logs for traceability of results.
# Dependencies
Python 3.6+

pandas, numpy, biopython, tqdm

System tools: mlst, spaTyper, blastn, makeblastdb
#Notes

High-throughput genome analysis is recommended with --threads set according to available CPU cores.

For SCCmec typing, the CGE SCCmecFinder database must be downloaded and configured separately.
