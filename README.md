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
conda create -n staphscope -c bioconda -c conda-forge blast python=3.8
conda activate staphscope
pip install staphscope
sudo apt-get staphscope
