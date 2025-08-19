# Staphscope Typing Tool Documentation

Welcome to the Staphscope Typing Tool documentation.

## Overview

Staphscope is a comprehensive typing tool for Staphylococcus aureus that performs:

- Multi-locus sequence typing (MLST)
- spa typing
- SCCmec typing

## Quick Start

```bash
# Install
pip install staphscope-typing-tool
github:bbeckley-hub/staphscope-typing-tool
conda install staphscope-typing-tool
sudo apt install staphscope-typing-tool

# Run on your samples
staphscope -i *.fna -o results
