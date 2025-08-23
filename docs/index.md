# Staphscope Typing Tool Documentation

Welcome to the Staphscope Typing Tool documentation.

## Overview

Staphscope is a comprehensive typing tool for Staphylococcus aureus that performs:

- Multi-locus sequence typing (MLST)
- spa typing
- SCCmec typing

## Quick Start

```bash
# Install using pip
pip install staphscope

# Or install from GitHub
pip install git+https://github.com/bbeckley-hub/staphscope-typing-tool.git

# Run on your samples
staphscope -i *.fna -o results
