# 🛠️ StaphScope Installation Guide

## System Requirements

- **Operating System**: Linux (Ubuntu 20.04+ recommended) or macOS
- **Python**: Python 3.8 or higher
- **Memory**: 4 GB RAM minimum (8+ GB recommended for large batches)
- **Storage**: 10 GB free space (for AMR database and results)

---

## Install via Conda (Recommended)

```bash
# Create a new conda environment
conda create -n staphscope -c conda-forge -c bioconda staphscope -y
conda activate staphscope

# Verify installation
staphscope --help
