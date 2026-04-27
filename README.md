
```markdown
███████╗████████╗ █████╗ ██████╗ ██╗  ██╗███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔════╝╚══██╔══╝██╔══██╗██╔══██╗██║  ██║██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
███████╗   ██║   ███████║██████╔╝███████║███████╗██║     ██║   ██║██████╔╝█████╗  
╚════██║   ██║   ██╔══██║██╔═══╝ ██╔══██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
███████║   ██║   ██║  ██║██║     ██║  ██║███████║╚██████╗╚██████╔╝██║     ███████╗
╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝     ╚═╝  ╚═╝╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝
```  

<div align="center">

# 🔬 StaphScope

### **A species-optimized computational pipeline for rapid, accessible *Staphylococcus aureus* genotyping and surveillance**

#### **Complete MRSA/MSSA genomic analysis in minutes — not hours**

![Version](https://anaconda.org/bbeckley-hub/staphscope/badges/version.svg)
![Latest Release Date](https://anaconda.org/bbeckley-hub/staphscope/badges/latest_release_date.svg)
![Platforms](https://anaconda.org/bbeckley-hub/staphscope/badges/platforms.svg)
![License](https://anaconda.org/bbeckley-hub/staphscope/badges/license.svg)
![Downloads](https://anaconda.org/bbeckley-hub/staphscope/badges/downloads.svg)
[![DOI](https://img.shields.io/badge/DOI-10.1186%2Fs12864--026--12609--x-blue)](https://doi.org/10.1186/s12864-026-12609-x)



[![Docker Pulls](https://img.shields.io/docker/pulls/bbeckleyhub/staphscope)](https://hub.docker.com/r/bbeckleyhub/staphscope)
[![Docker Image Size](https://img.shields.io/docker/image-size/bbeckleyhub/staphscope/latest)](https://hub.docker.com/r/bbeckleyhub/staphscope)
[![Docker Version](https://img.shields.io/docker/v/bbeckleyhub/staphscope?sort=semver)](https://hub.docker.com/r/bbeckleyhub/staphscope)
[![Contributions Welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg)](#)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-Profile-0A66C2?style=flat&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/brown-beckley-190315319)
[![Stage](https://img.shields.io/badge/status-active-brightgreen)](#)
![Conda Downloads](https://img.shields.io/conda/dn/bioconda/staphscope?label=Conda%20Downloads)

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Conda](https://img.shields.io/badge/conda-✓-green.svg)](https://docs.conda.io/en/latest/)
[![MIT License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![GitHub Issues](https://img.shields.io/github/issues/bbeckley-hub/staphscope-typing-tool)](https://github.com/bbeckley-hub/staphscope-typing-tool/issues)
[![GitHub Stars](https://img.shields.io/github/stars/bbeckley-hub/staphscope-typing-tool)](https://github.com/bbeckley-hub/staphscope-typing-tool/stargazers)
[![Sample Report](https://img.shields.io/badge/📊-View_Sample_Report-blue)](https://htmlpreview.github.io/?https://github.com/bbeckley-hub/staphscope-typing-tool/blob/main/staphscope_ultimate_report.html)
![Profile Views](https://komarev.com/ghpvc/?username=bbeckley-hub&label=Profile%20Views&color=0e75b6&style=flat)
[![Google Scholar](https://img.shields.io/badge/Google%20Scholar-Profile-4285F4?style=flat&logo=googlescholar&logoColor=white)](https://scholar.google.com/citations?user=CYNOsqIAAAAJ&hl=en)



![GitHub stats](https://github-readme-stats.vercel.app/api?username=bbeckley-hub&show_icons=true&theme=radical)
![Top Langs](https://github-readme-stats.vercel.app/api/top-langs/?username=bbeckley-hub&layout=compact&theme=radical)
[![GitHub Streak](https://streak-stats.demolab.com?user=bbeckley-hub&theme=radical&date_format=j%20M%5B%20Y%5D)](https://git.io/streak-stats)


**Two ways to use StaphScope:**  
🖥️ **Command-line tool** for high-throughput, local analysis  
🌐 **StaphScope Web** for non-bioinformaticians – [https://staphscope.dpdns.org](https://staphscope.dpdns.org)

</div>

---

## 📋 **Table of Contents**

- [🎯 Overview](#-overview)
- [✨ Key Features](#-key-features)
- [🌐 StaphScope Web Platform](#-staphscope-web-platform)
- [⚡ Quick Start (CLI)](#-quick-start-cli)
- [🔧 Installation (CLI)](#-installation-cli)
- [🐳 Staphscope Docker Usage](#-Staphscope-docker-usage)
- [🔗 Integrated External Tools & Dependencies](#-integrated-external-tools--dependencies)
- [🚀 Usage Guide (CLI)](#-usage-guide-cli)
- [📁 Output Structure](#-output-structure)
- [🔍 Analytical Modules](#-analytical-modules)
- [📈 Performance Benchmarks](#-performance-benchmarks)
- [🔬 Validation & Accuracy](#-validation--accuracy)
- [🆚 Tool Comparison](#-tool-comparison)
- [🤖 AI Integration Guide](#-ai-integration-guide)
- [🔮 Future Development](#-future-development)
- [❓ Frequently Asked Questions](#-frequently-asked-questions)
- [🐛 Troubleshooting](#-troubleshooting)
- [📚 Citation](#-citation)
- [🙏 Acknowledgements](#-acknowledgements)
- [👥 Authors & Contact](#-authors--contact)
- [📄 License](#-license)
- [📚 Third-Party Tool Citations](#-Third-Party-Tool-Citations)

---

## 🎯 **Overview**

**StaphScope** is an automated, locally-executable computational pipeline designed specifically for comprehensive *Staphylococcus aureus* genomic surveillance. It addresses the critical bottleneck in MRSA (Methicillin-Resistant *S. aureus*) research by integrating **six essential genotyping methods** into a single, cohesive workflow.

### 🌍 **The Problem**
- **Fragmented Bioinformatics**: Traditional MRSA analysis requires 5+ separate tools with conflicting dependencies
- **Resource Barriers**: Web-based services need constant internet and raise data privacy concerns
- **Time Constraints**: Generalist platforms take hours; outbreaks need answers in minutes
- **Interpretation Challenges**: Raw data without epidemiological context limits actionable insights

### 💡 **Our Solution**
StaphScope delivers:
- **✅ Single-command installation** via Conda
- **✅ 10-14 minute complete analysis** (24 samples, 16 cores)
- **✅ 100% local execution** with data privacy
- **✅ Intelligent resource management** using Python's psutil library
- **✅ Interactive HTML reports** with epidemiological context
- **✅ Automated MRSA/MSSA classification** with confidence scoring
- **✅ Web-based interface** for non-bioinformaticians

**Perfect for**: Clinical labs, outbreak investigations, research studies, and public health surveillance.

---

## ✨ **Key Features**

### 🔬 **Core Analytical Modules**

| Module | 🎯 Purpose | 📊 Key Outputs | ⚡ Speed |
|--------|------------|----------------|----------|
| **FASTA QC** | Comprehensive quality control (N50/N70/N90, GC%, contig stats) | HTML, TSV, JSON reports with visual summaries | <30 sec |
| **MLST Typing** | Phylogenetic classification via 7 housekeeping genes | ST, CC, allele profiles, epidemiological context | <1 min |
| ***spa* Typing** | Hypervariable region analysis of protein A gene | *spa* type, repeat patterns, alignment metrics | <1 min |
| **SCC*mec* Typing** | Methicillin resistance cassette characterization | SCC*mec* type (I-XIII), *mec*/*ccr* complexes, confidence scores | 1-2 min |
| **AMR Profiling** | Comprehensive resistance gene detection (AMRFinderPlus) | 5,000+ AMR genes, risk categorization, cross-sample patterns | 2-3 min |
| **ABRicate Screening** | Multi-database virulence/plasmid detection (9 databases) | Plasmid replicons, virulence factors, clinical flags | 3-4 min |
| **Visualization Suite** | Publication-ready graphics using seaborn, plotly, matplotlib | 14+ graph types in PDF, PNG, SVG, interactive HTML | 1-2 min |
| **Lineage Database** | Global epidemiological context | 44 major lineages, geographical distribution, outbreak potential | Instant |

### 🛡️ **MRSA-Specific Innovations**
- **Automated MRSA Classification**: Based on concurrent *mecA/mecC* + SCC*mec* detection
- **Clinical Gene Flagging**: Automatic highlighting of PVL, enterotoxins, *van* genes
- **Risk Assessment**: Categorizes genes as 'Critical Risk' (e.g., *mecA*, *vanA*) or 'High Risk'
- **Cross-Genome Pattern Discovery**: Summarizes gene frequencies across entire sample sets
- **Curated Lineage Database**: 44 major lineages with HA-MRSA, CA-MRSA, LA-MRSA classifications

### 🚀 **Performance Advantages**
- **8-10× faster** than Bactopia for *S. aureus*-specific analyses
- **Linear scaling** with sample numbers (R² = 0.931)
- **Dynamic resource allocation** using Python psutil
- **Low memory footprint**: Runs on 4GB RAM, scales to HPC clusters

## 📊 Sample Output

See a complete interactive report generated by StaphScope:

[![Sample Report](https://img.shields.io/badge/📊-View_Sample_Report-blue)](https://htmlpreview.github.io/?https://github.com/bbeckley-hub/staphscope-typing-tool/blob/main/staphscope_ultimate_report.html)

*The report includes AMR and virulence gene tables, filter buttons, combination tables, and FASTA QC metrics.*

---

## 🌐 **StaphScope Web Platform**

For researchers and clinicians who prefer a graphical interface, **StaphScope Web** provides all the power of the command-line tool in an easy-to-use web application.

### **Key Web Features**
- ✅ **Drag-and-drop file upload** (single, multiple, or ZIP archives)
- ✅ **Module selection** – choose which analyses to run
- ✅ **Real-time progress monitoring** with live logs
- ✅ **Beautiful HTML reports** with interactive visualizations
- ✅ **Download all results as a single ZIP** file
- ✅ **Responsive design** – works on desktop and tablet
- ✅ **No installation required** – works in any modern browser

### **Technology Stack**
- **Backend**: Flask (Python web framework)
- **Task Queue**: Celery with Redis broker
- **Bioinformatics Engine**: StaphScope CLI (via Conda)
- **Frontend**: Bootstrap 5, JavaScript
- **Deployment**: Gunicorn + Nginx

### **Quick Access**
> 🌐 **Try StaphScope Web today:** [https://staphscope.dpdns.org](https://staphscope.dpdns.org)  
> 📦 **Web Repository:** [https://github.com/bbeckley-hub/staphscope-web](https://github.com/bbeckley-hub/staphscope-web)

*Note: The web version limits uploads to 10 files per job for fair resource usage. For larger datasets, please use the command-line tool.*

*Note: Currently hosted on personal infrastructure; availability may vary as we work toward sustainable 24/7 hosting.*

---

## ⚡ **Quick Start (CLI)**

### **Install in 60 seconds**
```bash
# Method 1: Conda (Recommended)
conda create -n staphscope -c conda-forge -c bioconda staphscope -y
conda activate staphscope

# Method 2: Mamba (Faster installation)
mamba create -n staphscope -c conda-forge -c bioconda staphscope -y
mamba activate staphscope

# Method 3: From source (Needs external databases)
git clone https://github.com/bbeckley-hub/staphscope-typing-tool.git
cd staphscope-typing-tool
conda env create -f environment.yml
conda activate staphscope
pip install -e .
```

### **Run your first analysis**
```bash
# Single genome
staphscope -i genome.fasta -o results/

# Batch processing (24 genomes)
staphscope -i "*.fna" -o batch_results --threads 16
# Complete in ~14 minutes! 🎉
```

---

## STAPHSCOPE TERMINAL DISPLAY

```
usage: staphscope [-h] -i INPUT -o OUTPUT [-t THREADS] [--skip-amr]
                  [--skip-abricate] [--skip-mlst] [--skip-spa] [--skip-sccmec]
                  [--skip-lineage] [--skip-comprehensive]

StaphScope: Complete S. aureus Typing Pipeline

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input FASTA file(s) - can use glob patterns like
                        "*.fna" or "*.fasta"
  -o OUTPUT, --output OUTPUT
                        Output directory for all results
  -t THREADS, --threads THREADS
                        Number of threads (default: 2)
  --skip-amr            Skip AMR analysis (AMRfinderPlus)
  --skip-abricate       Skip ABRicate analysis
  --skip-mlst           Skip MLST analysis
  --skip-spa            Skip spa typing analysis
  --skip-sccmec         Skip SCCmec analysis
  --skip-lineage        Skip lineage reference generation
  --skip-comprehensive  Skip comprehensive report generation (MLST + spa +
                        SCCmec)

Examples:
  staphscope -i genome.fna -o results/
  staphscope -i "*.fna" -o batch_results --threads 8
  staphscope -i "*.fasta" -o analysis --threads 16 --skip-lineage
  staphscope -i "genome*.fa" -o results/ --threads 4 --skip-comprehensive

Supported FASTA formats: .fna, .fasta, .fa, .fn

Analysis Modules:
  • MLST (Multi-Locus Sequence Typing)
  • spa typing (Staphylococcal Protein A)  
  • SCCmec typing (Methicillin Resistance Cassette)
  • AMR profiling (Antimicrobial Resistance)
  • ABRicate (Comprehensive resistance/Plasmid/virulence)
  • Lineage reference database
  • Comprehensive report (MLST + spa + SCCmec summary)

Output: Comprehensive results for all analyses in organized directories
Please run abricate --setupdb for recent gene annotations!!!
⭐ Star us on GitHub if you find this tool useful!

Transforming fragmented genomic data into coherent biological narratives 🧬✨
```

---

## 🔧 **Installation (CLI)**

### **System Requirements**
| Resource | Minimum | Recommended | Production |
|----------|---------|-------------|------------|
| **CPU Cores** | 2 | 8+ | 16+ |
| **RAM** | 4 GB | 8 GB | 16 GB |
| **Storage** | 2 GB | 10 GB | 50 GB+ |
| **OS** | Linux, macOS, WSL2 | Linux | Linux Cluster |

### **Step-by-Step Installation**

#### **1. Install Miniconda (if needed)**
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

#### **2. Install StaphScope**
```bash
# Add channels in correct order
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels bbeckley-hub

# Create and activate environment
conda create -n staphscope python=3.9 staphscope -y
conda activate staphscope

# Verify installation
staphscope --help
```

#### **3. Update Databases (Recommended)**
```bash
abricate --setupdb
```

---

## 🐳 **Staphscope Docker Usage**


## 📦 Quick Start

```bash
# Pull the latest image
docker pull bbeckleyhub/staphscope:latest

# Test installation
docker run --rm bbeckleyhub/staphscope:latest --help

# Analyze your data
docker run --rm \
  -v $(pwd)/genomes:/data/input \
  -v $(pwd)/results:/data/output \
  bbeckleyhub/staphscope:latest \
  -i "*.fasta" -o /data/output -t 4

# Outputs
sudo chown -R $USER:$USER ./output

```

## 🖥️ Singularity for HPC (no `sudo`, correct ownership)

On HPC clusters that support [Singularity/Apptainer](https://sylabs.io/singularity/), you can run StaphScope **without `sudo`** and output files will be owned by your user automatically.

> **Important:** StaphScope writes temporary files inside its own installation directory (e.g., `/opt/staphscope/...`). Singularity mounts containers as read‑only by default, so you **must** add the `--writable-tmpfs` flag to allow these writes. The flag creates an ephemeral, writable overlay in memory – no permanent changes are made to the container.

### Option A: Direct pull (if network allows)

```bash
singularity pull staphscope.sif docker://bbeckleyhub/staphscope:latest
singularity run --writable-tmpfs -B $(pwd):/data staphscope.sif -i "/data/*.fasta" -o /data/output
```

### Option B: Convert from a local Docker image (when `singularity pull` fails)

If you encounter TLS timeouts or other network errors (common on some HPCs), convert an existing Docker image to a Singularity SIF file on a machine with Docker, then transfer the `.sif` file to the HPC.

**Step 1 – on a machine with Docker (e.g., your laptop):**

```bash
docker pull bbeckleyhub/staphscope:latest
docker save bbeckleyhub/staphscope:latest -o staphscope.tar
singularity build staphscope.sif docker-archive://staphscope.tar
```

Now copy `staphscope.sif` to your HPC home or project directory (e.g., using `scp`).

**Step 2 – on the HPC (no sudo needed):**

```bash
singularity run --writable-tmpfs -B $(pwd):/data staphscope.sif -i "/data/*.fasta" -o /data/output
```

### Explanation of flags

| Flag | Purpose |
|------|---------|
| `--writable-tmpfs` | Creates a temporary writable overlay – **required** for StaphScope to write intermediate files to `/opt/...` |
| `-B $(pwd):/data` | Binds your current directory to `/data` inside the container (input files are read from here, output is written here) |
| `-i "/data/*.fasta"` | Input pattern – use quotes to prevent shell expansion on the host |
| `-o /data/output` | Output directory (will appear as `./output` on your host) |

### Additional options

You can use any StaphScope flag, e.g.:

```bash
singularity run --writable-tmpfs -B $(pwd):/data staphscope.sif \
    -i "/data/*.fasta" -o /data/output --threads 8 --skip-amr
```

### Verify it works

After a successful run, you will see output indicating each module completed. All result files in `./output` will be owned by **your HPC user** – no `sudo chown` needed.

---

## 🔗 **Integrated External Tools & Dependencies**

StaphScope integrates several powerful open-source tools and databases. These are **not bundled directly in this repository**. Instead, they are automatically installed as **dependencies via Conda** (as defined in `environment.yml`). The MIT license that applies to the StaphScope pipeline code does not cover these external tools. Each tool is used under the terms of its own license, and we gratefully acknowledge their authors.

| Tool/Database | Purpose | Source | License |
|---------------|---------|--------|---------|
| **MLST** | Multi-locus sequence typing | [tseemann/mlst](https://github.com/tseemann/mlst) | GPL v2 |
| **ABRicate** | Mass screening for resistance/virulence | [tseemann/abricate](https://github.com/tseemann/abricate) | GPL v2 |
| **AMRFinderPlus** | Antimicrobial resistance gene detection | [ncbi/amr](https://github.com/ncbi/amr) | Public Domain |
| **SCCmecFinder** | SCCmec typing | [genomicepidemiology/Sccmecfinder](https://bitbucket.org/genomicepidemiology/Sccmecfinder) | Apache-2.0 |
| **spa typing** | *spa* gene typing | [spa.ridom.de](https://spa.ridom.de/) | Free for academic use |
| **PubMedST** | MLST allele database | [pubmlst.org](https://pubmlst.org/organisms/staphylococcus-aureus) | Open access for research |

---

## 🚀 **Usage Guide (CLI)**

### **Basic Commands**
```bash
# Single genome
staphscope -i genome.fasta -o results/

# Batch processing with wildcards
staphscope -i "*.fna" -o results_2025 --threads 8

# Skip specific modules
staphscope -i sample.fna -o results --skip-spa --skip-lineage
```

### **Input Formats**
- Accepted: `.fna`, `.fasta`, `.fa`, `.fn`
- Required: Assembled genomes (contigs or complete)
- Batch patterns: `*.fasta`, `sample_*.fna`, etc.

### **Real-World Examples**

#### **Clinical Laboratory Setting**
```bash
# Daily surveillance of 12 isolates
staphscope -i "daily_isolates/*.fasta" -o /mnt/shared/surveillance/$(date +%Y%m%d) --threads 12
# Complete in ~8 minutes
```

#### **Outbreak Response**
```bash
# Urgent investigation (8 suspected cases)
staphscope -i "outbreak/*.fasta" -o /tmp/urgent_analysis --skip-lineage
# Results in ~4 minutes
```

---

## 📁 **Output Structure**

StaphScope generates a comprehensive, organized output directory:

```
Staphscope/
├── abricate_results/          # Multi-database screening (9 DBs)
├── amr_results/               # AMR gene profiling (AMRFinder+)
├── mlst_results/              # MLST typing
├── sccmec_results/            # SCCmec typing
├── spa_results/               # spa typing
├── lineage_results/           # Phylogenetic lineage
├── qc_results/                # FASTA quality control
├── visualization_results/     # Publication-ready plots
└── Staphscope_final_report/   # Consolidated reports (HTML/JSON/TSV)
```

Each module contains:
- **Per-sample directories** with raw outputs
- **Summary files** (TSV/JSON) for cross-sample analysis
- **Interactive HTML reports** for visualization
- **Master reports** combining all results

---

## 🔍 **Analytical Modules**

### **1. FASTA QC**
- **Metrics**: N50/N70/N90, L50/L70/L90, GC content, total length, contig count
- **Outputs**: HTML reports with histograms, TSV/JSON for downstream analysis

### **2. MLST Typing**
- **Database**: PubMedST *S. aureus*
- **Method**: BLAST-based allele calling
- **Output**: ST, CC, 7-gene profile, epidemiological context

### **3. *spa* Typing**
- **Database**: Ridom *spa* repeat database
- **Method**: BLAST against repeat sequences
- **Output**: *spa* type, repeat pattern, alignment metrics

### **4. SCC*mec* Typing**
- **Method**: Hierarchical two-method system (gene-based + k-mer homology)
- **Output**: SCC*mec* type (I-XIII), confidence scores, *mec*/*ccr* complexes
- **Subtyping**: Types IV and V community-associated cassettes

### **5. AMR Profiling**
- **Tool**: NCBI-AMRFinderPlus v4.2.4
- **Coverage**: 5,000+ AMR genes
- **Risk Assessment**: Critical Risk (*mecA*, *vanA*, *cfr*), High Risk (*erm*, *tetM*)

### **6. ABRicate Screening**
- **Databases**: VFDB, ResFinder, CARD, PlasmidFinder, MegaRes, NCBI, ARG-ANNOT, ECOH, EcoLi_VF
- **Thresholds**: ≥80% identity and coverage
- **Clinical Flags**: PVL, enterotoxins, *van* genes

### **7. Visualization Suite**
- **Libraries**: seaborn, plotly, matplotlib
- **Plot Types**: Box plots, violin plots, bar charts, heatmaps, correlation matrices, pie charts, line graphs
- **Formats**: PNG, SVG, PDF, interactive HTML

### **8. Lineage Database**
- **Content**: 44 major *S. aureus* lineages (18 HA-MRSA, 19 CA-MRSA, 7 LA-MRSA)
- **Metadata**: Geographical distribution, clinical significance, outbreak potential

---

## 📈 **Performance Benchmarks**

| System | Samples | Time | Speed vs Bactopia |
|--------|---------|------|-------------------|
| Laptop (2 cores, 8GB) | 1 | 2m 33s | 5× faster |
| Laptop (2 cores, 8GB) | 24 | 28m 17s | 6× faster |
| Workstation (16 cores, 16GB) | 1 | 1m 31s | 8× faster |
| Workstation (16 cores, 16GB) | 24 | 14m 34s | 10× faster |
| Workstation (16 cores, 16GB) | 100 | ~60m | 12× faster |

### **Resource Efficiency**
- **Memory Usage**: 2-4 GB typical, scales linearly
- **Storage**: ~100 MB per sample
- **CPU**: Dynamic allocation via psutil

---

## 🔬 **Validation & Accuracy**

### **Reference Strain Validation**
**100% concordance** with gold-standard reference genomes:

| Reference Strain | Expected Type | StaphScope Result |
|------------------|---------------|-------------------|
| USA300 | ST8–t008–IV(2B) | ✅ ST8–t008–IV(2B) |
| N315 | ST5–t002–II(2A) | ✅ ST5–t002–II(2A) |
| MRSA252 | ST36–t018–II(2A) | ✅ ST36–t018–II(2A) |
| TW20 | ST239–t037–III(3A) | ✅ ST239–t037–III(3A) |
| NCTC8325 | ST8–t211–None | ✅ ST8–t211–Not Assigned |

### **Clinical Isolate Analysis (n=24)**
- **MRSA**: 21 isolates (87.5%)
- **MSSA**: 3 isolates (12.5%)
- **Dominant STs**: ST5 (9), ST8 (5), ST22 (2)
- **Critical Genes**: *mecA* (21), *mecC* (1), *fosB* (20)
- **PVL**: 7 isolates (29.2%), all ST8/ST59
- **Plasmids**: 14/24 genomes (58.3%) with plasmid replicons

---

## 🆚 **Tool Comparison**

| Feature | StaphScope | Bactopia | Nullarbor | Mykrobe |
|---------|------------|----------|-----------|---------|
| **Analysis Focus** | 🎯 *S. aureus*-optimized | Multi-species | Multi-species | Multi-species |
| **Input Format** | Assembled genomes | Raw reads | Raw reads | Raw reads |
| **Installation** | Single Conda package | Complex (Nextflow+Docker) | Conda + DB downloads | Single Conda |
| **Execution** | Local CLI + Web GUI | Local/Cluster | Local | CLI + Web GUI |
| **Parallelization** | Auto-resource detection | Pipeline-level | Sample-level | Single-threaded |
| **MRSA Features** | Integrated classification + lineage DB | General typing | General typing | Resistance only |
| **Critical Gene Flagging** | ✅ *mecA*, PVL, *van* genes | ❌ | ❌ | ❌ |
| **Resource Needs** | Low-moderate (2+ GB) | High (HPC recommended) | High (Cluster) | Low-moderate |
| **Web Interface** | ✅ StaphScope Web | ❌ | ❌ | ✅ Mykrobe web |

---

## 🤖 **AI Integration Guide**

StaphScope generates comprehensive HTML reports that are **perfect for AI analysis**. Here's how to use AI tools to get more from your data.

### 🚀 Quick Start
1. **Install any AI browser extension** (ChatGPT, Claude, Gemini)
2. **Open your report**: `staphscope_ultimate_report.html`
3. **Select text** in any section (AMR Genes, MLST Analysis, etc.)
4. **Right-click → Ask AI** with your question

### 💡 Example Questions

**For MLST Analysis:**
- "What is the clinical significance of ST5 vs ST8?"
- "Which samples are MRSA and what ST are they?"

**For AMR Genes:**
- "Explain the mecA gene and its importance"
- "Which samples have multiple resistance genes?"
- "What treatment implications do these genes have?"

**For Virulence Factors:**
- "Which samples carry PVL toxin?"
- "Are there any high-risk virulence combinations?"

**For Pattern Discovery:**
- "Are there correlations between ST and specific genes?"
- "Identify any concerning patterns in this dataset"

### 📊 Pro Tips
- **Provide context**: "I'm analyzing *S. aureus* genomics data..."
- **Be specific**: Instead of "tell me about this", ask "what does SCCmec type IV indicate?"
- **Ask for interpretations**: "What are the clinical implications of these findings?"
- **Request summaries**: "Summarize the resistance profile of sample XYZ"

### ⚡ Why This Works
StaphScope reports are structured with clear tables and organized data that AI can easily understand. Each gene is shown with all genomes that contain it, making pattern analysis straightforward.

> *"AI provides powerful insights but always verify critical findings with domain experts."*

---

## 🔮 **Future Development**

### **🚀 Upcoming Features (2025-2026)**
```python
# Planned machine learning module
staphscope --ml-predict --input results.json --model outbreak_risk

# Raw read support (in development)
staphscope --raw-reads sample_R1.fastq sample_R2.fastq --assembler shovill
```

### **Machine Learning Module**
- **Outbreak Prediction**: Identify emerging patterns and transmission networks
- **Phenotype Inference**: Predict virulence, transmissibility from genotype
- **Risk Scoring**: Automated risk assessment for clinical isolates
- **Anomaly Detection**: Flag novel or unexpected genetic combinations

### **Expansion Plans**
1. **Raw Read Support**: Direct FASTQ analysis with integrated assembly (Snippy)
2. **Real-Time Updates**: Live database synchronization
3. **Plugin System**: Community-contributed analysis modules
4. **Database Contributions**: User-submitted lineage updates
5. **Translation Support**: Help translate the interface

---

## ❓ **Frequently Asked Questions**

### **General Questions**

**Q: Is StaphScope free to use?**  
A: Yes! StaphScope is open-source under the MIT License. Free for academic, clinical, and commercial use.

**Q: What makes StaphScope different from other tools?**  
A: StaphScope is *S. aureus*-optimized, integrates 6 analysis types in one workflow, runs 8-10× faster than generalist tools, and includes a curated global lineage database.

**Q: Can I use StaphScope for clinical diagnosis?**  
A: StaphScope is a research tool. While highly accurate, results should be validated with orthogonal methods for clinical decision-making.

**Q: Which version should I use – CLI or Web?**  
A: Use the **Web version** for convenience, small batches (≤10 files), and if you prefer a graphical interface. Use the **CLI version** for large batches (100+ genomes), integration into pipelines, or when working with sensitive data locally.

### **Technical Questions**

**Q: Why only assembled genomes? When will raw read support be added?**  
A: We focused first on assembled genomes for speed and simplicity. Raw read support is our #1 priority for 2026 development.

**Q: How often are databases updated?**  
A: We have planned sequential releases when database updates are needed. The lineage database is manually curated every 6 months. Users can run `abricate --setupdb` anytime.

**Q: Can I run StaphScope on Windows?**  
A: Yes, via WSL2 (Windows Subsystem for Linux). Native Windows support is planned.

**Q: How do I handle very large batches (1000+ genomes)?**  
A: Use the CLI with glob patterns and appropriate threading. StaphScope scales linearly.

### **Analysis Questions**

**Q: What does "Not Assigned" mean for SCCmec typing?**  
A: This indicates insufficient evidence for cassette classification—usually MSSA or novel SCCmec types.

**Q: How is MRSA status determined?**  
A: MRSA = positive for both SCCmec element AND *mecA* or *mecC* gene. MSSA = lacks either criterion.

**Q: Are virulence factors from other species filtered out?**  
A: Yes! The ABRicate module uses *S. aureus*-optimized thresholds and databases to minimize cross-species false positives.

---

## 🐛 **Troubleshooting**

### **Common Issues & Solutions**

```bash
# Issue: Database errors
# Solution:
abricate --setupdb

# Issue: Missing dependencies
# Solution:
conda remove staphscope
conda clean --all
conda install -c bbeckley-hub staphscope  # Fresh install

# Issue: Web version not loading
# Solution: Check internet connection or try a different browser.
# The service may be temporarily down; check GitHub for updates.
```

### **Getting Help**
1. **Check existing issues**: [GitHub Issues](https://github.com/bbeckley-hub/staphscope-typing-tool/issues)
2. **Search closed issues**: Many problems already solved
3. **Create new issue**: Include:
   - Full error message
   - `staphscope --version`
   - Conda environment list (`conda list`)
   - Example command that failed
4. **Email support**: brownbeckley94@gmail.com (response within 48 hours)

---

## 📚 **Citation**

If you use StaphScope in your research, please cite:

> Beckley, B., Amarh, V. (2026). StaphScope: a species‑optimized computational pipeline for rapid and accessible *Staphylococcus aureus* genotyping and surveillance. *BMC Genomics*, 27:123.

**DOI**: [10.1186/s12864-026-12609-x](https://doi.org/10.1186/s12864-026-12609-x)

```bibtex
@article{beckley2026staphscope,
  title={StaphScope: a species‑optimized computational pipeline for rapid and accessible Staphylococcus aureus genotyping and surveillance},
  author={Beckley, Brown and Amarh, Vincent},
  journal={BMC Genomics},
  volume={27},
  pages={123},
  year={2026},
  doi={10.1186/s12864-026-12609-x}
}
```

### **Software Citation**
```bibtex
@software{staphscope2026,
  author = {Brown Beckley},
  title = {StaphScope: A species-optimized computational pipeline for Staphylococcus aureus genotyping},
  year = {2026},
  publisher = {GitHub},
  url = {https://github.com/bbeckley-hub/staphscope-typing-tool}
}
```

### **Integrated Tool Citations**
Please also cite the essential tools that make StaphScope possible (see BibTeX in the repository).

---

## 🙏 **Acknowledgements**

StaphScope stands on the shoulders of giants. We are deeply grateful to:

- **Torsten Seemann** for MLST, ABRicate, and countless foundational tools
- **NCBI team** for AMRFinderPlus
- **CGE team** for SCCmecFinder and database curation
- **PubMedST, Ridom, CARD, VFDB** for essential databases
- **Python community** for Biopython, pandas, plotly, seaborn, matplotlib
- **Early adopters and beta testers** for invaluable feedback
- **Peer reviewers & Editorial Team @BMC GENOMICS** for their constructive feedback, which significantly strengthened this tool and it manuscript. 
> *"If we ever meet in person, the drinks are on me!" – Brown Beckley*

---

## 👥 **Authors & Contact**

**Brown Beckley** (Primary Developer)  
- University of Ghana Medical School  
- 📧 brownbeckley94@gmail.com  
- 🐙 GitHub: [bbeckley-hub](https://github.com/bbeckley-hub)  
- LinkedIn: [@brownbeckley](https://www.linkedin.com/in/brown-beckley-190315319/)  
- 📞 +233 508820617

**Amarh Vincent** (Co-Author)  
- University of Ghana Medical School

### **Collaboration Opportunities**
We welcome collaborations on:
- MRSA epidemiology studies
- Clinical validation projects
- Bioinformatics tool development
- Global surveillance initiatives
- Public health applications

---

## 📄 **License**

### Core StaphScope Code
The StaphScope pipeline code (the workflow engine, report generation, HTML templates, and Python modules written by the authors) is licensed under the **MIT License** – see the [LICENSE](LICENSE) file for details.

### StaphScope Web Code
The web interface is also open-source and available under the MIT License in its [separate repository](https://github.com/bbeckley-hub/staphscope-web).

### Third-Party Tools
StaphScope executes several external bioinformatics tools, which are installed as Conda dependencies. Each tool is the property of its respective developers and is used under its own license:

| Tool | License |
|------|---------|
| `mlst` (Torsten Seemann) | GPL v2 |
| `ABRicate` (Torsten Seemann) | GPL v2 |
| `AMRFinderPlus` (NCBI) | Public Domain |
| `SCCmecFinder` (CGE) | Apache-2.0 |
| `spa typing` (Ridom) | Free for academic use |

By using StaphScope, you agree to comply with the licenses of these third-party tools.

---

### 📚 **Third-Party Tool Citations**

StaphScope integrates several powerful open-source tools and databases. If you use StaphScope in your research, please also cite the following essential tools:

#### **MLST (Torsten Seemann)**
```bibtex
@software{seemann_mlst_2018,
  author = {Seemann, T.},
  title = {MLST: Scan contig files against traditional PubMLST typing schemes},
  year = {2018},
  publisher = {GitHub},
  url = {https://github.com/tseemann/mlst}
}
```

#### **ABRicate (Torsten Seemann)**
```bibtex
@software{seemann_abricate_2018,
  author = {Seemann, T.},
  title = {ABRicate: Mass screening of contigs for antimicrobial resistance and virulence genes},
  year = {2018},
  publisher = {GitHub},
  url = {https://github.com/tseemann/abricate}
}
```

#### **AMRFinderPlus (NCBI)**
```bibtex
@article{feldgarden_amrfinderplus_2019,
  author = {Feldgarden, M. et al.},
  title = {AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence},
  journal = {Scientific Reports},
  volume = {11},
  pages = {12728},
  year = {2019},
  doi = {10.1038/s41598-021-91456-0}
}
```

#### **SCCmecFinder (CGE)**
```bibtex
@article{kaya_sccmecfinder_2018,
  author = {Kaya, H. et al.},
  title = {SCCmecFinder, a Web-Based Tool for Typing of Staphylococcal Cassette Chromosome mec in Staphylococcus aureus Using Whole-Genome Sequence Data},
  journal = {mSphere},
  volume = {3},
  number = {1},
  pages = {e00612-17},
  year = {2018},
  doi = {10.1128/mSphere.00612-17}
}
```

#### ***spa* Typing (Ridom)**
```bibtex
@article{mellmann_spa_typing_2005,
  author = {Mellmann, A. et al.},
  title = {Evidenzbasierte Hygienemassnahmen mittels spa-Typisierung bei MRSA-Häufungen im Krankenhaus},
  journal = {Deutsche Medizinische Wochenschrift},
  volume = {130},
  number = {22},
  pages = {1364-1368},
  year = {2005},
  doi = {10.1055/s-2005-868351},
  note = {Database: https://spa.ridom.de}
}
```

---

### **📊 Database Citations**

#### **CARD (Comprehensive Antibiotic Resistance Database)**
```bibtex
@article{alcock_card_2023,
  author = {Alcock, B. P. et al.},
  title = {CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database},
  journal = {Nucleic Acids Research},
  volume = {51},
  number = {D1},
  pages = {D690-D699},
  year = {2023},
  doi = {10.1093/nar/gkac920}
}
```

#### **ResFinder**
```bibtex
@article{bortolaia_resfinder_2020,
  author = {Bortolaia, V. et al.},
  title = {ResFinder 4.0 for predictions of phenotypes from genotypes},
  journal = {Journal of Antimicrobial Chemotherapy},
  volume = {75},
  number = {12},
  pages = {3491-3500},
  year = {2020},
  doi = {10.1093/jac/dkaa345}
}
```

#### **ARG-ANNOT**
```bibtex
@article{gupta_argannot_2014,
  author = {Gupta, S. K. et al.},
  title = {ARG-ANNOT, a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes},
  journal = {Antimicrobial Agents and Chemotherapy},
  volume = {58},
  number = {1},
  pages = {212-220},
  year = {2014},
  doi = {10.1128/AAC.01310-13}
}
```

#### **VFDB (Virulence Factor Database)**
```bibtex
@article{chen_vfdb_2016,
  author = {Chen, L. et al.},
  title = {VFDB 2016: hierarchical and refined dataset for big data analysis—10 years on},
  journal = {Nucleic Acids Research},
  volume = {44},
  number = {D1},
  pages = {D694-D697},
  year = {2016},
  doi = {10.1093/nar/gkv1239}
}
```

#### **PlasmidFinder**
```bibtex
@article{carattoli_plasmidfinder_2014,
  author = {Carattoli, A. et al.},
  title = {In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing},
  journal = {Antimicrobial Agents and Chemotherapy},
  volume = {58},
  number = {7},
  pages = {3895-3903},
  year = {2014},
  doi = {10.1128/AAC.02412-14}
}
```

#### **EcOH (E. coli O/H typing)**
```bibtex
@article{joensen_ecoh_2015,
  author = {Joensen, K. G. et al.},
  title = {Rapid and easy in silico serotyping of Escherichia coli isolates by use of whole-genome sequencing data},
  journal = {Journal of Clinical Microbiology},
  volume = {53},
  number = {8},
  pages = {2410-2426},
  year = {2015},
  doi = {10.1128/JCM.00008-15}
}
```

#### **MEGARes 3.0**
```bibtex
@article{bonin_megares_2023,
  author = {Bonin, N. et al.},
  title = {MEGARes and AMR++, v3.0: an updated comprehensive database of antimicrobial resistance determinants and an improved software pipeline for classification using high-throughput sequencing},
  journal = {Nucleic Acids Research},
  volume = {51},
  number = {D1},
  pages = {D744-D752},
  year = {2023},
  doi = {10.1093/nar/gkac1047}
}
```

---

### 📝 **Usage Note**

When citing StaphScope in your publications, please include the main StaphScope citation along with citations for the specific tools and databases you used:

> "Genomic analysis was performed using StaphScope [Beckley & Amarh, 2026], which integrates MLST [Seemann, 2018], ABRicate [Seemann, 2018], AMRFinderPlus [Feldgarden et al., 2019], and SCCmecFinder [Kaya et al., 2018] for comprehensive *S. aureus* characterization. Antimicrobial resistance genes were identified using the CARD [Alcock et al., 2023] and ResFinder [Bortolaia et al., 2020] databases."

---

<div align="center">

## **🚀 Ready to revolutionize your MRSA analysis?**

| **Choose Your Platform** | |
|--------------------------|-|
| 🖥️ **Command Line** | For high-throughput, local analysis |
| 🌐 **StaphScope Web** | For non-bioinformaticians – [https://staphscope.dpdns.org](https://staphscope.dpdns.org) |

[![Get Started CLI](https://img.shields.io/badge/GET_STARTED_CLI-Now-green?style=for-the-badge&logo=github)](https://github.com/bbeckley-hub/staphscope-typing-tool#-quick-start-cli)
[![Try Web Version](https://img.shields.io/badge/TRY_WEB_VERSION-Here-blue?style=for-the-badge&logo=html5)](https://staphscope.dpdns.org)
[![Report Issue](https://img.shields.io/badge/REPORT_ISSUE-Here-red?style=for-the-badge&logo=github)](https://github.com/bbeckley-hub/staphscope-typing-tool/issues)

**From days to minutes. From fragmented to integrated. From data to insights.**

*StaphScope: Precision surveillance for the antibiotic resistance era.*

⭐ **If you find this tool useful, please star the repository!** ⭐

*Join the Fight Against Antimicrobial Resistance*

Antimicrobial resistance (AMR) represents one of the most significant global health threats of our time. We invite researchers, clinicians, and public health professionals to collaborate with us in expanding and validating our database, sharing regional epidemiological data, and advancing AMR surveillance.

**Together, we can enhance global AMR monitoring and develop more effective treatment strategies.**

</div>

