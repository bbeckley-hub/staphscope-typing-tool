
```bash
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

![Version](https://anaconda.org/bbeckley-hub/staphscope/badges/version.svg)
![Latest Release Date](https://anaconda.org/bbeckley-hub/staphscope/badges/latest_release_date.svg)
![Platforms](https://anaconda.org/bbeckley-hub/staphscope/badges/platforms.svg)
![License](https://anaconda.org/bbeckley-hub/staphscope/badges/license.svg)
![Downloads](https://anaconda.org/bbeckley-hub/staphscope/badges/downloads.svg)

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Conda](https://img.shields.io/badge/conda-✓-green.svg)](https://docs.conda.io/en/latest/)
[![MIT License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![GitHub Issues](https://img.shields.io/github/issues/bbeckley-hub/staphscope-typing-tool)](https://github.com/bbeckley-hub/staphscope-typing-tool/issues)
[![GitHub Stars](https://img.shields.io/github/stars/bbeckley-hub/staphscope-typing-tool)](https://github.com/bbeckley-hub/staphscope-typing-tool/stargazers)

**Complete MRSA/MSSA genomic analysis in minutes — not hours**

[📖 Documentation](#-documentation) • [⚡ Quick Start](#-quick-start) • [✨ Features](#-features) • 
[🔧 Installation](#-installation) • [🚀 Usage](#-usage) • [📊 Output](#-output) • 
[📈 Performance](#-performance) • [🔮 Future Roadmap](#-future-roadmap) • [🤝 Contributing](#-contributing)

</div>

---

## 📋 **Table of Contents**
- [🎯 Overview](#-overview)
- [✨ Key Features](#-key-features)
- [⚡ Quick Start](#-quick-start)
- [🔧 Installation](#-installation)
- [🚀 Usage Guide](#-usage-guide)
- [📊 Output Structure](#-output-structure)
- [🔍 Analytical Modules](#-analytical-modules)
- [📈 Performance Benchmarks](#-performance-benchmarks)
- [🔬 Validation & Accuracy](#-validation--accuracy)
- [🆚 Tool Comparison](#-tool-comparison)
- [🔮 Future Development](#-future-development)
- [❓ FAQ](#-frequently-asked-questions)
- [🐛 Troubleshooting](#-troubleshooting)
- [📚 Citation](#-citation)
- [🙏 Acknowledgements](#-acknowledgements)
- [👥 Authors & Contact](#-authors--contact)
- [📄 License](#-license)

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

**Perfect for**: Clinical labs, outbreak investigations, research studies, and public health surveillance.

---

## ✨ **Key Features**

### 🔬 **Core Analytical Modules**
| Module | 🎯 Purpose | 📊 Key Outputs | ⚡ Speed |
|--------|------------|----------------|----------|
| **MLST Typing** | Phylogenetic classification via 7 housekeeping genes | ST, CC, allele profiles, epidemiological context | <1 min |
| ***spa* Typing** | Hypervariable region analysis of protein A gene | *spa* type, repeat patterns, alignment metrics | <1 min |
| **SCC*mec* Typing** | Methicillin resistance cassette characterization | SCC*mec* type (I-XIII), *mec*/*ccr* complexes, confidence scores | 1-2 min |
| **AMR Profiling** | Comprehensive resistance gene detection | 5,000+ AMR genes, risk categorization, cross-sample patterns | 2-3 min |
| **ABRicate Screening** | Multi-database virulence/plasmid detection | 9 databases, plasmid replicons, virulence factors, clinical flags | 3-4 min |
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

---

## ⚡ **Quick Start**

### **Install in 60 seconds**
```bash
# Method 1: Conda (Recommended - handles all dependencies)
conda install -c bbeckley-hub -c bioconda -c conda-forge staphscope

# Method 2: From source
git clone https://github.com/bbeckley-hub/staphscope-typing-tool.git
cd staphscope-typing-tool
conda env create -f environment.yml
conda activate staphscope
pip install -e .
```
Refer to [**Update Databases (Recommended)- Please refer to this resources for STAPHSCOPE's integrated databases.**]
```

### **Run your first analysis**
```bash
# Single genome analysis
staphscope -i genome.fasta -o results/

# Batch processing (24 genomes)
staphscope -i "*.fna" -o batch_results --threads 16

# Analysis complete in ~14 minutes! 🎉

```
----
STAPHSCOPE TERMINAL DISPLAY
----
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
`
Output: Comprehensive results for all analyses in organized directories
Please run abricate --setupdb for recent gene annotations!!!
⭐ Star us on GitHub if you find this tool useful!

Transforming fragmented genomic data into coherent biological narratives 🧬✨
```   
---
## 🔧 **Installation**

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
# Download Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Follow prompts, then:
source ~/.bashrc
```

#### **2. Install StaphScope**
```bash
# Create and activate environment
conda create -n staphscope python=3.10 or 3.8 or 3.9 or 3.11 or 3.12
conda activate staphscope

# Install from conda-forge channel
conda install -c bbeckley-hub -c bioconda -c conda-forge staphscope

# Verify installation
staphscope --version
```
---

NB: ALWAYS CREATE A PYTHON ENV BEFORE INSTALLING STAPHSCOPE TO AVOID CONDA ISSUES!!!!!!
--

#### **3. Update Databases (Recommended)**
```bash
# Update ABRicate databases(version ≥1.0.1)

abricate --setupdb
```
# AMRFinderPlus database
Bundles v4.24 (latest)

# All other Bundled Databases
[https://pubmlst.org/organisms/staphylococcus-aureus]
[https://github.com/tseemann/mlst]
[http://spa.ridom.de/dynamic/sparepeats.fasta]
[https://spa.ridom.de/dynamic/spatypes.txt]
[https://bitbucket.org/genomicepidemiology/Sccmecfinder]
[https://github.com/ncbi/amr]

### **Docker Installation (Alternative)**
---
```dockerfile
# Coming soon! Containerized version in development
# docker pull bbeckley/staphscope:latest
```

---

## 🚀 **Usage Guide**

### **Basic Commands**
```bash
# Single genome analysis
staphscope -i /path/to/genome.fasta -o /path/to/results

# Batch processing with wildcards
staphscope -i "*.fna" -o results_2025 --threads 8

# Specify custom number of threads
staphscope -i "*.fasta" -o analysis -t 16

# Skip specific modules (if already analyzed)
staphscope -i sample.fna -o results --skip-spa --skip-lineage
```
### **Input Formats**
- **Accepted**: `.fna`, `.fasta`, `.fa`, `.fn`
- **Required**: Assembled genomes (contigs or complete)
- **Not supported**: Raw reads (FASTQ) - see Future Roadmap
- **Batch patterns**: `*.fasta`, `sample_*.fna`, `[0-9].fa`

### **Real-World Examples**

#### **Clinical Laboratory Setting**
```bash
# Daily surveillance of 12 isolates
staphscope -i "daily_isolates/*.fasta" -o /mnt/shared/surveillance/$(date +%Y%m%d) --threads 12

# Expected: Complete analysis in ~8 minutes
# Output: Interactive HTML report for clinical team review
```

#### **Research Project Analysis**
```bash
# 100-genome phylogeny project
Supporst glob patterns

# Expected: All 100 genomes analyzed in ~90 minutes
# Output: Machine-readable JSON/TSV for downstream analysis
```

#### **Outbreak Response**
```bash
# Urgent outbreak investigation (8 suspected cases)
staphscope -i "outbreak/*.fasta" -o /tmp/urgent_analysis --skip-lineage

# Expected: Results in ~4 minutes
# Output: Immediate identification of shared SCCmec types and resistance profiles
```
---

## 📁 **Output Structure/Directory Layout**

Staphscope generates a comprehensive, organized output directory with results from each analysis module. Below is a typical directory tree (example from sample `MRSA252`):

```
Staphscope/
├── abricate_results/                    # Gene detection (ABRicate)
│   ├── MRSA252/                         # Per-sample detailed results
│   │   ├── abricate_*.txt               # Raw ABRicate outputs per database
│   │   ├── abricate_*_report.html       # HTML reports per database
│   │   └── MRSA252_comprehensive_abricate_report.html
│   ├── staph_*_abricate_summary.tsv     # Combined TSV summaries per database
│   ├── staph_*_summary.json             # JSON summaries per database
│   ├── staph_*_summary_report.html      # HTML summary reports per database
│   └── staph_abricate_master_summary.json  # Master summary across all DBs
│
├── amr_results/                         # AMR gene profiling (AMRFinder+)
│   ├── MRSA252/
│   │   ├── MRSA252_amrfinder.txt        # Raw AMRFinder output
│   │   └── MRSA252_amrfinder_report.html
│   ├── staph_amrfinder_summary.tsv      # Tabular summary
│   ├── staph_amrfinder_summary.json     # JSON summary
│   ├── staph_amrfinder_summary_report.html
│   ├── staph_amrfinder_statistics_summary.tsv  # Statistical summary
│   └── staph_amrfinder_master_summary.json     # Master JSON
│
├── mlst_results/                        # Multi-Locus Sequence Typing
│   ├── MRSA252/
│   │   ├── mlst_raw_output.txt          # Raw MLST output
│   │   ├── mlst_report.txt/.tsv/.html   # Formatted reports
│   ├── mlst_summary.tsv                 # Combined TSV summary
│   ├── mlst_summary.json                # Combined JSON summary
│   └── mlst_summary.html                # Combined HTML report
│
├── sccmec_results/                      # SCCmec typing (MyKmerFinder)
│   ├── s_MRSA252/                       # Per-sample SCCmec results
│   │   ├── results_MyKmerFinder.txt     # Kmer-based typing results
│   │   ├── results_tab_MyDbFinder.txt   # Database matching results
│   │   ├── sccmec_detailed_results.txt  # Detailed typing report
│   │   ├── sccmec_enhanced_report.json  # Enhanced JSON report
│   │   └── staphscope_comprehensive_report.html
│   ├── staphscope_summary.tsv           # Combined SCCmec summary
│   ├── staphscope_summary.html          # HTML summary
│   └── staphscope_detailed_results.csv  # Detailed combined results
│
├── spa_results/                         # spa typing
│   ├── MRSA252/
│   │   ├── spa_typing_raw.txt           # Raw spa typing output
│   │   ├── spa_typing_report.txt/.tsv/.html
│   ├── spa_summary.tsv                  # Combined TSV summary
│   ├── spa_summary.json                 # Combined JSON summary
│   └── spa_summary.html                 # Combined HTML report
│
├── lineage_results/                     # Phylogenetic lineage assignment
│   └── staphscope_lineage_reference.html  # Lineage reference report
│
└── Staphscope_final_report/             # Consolidated final reports
    ├── staphscope_comprehensive_report.html/.json/.tsv  # Master reports
    ├── staphscope_summary.html/.tsv                     # High-level summary
    ├── staphscope_detailed_results.csv                  # Detailed CSV
    ├── mlst_summary.*                                   # MLST summaries
    ├── spa_summary.*                                    # spa summaries
    ├── staph_amrfinder_*                                # AMR summaries
    ├── staph_*_abricate_summary.tsv                    # ABRicate summaries
    ├── staph_*_summary.json                            # JSON summaries
    └── staph_*_summary_report.html                     # HTML summaries
```

### **Key File Types**
- **`.txt` / `.tsv` / `.csv`**: Raw and tabulated data for downstream analysis
- **`.json`**: Structured data for programmatic access and integration
- **`.html`**: Interactive visual reports for manual inspection
- **Per-sample directories**: Contain raw and detailed outputs for individual isolates
- **Summary files**: Aggregated results across all processed samples

### **Quick Access**
- **Single-sample overview**: Check the sample-specific directory (e.g., `MRSA252/`) within each module
- **Cross-sample summaries**: Look for `*_summary.tsv` or `*_summary.json` in each module's root
- **Final consolidated report**: All key results merged in `Staphscope_final_report/`
- **Visualization-ready data**: Most `.tsv` and `.json` files are optimized for the Automatic Visualization Module

This organized structure ensures easy navigation, reproducibility, and integration with downstream bioinformatics workflows!!

### **Interactive HTML Report Features**
- **Dashboard Overview**: Summary statistics at a glance
- **Interactive Tables**: Sort, filter, search all results
- ** Cross genome pattern discovery
- **Clinical Alerts**: Color-coded risk indicators

### **Machine-Readable Outputs**
```json
// Example JSON output structure
{
  "sample": "USA300_FPR3757",
  "mlst": {"st": "ST8", "cc": "CC8", "alleles": ["1","1","1","1","1","1","1"]},
  "spa": {"type": "t008", "repeats": "11-19-12-21-17-34-24-34-22-25"},
  "sccmec": {"type": "IV(2B)", "confidence": "very-high", "mec_complex": "A", "ccr_complex": "2"},
  "mrsa_status": "MRSA",
  "amr_genes": [
    {"gene": "mecA", "risk": "CRITICAL", "coverage": 100, "identity": 99.8},
    {"gene": "fosB", "risk": "HIGH", "coverage": 100, "identity": 98.5}
  ],
  "virulence_factors": ["lukS-PV", "lukF-PV", "hlgA", "hlgB", "hlgC"],
  "lineage": {
    "name": "USA300",
    "classification": "CA-MRSA",
    "risk_level": "High",
    "geography": "North America",
    "pvl_status": "Positive"
  }
}
```

---

## 🔍 **Analytical Modules**

### **1. MLST Typing** 🧬
- **Database**: PubMedST *S. aureus* 
- **Method**: BLAST-based allele calling (100% coverage/identity default)
- **Output**: ST, CC, 7-gene profile (*arcC, aroE, glpF, gmk, pta, tpi, yqiL*)
- **Enhanced**: Automatic lineage database query for epidemiological context

### **2. *spa* Typing** 🧬
- **Database**: Ridom *spa* repeat database
- **Method**: BLAST against repeat sequences
- **Output**: *spa* type, repeat pattern, contig location, alignment metrics

### **3. SCC*mec* Typing** 🛡️
- **Tool**: SCCmecFinder (hierarchical two-method system)
- **Primary**: Gene-based (*ccr*/*mec* complexes, 90% ID, 60% coverage)
- **Secondary**: k-mer homology (types I-XIII, ≥50% template coverage)
- **Confidence Levels**: Very-high, High, Medium, Low, Not Assigned
- **Subtyping**: Types IV and V community-associated cassettes

### **4. AMR Profiling** 💊
- **Tool**: NCBI-AMRFinderPlus v4.2.4 (curated database 2025-12-03.1)-bundled 
- **Optimization**: *S. aureus*-specific database curation
- **Risk Assessment**: Critical Risk (*mecA*, *vanA*, *cfr*), High Risk (*erm*, *tetM*)
- **Pattern Discovery**: Cross-genome frequency analysis

### **5. ABRicate Screening** 🔍
- **Databases**: 
  - VFDB (Virulence factors)
  - ResFinder (Acquired resistance)
  - CARD (Comprehensive resistance)
  - PlasmidFinder (Replicon typing)
  - MegaRes, NCBI, ARG-ANNOT, ECOH, EcoLi_VF
- **Thresholds**: ≥80% identity and coverage
- **Clinical Flags**: Automatic highlighting of PVL, enterotoxins, *van* genes

### **6. Lineage Database** 🌍
- **Content**: 44 major *S. aureus* lineages (18 HA-MRSA, 19 CA-MRSA, 7 LA-MRSA)
- **Metadata**: Geographical distribution, clinical significance, virulence profiles
- **High-Risk**: 9 lineages flagged as high-risk, 5 PVL-positive
- **Manual Curation**: Updated via periodic literature review

---

## 📈 **Performance Benchmarks**

### **Speed Comparison**
| System | Samples | Time | Speed vs Bactopia |
|--------|---------|------|-------------------|
| 💻 Laptop (2 cores, 8GB RAM) | 1 | 2m 33s | 5× faster |
| 💻 Laptop (2 cores, 8GB RAM) | 24 | 28m 17s | 6× faster |
| 🖥️ Workstation (16 cores, 16GB RAM) | 1 | 1m 31s | 8× faster |
| 🖥️ Workstation (16 cores, 16GB RAM) | 24 | 14m 34s | 10× faster |
| 🖥️ Workstation (16 cores, 16GB RAM) | 100 | ~60m | 12× faster |

### **Resource Efficiency**
- **Memory Usage**: 2-4 GB typical, scales linearly with samples
- **CPU Utilization**: Dynamic allocation via psutil (no resource waste)
- **Storage**: ~100 MB per sample analysis
- **Parallelization**: Sample-level + intra-module threading

### **Validation Accuracy**
| Reference Strain | Expected Type | StaphScope Result | Concordance |
|------------------|---------------|-------------------|-------------|
| USA300 | ST8–t008–IV(2B) | ST8–t008–IV(2B) | ✅ 100% |
| N315 | ST5–t002–II(2A) | ST5–t002–II(2A) | ✅ 100% |
| MRSA252 | ST36–t018–II(2A) | ST36–t018–II(2A) | ✅ 100% |
| TW20 | ST239–t037–III(3A) | ST239–t037–III(3A) | ✅ 100% |
| NCTC8325 | ST8–t211–None | ST8–t211–Not Assigned | ✅ 100% |

### **Case Study: 24 Clinical Isolates**
- **MRSA**: 21 isolates (87.5%)
- **MSSA**: 3 isolates (12.5%)
- **Dominant STs**: ST5 (9), ST8 (5), ST22 (2)
- **Critical Genes**: *mecA* (21), *mecC* (1), *fosB* (20)
- **PVL**: 7 isolates (29.2%), all ST8/ST59
- **Plasmids**: 14/24 genomes (58.3%) with plasmid replicons

---

## 🔬 **Validation & Accuracy**

### **Reference Strain Validation**
StaphScope was validated against gold-standard reference genomes with **100% concordance** for:
- MLST types (PubMedST database)
- *spa* types (Ridom database)
- SCC*mec* types (CGE reference)
- AMR profiles (NCBI-AMRFinderPlus)

### **Clinical Isolate Analysis**
**24 diverse *S. aureus* genomes analyzed:**
```bash
# Detected lineages
ST5 (Healthcare-associated): 9 isolates (37.5%)
ST8 (USA300, CA-MRSA): 5 isolates (20.8%)
ST22 (EMRSA-15): 2 isolates (8.3%)
ST239 (Brazilian/Hungarian): 1 isolate (4.2%)
ST59 (Asian CA-MRSA): 1 isolate (4.2%)
ST398 (Livestock-associated): 1 isolate (4.2%)
ST9 (Livestock-associated): 2 isolates (8.3%)
ST36 (EMRSA-16): 1 isolate (4.2%)
ST425: 1 isolate (4.2%)
```

### **Resistance Gene Prevalence**
| Gene | Prevalence | Risk Level | Phenotype |
|------|------------|------------|-----------|
| *mecA* | 87.5% (21/24) | CRITICAL | Methicillin resistance |
| *fosB* | 83.3% (20/24) | HIGH | Fosfomycin resistance |
| *blaZ* | 37.5% (9/24) | HIGH | Beta-lactamase |
| *dfrG* | 16.7% (4/24) | HIGH | Trimethoprim resistance |
| *mecC* | 4.2% (1/24) | CRITICAL | Alternative methicillin resistance |

### **Core Virulence Factors**
- **100% prevalence**: *hlgA/B/C* (gamma-hemolysin), *hld* (delta-hemolysin), *aur* (aureolysin)
- **29.2% prevalence**: PVL genes (*lukS-PV*, *lukF-PV*) in ST8/ST59 lineages

---

## 🆚 **Tool Comparison**

### **Feature Comparison Table**
| Feature | StaphScope | Bactopia | Nullarbor | Mykrobe |
|---------|------------|----------|-----------|---------|
| **Analysis Focus** | 🎯 *S. aureus*-optimized | Multi-species | Multi-species | Multi-species |
| **Input Format** | Assembled genomes | Raw reads | Raw reads | Raw reads |
| **Installation** | Single Conda package | Complex (Nextflow+Docker) | Conda + DB downloads | Single Conda |
| **Execution** | Local CLI | Local/Cluster | Local | CLI + Web GUI |
| **Parallelization** | Auto-resource detection | Pipeline-level | Sample-level | Single-threaded |
| **MRSA Features** | Integrated classification + lineage DB | General typing | General typing | Resistance only |
| **Critical Gene Flagging** | ✅ *mecA*, PVL, *van* genes | ❌ Absent | ❌ Absent | ❌ Absent |
| **Resource Needs** | Low-moderate (2+ GB) | High (HPC recommended) | High (Cluster) | Low-moderate |
| **Setup Ease** | Single command | Multiple steps | Multiple steps | Single command |

### **When to Choose StaphScope**
- ✅ **Ideal for**: *S. aureus*-specific research, clinical MRSA surveillance, outbreak response
- ✅ **Best when**: You need integrated typing + resistance + virulence in one workflow
- ✅ **Perfect if**: You value speed (minutes vs hours) and data privacy (local execution)

### **When to Choose Other Tools**
- ⚠️ **Use Bactopia/Nullarbor**: Multi-species projects, raw read analysis, extensive QC
- ⚠️ **Use Mykrobe**: Quick resistance profiling only, web interface preferred

---

## 🔮 **Future Development**

### **🚀 Upcoming Features (2025-2026)**
```python
# Planned machine learning module
staphscope --ml-predict --input results.json --model outbreak_risk

# Raw read support (in development)
staphscope --raw-reads sample_R1.fastq sample_R2.fastq --assembler shovill


```
---
## 📊 **Automatic Visualization Module**
---
This module will automatically generates publication-quality visualizations from **Staphscope** analysis results using modern Python plotting libraries.

### **Features**
- **Multi-format Support**: Generate PNG, SVG, PDF, and interactive HTML visualizations
- **Comprehensive Plot Types**:
  - **Statistical Plots**: Box plots, violin plots, and distribution histograms
  - **Comparison Charts**: Bar charts, grouped bars, and stacked plots
  - **Trend Analysis**: Line graphs, scatter plots with regression lines
  - **Composition Views**: Pie charts, donut charts, and treemaps
  - **Correlation Insights**: Heatmaps, pair plots, and correlation matrices
- **Smart Defaults**: Automatically selects appropriate plot types based on data structure
- **Customizable Themes**: Built-in color palettes optimized for scientific publishing

### **Supported Libraries**
- **Seaborn**: Statistical visualizations with beautiful default styles
- **Matplotlib**: Foundation layer for complete customization
- **Plotly**: Interactive HTML plots for exploratory analysis
- **Pandas**: Built-in plotting for quick data exploration


### **Output Examples**
- `sample_distribution.png` - Diversity metrics across samples
- `variant_frequency_bar.svg` - Top mutations with confidence intervals
- `correlation_heatmap.html` - Interactive sample similarity matrix
- `time_series_trend.pdf` - Longitudinal tracking of key markers

------

### **Machine Learning Module**
- **Outbreak Prediction**: Identify emerging patterns and transmission networks
- **Phenotype Inference**: Predict virulence, transmissibility from genotype
- **Risk Scoring**: Automated risk assessment for clinical isolates
- **Anomaly Detection**: Flag novel or unexpected genetic combinations
  
### **Expansion Plans**
1. **Raw Read Support**: Direct FASTQ analysis with integrated assembly(Snippy)
2. **Real-Time Updates**: Live database synchronization
   
### **Community-Driven Development**
- **Plugin System**: Community-contributed analysis modules
- **Database Contributions**: User-submitted lineage updates
- **Benchmark Datasets**: Shared validation datasets
- **Translation Support**: Help translate the interface to your language

---

## ❓ **Frequently Asked Questions**

### **General Questions**
**Q: Is StaphScope free to use?**  
A: Yes! StaphScope is open-source under the MIT License. Free for academic, clinical, and commercial use.

**Q: What makes StaphScope different from other tools?**  
A: StaphScope is *S. aureus*-optimized, integrates 6 analysis types in one workflow, runs 8-10× faster than generalist tools, and includes a curated global lineage database.

**Q: Can I use StaphScope for clinical diagnosis?**  
A: StaphScope is a research tool. While highly accurate, results should be validated with orthogonal methods for clinical decision-making.

### **Technical Questions**
**Q: Why only assembled genomes? When will raw read support be added?**  
A: We focused first on assembled genomes for speed and simplicity. Raw read support is our #1 priority for 2026 development.

**Q: How often are databases updated?**  
A: We have planned sequential releases when databases updates are needed. The lineage database is manually curated every 6 months. Users can run `abricate --setupdb` anytime.

**Q: Can I run StaphScope on Windows?**  
A: Yes, via WSL2 (Windows Subsystem for Linux). Native Windows support is planned.

**Q: How do I handle very large batches (1000+ genomes)?**  
A: Just you use glob patterns and take a coffee break.

### **Analysis Questions**
**Q: What does "Not Assigned" mean for SCCmec typing?**  
A: This indicates insufficient evidence for cassette classification—usually MSSA or novel SCCmec types.

**Q: How is MRSA status determined?**  
A: MRSA = positive for both SCCmec element AND *mecA* or *mecC* gene. MSSA = lacks either criterion.

**Q: Are virulence factors from other species filtered out?**  
A: Yes! The ABRicate module uses *S. aureus*-optimized thresholds and databases minimize cross-species false positives.

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

### **Primary Citation**
If you use StaphScope in your research, please cite our manuscript:

```bibtex
@article{beckley2025staphscope,
  title={StaphScope: A species-optimized computational pipeline for rapid and accessible Staphylococcus aureus genotyping and surveillance},
  author={Beckley, Brown and Vincent, Amarh},
  journal={In preparation},
  year={2025},
  note={Manuscript submitted for publication}
}
```

### **Software Citation**
```bibtex
@software{staphscope2025,
  title = {StaphScope: Comprehensive Staphylococcus aureus Genotyping Pipeline},
  author = {Brown Beckley},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/bbeckley-hub/staphscope-typing-tool},
  version = {1.0.0}
}
```

### **Integrated Tool Citations**
Please also cite these essential tools that make StaphScope possible:

```bibtex
# MLST
@article{seemann2018mlst,
  title={mlst: Scan contig files against traditional PubMLST typing schemes},
  author={Seemann, Torsten},
  year={2018},
  publisher={GitHub}
}

# AMRFinderPlus
@article{feldgarden2021amrfinderplus,
  title={AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence},
  author={Feldgarden, Michael and others},
  journal={Scientific Reports},
  volume={11},
  number={1},
  pages={12728},
  year={2021}
}

# ABRicate
@software{seemann2024abricate,
  title={ABRicate: Mass screening of contigs for antimicrobial and virulence genes},
  author={Seemann, Torsten},
  year={2024},
  publisher={GitHub}
}

# SCCmecFinder
@article{kaya2018sccmecfinder,
  title={SCCmecFinder, a Web-Based Tool for Typing of Staphylococcal Cassette Chromosome mec in Staphylococcus aureus Using Whole-Genome Sequence Data},
  author={Kaya, H and others},
  journal={mSphere},
  volume={3},
  number={1},
  pages={e00612-17},
  year={2018}
}
```

---

## 🙏 **Acknowledgements**

### **Open Source Foundations**
StaphScope stands on the shoulders of giants. We are deeply grateful to:

- **Tool Developers**: Torsten Seemann (MLST, ABRicate), NCBI team (AMRFinderPlus), H. Kaya (SCCmecFinder)
- **Database Curators**: PubMedST, Ridom *spa*, CGE, CARD, VFDB teams
- **Python Ecosystem**: Biopython, psutil, pandas, plotly developers
- **Testing Community**: Early adopters who provided invaluable feedback

### **Special Thanks**
- **Reviewers & Editors**: For strengthening this manuscript
- **University of Ghana**: For institutional support
- **Open Science Community**: For making this work possible

> "If we ever meet in person, the drinks are on me!" - Brown Beckley

### **How to Contribute**
1. **Report Bugs**: GitHub Issues
2. **Suggest Features**: GitHub Discussions
3. **Improve Documentation**: Pull requests welcome
4. **Share Data**: Contribute to the lineage database
5. **Translate**: Help translate to your language

---

## 👥 **Authors & Contact**

### **Primary Developer**
**Brown Beckley**  
- 🎓 Department of Medical Biochemistry, University of Ghana Medical School  
- 🔬 Department of Biochemistry and Biotechnology, KNUST  
- 📧 brownbeckley94@gmail.com  
- 🐙 GitHub: [bbeckley-hub](https://github.com/bbeckley-hub)  
- 🐦 LinkedIn: [@brownbeckley](https://www.linkedin.com/in/brown-beckley-190315319/)  

### **Co-Author**
**Amarh Vincent**  
- 🎓 Department of Medical Biochemistry, University of Ghana Medical School  


### **Collaboration Opportunities**
We welcome collaborations on:
- 🧬 MRSA epidemiology studies
- 🏥 Clinical validation projects
- 💻 Bioinformatics tool development
- 🌍 Global surveillance initiatives
  

**Contact for collaboration**: brownbeckley94@gmail.com

### **Stay Updated**
- **GitHub Releases**: Star and watch the repository
- **LinkedIn**: Follow for announcements

---

## 📄 **License**

StaphScope is released under the **MIT License**:

### **Third-Party Licenses**
StaphScope integrates several open-source tools, each with their own licenses:
- MLST: GPL-3.0
- ABRicate: GPL-2.0
- AMRFinderPlus: Public Domain
- SCCmecFinder: Apache-2.0

All dependencies are properly credited and their licenses respected.

---

<div align="center">

## **🚀 Ready to revolutionize your MRSA analysis?**

[![Get Started](https://img.shields.io/badge/GET_STARTED-Now-green?style=for-the-badge&logo=github)](https://github.com/bbeckley-hub/staphscope-typing-tool#-quick-start)
[![View Demo](https://img.shields.io/badge/VIEW_DEMO-Report-blue?style=for-the-badge&logo=html5)](https://github.com/bbeckley-hub/staphscope-typing-tool/tree/main/examples)
[![Report Issue](https://img.shields.io/badge/REPORT_ISSUE-Here-red?style=for-the-badge&logo=github)](https://github.com/bbeckley-hub/staphscope-typing-tool/issues)

**From days to minutes. From fragmented to integrated. From data to insights.**

*StaphScope: Precision surveillance for the antibiotic resistance era.*

</div>


  Found this tool useful? Drop a star ⭐ and follow the page for more exciting updates on planned modules!!

---

