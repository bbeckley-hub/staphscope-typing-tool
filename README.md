````markdown
<details>
<summary><h2>📦 StaphScope: Advanced Staphylococcus aureus Typing & Lineage Analysis Platform</h2></summary>

<div align="center">

![Version](https://anaconda.org/bbeckley-hub/staphscope/badges/version.svg)
![Latest Release Date](https://anaconda.org/bbeckley-hub/staphscope/badges/latest_release_date.svg)
![Platforms](https://anaconda.org/bbeckley-hub/staphscope/badges/platforms.svg)
![License](https://anaconda.org/bbeckley-hub/staphscope/badges/license.svg)
![Downloads](https://anaconda.org/bbeckley-hub/staphscope/badges/downloads.svg)

**Comprehensive MRSA genomic analysis pipeline for typing, resistance profiling, and lineage analysis**  
**Supports Python 3.8 → 3.12**

[Quick Start](#-quick-start) • [Features](#-features) • [Installation](#-installation) • [Usage](#-usage)

</div>

---

## 🎯 **Purpose**
StaphScope is a unified bioinformatics workflow for **Methicillin-Resistant *Staphylococcus aureus* (MRSA)** analysis.  
It performs complete genomics-based characterization: typing, resistance profiling, virulence detection, and lineage assignment.

### **Key Applications**
- 🏥 **MRSA Surveillance**  
- 🔍 **Outbreak Investigation**  
- 🔬 **Research Genomics**  
- 🌍 **Public Health AMR Monitoring**

---

## ✨ **Features**

### 🔬 **Core Analysis Modules**
| Module | Description | Output |
|--------|-------------|--------|
| **MLST** | Sequence typing | ST, CC |
| **spa Typing** | Protein A typing | spa type, repeats |
| **SCCmec** | MRSA cassette typing | SCCmec type, mec/ccr complexes |
| **AMR Profiling** | AMRFinderPlus | AMR genes, mechanisms |
| **ABRicate** | Comprehensive screening | Resistance, virulence, plasmids |
| **Lineage Analysis** | Phylo-reference | Interactive HTML report |

### 🛡️ **MRSA-Specific Capabilities**
- SCCmec I–XIII detection  
- mecA/mecC detection  
- PVL toxin screening  
- Epidemic clone identification  
- Complete AMR pattern prediction  

---

## 🚀 **Quick Start**

### **Installation**

#### **Option 1 — Conda (Recommended)**
```bash
conda install -c bbeckley-hub -c bioconda -c conda-forge staphscope
````

#### **Option 2 — From Source**

```bash
git clone https://github.com/bbeckley-hub/staphscope-typing-tool.git
cd staphscope-typing-tool
conda env create -f environment.yml
conda activate staphscope
pip install -e .
```

---

## 📌 **Basic Usage**

```bash
# Single genome
staphscope -i genome.fasta -o results/

# Batch analysis
staphscope -i "*.fna" -o batch_results --threads 8

# Skip modules
staphscope -i "MRSA_*.fasta" -o analysis --threads 16 --skip-lineage
```

---

## 📋 **Complete Usage**

```
usage: staphscope [-h] -i INPUT -o OUTPUT [-t THREADS] [--skip-amr] [--skip-abricate] 
                  [--skip-mlst] [--skip-spa] [--skip-sccmec] [--skip-lineage]
```

### **Input Formats**

* `.fna`, `.fasta`, `.fa`, `.fn`
* Single genomes or glob patterns
* Contigs or complete assemblies

---

## 🔧 **Analysis Modules**

### **1. MLST**

* Database: PubMLST
* Outputs: ST, CC, allele profiles

### **2. spa Typing**

* Tool: spaTyper
* Outputs: spa type, repeats

### **3. SCCmec Finder**

* Types I–XI + subtypes
* Determines mec & ccr complexes

### **4. AMR Profiling**

* Tool: AMRFinderPlus
* 5,000+ resistance genes

### **5. ABRicate**

* Databases: CARD, ResFinder, NCBI, VFDB, MEGARES, PlasmidFinder

### **6. Lineage Reference**

* Interactive HTML report
* Global MRSA context

---

## 📊 **Output Structure**

```
output/
├── mlst_results/
├── spa_results/
├── sccmec_results/
├── amr_results/
├── abricate_results/
└── lineage_results/
```

---

## 🦠 **Key MRSA Markers**

| Category       | Markers              | Significance        |
| -------------- | -------------------- | ------------------- |
| **Resistance** | mecA, mecC, blaZ     | β-lactam resistance |
| **Virulence**  | PVL                  | Severe infections   |
| **Toxins**     | TSST-1, enterotoxins | Toxic shock         |
| **Adhesion**   | fnbA/B, clfA/B       | Biofilms            |

Major clones identified: **CC5, CC8, CC22, CC30, CC45**

---

## 🔬 **Use Cases**

### **Hospital Outbreak**

```bash
staphscope -i "outbreak_*.fasta" -o outbreak_analysis --threads 16
```

### **Surveillance**

```bash
staphscope -i "surveillance_*.fna" -o yearly_surveillance
```

### **Research**

```bash
staphscope -i isolate.fasta -o complete_analysis
```

---

## 💾 **System Requirements**

### Minimum

* 4 cores
* 4 GB RAM
* 4 GB disk

### Recommended

* 8+ cores
* 8+ GB RAM
* 10+ GB disk

### Dependencies (Auto-Installed)

* **Python 3.8 → 3.12**
* AMRFinderPlus
* ABRicate
* MLST
* BLAST+

---

## 🐛 **Troubleshooting**

### Database Fix

```bash
amrfinder --update
abricate --setupdb
```

### Memory Issues

```bash
staphscope -i "batch1_*.fna" -o results1 --threads 4
```

---

## 🧩 **Support**

* Issues: [https://github.com/bbeckley-hub/staphscope-typing-tool/issues](https://github.com/bbeckley-hub/staphscope-typing-tool/issues)
* Email: [brownbeckley94@gmail.com](mailto:brownbeckley94@gmail.com)
* Docs: `docs/` directory

---

## 📚 **Citation**

```bibtex
@software{staphscope2024,
  title = {StaphScope: Advanced Staphylococcus aureus Typing and Lineage Analysis Platform},
  author = {Brown Beckley},
  year = {2025},
  url = {https://github.com/bbeckley-hub/staphscope-typing-tool},
  note = {Comprehensive MRSA genomic analysis tool}
}
```

---

## 👨‍💻 **Author**

**Brown Beckley**
University of Ghana Medical School
Department of Medical Biochemistry
Email: [brownbeckley94@gmail.com](mailto:brownbeckley94@gmail.com)
GitHub: [https://github.com/bbeckley-hub](https://github.com/bbeckley-hub)

---

## 📄 **License**

MIT License – see `LICENSE` file.

---

</details>
```

---

