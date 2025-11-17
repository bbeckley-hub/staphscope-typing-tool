<div align="center"> 
  
  ![Version](https://anaconda.org/bbeckley-hub/staphscope/badges/version.svg)![Latest Release Date](https://anaconda.org/bbeckley-hub/staphscope/badges/latest_release_date.svg)![Platforms](https://anaconda.org/bbeckley-hub/staphscope/badges/platforms.svg)]![License](https://anaconda.org/bbeckley-hub/staphscope/badges/license.svg)]![Downloads](https://anaconda.org/bbeckley-hub/staphscope/badges/downloads.svg)]
  
  **Comprehensive MRSA genomic analysis pipeline for typing, resistance profiling, and lineage analysis** 
  
  [Quick Start](#-quick-start) • [Features](#-features) • [Installation](#-installation) • [Usage](#-usage) 
  
  </div> 
  
  ## 🎯 Purpose StaphScope is a comprehensive bioinformatics pipeline specifically designed for Methicillin-Resistant Staphylococcus aureus (MRSA)  genomic analysis. This all-in-one tool provides complete characterization of S. aureus isolates through multiple typing methods, antimicrobial resistance profiling, virulence factor detection, and lineage analysis.
  
  ### Key Applications 
  - **🏥 MRSA Surveillance**: Track and characterize MRSA strains in clinical and research settings
  - **🔍 Outbreak Investigation**: Identify related strains and transmission patterns
  - **🔬 Research Analysis**: Comprehensive genomic profiling for academic studies
  - **🌍 Public Health**: Support antimicrobial resistance monitoring programs

## ✨ Features ### 🔬 Core Analysis Modules 

| Module | Description | Key Outputs |  
| **MLST** | Multi-Locus Sequence Typing | Sequence Type (ST), Clonal Complex (CC) | 
| **spa Typing** | Staphylococcal Protein A typing | spa type, repeat sequence | 
| **SCCmec** | Staphylococcal Cassette Chromosome mec | SCCmec type, mec gene complex, ccr complex | 
| **AMR Profiling** | Antimicrobial Resistance genes | Resistance genes, drug classes, mechanisms | 
| **ABRicate** | Comprehensive resistance & virulence | Plasmid markers, virulence factors, resistance databases | 
| **Lineage Analysis** | Strain lineage reference | HTML report with strain classification | 

### 🛡️ MRSA-Specific Capabilities 

- **SCCmec Typing**: Accurate identification of SCCmec types I-XIII
- **mecA/mecC Detection**: Methicillin resistance determinant detection 
- **PVL Toxin Screening**: Panton-Valentine Leukocidin gene detection 
- **AMR Profile**: Comprehensive antimicrobial resistance pattern 
- **Epidemic Clones**: Identification of major MRSA clonal complexes (CC5, CC8, CC22, CC30, CC45)

## 🚀 Quick Start ### Installation **Option 1: Conda Installation (Recommended)**

conda install -c bbeckley-hub -c bioconda -c conda-forge staphscope

Option 2: From Source

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

