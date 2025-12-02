<div align="center"> 
  
  ![Version](https://anaconda.org/bbeckley-hub/staphscope/badges/version.svg)![Latest Release Date](https://anaconda.org/bbeckley-hub/staphscope/badges/latest_release_date.svg)![Platforms](https://anaconda.org/bbeckley-hub/staphscope/badges/platforms.svg)]![License](https://anaconda.org/bbeckley-hub/staphscope/badges/license.svg)]![Downloads](https://anaconda.org/bbeckley-hub/staphscope/badges/downloads.svg)]
  
  **Comprehensive MRSA genomic analysis pipeline for typing, resistance profiling, and lineage analysis** 
  
  [Quick Start](#-quick-start) • [Features](#-features) • [Installation](#-installation) • [Usage](#-usage) 
  
  </div> 

```bash
███████╗████████╗ █████╗ ██████╗ ██╗  ██╗███████╗ ██████╗ ██████╗ ██████╗ ███████╗
██╔════╝╚══██╔══╝██╔══██╗██╔══██╗██║  ██║██╔════╝██╔════╝██╔═══██╗██╔══██╗██╔════╝
███████╗   ██║   ███████║██████╔╝███████║███████╗██║     ██║   ██║██████╔╝█████╗  
╚════██║   ██║   ██╔══██║██╔═══╝ ██╔══██║╚════██║██║     ██║   ██║██╔═══╝ ██╔══╝  
███████║   ██║   ██║  ██║██║     ██║  ██║███████║╚██████╗╚██████╔╝██║     ███████╗
╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝     ╚═╝  ╚═╝╚══════╝ ╚═════╝ ╚═════╝ ╚═╝     ╚══════╝
```  
---
  
  ## 🎯 Purpose 
  StaphScope is a comprehensive bioinformatics pipeline specifically designed for Methicillin-Resistant Staphylococcus aureus (MRSA)  genomic analysis. This all-in-one tool provides complete characterization of S. aureus isolates through multiple typing methods, antimicrobial resistance profiling, virulence factor detection, and lineage analysis.

---
  
  ### Key Applications 
  - **🏥 MRSA Surveillance**: Track and characterize MRSA strains in clinical and research settings
  - **🔍 Outbreak Investigation**: Identify related strains and transmission patterns
  - **🔬 Research Analysis**: Comprehensive genomic profiling for academic studies
  - **🌍 Public Health**: Support antimicrobial resistance monitoring programs
---

## ✨ Features 
### 🔬 Core Analysis Modules 

| Module | Description | Key Outputs |  
| **MLST** | Multi-Locus Sequence Typing | Sequence Type (ST), Clonal Complex (CC) | 

| **spa Typing** | Staphylococcal Protein A typing | spa type, repeat sequence | 

| **SCCmec** | Staphylococcal Cassette Chromosome mec | SCCmec type, mec gene complex, ccr complex | 

| **AMR Profiling** | Antimicrobial Resistance genes | Resistance genes, drug classes, mechanisms | 

| **ABRicate** | Comprehensive resistance & virulence | Plasmid markers, virulence factors, resistance databases | 

| **Lineage Analysis** | Strain lineage reference | HTML report with strain classification | 

---

### 🛡️ MRSA-Specific Capabilities 

- **SCCmec Typing**: Accurate identification of SCCmec types I-XIII
- **mecA/mecC Detection**: Methicillin resistance determinant detection 
- **PVL Toxin Screening**: Panton-Valentine Leukocidin gene detection 
- **AMR Profile**: Comprehensive antimicrobial resistance pattern 
- **Epidemic Clones**: Identification of major MRSA clonal complexes (CC5, CC8, CC22, CC30, CC45)

---

## 🚀 Quick Start 
### Installation 
**Option 1: Conda Installation (Recommended)**
```bash
conda install -c bbeckley-hub -c bioconda -c conda-forge staphscope
```
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

* 2 cores
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
* CGECORE
* KMA
  

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
Acknowledgements
---
The creation of Staphscope is a testament to the power of open-source collaboration. It is, in every sense, a synthesis of the incredible work done by others who generously shared their tools and data with the world.

From the bottom of my heart, I thank the authors and maintainers of these integrated software packages and public datasets. Your contributions are the foundation upon which Staphscope is built.

And if you're reading this—thank you. If we ever meet in person, the drinks are on me! My greatest thanks to you all.

---

## 👨‍💻 **Author**

**Brown Beckley**
University of Ghana Medical School-
Department of Medical Biochemistry
Email: [brownbeckley94@gmail.com](mailto:brownbeckley94@gmail.com)
GitHub: [https://github.com/bbeckley-hub](https://github.com/bbeckley-hub)

---

## 📄 **License**

MIT License – see `LICENSE` file.

---
## 📄 **OTHER CITATIONS**

Please cite the following integrated tools:

>Larsen, M., Cosentino, S., Rasmussen, S., Rundsten, C., Hasman, H., Marvig, R., Jelsbak, L., Sicheritz-PontÃ©n, T., Ussery, D., Aarestrup, F., & Lund, O. (2012). Multilocus Sequence Typing of Total Genome Sequenced Bacteria.
Journal of Clinical Microbiology, 50(4), 1355-1361. doi: 10.12.0/JCM.06094-11

>Clausen, P., Aarestrup, F., & Lund, O. (2018). Rapid and precise alignment of raw reads against redundant databases with KMA.
Bmc Bioinformatics,19(1), 307

> 
    Seemann T, Abricate, Github https://github.com/tseemann/abricate
    NCBI AMRFinderPlus - doi: 10.1128/AAC.00483-19
    CARD - doi:10.1093/nar/gkw1004
    Resfinder - doi:10.1093/jac/dks261
    ARG-ANNOT - doi:10.1128/AAC.01310-13
    VFDB - doi:10.1093/nar/gkv1239
    PlasmidFinder - doi:10.1128/AAC.02412-14
    EcOH - doi:10.1099/mgen.0.000064
    MEGARES 2.00 - doi:10.1093/nar/gkz1010

> Feldgarden M, Brover V, Gonzalez-Escalona N, Frye JG, Haendiges J, Haft DH, Hoffmann M, Pettengill JB, Prasad AB, Tillman GE, Tyson GH, Klimke W. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Sci Rep. 2021 Jun 16;11(1):12728. doi: 10.1038/s41598-021-91456-0. PMID: 34135355; PMCID: PMC8208984.

---
UPCOMING FEATURES
---

Machine Learning analysis pattern discovery,
Regular database updates to strength MRSA surveillance

---
FOR COLLABORATION AND FEATURE SUGGESTIONS, DO NOT HESISTATE TO REACH OUT BY MAIL.
---
