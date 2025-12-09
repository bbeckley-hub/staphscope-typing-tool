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
import React from 'react';
import { Database, GitBranch, FileText, Activity, Shield, Dna, Settings, CheckCircle, FileSpreadsheet, Globe } from 'lucide-react';

const StaphScopeWorkflow = () => {
  const modules = [
    { name: 'MLST', icon: Dna, color: '#3b82f6', desc: 'Multi-Locus Sequence Typing' },
    { name: 'spa typing', icon: Activity, color: '#8b5cf6', desc: 'Protein A Typing' },
    { name: 'SCCmec', icon: Shield, color: '#ec4899', desc: 'Methicillin Resistance Cassette' },
    { name: 'AMRFinderPlus', icon: Settings, color: '#f59e0b', desc: 'AMR Gene Detection' },
    { name: 'ABRicate', icon: Database, color: '#10b981', desc: 'Resistance & Virulence Screening' },
    { name: 'Lineage DB', icon: GitBranch, color: '#06b6d4', desc: 'Lineage Reference Generation' }
  ];

  const outputs = [
    { format: 'HTML', icon: Globe, color: '#3b82f6' },
    { format: 'JSON', icon: FileText, color: '#8b5cf6' },
    { format: 'TSV', icon: FileSpreadsheet, color: '#10b981' }
  ];

  return (
    <div className="w-full h-full bg-gradient-to-br from-slate-50 to-blue-50 p-8 font-sans">
      <div className="max-w-7xl mx-auto">
        {/* Title */}
        <div className="text-center mb-8">
          <h1 className="text-3xl font-bold text-slate-800 mb-2">
            StaphScope Workflow and Analysis Architecture
          </h1>
          <p className="text-slate-600 text-sm max-w-3xl mx-auto">
            Schematic overview of StaphScope's modular parallel execution architecture, showing 
            simultaneous analysis across six specialized modules with comprehensive multi-format reporting
          </p>
        </div>

        {/* Main Workflow Container */}
        <div className="bg-white rounded-2xl shadow-xl p-8 border-2 border-slate-200">
          
          {/* Input Stage */}
          <div className="flex flex-col items-center mb-8">
            <div className="bg-gradient-to-r from-blue-500 to-blue-600 text-white px-8 py-4 rounded-xl shadow-lg mb-4 flex items-center gap-3">
              <FileText size={28} />
              <div>
                <div className="font-bold text-lg">Input Genome Assembly</div>
                <div className="text-sm text-blue-100">FASTA format (.fna, .fasta, .fa)</div>
              </div>
            </div>
            <div className="h-8 w-0.5 bg-slate-300"></div>
          </div>

          {/* Initialization */}
          <div className="flex flex-col items-center mb-8">
            <div className="bg-gradient-to-r from-purple-500 to-purple-600 text-white px-6 py-3 rounded-lg shadow-md flex items-center gap-2">
              <Settings className="animate-spin" size={20} />
              <span className="font-semibold">Initialize StaphScope Platform</span>
            </div>
            <div className="text-xs text-slate-500 mt-2">Database loading • Engine configuration • Thread optimization</div>
            <div className="h-8 w-0.5 bg-slate-300 mt-4"></div>
          </div>

          {/* Parallel Processing Header */}
          <div className="text-center mb-6">
            <div className="inline-block bg-gradient-to-r from-emerald-500 to-teal-500 text-white px-6 py-2 rounded-full font-bold text-sm shadow-md">
              ⚡ PARALLEL EXECUTION ENGINE
            </div>
          </div>

          {/* Module Grid */}
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4 mb-8 relative">
            {/* Connecting lines background */}
            <div className="absolute inset-0 flex items-center justify-center pointer-events-none">
              <div className="w-full h-px bg-gradient-to-r from-transparent via-slate-300 to-transparent"></div>
            </div>

            {modules.map((module, idx) => {
              const Icon = module.icon;
              return (
                <div 
                  key={idx}
                  className="bg-white border-2 rounded-xl p-4 shadow-md hover:shadow-xl transition-all duration-300 hover:-translate-y-1 relative z-10"
                  style={{ borderColor: module.color }}
                >
                  <div className="flex items-start gap-3 mb-3">
                    <div 
                      className="p-2 rounded-lg"
                      style={{ backgroundColor: `${module.color}20` }}
                    >
                      <Icon size={24} style={{ color: module.color }} />
                    </div>
                    <div className="flex-1">
                      <div className="font-bold text-slate-800">{module.name}</div>
                      <div className="text-xs text-slate-500 leading-tight">{module.desc}</div>
                    </div>
                  </div>
                  
                  {/* Progress indicator */}
                  <div className="space-y-1.5">
                    <div className="flex items-center gap-2 text-xs">
                      <CheckCircle size={14} className="text-green-500" />
                      <span className="text-slate-600">Analysis Complete</span>
                    </div>
                    <div className="h-1.5 bg-slate-100 rounded-full overflow-hidden">
                      <div 
                        className="h-full rounded-full transition-all duration-1000"
                        style={{ 
                          width: '100%',
                          backgroundColor: module.color
                        }}
                      ></div>
                    </div>
                  </div>
                </div>
              );
            })}
          </div>

          {/* Convergence Arrow */}
          <div className="flex flex-col items-center my-6">
            <div className="flex items-center gap-4">
              <div className="h-px w-24 bg-gradient-to-r from-transparent to-slate-300"></div>
              <div className="text-2xl">▼</div>
              <div className="h-px w-24 bg-gradient-to-l from-transparent to-slate-300"></div>
            </div>
          </div>

          {/* Integration Layer */}
          <div className="bg-gradient-to-r from-indigo-500 to-purple-500 rounded-xl p-6 shadow-lg mb-8">
            <div className="flex items-center justify-center gap-3 mb-4">
              <Database className="text-white" size={28} />
              <span className="text-white font-bold text-xl">Data Integration & Processing</span>
            </div>
            <div className="grid grid-cols-1 md:grid-cols-3 gap-3 text-sm text-white">
              <div className="bg-white/20 rounded-lg p-3 text-center backdrop-blur">
                <div className="font-semibold mb-1">Results Consolidation</div>
                <div className="text-xs text-indigo-100">Merging module outputs</div>
              </div>
              <div className="bg-white/20 rounded-lg p-3 text-center backdrop-blur">
                <div className="font-semibold mb-1">Lineage Mapping</div>
                <div className="text-xs text-indigo-100">Reference database matching</div>
              </div>
              <div className="bg-white/20 rounded-lg p-3 text-center backdrop-blur">
                <div className="font-semibold mb-1">Quality Control</div>
                <div className="text-xs text-indigo-100">Data validation & cleanup</div>
              </div>
            </div>
          </div>

          {/* Output Stage */}
          <div className="flex flex-col items-center">
            <div className="h-8 w-0.5 bg-slate-300 mb-4"></div>
            <div className="bg-gradient-to-r from-green-500 to-emerald-600 text-white px-8 py-4 rounded-xl shadow-lg mb-4 w-full max-w-2xl">
              <div className="flex items-center justify-center gap-3 mb-3">
                <FileText size={28} />
                <span className="font-bold text-xl">Comprehensive Report Generation</span>
              </div>
              
              {/* Output Formats */}
              <div className="grid grid-cols-3 gap-3 mt-4">
                {outputs.map((output, idx) => {
                  const Icon = output.icon;
                  return (
                    <div 
                      key={idx}
                      className="bg-white/20 backdrop-blur rounded-lg p-3 text-center"
                    >
                      <Icon className="mx-auto mb-2" size={24} />
                      <div className="font-semibold text-sm">{output.format}</div>
                      <div className="text-xs text-green-100 mt-1">Report Format</div>
                    </div>
                  );
                })}
              </div>
            </div>

            {/* Final Output Directory */}
            <div className="bg-slate-100 border-2 border-slate-300 rounded-lg px-6 py-3 text-center">
              <div className="font-mono text-sm text-slate-700 font-bold">
                📁 Staphscope/Staphscope_final_report/
              </div>
              <div className="text-xs text-slate-500 mt-1">
                Complete analysis results with typing data, resistance profiles, and virulence factors
              </div>
            </div>
          </div>
        </div>

        {/* Footer Statistics */}
        <div className="mt-6 grid grid-cols-3 gap-4 text-center">
          <div className="bg-white rounded-lg p-4 shadow border border-slate-200">
            <div className="text-2xl font-bold text-blue-600">~3 min</div>
            <div className="text-xs text-slate-600">Average Runtime</div>
          </div>
          <div className="bg-white rounded-lg p-4 shadow border border-slate-200">
            <div className="text-2xl font-bold text-purple-600">6 Modules</div>
            <div className="text-xs text-slate-600">Parallel Analyses</div>
          </div>
          <div className="bg-white rounded-lg p-4 shadow border border-slate-200">
            <div className="text-2xl font-bold text-green-600">3 Formats</div>
            <div className="text-xs text-slate-600">Report Outputs</div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default StaphScopeWorkflow;
--
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

* Types I–XIII + subtypes
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
