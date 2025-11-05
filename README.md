https://anaconda.org/bbeckley-hub/staphscope/badges/version.svg
badge
https://anaconda.org/bbeckley-hub/staphscope/badges/latest_release_date.svg
badge
https://anaconda.org/bbeckley-hub/staphscope/badges/latest_release_relative_date.svg
badge
https://anaconda.org/bbeckley-hub/staphscope/badges/platforms.svg
badge
https://anaconda.org/bbeckley-hub/staphscope/badges/license.svg
badge
https://anaconda.org/bbeckley-hub/staphscope/badges/downloads.svg


**StaphScope: Advanced Staphylococcus aureus Typing & Lineage Analysis Platform**



StaphScope is a comprehensive bioinformatics pipeline specifically designed for Methicillin-Resistant Staphylococcus aureus (MRSA) genomic analysis. This all-in-one tool provides complete characterization of S. aureus isolates through multiple typing methods, antimicrobial resistance profiling, virulence factor detection, and lineage analysis.
🎯 Purpose

   ** MRSA Surveillance: Track and characterize MRSA strains in clinical and research settings

   ** Outbreak Investigation: Identify related strains and transmission patterns
**
    Research Analysis: Comprehensive genomic profiling for academic studies
**
    Public Health: Support antimicrobial resistance monitoring programs**
**
✨ Features
🔬 Core Analysis Modules
Module	Description	Key Outputs
MLST	Multi-Locus Sequence Typing	Sequence Type (ST), Clonal Complex (CC)
spa Typing	Staphylococcal Protein A typing	spa type, repeat sequence
SCCmec	Staphylococcal Cassette Chromosome mec	SCCmec type, mec gene complex, ccr complex
AMR Profiling	Antimicrobial Resistance genes	Resistance genes, drug classes, mechanisms
ABRicate	Comprehensive resistance & virulence	Plasmid markers, virulence factors, resistance databases
Lineage Analysis	Strain lineage reference	HTML report with strain classification
🛡️ MRSA-Specific Capabilities
**
    SCCmec Typing: Accurate identification of SCCmec types I-XIII

   ** mecA/mecC Detection: Methicillin resistance determinant detection
**
    PVL Toxin Screening: Panton-Valentine Leukocidin gene detection

   ** AMR Profile: Comprehensive antimicrobial resistance pattern**
**
    Epidemic Clones: Identification of major MRSA clonal complexes (CC5, CC8, CC22, CC30, CC45)******
**

🚀 Quick Start
Installation
Option 1: Conda Installation (Recommended)

conda install -c bbeckley-hub -c bioconda -c conda-forge staphscope

Option 2: From Source

git clone https://github.com/bbeckley-hub/staphscope-typing-tool.git

cd staphscope-typing-tool

conda env create -f environment.yml

conda activate staphscope

pip install -e .

Basic Usage

# Single genome analysis
staphscope -i genome.fasta -o results/

# Batch analysis of multiple genomes
staphscope -i "*.fna" -o batch_results --threads 8

# Custom analysis (skip specific modules)
staphscope -i "MRSA_*.fasta" -o analysis --threads 16 --skip-lineage

**📋 Complete Usage
Command Line Options**
bash

usage: staphscope [-h] -i INPUT -o OUTPUT [-t THREADS] [--skip-amr] [--skip-abricate] [--skip-mlst]
                  [--skip-spa] [--skip-sccmec] [--skip-lineage]

StaphScope: Complete S. aureus Typing Pipeline

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input FASTA file(s) - can use glob patterns like "*.fna" or "*.fasta"
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

Supported Input Formats

   ** .fna, .fasta, .fa, .fn (standard FASTA formats)
**
    Single files or batch processing with glob patterns

   ** **Assembled genomes or contigs**

**🔧 Analysis Modules Details
1. MLST Analysis

    Tool: MLST (https://github.com/tseemann/mlst)

    Database: PubMedST S. aureus scheme

    Output: Sequence Type (ST), Clonal Complex (CC), allele profiles

    MRSA Relevance: Identifies major MRSA clonal complexes

2. spa Typing

    Tool: spaTyper

    Database: Ridom StaphType scheme

    Output: spa type, repeat sequence, Ridom classification

    MRSA Relevance: High-resolution strain discrimination

3. SCCmec Analysis

    Tool: SCCmecFinder with custom StaphScope enhancements

    Coverage: Types I-XI and subtypes

    Output: SCCmec type, mec complex, ccr complex, subtypes

    MRSA Relevance: Core MRSA characterization - identifies resistance cassette

4. AMR Profiling

    Tool: NCBI AMRFinderPlus v3.12.8

    Coverage: 5,000+ resistance genes across all drug classes

    Output: Resistance genes, drug classes, mechanisms, point mutations

    MRSA Relevance: Comprehensive resistance profile including:

        β-lactams (mecA, mecC, blaZ)

        Aminoglycosides

        Macrolides

        Tetracyclines

        Fluoroquinolones

5. ABRicate Analysis

    Tool: ABRicate v1.0.1

    Databases**:

        CARD: Comprehensive Antibiotic Resistance Database

        ResFinder: Acquired resistance genes

        NCBI: National Center for Biotechnology Information

        VFDB: Virulence Factors Database

        PlasmidFinder: Plasmid replicon types

        ARG-ANNOT: Antibiotic Resistance Gene Annotation

        MEGARES: Comprehensive resistance database

    Output: Plasmid markers, virulence factors, comprehensive resistance profile

6. Lineage Reference

    Output: Interactive HTML report

    Content: Strain classification, typing results summary, epidemiological data

    MRSA Relevance: Contextualizes isolates within global MRSA populations

📊 Output Structure
text

output_directory/
├── mlst_results/
│   ├── mlst_summary.csv
│   └── individual_sample_results/
├── spa_results/
│   ├── spa_typing_summary.csv
│   └── detailed_reports/
├── sccmec_results/
│   ├── sccmec_summary.csv
│   ├── sccmec_detailed.csv
│   └── sample_directories/
├── amr_results/
│   ├── amr_summary.csv
│   └── amr_detailed_results/
├── abricate_results/
│   ├── summary/
│   │   ├── card_summary.txt
│   │   ├── resfinder_summary.txt
│   │   ├── ncbi_summary.txt
│   │   ├── vfdb_summary.txt
│   │   └── plasmidfinder_summary.txt
│   └── individual_reports/
└── lineage_results/
    └── staphscope_lineage_reference.html

🦠 MRSA-Specific Analysis
Key MRSA Markers Detected
Category	Markers	Clinical Significance
Resistance	mecA, mecC, blaZ	β-lactam resistance
Virulence	PVL (lukS-PV, lukF-PV)	Necrotizing infections
Toxins	TSST-1, enterotoxins	Toxic shock syndrome
Adhesion	fnbA, fnbB, clfA, clfB	Biofilm formation
Major MRSA Clonal Complexes Identified

    CC5: USA100, NY/Japan clone

    CC8: USA300, Brazilian/Hungarian clone

    CC22: EMRSA-15, UK hospital clone

    CC30: EMRSA-16, Southwest Pacific clone

    CC45: Berlin clone, community-associated MRSA

🔬 Example Use Cases
Hospital Outbreak Investigation
bash

# Analyze outbreak isolates
staphscope -i "outbreak_*.fasta" -o outbreak_analysis --threads 16

# Key outputs:
# - Relatedness via MLST/spa typing
# - SCCmec type consistency
# - Resistance gene profile comparison
# - Virulence factor patterns

Surveillance Studies
bash

# Process surveillance isolates
staphscope -i "surveillance_*.fna" -o yearly_surveillance --threads 8

# Surveillance insights:
# - Predominant SCCmec types
# - Emerging resistance patterns
# - Clonal complex distribution
# - Temporal trends analysis

Research Characterization
bash

# Comprehensive isolate characterization
staphscope -i research_isolate.fasta -o complete_analysis

# Research outputs:
# - Complete typing profile
# - Resistance mechanism details
# - Virulence potential assessment
# - Epidemiological context

💾 System Requirements
Minimum

    CPU: 4 cores

    RAM: 4 GB

    Storage: 4 GB free space

    OS: Linux (Ubuntu/CentOS) or macOS

Recommended

    CPU: 2+ cores

    RAM: 4+ GB

    Storage: 10GB+ free space(NOT REQUIRED)- GOOD SPACE IS ENOUGH!!

    OS: Linux with Conda/mamba

Dependencies (Automatically Installed)

    Python: 3.8, 3.9, 3.10, 3.11 or 3.12

    Bioinformatics Tools: ABRicate, AMRFinderPlus, MLST, BLAST+

    Perl: Required for several analysis tools

    Databases: All required databases downloaded automatically

🐛 Troubleshooting
Common Issues

Database download failures:
bash

# Manual database update
amrfinder --update
abricate --setupdb card

Memory issues with large batches:
bash

# Process in smaller batches
staphscope -i "batch1_*.fna" -o results_batch1 --threads 4
staphscope -i "batch2_*.fna" -o results_batch2 --threads 4

Support

    Issues: https://github.com/bbeckley-hub/staphscope-typing-tool/issues

    Email: brownbeckley94@gmail.com

    Documentation: See docs/ directory for detailed documentation

📚 Citation

If you use StaphScope in your research, please cite:
bibtex

@software{staphscope2024,
  title = {StaphScope: Advanced Staphylococcus aureus Typing and Lineage Analysis Platform},
  author = {Brown Beckley},
  year = {2025},
  url = {https://github.com/bbeckley-hub/staphscope-typing-tool},
  note = {Comprehensive MRSA genomic analysis tool}
}

👨‍💻 Author

Brown Beckley

    University of Ghana Medical School

    Department of Medical Biochemistry

    Email: brownbeckley94@gmail.com

    GitHub: bbeckley-hub

📄 License

This project is licensed under the MIT License - see the LICENSE file for details.
🔗 Related Resources

    PubMedST: https://pubmlst.org/organisms/staphylococcus-aureus

    SCCmec Database: https://www.sccmec.org/

    NCBI AMRFinderPlus: https://github.com/ncbi/amr

    CARD Database: https://card.mcmaster.ca/

StaphScope: Empowering MRSA research through comprehensive genomic analysis 🧫🔬
