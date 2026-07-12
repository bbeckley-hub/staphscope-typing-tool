
# 🧬 StaphScope Module Descriptions

StaphScope runs a series of modules in sequence. Here's what each module does.

---

## Analysis Modules

### 1. FASTA QC
**Description**: Quality control for input assemblies.
**Outputs**:
- `FASTA_QC_summary.tsv`, `.html`, `.json`
- Individual per‑sample QC reports
**Key metrics**: Contig count, N50, GC%, total length

---

### 2. MLST (Multi‑Locus Sequence Typing)
**Description**: Determines sequence type (ST) using the PubMLST scheme for *S. aureus*.
**Implementation**: [Torsten Seemann's MLST](https://github.com/tseemann/mlst)
**Outputs**:
- `mlst_summary.tsv`, `.html`, `.json`
- Per‑sample MLST reports

---

### 3. spa Typing
**Description**: Determines spa types based on the repeat region of protein A.
**Implementation**: Original code by [mjsull](https://github.com/mjsull/spa_typing), modified by JFSanchezHerrero, using [Ridom SpaServer](https://spa.ridom.de/)
**Outputs**:
- `spa_summary.tsv`, `.html`, `.json`
- Per‑sample spa reports

---

### 4. SCCmec Typing
**Description**: Determines staphylococcal cassette chromosome *mec* (SCCmec) type.
**Implementation**: [SCCmecFinder](https://cge.cbs.dtu.dk/services/SCCmecFinder/) (Kaya et al.)
**Outputs**:
- `staphscope_summary.tsv`, `.html`
- Per‑sample SCCmec reports

---

### 5. Agr Typing (NEW in v1.3.0)
**Description**: Determines accessory gene regulator (agr) type (I‑IV) using AgrV.
**Implementation**: [AgrV](https://github.com/VishnuRaghuram94/AgrV) (Raghuram et al.)
**Outputs**:
- `agr_summary.tsv`, `.html`, `.json`
- Per‑sample agr reports

---

### 6. AMRFinderPlus
**Description**: Detects antimicrobial resistance genes and point mutations.
**Implementation**: [NCBI AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/)
**Outputs**:
- `staph_amrfinder_summary.tsv`, `.html`, `.json`
- `mutation_summary.tsv`, `.html`, `.json`
- Per‑sample AMR reports

---

### 7. ABRicate
**Description**: Comprehensive screening of resistance, virulence, and plasmid genes.
**Implementation**: [ABRicate](https://github.com/tseemann/abricate) (Torsten Seemann)
**Databases**:
- CARD
- ResFinder
- VFDB
- NCBI
- MEGARes
- ARG-ANNOT
- BacMet2
- PlasmidFinder
**Outputs**:
- `staph_*_abricate_summary.tsv` (one per database)
- Summary HTML and JSON reports

---

### 8. Lineage Reference
**Description**: Generates a reference database for S. aureus lineages.
**Outputs**:
- `staphscope_lineage_reference.html`

---

### 9. Comprehensive Report
**Description**: Combines MLST, spa, SCCmec, and agr results into a unified report.
**Outputs**:
- `staphscope_comprehensive_report.html`, `.json`, `.tsv`

---

### 10. Ultimate Reporter (Gene‑centric)
**Description**: Interactive HTML report showing all genes with their genome lists. Includes dynamic grouping by typing.
**Outputs**:
- `staphscope_ultimate_report.html`, `.json`
- CSV files for all gene tables

---

### 11. Sample‑Centric Reporter (NEW in v1.3.0)
**Description**: Interactive HTML report showing all genomes with their gene profiles. Each genome is displayed as an interactive box.
**Outputs**:
- `staphscope_ultimate_report.html`, `.json`
- Interactive per‑sample tables

---

### 12. Visualization (Optional)
**Description**: Generates plots and dashboards for exploratory data analysis.
**Outputs**:
- `STAPHSCOPE_VISUALIZATIONS/` directory
- PNG, SVG, PDF, and CSV files

---

## Skip Flags

| Module | Skip Flag |
|--------|-----------|
| FASTA QC | `--skip-fasta-qc` |
| MLST | `--skip-mlst` |
| spa | `--skip-spa` |
| SCCmec | `--skip-sccmec` |
| Agr | `--skip-agr` |
| AMR | `--skip-amr` |
| ABRicate | `--skip-abricate` |
| Lineage | `--skip-lineage` |
| Comprehensive | `--skip-comprehensive` |
| Sample‑centric | `--skip-sample-centric` |
| Visualization | `--skip-visualization` |
