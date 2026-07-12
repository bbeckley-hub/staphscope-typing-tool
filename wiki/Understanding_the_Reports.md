# 📊 Understanding StaphScope Reports

StaphScope generates two main interactive HTML reports, plus supporting data files.

---

## 1. Gene‑Centric Report (Ultimate Reporter)

**Location**: `Staphscope_final_report/STAPHSCOPE_ULTIMATE_REPORTS/staphscope_ultimate_report.html`

**Philosophy**: **Each gene** is displayed with **all genomes** that carry it.

### Tabs

| Tab | Description |
|-----|-------------|
| Summary | Overview of the dataset |
| Sample Overview | All samples with typing results |
| FASTA QC | Assembly quality metrics |
| MLST | ST distribution, ST‑spa, ST‑SCCmec, ST‑agr combinations |
| spa | spa type distribution, spa‑ST, spa‑SCCmec, spa‑agr combinations |
| SCCmec | SCCmec distribution, SCCmec‑ST, SCCmec‑spa, SCCmec‑agr combinations |
| MRSA Analysis | MRSA‑specific combinations with agr |
| agr Typing | Agr type distribution and all agr combinations |
| AMR | Gene‑centric AMR table with grouping |
| Virulence | Gene‑centric virulence table with grouping |
| BACMET | Biocide and heavy metal resistance genes |
| Plasmids | Plasmid replicon table |
| Mutations | Point mutations with grouping |
| Patterns | Triple/four‑way typing, gene co‑occurrence (top 500) |
| AI Guide | How to use with ChatGPT, Claude, Gemini |
| Call to Action | Global AMR context |
| Citation | How to cite StaphScope and its dependencies |
| Funding | Support information |
| Export | Download CSV files and JSON data |

### Dynamic Grouping

In the AMR, Virulence, BACMET, Plasmids, and Mutations tabs, you can group genomes by:
- MLST
- spa
- SCCmec
- agr
- Any combination (e.g., ST‑spa‑SCCmec‑agr)

---

## 2. Sample‑Centric Report

**Location**: `Staphscope_final_report/STAPHSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS/staphscope_ultimate_report.html`

**Philosophy**: **Each genome** is displayed as an interactive box with all its genes.

### Features
- Each sample shows MLST, spa, SCCmec, MRSA/MSSA status, and agr type as badges.
- Horizontally scrollable tables – no truncation.
- Filter by sample name or database.
- Per‑sample gene lists for AMR, Virulence, BACMET, Plasmids, and Mutations.

---

## 3. Comprehensive Report

**Location**: `Staphscope_final_report/staphscope_comprehensive_report.html`

A single‑page HTML report combining MLST, spa, SCCmec, and agr results.

---

## 4. CSV and JSON Data

All data is available as CSV and JSON for downstream analysis:

- `sample_overview.csv`
- `amr_genes.csv`
- `virulence_genes.csv`
- `bacmet_genes.csv`
- `plasmid_replicons.csv`
- `mutations.csv`
- `fasta_qc.csv`
- `pattern_discovery.csv`
- `staphscope_ultimate_report.json` (complete dataset)

---

## 5. Visualization

**Location**: `STAPHSCOPE_VISUALIZATIONS/`

Plots and dashboards:
- MLST distribution
- spa distribution
- SCCmec distribution
- Database comparison
- Virulence gene analysis
- And more (PNG, SVG, PDF formats)

---

## 6. Log File

**Location**: `staphscope_run.log`

Detailed logs of the entire analysis – useful for debugging.
