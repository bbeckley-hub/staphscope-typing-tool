# 🧬 Agr Typing Guide

The **accessory gene regulator (agr)** system is a quorum‑sensing circuit that controls virulence gene expression in *S. aureus*. Four agr types (I‑IV) are recognised.

---

## What is Agr Typing?

- **Agr** regulates the production of many virulence factors, including toxins, exoenzymes, and surface proteins.
- Different agr types are associated with different epidemiological profiles.
- Agr dysfunction (e.g., frameshift mutations) is common in clinical isolates and has been linked to persistent infections.

---

## How StaphScope Performs Agr Typing

StaphScope uses **[AgrV](https://github.com/VishnuRaghuram94/AgrV)** (Raghuram et al., 2022) to perform agr typing.

### Input
- FASTA files (same as other modules)

### Output
- `agr_summary.tsv` – TSV summary with sample, agr type, group, and status
- `agr_summary.html` – HTML summary report
- `agr_summary.json` – JSON summary
- Per‑sample reports in `agr_results/` subdirectories

---

## Agr Types and Their Significance

| Agr Type | Characteristics |
|----------|-----------------|
| **I** | Common in community‑acquired MRSA (CA‑MRSA) |
| **II** | Frequently associated with hospital‑acquired MRSA (HA‑MRSA) |
| **III** | Less common; often associated with specific lineages |
| **IV** | Emerging; associated with some livestock‑associated MRSA (LA‑MRSA) |

---

## Agr in the Reports

### Gene‑centric Report
- **Agr Tab**: Shows agr type distribution, samples by agr type, and all agr combinations (agr‑MLST, agr‑spa, agr‑SCCmec, agr‑MLST‑spa, agr‑MLST‑SCCmec, agr‑spa‑SCCmec, and four‑way).

### Sample‑centric Report
- Each sample box displays the agr type as a color‑coded badge.

### Grouping Feature
- You can group genomes by agr type in the AMR, Virulence, BACMET, Plasmids, and Mutations tabs.

---

## Agr in the Orchestrator

Agr typing is integrated into the orchestrator as `--skip-agr`:

```bash
# Run all modules including agr
staphscope -i "*.fna" -o results

# Skip agr typing
staphscope -i "*.fna" -o results --skip-agr
