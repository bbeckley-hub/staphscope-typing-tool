# 🤖 AI Integration Guide

You can upload your StaphScope reports to AI tools like **ChatGPT, Claude, or Gemini** for interactive analysis.

---

## Why Use AI with StaphScope?

- **Natural language queries** – ask questions in plain English.
- **Pattern recognition** – AI can spot epidemiological trends and clone associations.
- **Hypothesis generation** – uncover unexpected correlations for follow‑up.
- **Literature synthesis** – connect your findings with published literature.

---

## How to Upload Data

### Option 1: Upload the JSON File (Recommended)

File: `staphscope_ultimate_report.json` (located in `STAPHSCOPE_ULTIMATE_REPORTS/` or `STAPHSCOPE_ULTIMATE_SAMPLE_CENTRIC_REPORTS/`)

Most AI tools support file uploads.

### Option 2: Upload the HTML Report

The HTML report is self‑contained – upload it directly to the AI.

### Option 3: Copy‑Paste Tables

Copy specific tables from the report and paste them into the chat.

---

## Example Questions to Ask

### Epidemiology & Clonality
- "What are the most common MLST sequence types in this dataset?"
- "Which spa types are dominant in MRSA vs MSSA?"
- "Are there any ST‑spa‑SCCmec‑agr combinations with >2 isolates?"
- "What is the agr type distribution, and does it correlate with MRSA status?"

### Antimicrobial Resistance
- "How many samples carry mecA? What are their STs, spa types, and agr types?"
- "Are there any vanA/vanB positive samples? What is their SCCmec type?"
- "Which AMR genes co‑occur most frequently?"
- "Do any samples have combined β‑lactam + macrolide resistance?"

### Virulence & Toxins
- "Which samples carry PVL (lukF/S-PV)? Are they associated with specific STs or agr types?"
- "List all samples with TSST‑1 (tsst)."
- "Which enterotoxin genes are most prevalent?"
- "Is there a correlation between biofilm (ica) genes and MRSA?"

### Mutations & Biocides
- "What are the most frequent point mutations in gyrA or parC?"
- "Are there any linezolid‑related mutations (23S rRNA)?"
- "Which samples carry qac genes (disinfectant resistance)?"
- "Is mer (mercury) resistance linked to specific STs?"

### Advanced Queries
- "Write a summary of the resistance profile for ST5."
- "Compare virulence gene carriage between MRSA and MSSA."
- "Which clones have both AMR and virulence genes?"

---

## Ethical Guidelines

| Guideline | Description |
|-----------|-------------|
| **AI assists, experts decide** | Always interpret AI‑generated insights in your clinical/epidemiological context. |
| **No patient data** | Only upload aggregated, de‑identified genomic data. |
| **Verify critical findings** | Important resistance/virulence calls should be confirmed by human review or secondary tools. |
| **Be transparent** | When publishing, disclose that AI tools were used for exploratory analysis. |

---

## A (Mostly Serious) AI Survival Guide

- **If the AI says "I don't know"** – trust it. It's being honest.
- **If the AI says "It is widely known"** – ask for a reference. It may have made it up.
- **If the AI offers a completely novel evolutionary theory** – check if your coffee is spiked.
- **Remember:** AI won't take your job – but a microbiologist who knows how to use AI might!

---

## Suggested Prompt

When uploading, start with:

> "You are a bioinformatician analysing *Staphylococcus aureus* genomes. The attached data contains typing, AMR, virulence, BACMET, and mutation information. Answer my questions with references to the data."

---

## Feedback

If you find a great use case for AI with StaphScope, let us know! Email: brownbeckley94@gmail.com
