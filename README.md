<p align="center">
  <img src="https://raw.githubusercontent.com/bbeckley-hub/staphscope-typing-tool/main/images/staphscope.png" alt="StaphScope Banner" width="100%">
</p>

<div align="center">

# 🔬 StaphScope

### **A species-optimized computational pipeline for rapid, accessible *Staphylococcus aureus* genotyping and surveillance**

#### **Complete MRSA/MSSA genomic analysis in minutes — not hours**

[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)
![Version](https://anaconda.org/bioconda/staphscope/badges/version.svg)
![Last Updated](https://anaconda.org/bioconda/staphscope/badges/latest_release_relative_date.svg)
![Platforms](https://anaconda.org/bbeckley-hub/staphscope/badges/platforms.svg)
![License](https://anaconda.org/bbeckley-hub/staphscope/badges/license.svg)
[![DOI](https://img.shields.io/badge/DOI-10.1186%2Fs12864--026--12609--x-blue)](https://doi.org/10.1186/s12864-026-12609-x)

[![Docker Pulls](https://img.shields.io/docker/pulls/bbeckleyhub/staphscope)](https://hub.docker.com/r/bbeckleyhub/staphscope)
[![Docker Image Size](https://img.shields.io/docker/image-size/bbeckleyhub/staphscope/latest)](https://hub.docker.com/r/bbeckleyhub/staphscope)
[![Docker Version](https://img.shields.io/docker/v/bbeckleyhub/staphscope?sort=semver)](https://hub.docker.com/r/bbeckleyhub/staphscope)
[![Contributions Welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg)](#)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-Profile-0A66C2?style=flat&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/brown-beckley-190315319)
[![Stage](https://img.shields.io/badge/status-active-brightgreen)](#)
![Conda Downloads](https://img.shields.io/conda/dn/bioconda/staphscope?label=Conda%20Downloads)

[![Powered by 🧠](https://img.shields.io/badge/powered%20by-science%20🔬-purple)](https://github.com/bbeckley-hub/staphscope-typing-tool)
[![Coffee](https://img.shields.io/badge/built%20with-%E2%98%95%20coffee-orange)](https://github.com/bbeckley-hub/staphscope-typing-tool)
[![Made with ❤️](https://img.shields.io/badge/made%20with-%E2%9D%A4%EF%B8%8F-red)](https://github.com/bbeckley-hub/staphscope-typing-tool)
[![Open Source Love](https://badges.frapsoft.com/os/v1/open-source.svg?v=103)](https://github.com/ellerbrock/open-source-badges/)
[![Made for Research](https://img.shields.io/badge/made%20for-Research-0066cc.svg)](https://github.com/bbeckley-hub/staphscope-typing-tool)

[![Documentation](https://img.shields.io/badge/docs-mkdocs-526CFE?logo=materialformkdocs)](https://bbeckley-hub.github.io/staphscope-typing-tool)
[![RST Badge](https://img.shields.io/badge/documentation-RST-4CAF50.svg)](https://www.sphinx-doc.org/)
[![Last Commit](https://img.shields.io/github/last-commit/bbeckley-hub/staphscope-typing-tool)](https://github.com/bbeckley-hub/staphscope-typing-tool/commits)
[![Contributors](https://img.shields.io/github/contributors/bbeckley-hub/staphscope-typing-tool)](https://github.com/bbeckley-hub/staphscope-typing-tool/graphs/contributors)
[![Security: bandit](https://img.shields.io/badge/security-bandit-yellow.svg)](https://github.com/PyCQA/bandit)

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336)](https://pycqa.github.io/isort/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![CI](https://img.shields.io/github/actions/workflow/status/bbeckley-hub/staphscope-typing-tool/ci.yml?branch=main&label=CI)](https://github.com/bbeckley-hub/staphscope-typing-tool/actions)
[![Tests](https://img.shields.io/badge/tests-passing-brightgreen.svg)](https://github.com/bbeckley-hub/staphscope-typing-tool/tests)
[![Open in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/bbeckley-hub/staphscope-typing-tool)

[![Speed](https://img.shields.io/badge/Speed-14%20min%2F24%20samples-FF6D00.svg)](https://github.com/bbeckley-hub/staphscope-typing-tool#performance-benchmarks)
[![StaphScope Web](https://img.shields.io/badge/🌐-Try%20StaphScope%20Web-4285F4?logo=googlechrome&logoColor=white)](https://eskape.bio)
[![MRSA|MSSA](https://img.shields.io/badge/Classifies-MRSA%20%7C%20MSSA-00BCD4.svg)](https://github.com/bbeckley-hub/staphscope-typing-tool)
[![Lineages](https://img.shields.io/badge/Lineages-50%20major-673AB7.svg)](https://github.com/bbeckley-hub/staphscope-typing-tool)
[![AI Ready](https://img.shields.io/badge/AI-Ready%20Reports-00C853.svg)](https://github.com/bbeckley-hub/staphscope-typing-tool#ai-integration-guide)

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Conda](https://img.shields.io/badge/conda-✓-green.svg)](https://docs.conda.io/en/latest/)
[![MIT License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![GitHub Issues](https://img.shields.io/badge/issues-Open-orange)](https://github.com/bbeckley-hub/staphscope-typing-tool/issues)
[![GitHub Stars](https://img.shields.io/github/stars/bbeckley-hub/staphscope-typing-tool)](https://github.com/bbeckley-hub/staphscope-typing-tool/stargazers)
[![Sample Report](https://img.shields.io/badge/📊-View_Sample_Report-blue)](https://htmlpreview.github.io/?https://bbeckley-hub.github.io/staphscope-typing-tool/#summary)
![Profile Views](https://komarev.com/ghpvc/?username=bbeckley-hub&label=Profile%20Views&color=0e75b6&style=flat)
[![Google Scholar](https://img.shields.io/badge/Google%20Scholar-Profile-4285F4?style=flat&logo=googlescholar&logoColor=white)](https://scholar.google.com/citations?user=CYNOsqIAAAAJ&hl=en)

![GitHub stats](https://github-readme-stats.vercel.app/api?username=bbeckley-hub&show_icons=true&theme=radical)
![Top Langs](https://github-readme-stats.vercel.app/api/top-langs/?username=bbeckley-hub&layout=compact&theme=radical)
[![GitHub Streak](https://streak-stats.demolab.com?user=bbeckley-hub&theme=radical&date_format=j%20M%5B%20Y%5D)](https://git.io/streak-stats)

**Two ways to use StaphScope:**  
🖥️ **Command-line tool** for high-throughput, local analysis  
🌐 **StaphScope Web** for non-bioinformaticians – [https://eskape.bio](https://eskape.bio)

</div>

---

## 🎉 **What's New in v2.1.0 — SAURUS Integration & Database Refresh**

**SAURUS is live.** We are proud to introduce [SAURUS (Staphylococcus aureus Unified Resource for Universal Surveillance)](https://bbeckley-hub.github.io/saurus/), a free, interactive, browser-based database that solves the biggest bottleneck in *S. aureus* genomics: **retrieving the exact subset of genomes you need.**

### The Problem
You want all ST5 MRSA bloodstream isolates from the United States carrying with spa type t002, SCCmec type II and agr type II. In existing repositories, this means downloading multiple compressed tables, writing custom scripts to cross-reference typing, gene presence, and metadata, and manually merging files. It is slow, error-prone, and requires coding skills that many clinical and public health researchers do not have.

### The Solution
SAURUS puts the entire public *S. aureus* repository in your browser. With **134,299 unique genomes** integrated across **eleven typing, gene-annotation, and epidemiological modules**, you can filter by MLST, spa, SCCmec (two independent callers), agr, capsule, MRSA status, AMR genes, virulence factors, plasmid replicons, country, isolation source, host, and collection year — all in real time.

**Click. Filter. Export.** The result is a plain-text accession list, ready for:
`datasets download genome accession --inputfile list.txt`

### SAURUS by the Numbers
| Metric | Value |
|---|---|
| **Genomes** | **134,299** unique *S. aureus* accessions |
| **Typed Genomes** | **98,138** with de novo annotations |
| **Typing Modules** | 11 (MLST, spa, SCCmec-CGE, SCCmec-RPet, agr, capsule, MRSA, AMR, Virulence, Plasmid, Geography) |
| **Countries** | 139 |
| **MRSA Rate** | 59.2% |
| **Unique STs** | 2,116 |
| **Unique spa Types** | 3,932 |

> **No installation. No account. No API key.** SAURUS runs entirely over HTTPS on GitHub Pages, making it immune to server crashes, API rate limits, and funding deprecation.

**[🚀 Explore SAURUS Now →](https://bbeckley-hub.github.io/saurus/)**

### StaphScope v2.1.0 Improvements
- **Refreshed spa Database** – Bundled Ridom SpaServer snapshot updated to **22,848 spa types** (865 repeats, 472,643 strains, 185 countries).
- **Report Level Control** – New `--report-level` flag (`final`, `summary`, `full`) prevents disk exhaustion on 300+ sample cohorts by managing intermediate files.
- **MLST Database Management** – New `--pull-mlst-db` and `--update-mlst-db` commands for streamlined database setup.
- **Performance Fixes** – Resolved AMR percentage display anomalies and accelerated the Compare/Cluster tab for large cohorts via precomputed similarity matrices and canvas-based heatmaps.

---

## 🎯 **Overview**

**StaphScope** is an automated, locally-executable computational pipeline designed specifically for comprehensive *Staphylococcus aureus* genomic surveillance. It addresses the critical bottleneck in MRSA research by integrating **ten essential genotyping and profiling methods** into a single, cohesive workflow.

### The Problem
- **Fragmented Bioinformatics**: Traditional MRSA analysis requires 5+ separate tools with conflicting dependencies.
- **Resource Barriers**: Web-based services need constant internet and raise data privacy concerns.
- **Time Constraints**: Generalist platforms take hours; outbreaks need answers in minutes.
- **Interpretation Challenges**: Raw data without epidemiological context limits actionable insights.

### Our Solution
StaphScope delivers:
- ✅ **Single-command installation** via Conda.
- ✅ **10–14 minute complete analysis** (24 samples, 16 cores).
- ✅ **100% local execution** with data privacy.
- ✅ **Interactive HTML reports** with epidemiological context.
- ✅ **Automated MRSA/MSSA classification** with confidence scoring.
- ✅ **Web-based interface** for non-bioinformaticians.

**Perfect for**: Clinical labs, outbreak investigations, research studies, and public health surveillance.

---

## ✨ **Key Features**

### 🔬 **Core Analytical Modules**

| Module | 🎯 Purpose | 📊 Key Outputs | ⚡ Speed |
|--------|------------|----------------|----------|
| **FASTA QC** | Quality control + fastANI species confirmation | N50, GC%, contig stats, ANI %, species | <30 sec |
| **MLST Typing** | Phylogenetic classification via 7 housekeeping genes | ST, CC, allele profiles | <1 min |
| ***spa* Typing** | Hypervariable region analysis of protein A gene | *spa* type, repeat patterns | <1 min |
| **SCC*mec* Typing (CGE)** | Methicillin resistance cassette characterization | SCC*mec* type (I-XIII), confidence scores | 1-2 min |
| **SCC*mec* Typing (RPet)** | Independent second-opinion caller | SCC*mec* type + subtype (IVa, IIb, ...) | 1-2 min |
| **Capsule Typing** | cap5/cap8 serotype determination | Type 5 / Type 8, completeness score | <30 sec |
| **Agr Typing** | Accessory gene regulator (agr) type determination | agr type I-IV, group, status | <1 min |
| **AMR Profiling** | Comprehensive resistance gene detection | 5,000+ AMR genes, risk categorization | 2-3 min |
| **ABRicate Screening** | Multi-database virulence/plasmid detection (10 DBs) | Plasmid replicons, virulence factors | 3-4 min |
| **MGE Profiling** | mobileOG-db per-sample category counts | 10 categories, key MGE signatures | 2-3 min |
| **Visualization Suite** | Interactive dashboard + publication plots | HTML dashboard, PNG, PDF, SVG, ZIP bundle | 1-2 min |
| **Lineage Database** | Global epidemiological context | 50 major lineages, geographical distribution | Instant |

---

### 📊 **Sample Integrated Reports**

Every report in one interactive hub — browse, search, preview, and open any StaphScope output. **[Open the Live Report Hub →](https://bbeckley-hub.github.io/staphscope-typing-tool/)**

| Report Type | What it shows | Direct link |
|---|---|---|
| **📄 Comprehensive Report** | Master typing table — MLST, spa, SCCmec (CGE + RPet + subtype), agr, capsule, MRSA/MSSA | [![View](https://img.shields.io/badge/📄-Comprehensive_Report-FF9800)](https://bbeckley-hub.github.io/staphscope-typing-tool/staphscope_comprehensive_report.html) |
| **🧬 Gene-Centric Report** | Every AMR / virulence / BACMET / plasmid / mutation gene, with all genomes that carry it | [![View](https://img.shields.io/badge/🧬-Gene_Centric_Report-red)](https://bbeckley-hub.github.io/staphscope-typing-tool/staphscope_ultimate_gene_centric_report.html) |
| **📋 Sample-Centric Report** | Per-isolate lazy-loaded boxes with a full Compare &amp; Cluster tab | [![View](https://img.shields.io/badge/📋-Sample_Centric_Report-purple)](https://bbeckley-hub.github.io/staphscope-typing-tool/staphscope_ultimate_sample_centric_report.html) |
| **📈 Interactive Dashboard** | Plotly dashboard — Overview, Typing, QC, AMR, Virulence, MGE, Alerts, Story, Compare | [![View](https://img.shields.io/badge/📈-Interactive_Dashboard-16A085)](https://bbeckley-hub.github.io/staphscope-typing-tool/staphscope_dashboard.html) |
| **🔬 MLST Summary** | Sequence types + allele profiles across your cohort | [![View](https://img.shields.io/badge/🔬-MLST_Summary-10b981)](https://bbeckley-hub.github.io/staphscope-typing-tool/mlst_summary.html) |
| **🧬 spa Summary** | spa types + repeat patterns | [![View](https://img.shields.io/badge/🧬-spa_Summary-10b981)](https://bbeckley-hub.github.io/staphscope-typing-tool/spa_summary.html) |
| **🛡️ SCCmec (CGE)** | Cassette type, *mec* class, *ccr* complex | [![View](https://img.shields.io/badge/🛡️-SCCmec_CGE-10b981)](https://bbeckley-hub.github.io/staphscope-typing-tool/staphscope_sccmec_cge_summary.html) |
| **🛡️ SCCmec (RPet)** | Independent caller + fine-grained subtype | [![View](https://img.shields.io/badge/🛡️-SCCmec_RPet-10b981)](https://bbeckley-hub.github.io/staphscope-typing-tool/staphscope_sccmec_rpet_summary.html) |
| **🧬 agr Summary** | agr types I–IV (AgrVATE) | [![View](https://img.shields.io/badge/🧬-agr_Summary-10b981)](https://bbeckley-hub.github.io/staphscope-typing-tool/agr_summary.html) |
| **📏 FASTA QC** | N50 / GC / contigs + fastANI species confirmation | [![View](https://img.shields.io/badge/📏-FASTA_QC-06b6d4)](https://bbeckley-hub.github.io/staphscope-typing-tool/FASTA_QC_summary.html) |
| **💊 Capsule Typing** | cap5 / cap8 serotype with completeness score | [![View](https://img.shields.io/badge/💊-Capsule_Typing-06b6d4)](https://bbeckley-hub.github.io/staphscope-typing-tool/staphscope_capsule_summary.html) |
| **📱 MGE Profile** | mobileOG-db counts across 10 functional categories | [![View](https://img.shields.io/badge/📱-MGE_Profile-14b8a6)](https://bbeckley-hub.github.io/staphscope-typing-tool/staphscope_mge_summary.html) |

> **Note:** Reports may take a few seconds to load. For the best experience, download the HTML files and open them locally — every report is fully self-contained and offline-capable.

---

### 🛡️ **MRSA-Specific Innovations**
- **Automated MRSA Classification**: Based on concurrent *mecA/mecC* + SCC*mec* detection.
- **Dual SCCmec Caller Consensus**: CGE + RPet, with automatic naming-format normalization.
- **Clinical Gene Flagging**: Automatic highlighting of PVL, enterotoxins, *van* genes.
- **Risk Assessment**: Categorizes genes as 'Critical Risk' (e.g., *mecA*, *vanA*) or 'High Risk'.
- **Cross-Genome Pattern Discovery**: Summarizes gene frequencies across entire sample sets.
- **Cluster Detection**: Union-find algorithm on typing data to auto-detect transmission clusters.
- **Curated Lineage Database**: 50 major lineages with HA-MRSA, CA-MRSA, LA-MRSA classifications.

---

## 🔥 **Feature Deep-Dives**

<details>
<summary><b>⚖️ The Compare Tab — Because "I Think They're the Same Bug" Isn't Good Enough (Click to expand)</b></summary>
<br>

Let's be honest about how most outbreak investigations actually happen. Someone has two genomes. They run typing on both. Then they open two HTML files, scroll up and down, mentally diff the MLST, the spa type, the SCCmec, the agr, the capsule, the AMR genes, the virulence genes, the plasmid replicons… and at some point they say *"yeah, I think they're the same."*

That is not a methodology. That's an eyeball.

**The Compare tab replaces the eyeball.** Pick any two isolates, click a button, and get:
1. **A verdict banner** — 🔴 near-identical (≥95% similarity), 🟠 high (≥85%), 🟡 moderate (≥60%), 🔵 distinct
2. **A similarity gauge** — an animated colored ring showing overall match percentage
3. **A typing table** — every field side by side, green pill for match, red pill for diff
4. **Gene content split three ways** — 🟢 shared, 🔵 only in A, 🟣 only in B — per category (AMR, Virulence, BACMET, Plasmids, Mutations)
5. **Full profile expanders** — click any category to see both samples' complete gene lists side by side

**Why this matters for research:**
- **Outbreak confirmation:** Two isolates from the same ward with identical MLST + spa + SCCmec + agr are very likely a transmission pair. Your infection control team will actually do something about that.
- **Discordance analysis:** An ancestor and a resistant descendant look identical in *every* field except one — and that one field tells you which mobile element moved.
- **QC checks:** Two genomes from the same sequencing run that come back as "identical" — one is probably a duplicate. Better to catch that now than publish it.
- **Teaching:** If you have ever tried to explain "why isn't this MSSA?" to a clinician, showing them the Compare tab takes 30 seconds and produces actual agreement.

And when you have more than two samples? **Switch to Cluster Mode.** Automatic pairwise similarity matrix, colour-coded heatmap, automatic cluster detection (union-find algorithm on your own typing data), and a per-cluster summary card showing mean intra-cluster similarity. It's the difference between "we have 40 isolates" and "we have 40 isolates in 6 clusters, 3 of which are outbreak-suspicious."

**A caveat we're proud to state openly:** "Identical typing = transmission pair" is a **hypothesis generator, not a confirmation**. Two isolates can share MLST + spa + SCCmec + agr + capsule and still be unrelated — those markers evolve too slowly for hospital-level resolution. The gold standard remains core-genome SNP distance, which StaphScope doesn't compute (and shouldn't, at this scope). We use the phrase *"possible transmission pair"* deliberately. If we ever lose that hedge, shout at us.
</details>

<details>
<summary><b>🦠 The Virulence Expander — Because 400 Gene Names Won't Fit in a Cell (Click to expand)</b></summary>
<br>

Before v2.0.0, if you wanted to see which virulence genes a specific isolate carried, you had to open the Virulence tab, find the row for each gene, and mentally tick off which genomes were listed. For a cohort of 40 samples, this was a job you did once and never wanted to do again.

Now: **Sample Overview** has a **Virulence column with a click-to-expand count**. Click the number, and a panel slides open showing every VFDB gene that isolate carries — colour-coded, in the same visual language as everything else in the report.

**Why this is more than cosmetic:**
- **Case-by-case drill-down:** When a clinician asks *"what does this isolate actually carry?"*, you have one answer in one place.
- **Teaching tool:** Because explaining PVL + TSST-1 + enterotoxins is a lot easier when you can expand and show the gene list in real time.
- **Screening:** Sort by the count column and the top five highest-virulence isolates float to the top. Useful for prioritising further work.
- **No page reload:** This is a `details/summary` element. The whole report stays one file. Offline. No server. No JavaScript framework. Just HTML.
</details>

<details>
<summary><b>📱 The MGE Tab — Where Resistance Gets a Postal Address (Click to expand)</b></summary>
<br>

You have AMR genes. You have virulence genes. You know *what* a genome carries. But *how did it get there?*

**Mobile genetic elements are the answer.** Plasmids, prophages, transposons, integrative conjugative elements — these are the physical vehicles that shuttle resistance and virulence genes between unrelated bacteria. If you don't understand the MGE landscape, you're looking at the mail without ever seeing the postal system.

**What the MGE tab gives you:**
- **Per-sample counts** across ten functional categories (Integrase, Transfer, Stability, Phage, Replication, IS-associated, ICE-associated, Plasmid-associated, Phage-associated, and Key MGE-signatures).
- **Aggregate profiles by typing** — mean MGE load per ST, per SCCmec type, per agr type, per capsule type. High MGE burden in one lineage but not another is a real signal.
- **Top-N views** — the five isolates with the most IS-associated proteins, plasmid proteins, or high-confidence MGE signatures.
- **Education cards** — because "what is an IS element and why do I care?" is a question that deserves a proper answer.

**What it helps with:**
- **Resistance dissemination:** A plasmid-associated gene cluster in two unrelated STs is a smoking gun for horizontal transfer.
- **Virulence acquisition:** PVL and TSST-1 are often prophage-borne.
- **Persistent infections:** Strong stability systems let plasmids survive without selection.
- **Outlier detection:** An isolate with 2× the cohort's MGE load is either something new, or something broken. Either way you should look.

**A word of honesty:** mobileOG-db counts protein *signatures*, not reconstructed elements. Multiple proteins from the same element produce multiple hits, and some categories also fire on chromosomal genes. That's why we use the word *"signatures"* everywhere. If you want element-level claims, you need co-localization analysis — which is a bioinformatics project, not a column in a dashboard.
</details>

<details>
<summary><b>💊 The Capsule Tab — Because Phagocytes Have Feelings Too (Click to expand)</b></summary>
<br>

*S. aureus* wears a sugar coat. Specifically, it produces one of two dominant capsular polysaccharides: **Type 5** or **Type 8**. This coat is what stops your neutrophils from phagocytosing it efficiently. It's also one of the primary targets for investigational *S. aureus* vaccines.

**Why typing the capsule matters:**
- **Vaccine coverage:** If a candidate vaccine targets Type 5 and your entire cohort is Type 8, that vaccine is going to have a bad day.
- **Lineage correlation:** Type 5 is common in MSSA. Type 8 dominates MRSA ST8 and ST239 lineages.
- **Immune-evasion potential:** Capsule-positive isolates resist phagocytosis better than capsule-negative ones.
- **Cross-tabulation:** The Capsule tab cross-tabs with MLST, spa, agr, SCCmec-CGE, and SCCmec-Subtype.

**The visual touch:** In the Sample Overview table, **Type 5 is green** and **Type 8 is red** everywhere they appear.
</details>

<details>
<summary><b>🛡️ The SCCmec Tab — Because One Opinion Is Never Enough (Click to expand)</b></summary>
<br>

SCCmec typing tells you which *mec* cassette a MRSA isolate carries — Types I through XIII, with a fine-grained **subtype** on top (IIa, IVc, etc.). It's the single most important structural marker of MRSA lineage evolution.

**The tab shows you CGE SCCmecFinder and RPet's `sccmec` side by side**, with a normalizer that understands they use different naming conventions (`SCCmec_type_II(2A)` from CGE and `Type II(2A)` from RPet are the same cassette). You get:
- **CGE SCCmec type** distribution, cross-tabbed with ST, spa, agr, and capsule.
- **RPet SCCmec type** distribution, same cross-tabs.
- **SCCmec subtype** distribution — the finest-grained layer available.

**Why two callers matter:**
- **Cross-validation:** When both callers agree on a cassette type, you have high confidence. When they disagree, you have a genuine signal.
- **Subtype resolution:** RPet's subtype call (IVa vs IVc, IIa vs IIb) is often finer than CGE's.
- **Publication-ready:** Reviewers ask "which SCCmec caller did you use?" Having a two-caller consensus answer is a much stronger response.
</details>

---

## 🧰 **StaphScope Toolkit – iTOL Metadata Made Easy**

Data wrangling got you down? We built a companion toolkit that automatically generates iTOL-ready annotation files from your StaphScope outputs — **no Excel gymnastics required.**

- ✅ Converts your StaphScope CSVs into **iTOL colour strips** and **binary matrices**
- ✅ Handles **MLST, spa, SCCmec, agr, MRSA status, AMR genes, virulence factors, mutations, BacMet, and plasmid replicons**
- ✅ **Auto-selects top genes/mutations** so you never get empty files

```bash
# Clone the toolkit
git clone https://github.com/bbeckley-hub/staphscope-toolkit.git
cd staphscope-toolkit

# Install dependencies
pip install -r requirements.txt

# Generate iTOL metadata from your StaphScope run
python staphscope_itol.py --input_dir /path/to/Staphscope_final_report
```

> 📖 **[Full documentation & examples →](https://github.com/bbeckley-hub/staphscope-toolkit)**

---

## 🌐 **StaphScope Web Platform**

StaphScope is a core component of the **ESKAPE AMR Platform** – a unified web suite for species‑optimized genomic surveillance. 

- ✅ **Drag‑and‑drop file upload** (single, multiple, or ZIP archives)
- ✅ **Module selection** – choose which analyses to run
- ✅ **Real‑time progress monitoring** with live logs
- ✅ **Download all results as a single ZIP**
- ✅ **No installation required** – works in any modern browser

> 🌐 **Try StaphScope Web today:** [https://eskape.bio](https://eskape.bio)  
> 📦 **Web Repository:** [https://github.com/bbeckley-hub/eskape-web-platform](https://github.com/bbeckley-hub/eskape-web-platform)

*Note: The web version limits uploads to 10 files per job. For larger datasets, use the CLI.*

---

## ⚡ **Quick Start (CLI)**

### Install in 60 seconds
```bash
# Method 1: Conda (Recommended)
conda create -n staphscope -c conda-forge -c bioconda staphscope -y
conda activate staphscope

# Method 2: From source
git clone https://github.com/bbeckley-hub/staphscope-typing-tool.git
cd staphscope-typing-tool
conda env create -f environment.yml
conda activate staphscope
pip install -e .
```

### Run your first analysis
```bash
# Single genome
staphscope -i genome.fasta -o results/

# Batch processing (24 genomes)
staphscope -i "*.fna" -o batch_results --threads 16
# Complete in ~14 minutes! 🎉
```

### Update Databases (Recommended)
```bash
# ABRicate databases
abricate --setupdb

# AMR database
staphscope --update-amr-db

# MLST database
staphscope --pull-mlst-db
```
---
## 🐳 **Docker & Singularity Usage – avoid the padlock 🔓**

Docker and Singularity users get **everything pre-installed** — no first-run database setup, no `--pull-mlst-db`, no `--update-amr-db`, no `abricate --setupdb`. Every database (MLST scheme, AMRFinderPlus, ABRicate DBs, mobileOG-db) and every tool (BLAST, DIAMOND, Prodigal, fastANI, agrVATE) is baked into the image. **Pull and run.**

> **Conda users:** you must run `staphscope --pull-mlst-db` and `staphscope --update-amr-db` once after install. Docker and Singularity users skip that entirely.

---

### 🐳 **Docker (laptops, workstations, cloud VMs)**

Docker runs containers as `root` by default, so files written to bind-mounted directories come out owned by `root:root`. **Fix: add `-u $(id -u):$(id -g)`.**

```bash
# Pull the image
docker pull bbeckleyhub/staphscope:2.0.0

# Run (no padlock, no sudo chown)
docker run --rm \
  -u $(id -u):$(id -g) \
  -v $(pwd):/data \
  bbeckleyhub/staphscope:2.0.0 \
  -i "/data/*.fna" -o /data/output -t 4
```

- `-u $(id -u):$(id -g)` → files owned by **you**
- `-v $(pwd):/data` → mounts your current dir as `/data` inside the container
- `-i "/data/*.fna"` → uses the **container path**, not your host path
- `-o /data/output` → lands at `./output/` on your host

<details>
<summary><b>Verify the image and open a debugging shell (Click to expand)</b></summary>

```bash
# Verify the container runs
docker run --rm bbeckleyhub/staphscope:2.0.0 --help

# Interactive shell (for debugging)
docker run --rm -it bbeckleyhub/staphscope:2.0.0 bash
```

**Image size:** ~1.7 GB compressed / ~8 GB on disk.
</details>

---

### 🔬 **Singularity / Apptainer (HPC clusters) — recommended for HPC**

**HPC users should use Singularity, not Docker.** HPC login and compute nodes don't run the Docker daemon; Singularity (or its successor Apptainer) is the standard. Same SIF format, rootless by default, and no UID juggling — it runs as your real user automatically.

<details>
<summary><b>1. Load the runtime (Click to expand)</b></summary>

Module names vary by cluster — check with `module avail`:

```bash
# Newer clusters — Apptainer
module load apptainer

# Older clusters — Singularity
module load containers/singularity/3.6.4

# Pulling the image ALSO needs mksquashfs from squashfs-tools
module load squashfs-tools/4.4

singularity --version   # or: apptainer --version
```
</details>

<details>
<summary><b>2. Pull the image (one-time, needs internet on login node) (Click to expand)</b></summary>

```bash
singularity pull staphscope_2.1.0.sif docker://bbeckleyhub/staphscope:2.1.0
```

Creates `staphscope_2.1.0.sif` (~1.7 GB) in the current directory. **`squashfs-tools` is only needed for the pull, not for running.**
</details>

<details>
<summary><b>3. Run the pipeline (Click to expand)</b></summary>

```bash
singularity run --bind $(pwd):/data \
  staphscope_2.1.0.sif \
  -i "/data/*.fna" -o /data/output -t 16
```

- `--bind $(pwd):/data` — HPC equivalent of Docker's `-v`
- **No `-u` needed** — Singularity runs as your real UID
- **No `-e HOME=` needed** — Singularity uses your real `$HOME`, so MLST and AMR caches persist between runs

> ⚠️ **Use `singularity run`, not `singularity exec`.** The image's entrypoint sets up `PATH` and `LD_LIBRARY_PATH` for the conda environment — `run` invokes it, `exec` skips it, which will break AMRFinderPlus.
</details>

<details>
<summary><b>4. SLURM job script (Click to expand)</b></summary>

```bash
#!/bin/bash
#SBATCH --job-name=staphscope
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=staphscope_%j.log

module load containers/singularity/3.6.4

cd $SLURM_SUBMIT_DIR

singularity run --bind $(pwd):/data \
  /path/to/staphscope_2.1.0.sif \
  -i "/data/*.fna" \
  -o /data/output \
  -t ${SLURM_CPUS_PER_TASK}
```
</details>

<details>
<summary><b>5. PBS / Torque job script (Click to expand)</b></summary>

```bash
#!/bin/bash
#PBS -N staphscope
#PBS -l nodes=1:ppn=16
#PBS -l walltime=04:00:00
#PBS -j oe

module load containers/singularity/3.6.4

cd $PBS_O_WORKDIR

singularity run --bind $(pwd):/data \
  /path/to/staphscope_2.1.0.sif \
  -i "/data/*.fna" \
  -o /data/output \
  -t 16
```
</details>

<details>
<summary><b>6. HPC notes (Click to expand)</b></summary>

- **No `module load` needed at runtime** — only the login-node pull requires `squashfs-tools`. Compute nodes just need Singularity.
- **Pull once, run everywhere** — the SIF is a single self-contained file. Share it on your scratch filesystem and every compute node can use it.
- **Small home quota?** Redirect the Singularity build cache to scratch:
  ```bash
  export SINGULARITY_CACHEDIR=/scratch/$USER/singularity_cache
  ```
- **MLST / AMR caches** live under `$HOME/.local/share/staphscope/` by default. On tight home quotas:
  ```bash
  export STAPHSCOPE_MLST_DB=/scratch/$USER/staphscope/mlst_db
  ```
</details>

<details>
<summary><b>7. Troubleshooting (Click to expand)</b></summary>

| Error | Fix |
|---|---|
| `mksquashfs: executable file not found` | `module load squashfs-tools/4.4` before `singularity pull` |
| `manifest unknown` | The Docker tag hasn't been pushed, or the name is wrong — check https://hub.docker.com/r/bbeckleyhub/staphscope/tags |
| `singularity: command not found` | `module avail 2>&1 \| grep -i singularity` to find the module name |
| Pull hangs (blocked registry) | Download the SIF from the GitHub release page instead |
| `BUNDLED AMRfinderPlus not properly installed` | You used `singularity exec` — use `singularity run` instead |
</details>

---

### 🤔 **Which should you use?**

| Platform | Recommendation | Command pattern | Output ownership |
|---|---|---|---|
| **Laptop / workstation** | Docker | `docker run --rm -u $(id -u):$(id -g) -v "$PWD:/data" …` | You |
| **HPC login + compute nodes** | **Singularity (or Apptainer)** | `singularity run --bind "$PWD:/data" staphscope.sif …` | You (automatic) |
| **Cloud VM** | Either — Docker is simpler | see above | You |
| **Shared cluster without Docker** | Singularity | see above | You (automatic) |

**The rule:** Docker on machines where you control the daemon; Singularity on shared HPC systems where you don't.

---

## 🔗 **Integrated External Tools & Dependencies**

<details>
<summary><b>Click to expand the full list of integrated tools and their licenses</b></summary>

| Tool/Database | Purpose | Source | License |
|---------------|---------|--------|---------|
| **MLST** | Multi-locus sequence typing | [tseemann/mlst](https://github.com/tseemann/mlst) | GPL v2 |
| **ABRicate** | Mass screening for resistance/virulence | [tseemann/abricate](https://github.com/tseemann/abricate) | GPL v2 |
| **AMRFinderPlus** | Antimicrobial resistance gene detection | [ncbi/amr](https://github.com/ncbi/amr) | Public Domain |
| **SCCmecFinder** | SCCmec typing (CGE caller) | [genomicepidemiology/Sccmecfinder](https://bitbucket.org/genomicepidemiology/Sccmecfinder) | Apache-2.0 |
| **sccmec (RPet)** | SCCmec typing (RPet caller) | [rpetit3/sccmec](https://github.com/rpetit3/sccmec) | MIT |
| **AgrVATE** | Agr typing | [VishnuRaghuram94/AgrVATE](https://github.com/VishnuRaghuram94/AgrVATE) | MIT |
| **spa typing** | *spa* gene typing | [spa.ridom.de](https://spa.ridom.de/) | Free for academic use |
| **PubMLST** | MLST allele database | [pubmlst.org](https://pubmlst.org/organisms/staphylococcus-aureus) | Open access for research |
| **fastANI** | Species-level ANI confirmation | [ParBLiSS/FastANI](https://github.com/ParBLiSS/FastANI) | Apache-2.0 |
| **mobileOG-db** | Mobile genetic element DB | [clb21565/mobileOG-db](https://github.com/clb21565/mobileOG-db) | MIT |
| **Prodigal** | Prokaryotic gene prediction | [hyattpd/Prodigal](https://github.com/hyattpd/Prodigal) | GPL v3 |
| **DIAMOND** | Fast protein alignment | [bbuchfink/diamond](https://github.com/bbuchfink/diamond) | BSD-3 |
| **Chart.js** | Interactive visualization | [chartjs.org](https://www.chartjs.org/) | MIT |

</details>

---

## 📈 **Performance Benchmarks**

| System | Samples | Time | Speed vs Bactopia |
|--------|---------|------|-------------------|
| Laptop (2 cores, 8GB) | 1 | 2m 33s | 5× faster |
| Laptop (2 cores, 8GB) | 24 | 28m 17s | 6× faster |
| Workstation (16 cores, 16GB) | 1 | 1m 31s | 8× faster |
| Workstation (16 cores, 16GB) | 24 | 14m 34s | 10× faster |
| Workstation (16 cores, 16GB) | 100 | ~60m | 12× faster |

---

## 🔬 **Validation & Accuracy**

### Reference Strain Validation
**100% concordance** with gold-standard reference genomes:

| Reference Strain | Expected Type | StaphScope Result |
|------------------|---------------|-------------------|
| USA300 | ST8–t008–IV(2B) | ✅ ST8–t008–IV(2B) |
| N315 | ST5–t002–II(2A) | ✅ ST5–t002–II(2A) |
| MRSA252 | ST36–t018–II(2A) | ✅ ST36–t018–II(2A) |
| TW20 | ST239–t037–III(3A) | ✅ ST239–t037–III(3A) |
| NCTC8325 | ST8–t211–None | ✅ ST8–t211–Not Assigned |

---

## 🤖 **AI Integration Guide**

StaphScope reports are structured with clear tables and organized data that AI can easily understand.

1. **Install any AI browser extension** (ChatGPT, Claude, Gemini).
2. **Open your report**: `staphscope_dashboard.html` or `staphscope_ultimate_gene_centric_report.html`.
3. **Select text** in any section (AMR Genes, MLST Analysis, etc.).
4. **Right-click → Ask AI** with your question.

**Example Questions:**
- *"These two samples share MLST, spa, and SCCmec. What is the probability they are a transmission pair?"*
- *"Which lineage carries the highest mean plasmid-associated protein count?"*
- *"List all samples with vancomycin resistance genes and their SCCmec types."*

---

## 🔮 **Future Development**

### 🚀 Upcoming Features (2026-2027)
- **Raw read support** – Direct FASTQ analysis with integrated assembly (Shovill).
- **Machine learning module** – Outbreak prediction, phenotype inference, risk scoring.
- **Real‑time database updates** – Live synchronization of lineage and AMR databases.
- **Plugin system** – Community‑contributed analysis modules.
- **Expanded ESKAPE coverage** – Porting StaphScope's architecture to other ESKAPE pathogens.

---

## ❓ **Frequently Asked Questions**

**Q: Is StaphScope free to use?**  
A: Yes! Open-source under MIT License. Free for academic, clinical, and commercial use.

**Q: What makes StaphScope different from other tools?**  
A: *S. aureus*-optimized, integrates 10 analysis types including dual SCCmec callers, and runs 8–10× faster than generalist platforms.

**Q: Can I use StaphScope for clinical diagnosis?**  
A: StaphScope is a research tool. While highly accurate, results should be validated with orthogonal methods for clinical decision-making.

**Q: Which version should I use – CLI or Web?**  
A: Use the **Web version** for convenience, small batches, and graphical interface. Use the **CLI** for large batches, pipeline integration, or when working with sensitive data locally.

---

## 🐛 **Troubleshooting**

### Common Issues & Solutions

```bash
# Issue: AMR database missing or outdated
staphscope --force-update-amr-db

# Issue: ABRicate database not found
abricate --setupdb

# Issue: MLST database missing
staphscope --pull-mlst-db

# Issue: Permission errors in Docker
docker run --rm -u $(id -u):$(id -g) -v ... bbeckleyhub/staphscope ...

# Issue: Visualization dashboard fails with "master TSV not found"
# Solution: Ensure the pipeline ran to completion (Staphscope_final_report/ exists)
```

### Getting Help
1. **Check existing issues**: [GitHub Issues](https://github.com/bbeckley-hub/staphscope-typing-tool/issues)
2. **Create new issue** with: full error, `conda list` output, the failing command, and `staphscope_run.log`
3. **Email support**: brownbeckley94@gmail.com

---

## 📚 **Citation**

If you use StaphScope or SAURUS in your research, please cite:

> Beckley, B., Amarh, V. (2026). StaphScope: a species-optimized computational pipeline for rapid and accessible *Staphylococcus aureus* genotyping and surveillance. *BMC Genomics*, 27:261.

**DOI**: [10.1186/s12864-026-12609-x](https://doi.org/10.1186/s12864-026-12609-x)

```bibtex
@article{beckley2026staphscope,
  title={StaphScope: a species-optimized computational pipeline for rapid and accessible Staphylococcus aureus genotyping and surveillance},
  author={Beckley, Brown and Amarh, Vincent},
  journal={BMC Genomics},
  volume={27},
  pages={261},
  year={2026},
  doi={10.1186/s12864-026-12609-x}
}
```

---

## 🙏 **Acknowledgements**

StaphScope stands on the shoulders of giants. We are deeply grateful to:

- **Torsten Seemann** for MLST, ABRicate, and countless foundational tools.
- **NCBI team** for AMRFinderPlus.
- **CGE team** for SCCmecFinder and database curation.
- **Robert A. Petit III** for the `sccmec` RPet caller.
- **Vishnu Raghuram & Robert A. Petit III** for AgrVATE.
- **PubMLST, Ridom, CARD, VFDB, BacMet, mobileOG-db** for essential databases.
- **ParBLiSS** team for fastANI.
- **The Shelenkov et al. team** for openly publishing their gene dataset.
- **The Chart.js project** for interactive visualization support.
- **[Alyssa-Kent](https://github.com/Alyssa-Kent)** for catching the *mecC* misclassification in v1.3.2.
- **Peer reviewers & Editorial Team @ BMC Genomics** for their constructive feedback.

> *"If we ever meet in person, the drinks are on me!" – Brown Beckley*

---

## 👥 **Authors & Contact**

**Brown Beckley** (Primary Developer)
- University of Ghana Medical School
- 📧 brownbeckley94@gmail.com
- 🐙 GitHub: [bbeckley-hub](https://github.com/bbeckley-hub)
- LinkedIn: [@brownbeckley](https://www.linkedin.com/in/brown-beckley-190315319/)

**Amarh Vincent** (Co-Author)
- University of Ghana Medical School

---

## 📄 **License**

### Core StaphScope Code
The StaphScope pipeline code (workflow engine, report generation, HTML templates, and Python modules) is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.

### Third-Party Tools
StaphScope executes several external bioinformatics tools installed as Conda dependencies. Each tool is the property of its respective developers and is used under its own license.

---

<details>
<summary><b>📚 Third-Party Tool Citations (Click to expand)</b></summary>

#### **AgrVATE**
```bibtex
@article{raghuram_agrv_2022,
  author = {Raghuram, V. and Alexander, A. M. and Loo, H. Q. and Petit, R. A. 3rd and Goldberg, J. B. and Read, T. D.},
  title = {Species-Wide Phylogenomics of the Staphylococcus aureus Agr Operon Revealed Convergent Evolution of Frameshift Mutations},
  journal = {Microbiology Spectrum},
  volume = {10},
  number = {1},
  pages = {e0133421},
  year = {2022},
  doi = {10.1128/spectrum.01334-21}
}
```

#### **sccmec (RPet)**
```bibtex
@article{petit_staphopia_2018,
  author = {Petit, R. A. 3rd and Read, T. D.},
  title = {Staphylococcus aureus viewed from the perspective of 40,000+ genomes},
  journal = {PeerJ},
  volume = {6},
  pages = {e5261},
  year = {2018},
  doi = {10.7717/peerj.5261}
}
```

#### **fastANI**
```bibtex
@article{jain_fastani_2018,
  author = {Jain, C. and Rodriguez-R, L. M. and Phillippy, A. M. and Konstantinidis, K. T. and Aluru, S.},
  title = {High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries},
  journal = {Nature Communications},
  volume = {9},
  number = {1},
  pages = {5114},
  year = {2018},
  doi = {10.1038/s41467-018-07641-9}
}
```

#### **mobileOG-db**
```bibtex
@article{brown_mobileogdb_2022,
  author = {Brown, C. L. and Mullet, J. and Hindi, F. and Stoll, J. E. and Gupta, S. and Choi, M. and Keenum, I. and Vikesland, P. and Pruden, A. and Zhang, L.},
  title = {mobileOG-db: a Manually Curated Database of Protein Families Mediating the Life Cycle of Bacterial Mobile Genetic Elements},
  journal = {Applied and Environmental Microbiology},
  volume = {88},
  number = {18},
  pages = {e00991-22},
  year = {2022},
  doi = {10.1128/aem.00991-22}
}
```

#### **MLST**
```bibtex
@software{seemann_mlst_2018,
  author = {Seemann, T.},
  title = {MLST: Scan contig files against traditional PubMLST typing schemes},
  year = {2018},
  publisher = {GitHub},
  url = {https://github.com/tseemann/mlst}
}
```

#### **PubMLST**
```bibtex
@article{jolley_pubmlst_2018,
  author = {Jolley, K. A. and Bray, J. E. and Maiden, M. C. J.},
  title = {Open-access bacterial population genomics: {BIGSdb} software, the {PubMLST.org} website and their applications},
  journal = {Wellcome Open Research},
  volume = {3},
  pages = {124},
  year = {2018},
  doi = {10.12688/wellcomeopenres.14826.1}
}
```

#### **ABRicate**
```bibtex
@software{seemann_abricate_2018,
  author = {Seemann, T.},
  title = {ABRicate: Mass screening of contigs for antimicrobial resistance and virulence genes},
  year = {2018},
  publisher = {GitHub},
  url = {https://github.com/tseemann/abricate}
}
```

#### **AMRFinderPlus**
```bibtex
@article{feldgarden_amrfinderplus_2021,
  author = {Feldgarden, M. et al.},
  title = {AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence},
  journal = {Scientific Reports},
  volume = {11},
  pages = {12728},
  year = {2021},
  doi = {10.1038/s41598-021-91456-0}
}
```

#### **SCCmecFinder (CGE)**
```bibtex
@article{kaya_sccmecfinder_2018,
  author = {Kaya, H. et al.},
  title = {SCCmecFinder, a Web-Based Tool for Typing of Staphylococcal Cassette Chromosome mec in Staphylococcus aureus Using Whole-Genome Sequence Data},
  journal = {mSphere},
  volume = {3},
  number = {1},
  pages = {e00612-17},
  year = {2018},
  doi = {10.1128/mSphere.00612-17}
}
```

#### ***spa* Typing (Ridom)**
```bibtex
@article{mellmann_spa_typing_2005,
  author = {Mellmann, A. et al.},
  title = {Evidenzbasierte Hygienemassnahmen mittels spa-Typisierung bei MRSA-Häufungen im Krankenhaus},
  journal = {Deutsche Medizinische Wochenschrift},
  volume = {130},
  number = {22},
  pages = {1364-1368},
  year = {2005},
  doi = {10.1055/s-2005-868351}
}
```

#### **Biopython**
```bibtex
@article{cock_biopython_2009,
  author = {Cock, P. J. A. et al.},
  title = {Biopython: freely available Python tools for computational molecular biology and bioinformatics},
  journal = {Bioinformatics},
  volume = {25},
  number = {11},
  pages = {1422-1423},
  year = {2009},
  doi = {10.1093/bioinformatics/btp163}
}
```

</details>

<details>
<summary><b>🤖 AI Disclosure (Click to expand)</b></summary>

StaphScope is developed with the assistance of AI tools, and we believe in being transparent about it.

**How we use AI:**
- **Code drafting and refactoring.** AI assistants help generate initial implementations, suggest refactors, and catch patterns we might miss.
- **Documentation.** Some prose in this README, module docstrings, and the HTML report text has been drafted with AI assistance and subsequently edited by human authors.
- **Testing scaffolds.** AI helps write boilerplate tests and mock fixtures.

**What we do NOT do:**
- **We do not ship code we don't understand.** Every line that lands in the repository has been reviewed by at least one human author.
- **We do not ship code we can't debug.** If a bug appears in code you can't reason about, that's a bug in the process, not a feature.
- **We do not use AI for scientific claims.** Biological interpretations, thresholds, and validation calls are made by the human authors based on peer-reviewed literature and reference strains.
- **We do not train models on your data.** StaphScope is a local tool. Nothing runs through a third-party service. Your genomes stay on your machine.

**Responsibility:**
The authors of StaphScope take **full responsibility** for the correctness, behaviour, and safety of the released code. AI is a productivity tool in our hands, not an author on our paper. If something in StaphScope is wrong, the fault is ours — not the tools we used to write it.

</details>

---

<div align="center">

## **🚀 Ready to revolutionize your MRSA analysis?**

| **Choose Your Platform** | |
|--------------------------|-|
| 🖥️ **Command Line** | For high-throughput, local analysis |
| 🌐 **StaphScope Web** | For non-bioinformaticians – [https://eskape.bio](https://eskape.bio) |
| 🧬 **SAURUS Database** | Browse, filter, and download genomes – [https://bbeckley-hub.github.io/saurus/](https://bbeckley-hub.github.io/saurus/) |

[![Get Started CLI](https://img.shields.io/badge/GET_STARTED_CLI-Now-green?style=for-the-badge&logo=github)](https://github.com/bbeckley-hub/staphscope-typing-tool#-quick-start-cli)
[![Try Web Version](https://img.shields.io/badge/TRY_WEB_VERSION-Here-blue?style=for-the-badge&logo=html5)](https://eskape.bio)
[![Report Issue](https://img.shields.io/badge/REPORT_ISSUE-Here-red?style=for-the-badge&logo=github)](https://github.com/bbeckley-hub/staphscope-typing-tool/issues)

**From days to minutes. From fragmented to integrated. From data to insights.**

*StaphScope: Precision surveillance for the antibiotic resistance era.*

⭐ **If you find this tool useful, please star the repository!** ⭐

*Join the Fight Against Antimicrobial Resistance*

Antimicrobial resistance (AMR) represents one of the most significant global health threats of our time. We invite researchers, clinicians, and public health professionals to collaborate with us in expanding and validating our database, sharing regional epidemiological data, and advancing AMR surveillance.

**Together, we can enhance global AMR monitoring and develop more effective treatment strategies.**

</div>
