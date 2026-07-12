# 🔍 Using the Grouping Feature

One of the most powerful features in the StaphScope reports is **dynamic grouping by typing**.

---

## Where to Find It

In the **gene‑centric report** (`staphscope_ultimate_report.html`), go to any of these tabs:
- AMR
- Virulence
- BACMET
- Plasmids
- Mutations

Above each table, you'll see a bar labeled:

> **Group genomes by:**

---

## Available Grouping Options

| Grouping Option | Example |
|-----------------|---------|
| MLST | ST5, ST8, ST22, etc. |
| spa | t008, t011, t1234, etc. |
| SCCmec | I, II, III, IV, etc. |
| agr | I, II, III, IV |
| ST‑spa | ST5 - t011 |
| ST‑SCCmec | ST8 - IV |
| spa‑SCCmec | t008 - IV |
| ST‑agr | ST5 - II |
| spa‑agr | t011 - II |
| SCCmec‑agr | IV - II |
| ST‑spa‑agr | ST5 - t011 - II |
| ST‑SCCmec‑agr | ST5 - IV - II |
| spa‑SCCmec‑agr | t011 - IV - II |
| Triple (ST‑spa‑SCCmec) | ST5 - t011 - IV |
| Four‑way (ST‑spa‑SCCmec‑agr) | ST5 - t011 - IV - II |

---

## How to Use It

1. Click any button in the grouping bar.
2. The genome list in that table will instantly reorganise – genomes are grouped under sub‑headers for each typing value.
3. You can still use the search and highlight boxes to find specific genomes across all groups.
4. Click **"Reset"** to return to the flat, scrollable list.

---

## Why This Is Powerful

Instead of manually scanning a long list of genome names, you can instantly see:

- **"Does *mecA* appear only in ST5 and ST8?"**
- **"Does PVL appear only in spa type t008?"**
- **"Which agr types carry the most resistance genes?"**
- **"Is there a clone (ST‑spa‑SCCmec‑agr) that has both AMR and virulence genes?"**

---

## Example Questions

- *"Group AMR genes by MLST – which STs carry the most resistance genes?"*
- *"Group virulence genes by agr – which agr types carry PVL?"*
- *"Group BACMET genes by SCCmec – is biocide resistance linked to MRSA cassettes?"*
- *"Group mutations by four‑way typing – which clones have linezolid‑related mutations?"*

---

## Pro Tips

- Combine grouping with the **search boxes** to zoom in on specific genes or samples.
- Use the **highlight** box to find all occurrences of a specific genome across groups.
- The grouping works on the **fly** – no need to regenerate the report.
