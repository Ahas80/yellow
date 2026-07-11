# Documentation for Lecturer: rUTIs Pipeline Overview

**Author:** Aamir  
**Project:** Yellow RoUTIne (rUTIs) - Recurrent UTI Genomic Analysis  
**Date:** December 2025  

---

## 🎯 Quick Start (30-Second Summary)

**What:** Longitudinal genomic study of recurrent UTIs in nursing home residents  
**Question:** Do the same bacteria cause both ASB (no symptoms) and UTI (symptoms)?  
**Answer:** YES! We found phenotype switching driven by host factors, not bacterial evolution  

**Critical Finding:** 2 participants had the **same bacterial strain** (≤12 SNPs) cause ASB at one visit and UTI at another—**without acquiring new genes**. This is rare but scientifically important.

**Data:** 274 episodes, ~380 E. coli genomes, 100 participants, 6-year longitudinal follow-up

---

## 📚 Essential Reading (In Priority Order)

### 1️⃣ **Script & Figure Guide** - START HERE! 📊
**File:** `script_and_figure_guide.md` (201 lines, ~15 min read)

Complete overview of all 40+ scripts organized by analysis phase:
- Phase 0: Clinical data (00a-00c)
- Phase 1: WGS typing (02-13)
- Phase 2: GWAS & lineage (14-17)
- Phase 3: Longitudinal evolution (15-21)

**What you'll learn:** Purpose, methodology, inputs/outputs for every script

---

### 2️⃣ **Figure Explanations** - FOR PRESENTATIONS! 🎨
**File:** `figures/figure_legends_master.md` (120 lines, ~10 min read)

Clear "What/Why/How/Takeaway" explanations for 12 key figures:
- Clinical: Heatmaps, transitions, waterfall plots
- Genomics: VF prevalence, plasmid heatmaps, phylogeny trees
- Longitudinal: Swimmer plots showing patient timelines
- GWAS: Volcano/forest plots of gene associations

**Includes figure manifest table at end**

---

### 3️⃣ **Color Coding Reference** 🌈
**File:** `figures/colour_scheme.md` (30 lines, ~2 min read)

**CRITICAL for interpretation:**
- **UTI** = `#D55E00` (Vermilion/Orange)
- **ASB** = `#0072B2` (Blue)
- **Negative** = `#909090` (Grey)

All infection status plots use these canonical colors consistently.

---

### 4️⃣ **Terminology Clarification** ⚠️
**File:** `figures/timepoint_vs_isolate_clarification.md` (40 lines, ~3 min read)

**IMPORTANT:** Explains why Panaroo plots show 6-10 "isolates" not "timepoints"
- Each biological sample has 2 assemblers (Flye + Unicycler)
- Panaroo receives 2 assembly files per sample
- Example: 4 timepoints × 2 assemblers = 8 isolates

---

### 5️⃣ **Research Findings** 📈
**File:** `research_outcomes.md` (107 lines, ~10 min read)

Key discoveries:
- **Phenotype switching:** Same strain causes ASB then UTI
- **Genomic stasis:** No virulence gene acquisition during switches
- **Host-driven:** Symptoms driven by host factors (immune, catheter)
- **No lineage effect:** ST131 not more "symptomatic" in this cohort

Includes detailed case studies of 2 switch events with specific mutations.

---

## 🔧 Technical Documentation (Optional)

### Analysis Outputs
**File:** `analysis_outputs_guide.md` (220 lines)  
Detailed explanation of output file formats and result interpretation

### Pipeline Architecture
**File:** `pipeline_architecture.md` (185 lines)  
Technical workflow, computational requirements, dependencies

---

## 📋 Complete File Index

**Total: 8 essential documents**

| Priority | File | Lines | Topic |
|:---------|:-----|------:|:------|
| ⭐⭐⭐ | script_and_figure_guide.md | 201 | All scripts explained |
| ⭐⭐⭐ | figures/figure_legends_master.md | 120 | Figure explanations + manifest |
| ⭐⭐⭐ | figures/colour_scheme.md | 30 | Color reference |
| ⭐⭐⭐ | figures/timepoint_vs_isolate_clarification.md | 40 | Terminology |
| ⭐⭐ | research_outcomes.md | 107 | Key findings |
| ⭐ | analysis_outputs_guide.md | 220 | Output formats |
| ⭐ | pipeline_architecture.md | 185 | Technical details |
| 📝 | methods_draft.md | 45 | Manuscript methods |

**Estimated reading time:**
- Essential (★★★): 1 hour
- Core (★★): +40 min
- Complete: +1 hour

---

## 🔬 Research Context

### The Scientific Question
*Can the same E. coli strain cause both asymptomatic bacteriuria (ASB) and symptomatic UTI in the same patient?*

### Approach
1. Longitudinal sampling (T0-T4 + Uricult symptom events)
2. Whole-genome sequencing (WGS)
3. SNP-based strain tracking (0-25 SNPs = strong same strain)
4. Gene content comparison (VF, plasmids, AMR)
5. Clinical phenotype association (GWAS)

### Key Finding
Under the current primary UTI-vs-Not_UTI definition, switch events should be read from the freshly generated `results/longitudinal/phenotype_switch_candidates.csv`.

Candidate switch events are those where:
- Strong same-strain evidence under the current 0-25 SNP threshold
- Different primary clinical status (`Not_UTI` → `UTI` or reverse)
- Gene-content and SNP evidence can then be inspected descriptively
- ∴ Interpret host-versus-bacterial mechanisms cautiously because the UTI denominator is small

---

## 📊 Data Overview

| Metric | Count |
|:-------|------:|
| **Participants** | ~100 |
| **Clinical episodes** | 274 |
| **E. coli genomes** | 382 assemblies |
| **UTI episodes** | 42 |
| **ASB episodes** | 282 |
| **Negative episodes** | 66 |
| **Phenotype switches** | 6 confirmed |
| **Assemblers** | 2 per sample (Flye + Unicycler) |

---

## 🎨 What Has Been Standardized

✅ **Color scheme:** All infection plots use Orange/Blue/Grey consistently  
✅ **Terminology:** "Isolates" vs "Timepoints" clarified  
✅ **Documentation:** All scripts have purpose/methodology/outputs  
✅ **Figure legends:** 12 key figures explained for lecturers  
✅ **Code quality:** GLMM models, QC thresholds, statistical rigor applied

---

## 💡 Discussion Points for Your Lecturer

1. **Clinical Impact:** Should nursing homes treat ASB differently knowing it can switch to UTI without gene changes?
2. **Statistical Power:** N=6 switches - sufficient for strong claims?
3. **Causality:** Can we infer host-driven symptoms from SNP distance alone?
4. **Publication Strategy:** Which analyses are manuscript-ready vs. supplementary?
5. **Future Work:** Transcriptomics to validate "same genotype, different phenotype"?

---

## 📞 Project Details

**Student:** Aamir  
**Directory:** `/Users/Aamir/Desktop/rUTIs/`  
**Documentation:** `/Users/Aamir/Desktop/rUTIs/docs/`  
**Updated:** December 2, 2025

---

## ✅ Files You Can Delete (Merged Content)

~~`DOCUMENTATION_INDEX.md`~~ (merged into this file)  
~~`figure_manifest.md`~~ (merged into figure_legends_master.md)  
~~`figure_audit_summary.md` + `script_audit_infection_colors.md`~~ (merged, see note below)

**Note:** Technical audit details about color standardization are preserved in `script_and_figure_guide.md` and this README.

---

## Research Context

### What is this project about?

**Scientific Question:**  
*Do the same bacterial strains cause both asymptomatic bacteriuria (ASB) and symptomatic urinary tract infections (UTI) in nursing home residents?*

**Approach:**
1. Longitudinal clinical sampling (T0, T1, T2, T3, T4, Uricult)
2. Whole-genome sequencing of E. coli isolates
3. Comparative genomics to track strain persistence
4. Genotype-phenotype association (GWAS)
5. Within-host evolution analysis

**Key Finding:**  
Yes! We identified "phenotype switch events" where the **same bacterial strain** (0-25 SNPs under the current same-strain rule) caused ASB at one timepoint and UTI at another in the same patient.

**Implication:**  
Symptoms may be driven by **host factors** (catheter status, immune state) rather than bacterial virulence alone.

---

## Key Terminology

| Term | Definition |
|:-----|:-----------|
| **ASB** | Asymptomatic Bacteriuria (bacteria in urine, no symptoms) |
| **UTI** | Urinary Tract Infection (bacteria + symptoms: dysuria, urgency, fever) |
| **ST** | Sequence Type (MLST-based lineage classification, e.g., ST131) |
| **VF** | Virulence Factor (bacterial genes that enable pathogenesis) |
| **Panaroo** | Pangenome tool that identifies core vs. accessory genes |
| **GWAS** | Genome-Wide Association Study (linking genes to phenotypes) |
| **Swimmer Plot** | Timeline visualization showing patient infection history |

---

## Data Overview

- **Cohort:** ~274 episodes from nursing home residents
- **Participants:** ~100 individuals (some with 6+ timepoints)
- **Isolates:** ~380 E. coli genome assemblies
- **Phenotypes:**
  - UTI: 42 episodes
  - ASB: 282 episodes  
  - Negative: 66 episodes

---

## Quick Reference: Where to Find Things

| What you need | Where to look |
|:--------------|:--------------|
| **Script descriptions** | `docs/script_and_figure_guide.md` |
| **Figure explanations** | `docs/figures/figure_legends_master.md` |
| **Color meanings** | `docs/figures/colour_scheme.md` |
| **Complete figure list** | `docs/figures/figure_manifest.md` |
| **Results interpretation** | `docs/analysis_outputs_guide.md` |
| **Research findings** | `docs/research_outcomes.md` |

---

## Questions for Discussion

When reviewing this work, consider:

1. **Clinical Relevance:** How do these findings change ASB management in nursing homes?
2. **Statistical Rigor:** Are the GLMM models appropriate for this nested data structure?
3. **Sample Size:** Is N=6 phenotype switches sufficient for strong claims?
4. **Causality:** Can we infer host-driven symptoms from SNP distances alone?
5. **Publication Strategy:** Which analyses are ready for manuscript vs. supplementary?

---

## Contact & Support

**Student:** Aamir  
**Project Repository:** `/Users/Aamir/Desktop/rUTIs/`  
**Documentation:** `/Users/Aamir/Desktop/rUTIs/docs/`  

For questions about specific scripts, see the **Script & Figure Guide** first!
