# VF Abstract — Independent Verification Report

> Full audit produced by independent recomputation from anchor files. No project R scripts sourced — only raw CSV reads.
> Verification script: `/tmp/verify_vf_abstract.R`

---

## 1) Executive Verdict

### **PASS — with 3 minor issues requiring wording changes**

No major numerical errors were found. Every statistic in the abstract was independently recomputed from the anchor files (`vf_pa_all.csv`, `status_map.csv`, `gene_map.csv`) and all values match exactly or to rounding precision.

#### Top Issues (ranked by risk)

| # | Issue | Severity | Impact |
|---|-------|----------|--------|
| 1 | **Species "confirmed" is misleading** | Minor | The Organism column comes from the clinical microbiology lab database, not bioinformatic tools (Kraken/Mash/Prokka). `00_make_assembly_metadata.r` does not assign species. Recommend: "designated as *E. coli* in clinical microbiology records" |
| 2 | **"Negative" status label needs clarification** | Minor | 31 episodes labelled "Negative" still have VF data (min 43 VFDB genes). These are isolates that were sequenced but clinically classified as Negative. Abstract should clarify what "culture-negative" means in this context. |
| 3 | **Stratified abstract (≥2tp) has 13 UTI, not 16** | Minor | The all-participants abstract correctly reports 16 UTI. The ≥2tp stratification drops 3 UTI participants who have only 1 timepoint. This is correct behaviour but should be explicitly stated. |
| 4 | **status_map duplicate keys** | Cosmetic | 2 duplicate rows in status_map.csv after tp normalization; handled by taking first row per key — no impact on results |
| 5 | **Repeated-measures limitation acknowledged** | OK | Fisher tests correctly labelled as "exploratory/unadjusted" ✅ |
| 6 | **BH-adjusted p = 1.0** | OK | Mathematically correct — all nominal p-values are too large relative to the number of tests ✅ |
| 7 | **Uricult ordering** | Potential risk | Uricult placed AFTER T4 in ordering. If a Uricult event occurred between T0 and T1, the transition would be misordered. No collection dates are used. |
| 8 | **iro delta rounding** | Cosmetic | Abstract says 15.0pp, verified as 15.1pp (56.25% − 41.18% = 15.07pp). Acceptable rounding. |

---

## 2) Project File Map and Scripts Used

### Anchor Files (primary data sources)

| File | Path | Rows × Cols | Role |
|------|------|-------------|------|
| VF P/A matrix | `results/vf/vf_pa_all.csv` | 183 × 166 (164 genes) | Binary VF presence per Participant×Timepoint |
| Clinical status | `status_map.csv` | 276 × 13 | Infection_Status per Participant×Timepoint |
| Gene categories | `results/vf/gene_map.csv` | 209 × 3 | Gene → Category mapping |
| Assembly metadata | `assembly_metadata.csv` | 382 × 21 | Organism column (species claim) |

### Scripts in Pipeline

| Script | Role | Key Lines |
|--------|------|-----------|
| `02_gene_presence_analysis.R` | Runs Abricate, builds vf_pa_all.csv | L45-50 (thresholds), L134 (VFDB), L185-188 (pivot) |
| `04_gene_breakdown.R` | Annotates genes, creates gene_map.csv | L~50-80 (regex categories) |
| `00_make_assembly_metadata.r` | Creates assembly_metadata.csv | Does **NOT** assign species |
| `compute_vf_abstract_stats.R` | Produces abstract CSVs | All deliverables B-D |
| `compute_vf_stratified_by_depth.R` | Stratified analysis | ≥2/≥3/≥4tp cohorts |

---

## 3) Traceability Table

| # | Claim | Verified? | Recomputed | Source | Notes |
|---|-------|-----------|------------|--------|-------|
| 1 | "183 *E. coli* episodes from 87 participants" | ✅ Yes | 183 rows, 87 participants | `vf_pa_all.csv` direct count | Exact match |
| 2 | "All 382 assemblies confirmed as E. coli" | ⚠️ Partially | 382/382 say "Escherichia coli" | `assembly_metadata.csv` col 13 | "Confirmed" overstates; it's a lab designation, not bioinformatic verification |
| 3 | "136 ASB, 16 UTI, 31 Negative" | ✅ Yes | 136/16/31 | `vf_pa_all.csv` ⋈ `status_map.csv` | Exact match |
| 4 | "75/16/21 participants" | ✅ Yes | 75/16/21 | Same join, `n_distinct(Participant_id)` | Exact match |
| 5 | "ASB median 80 (IQR 70–92)" | ✅ Yes | 80.0 (69.75–92.0) | `vf_analysis_ready.csv` | IQR boundaries are 69.75 and 92.0 — abstract rounds to 70 and 92 ✅ |
| 6 | "UTI median 80.5 (IQR 70–98)" | ✅ Yes | 80.5 (70.0–97.5) | Same | IQR upper is 97.5, rounded to 98 in abstract ⚠️ (should be 97.5 or 98) |
| 7 | "Negative median 88 (IQR 77–100)" | ✅ Yes | 88.0 (77.0–99.5) | Same | IQR upper is 99.5, rounded to 100 ⚠️ |
| 8 | "astA 25.0% UTI vs 6.6% ASB; OR=4.63, p=0.033" | ✅ Yes | 4/16=25.0%, 9/136=6.6%; OR=4.628, p=0.0331 | Fisher's exact on 2×2 table | Exact match to 3 significant figures |
| 9 | "pic 0.0% UTI vs 18.4% ASB; p=0.075" | ✅ Yes | 0/16=0.0%, 25/136=18.4%; p=0.0752 | Fisher exact | Exact match |
| 10 | "iro operon 56.2% UTI vs 41.2% ASB" | ✅ Yes | 9/16=56.25%, 56/136=41.18% | Per-gene count | Δ=15.07pp (abstract says 15.0pp — acceptable rounding) |
| 11 | "No genes significant after BH" | ✅ Yes | All adjusted p ≥ 1.0 | `p.adjust(method="BH")` | Correct |
| 12 | "96 consecutive transitions from 83 participants" | ✅ Yes | 96 transitions, 83 participants | Consecutive-pair loop, `tp_order` | Exact match |
| 13 | "72.9% zero VF change" | ✅ Yes | 70/96 = 72.917% | `mean(!any_change)` | Exact match |
| 14 | "Median Jaccard = 1.000" | ✅ Yes | 1.000 | `median(jaccard)` | Exact match |
| 15 | "12 ASB→UTI transitions" | ✅ Yes | 12 | `filter(tt == "ASB→UTI")` | Exact match |
| 16 | "iroB–N gained 4/12 in ASB→UTI" | ✅ Yes | iroB/C/D/E/N each = 4 | Parsed `genes_gained` | Exact match |
| 17 | "kpsD/M gained 3/12" | ✅ Yes | kpsD=3, kpsM=3 | Same | Exact match |
| 18 | "≥2tp: 83 participants, 179 episodes" | ✅ Yes | 83 participants, 179 episodes | Stratified filter | Exact match |
| 19 | "≥3tp: 8 participants, 29 episodes" | ✅ Yes | 8/29 | Same | Exact match |
| 20 | "≥4tp: 5 participants, 20 episodes" | ✅ Yes | 5/20 | Same | Exact match |
| 21 | "iro gap: +20pp / +26pp / +35pp" | ✅ Yes | +20.1 / +26.2 / +34.6 | Per-stratum prevalence | Abstract rounds; verified values are consistent |
| 22 | "All 16 UTI at Uricult" | ✅ Yes | Cross-tab: UTI×Uricult=16, all others=0 | `table(status, tp_lab)` | Exact match |
| 23 | "164 VF genes screened" | ✅ Yes | 164 gene columns | `ncol(vf_pa) - 2` | Exact match |
| 24 | "Abricate VFDB ≥80% identity, ≥80% coverage" | ✅ Yes | Lines 45-50, 134, 148 | `02_gene_presence_analysis.R` | Exact match |

---

## 4) Recomputed Results (Tables)

### Table A: Cohort Counts

| | All | ≥2tp | ≥3tp | ≥4tp |
|---|---|---|---|---|
| Participants | 87 | 83 | 8 | 5 |
| Episodes | 183 | 179 | 29 | 20 |
| ASB | 136 | 135 | 18 | 13 |
| UTI | 16 | 13 | 7 | 4 |
| Negative | 31 | 31 | 4 | 3 |

**Status × Timepoint (all participants):**

| | T0 | T1 | T2 | Uricult |
|---|---|---|---|---|
| ASB | 64 | 64 | 8 | 0 |
| Negative | 16 | 13 | 1 | 1 |
| UTI | 0 | 0 | 0 | 16 |

### Table B: VF Burden by Status (all participants, no filter)

| Status | n | Median (IQR) | Mean ± SD | Range |
|--------|---|-------------|-----------|-------|
| ASB | 136 | 80.0 (69.75–92.0) | 79.3 ± 16.5 | 32–109 |
| UTI | 16 | 80.5 (70.0–97.5) | 81.3 ± 19.6 | 32–112 |
| Negative | 31 | 88.0 (77.0–99.5) | 85.7 ± 20.1 | 43–121 |

### Table C: Top Genes by |Δ(UTI−ASB)|

| Gene | ASB (n/136) | UTI (n/16) | Δ pp | Fisher OR | p | p_adj |
|------|-------------|-----------|------|-----------|---|-------|
| astA | 9 (6.6%) | 4 (25.0%) | +18.4 | 4.628 | 0.033 | 1.0 |
| pic | 25 (18.4%) | 0 (0.0%) | −18.4 | 0.000 | 0.075 | 1.0 |
| iroB–N | 56 (41.2%) | 9 (56.2%) | +15.1 | 1.829 | 0.291 | 1.0 |
| papE | 28 (20.6%) | 1 (6.2%) | −14.3 | 0.259 | 0.310 | 1.0 |
| vat | 76 (55.9%) | 11 (68.8%) | +12.9 | 1.731 | 0.426 | 1.0 |

### Table D: Longitudinal Transitions (all participants with ≥2tp)

| Transition | n | Median Jaccard | % Zero Change | Mean Gained | Mean Lost |
|-----------|---|---------------|---------------|-------------|----------|
| ASB→ASB | 61 | 1.000 | 82.0% | 2.3 | 0.9 |
| ASB→UTI | 12 | 1.000 | 58.3% | 9.6 | 3.2 |
| Neg→Neg | 9 | 0.989 | 44.4% | 12.8 | 3.8 |
| Neg→ASB | 7 | 1.000 | 57.1% | 9.0 | 11.9 |
| ASB→Neg | 6 | 1.000 | 66.7% | 13.0 | 1.3 |
| Neg→UTI | 1 | 1.000 | 100.0% | 0.0 | 0.0 |
| **Total** | **96** | **1.000** | **72.9%** | | |

### Table E: Top Gained Genes in ASB→UTI (n=12 transitions)

| Gene | Times Gained | Category |
|------|-------------|----------|
| iroB, iroC, iroD, iroE, iroN | 4 each | Iron acquisition |
| kpsD, kpsM | 3 each | Capsule/Surface |
| chuA, chuT, chuX | 2 each | Iron acquisition |

---

## 5) Discrepancy Report

| # | Claim in Abstract | Verified Value | Discrepancy? | Action |
|---|------------------|---------------|-------------|--------|
| 1 | "confirmed as E. coli" | Organism column = "Escherichia coli" in all 382 rows, but source is clinical lab not bioinformatic tool | **Yes — wording issue** | Change to "designated as *E. coli* in clinical microbiology records" |
| 2 | "IQR 70–92" (ASB burden) | 69.75–92.0 | **Minor rounding** | Change to "IQR 70–92" (acceptable) OR "IQR 69.8–92.0" for precision |
| 3 | "IQR 70–98" (UTI burden) | 70.0–97.5 | **Minor** | Change to "IQR 70–97.5" or "IQR 70–98" |
| 4 | "IQR 77–100" (Neg burden) | 77.0–99.5 | **Minor** | Change to "IQR 77–99.5" or "IQR 77–100" |
| 5 | "Δ=+15.0 pp" (iro operon) | 56.25% − 41.18% = 15.07pp | **Minor rounding** | Change to "+15.1 pp" or keep "+15.0 pp" |
| 6 | ≥2tp stratified table says "179 episodes" but status breakdown sums to 135+31+13=179 | 179 ✅ | **No discrepancy** | — |

**All numerical claims match to within acceptable rounding precision.** No major discrepancies found.

---

## 6) Corrected Abstract

### Virulence Factor Profiles of Urinary *Escherichia coli* Are Largely Conserved Across Clinical States and Within Individuals in a Prospective Cohort of Elderly Residents

**Background**
Asymptomatic bacteriuria (ASB) is highly prevalent among elderly residents of long-term care facilities, yet most episodes do not progress to symptomatic urinary tract infection (UTI). Whether the virulence factor (VF) repertoire of colonising *Escherichia coli* strains differs between ASB and UTI — and whether VF content changes within individuals over time — remains poorly characterised in this population. We aimed to systematically compare VF profiles between clinical states and quantify within-host VF dynamics across serial timepoints.

**Methods**
*E. coli* isolates from the YELLOW RoUTIne rUTI cohort, a prospective surveillance study of elderly residents, were subjected to whole-genome sequencing (Oxford Nanopore Technologies) and assembled using Flye and/or LongCycler. Virulence factor genes were detected using Abricate against the VFDB database (minimum 80% nucleotide identity, minimum 80% query coverage). If a gene was detected by either assembler for a given participant-timepoint, it was marked as present. Clinical status (ASB, UTI, or culture-negative) was assigned based on microbiological and clinical criteria per the study protocol; all isolates were designated as *E. coli* in the clinical microbiology records.

A total of 183 participant×timepoint observations from 87 participants were analysed (136 ASB, 16 UTI, 31 negative). VF burden was defined as the count of distinct VFDB genes detected per sample (164 genes screened). Between-group prevalence comparisons used Fisher's exact test (exploratory and unadjusted; the same participant can contribute multiple episodes, violating independence). Within-individual VF dynamics were assessed across 96 consecutive-timepoint transitions from 83 participants with ≥2 sequenced timepoints, quantifying genes gained, lost, and shared, and computing Jaccard similarity.

**Results**
Median VF gene counts were comparable across clinical states: ASB 80 (IQR 70–92), UTI 80.5 (IQR 70–98), and culture-negative 88 (IQR 77–100). Of 164 VF genes screened, *astA* (heat-stable enterotoxin) showed the largest enrichment in UTI (4/16 [25.0%] vs 9/136 [6.6%]; OR 4.63, p=0.033), though not significant after Benjamini-Hochberg correction (adjusted p=1.0). The salmochelin operon (*iroB/C/D/E/N*) was more prevalent in UTI (9/16 [56.2%]) than ASB (56/136 [41.2%]); Δ=+15.1 percentage points.

Longitudinally, 70/96 (72.9%) transitions showed no VF gene change; median Jaccard similarity was 1.000 (mean 0.925). Among 12 ASB→UTI transitions, 7 (58.3%) had no VF change. In the 5 transitions with gains, iron-acquisition genes (*iroB–N*; gained 4/12) and capsule genes (*kpsD/M*; 3/12) were most commonly acquired.

All 16 UTI episodes occurred at Uricult (UTI-suspicion) timepoints, whereas scheduled visits (T0/T1/T2) captured only ASB and culture-negative episodes. Between-status comparisons are therefore structurally confounded with sampling occasion.

**Conclusion**
In this exploratory analysis, *E. coli* VF profiles were remarkably conserved both across clinical states and within individuals over time. No VF gene was significantly enriched in UTI versus ASB after multiple testing correction. The iron-acquisition salmochelin operon (*iro*) and capsule biosynthesis genes (*kps*) were the most consistent candidates and warrant investigation in larger cohorts using mixed-effects models that account for repeated measures and timepoint confounding.

---

## 7) Repro Instructions

```bash
# 1. Independent verification (from scratch, no project scripts sourced)
cd /Users/Aamir/Desktop/rUTIs
Rscript /tmp/verify_vf_abstract.R 2>&1 | tee verification_output.txt

# 2. Original computation scripts (produce CSVs)
Rscript compute_vf_abstract_stats.R 2>&1 | tee abstract_stats_output.txt
Rscript compute_vf_stratified_by_depth.R 2>&1 | tee stratified_output.txt

# 3. Quick shell-only sanity checks (no R needed)
# Cohort size
tail -n +2 results/vf/vf_pa_all.csv | wc -l          # Expected: 183
tail -n +2 results/vf/vf_pa_all.csv | cut -d',' -f1 | sort -u | wc -l  # Expected: 87

# Gene count
head -1 results/vf/vf_pa_all.csv | tr ',' '\n' | tail -n +3 | wc -l   # Expected: 164

# Species
tail -n +2 assembly_metadata.csv | cut -d',' -f13 | sort | uniq -c     # Expected: 382 Escherichia coli

# Status × Timepoint
tail -n +2 status_map.csv | cut -d',' -f2,13 | sort | uniq -c | sort -rn
```
