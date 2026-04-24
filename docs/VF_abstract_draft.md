# VF Abstract Draft

## Title Options

1. **Virulence Factor Profiles of Urinary *Escherichia coli* Remain Largely Stable Within Elderly Residents Despite Transitions Between Asymptomatic Bacteriuria and Symptomatic UTI**
2. **Cross-Sectional and Longitudinal Analysis of *E. coli* Virulence Factor Repertoires in a Prospective Cohort of Elderly Residents With Recurrent Urinary Tract Infections**

---

## Background

Asymptomatic bacteriuria (ASB) is highly prevalent among elderly residents of long-term care facilities, yet most episodes do not progress to symptomatic urinary tract infection (UTI). Whether the virulence factor (VF) repertoire of colonizing *Escherichia coli* strains differs between ASB and UTI — and whether VF content changes within individuals over time — remains poorly characterized in this population. We aimed to systematically compare VF profiles between clinical states and to quantify within-host VF dynamics across serial timepoints.

## Methods

We analysed *E. coli* isolates from the YELLOW RoUTIne rUTI cohort, a prospective surveillance study of elderly residents. Urine samples were collected at scheduled visits (T0, T1, T2) and upon clinical suspicion of UTI (Uricult). Bacterial isolates were subjected to whole-genome sequencing (Oxford Nanopore Technologies), assembled using Flye and/or LongCycler, and screened for virulence factor genes using Abricate with the VFDB database (minimum 80% nucleotide identity, minimum 80% query coverage). Clinical status (ASB, UTI, or culture-negative) was assigned based on clinical and microbiological criteria recorded in the study protocol.

A total of 183 participant×timepoint observations were available from 87 participants (136 ASB, 16 UTI, 31 negative; all *E. coli*). Virulence factor burden was defined as the number of distinct VFDB genes detected per sample (164 genes total). Between-group comparisons used Fisher's exact test (exploratory, unadjusted). Longitudinal VF dynamics were assessed by computing consecutive-timepoint transitions for participants with ≥2 sequenced timepoints (n=83 participants, 96 transitions), quantifying VF genes gained, lost, and shared, and calculating Jaccard similarity of VF gene sets.

## Results

### VF Burden

Median VF gene counts were comparable across clinical states: ASB 80 (IQR 70–92), UTI 80.5 (IQR 70–98), and culture-negative 88 (IQR 77–100). No statistically significant difference in overall VF burden was observed between ASB and UTI episodes.

### Gene-Level Prevalence

Of 164 VF genes screened, the largest prevalence differences between UTI and ASB were observed for:
- **astA** (heat-stable enterotoxin): 25.0% UTI vs 6.6% ASB (Δ=+18.4 pp; Fisher exact OR=4.63, p=0.033, not significant after BH correction)
- **pic** (serine protease autotransporter): 0.0% UTI vs 18.4% ASB (Δ=−18.4 pp; p=0.075)
- **iro operon** (salmochelin siderophore: iroB/C/D/E/N): 56.2% UTI vs 41.2% ASB (Δ=+15.0 pp each)

No genes reached significance after Benjamini-Hochberg correction (all adjusted p>0.05), consistent with limited UTI sample size (n=16).

### Longitudinal VF Dynamics

Among 96 consecutive transitions from 83 participants, VF profiles were highly conserved. The overall median Jaccard similarity was 1.000 (mean 0.925), and 72.9% (70/96) of transitions showed no VF gene content change whatsoever. In the 12 ASB→UTI transitions observed, 7/12 (58.3%) had no VF change. Among the 5 transitions with VF gains, iron-acquisition genes (iroB/C/D/E/N, gained in 4/12) and capsule genes (kpsD/M, gained in 3/12) were the most commonly acquired.

## Conclusion

In this exploratory analysis of elderly residents with recurrent UTIs, *E. coli* VF repertoires were remarkably stable both between clinical states and within individuals over time. No individual VF gene was significantly enriched in UTI versus ASB after multiple testing correction, consistent with the hypothesis that host factors rather than bacterial virulence gene content may primarily determine clinical outcome. The iron acquisition (iro) and capsule (kps) gene clusters were the most commonly gained in ASB→UTI transitions, warranting further investigation in larger cohorts with appropriate longitudinal mixed models.

---

> [!NOTE]
> **Limitations acknowledged in Methods/Discussion**:
> - Small UTI sample (n=16), all from Uricult timepoint
> - Fisher tests are exploratory (repeated measures violate independence)  
> - VF detection by Abricate/VFDB at ≥80% identity — threshold-dependent
> - Timepoint ordering: Uricult is event-triggered, not at fixed intervals
> - All numbers from [compute_vf_abstract_stats.R](file:///Users/Aamir/Desktop/rUTIs/compute_vf_abstract_stats.R) using anchor files only
