# VF Next Steps Checklist

## ✅ What Is Verified and Reproducible Now

- [x] **VF detection pipeline**: Abricate+VFDB (≥80% id, ≥80% cov) traced from assemblies to P/A matrix
- [x] **Analysis-ready dataset**: 183 rows (87 participants), 100% match to clinical status, all *E. coli*
- [x] **VF burden by status**: No significant difference (ASB median=80, UTI median=80.5, Neg median=88)
- [x] **Gene-level prevalence**: 164 genes tested; top candidate astA (OR=4.63, p=0.033, BH-adjusted p>0.05)
- [x] **Longitudinal dynamics**: 96 transitions from 83 participants; 72.9% with zero VF change; median Jaccard=1.0
- [x] **Cross-checks**: Discrepancies with existing outputs fully explained by definitional differences
- [x] **Traceability**: Every number linked to anchor files, joins, and code

## ⚠️ Known Limitations

### Sample Size
- Only **16 UTI episodes** (all from Uricult timepoints) — severely limits statistical power
- **No UTIs at scheduled timepoints** (T0/T1/T2) — UTI vs ASB comparisons are confounded with timepoint

### Methodological
- **Fisher tests violate independence**: Same participant contributes multiple episodes → pseudoreplication. Fisher p-values are exploratory only.
- **GLMM degenerate**: Existing GLMM results (OR≈1.0, p=1) are artifacts of timepoint confounding (Infection_Status near-perfectly predicted by tp_lab=Uricult)
- **VF detection threshold**: Gene presence is binary (≥80% id & cov). Copy number variation, expression levels, and phase variation are not captured.
- **Category annotation heuristic**: `gene_map.csv` categories are regex-based, not manually curated. ~26 genes per sample are "Unassigned."
- **Timepoint ordering**: Uricult is event-triggered (UTI suspicion), not a fixed time interval from T0/T1/T2. This creates a structural confound.

### Coverage
- **93 clinical episodes lack VF data**: 77 ASB, 9 Negative, 6 UTI, 1 Culture-positive — these are episodes without sequenced isolates
- **Sparse deep follow-up**: Only 8 participants have ≥3 timepoints, 5 have ≥4

## 🔲 Concrete Next Steps for Publication-Grade Inference

### Priority 1: Address Confounding (Critical)
- [ ] **Acknowledge in abstract/methods**: UTI status and Uricult timepoint are confounded — any VF difference could reflect sampling occasion rather than true virulence
- [ ] **Sensitivity analysis**: For participants who have BOTH scheduled AND Uricult samples, compare within-person VF profiles at scheduled vs Uricult visits
- [ ] **Alternative comparisons**: Compare VF profiles of participants who EVER had UTI vs those who never did (participant-level, not episode-level) to avoid pseudoreplication

### Priority 2: Appropriate Longitudinal Models
- [ ] **GEE or mixed models**: Replace Fisher tests with GEE (working correlation = exchangeable) or GLMM with participant random intercept, using only participants with both ASB and UTI episodes
- [ ] **Longitudinal VF trajectory modelling**: Linear mixed model of VF count over time, with participant random slope, testing interaction between time and episode type

### Priority 3: VF Annotation Quality
- [ ] **Manual curation of gene_map.csv**: Replace regex heuristic with VFDB's own category annotations (download gene-to-category mapping from VFDB website)
- [ ] **Remove/reclassify AMR genes**: Some genes classified as AMR by the heuristic are not true virulence factors

### Priority 4: Expand Coverage
- [ ] **Sequence missing UTI isolates**: 6 UTI episodes lack VF data — these should be priority for sequencing
- [ ] **Increase temporal coverage**: More participants with ≥3 timepoints would strengthen longitudinal claims

### Priority 5: Complementary Analyses
- [ ] **Plasmid-linked VF analysis**: The iro operon (most commonly gained in ASB→UTI) is typically plasmid-borne — check correlation with PlasmidFinder results
- [ ] **ST-stratified VF profiles**: VF repertoire varies by ST; controlling for ST would clarify whether VF differences are strain-level or status-level
- [ ] **Phylogenetic correction**: VF carriage is phylogenetically structured — phylogenetic logistic regression would account for population structure
