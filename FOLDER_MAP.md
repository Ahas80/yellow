# Current Longcycler-only folder map

The paths below identify the active release inputs, analysis products, quality gates and communication artifacts.

## Active paths

| Area | Current path | Role |
|---|---|---|
| Clinical cohort | `results/clinical/analysis_cohort_longcycler.csv` | Exact analytical episode set |
| FASTA provenance | `results/qc/analysis_assembly_manifest.csv` | Selected paths and content hashes |
| VF | `results/vf/` | Selected-cohort VF calls and matrices |
| MLST | `results/mlst/` | Lineage context |
| Direct comparisons | `results/strain_compare/` | Pair-specific genomic distances |
| Longitudinal | `results/longitudinal/` | Canonical adjacent transition table |
| Genomic AMR | `results/amr/` | Validated determinant, prevalence, longitudinal, prediction and provenance tables |
| Genomic AMR results | `results/summary/table_13_genomic_amr_summary.csv` | Numeric supporting results used in the final narrative |
| Mechanism | `results/mechanism/` | Linked focused-transition casebook |
| Research questions | `results/research_questions/RQ01` to `RQ10` | Audited analyses |
| Registry and markers | `results/pipeline/` | Release contract and completion state |
| QC reports | `results/qc/` | Cohort, provenance and delivery gates |
| Figures | `plots/` | Generated analytical figures |
| Communication | `outputs/` | Registry-bound decks, handouts and codebook |

## Audited release contract

- Scope: Longcycler-only selected QC-passing assemblies; one selected assembly per analytical episode.
- Clinical definition: operational UTI phenotype. It is a versioned culture-plus-compatible-symptom rule, not a reconstruction of the full published protocol.
- Analytical cohort: 532/161/16/516 (episodes/residents/operational UTI/operational Not_UTI).
- Direct evidence: 893 all within-resident pairs.
- Adjacent evidence: 371/139/140 (pairs/residents/pairs at or below 25 SNP).
- Focused transitions: 9 Not_UTI -> UTI; 5/9 at or below 25 SNP.
- Mechanism casebook: 9/9/0 (cases/linked/missing).
- Near-miss audit: 17 rows; these are not operational UTI cases.
- Attrition/QC context only: 583/166/18/565 (full-source episodes/residents/operational UTI/operational Not_UTI); these are not analytical denominators.
- Research-question boundary: RQ01-RQ10 only.
- Interpretation: exploratory, observational and non-causal.

## Methods fixed by the registry

- Operational phenotype: culture lower bound >=1,000 CFU/mL plus compatible symptoms under the versioned rule.
- Assembly QC: <=200 contigs, N50 >=20,000 bp and genome size 4,000,000-6,000,000 bp; read coverage, completeness and contamination are excluded metrics.
- VF calls: ABRicate with VFDB at >=80% identity and >=80% coverage, SHA-bound to the selected Longcycler FASTA manifest.
- Genomic AMR: AMRFinderPlus 4.2.7 acquired genes and known resistance mutations are primary; ResFinder 4.7.2 with PointFinder is complementary and ABRicate-ResFinder at 80/80 is the legacy comparison. These are genomic predictions, not phenotypic AST.
- MLST: lineage context only; provider calls require >=95% good targets. Policy: provider_qc95 call key/path-linked to the selected Longcycler episode; local fallback excluded. When required, local calls are labelled and use the same selected FASTA.
- Pair-specific evidence: dnadiff is primary, with an operational threshold of <=25 SNP; MLST or graph context cannot overrule a conflicting direct pair.
- Population context: Parsnp core-genome and Panaroo pangenome outputs provide context, not pair-specific continuity proof.

## Provenance

- Claim registry: `results/pipeline/longcycler_release_claim_registry.json`
- Registry SHA-256: `87f65979259a4c92545afbb74b5d3a8efffb8a3d18cd2481df165b97e5c6d30e`
- Registry generated: 2026-07-15 01:36:53 CEST
- This file is generated; edit the registry-producing analysis or this writer rather than hand-editing release claims.
