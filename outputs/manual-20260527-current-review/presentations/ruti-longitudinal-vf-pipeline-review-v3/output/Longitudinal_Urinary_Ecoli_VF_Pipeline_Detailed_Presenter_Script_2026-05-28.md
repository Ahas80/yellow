# Longitudinal urinary E. coli virulence-factor review

Canonical Longcycler-only presenter guide for `v3`.

## Release anchors

- Analytical cohort: 532 episodes; 161 residents; 16 operational UTI; 516 operational Not_UTI.
- Longitudinal evidence: 893 direct pairs; 371 adjacent transitions from 139 residents; 140 at or below 25 SNPs.
- Focused transition view: 9 Not_UTI→UTI; 5 at or below 25 SNPs; casebook 9/9/0 cases/linked/missing.
- Near-miss audit: 17 rows, separate from operational UTI cases.
- Research questions: RQ01–RQ10.
- Interpretation: exploratory observational analysis; no causal claim.

## Release provenance

- Registry: `/Users/Aamir/Desktop/rUTIs/results/pipeline/longcycler_release_claim_registry.json`
- Registry SHA-256: `87f65979259a4c92545afbb74b5d3a8efffb8a3d18cd2481df165b97e5c6d30e`

## Slide-by-slide script

### Slide 01 — Title and release anchors

**Say:** Open with the selected analytical denominator: 532 episodes from 161 residents, with 371 adjacent transitions and 16 operational UTI episodes.

**Emphasise:** The whole review is selected QC-passing Longcycler only.

**Boundary:** Operational UTI is an exploratory annotation, not the organising genome source.

### Slide 02 — Scope lock

**Say:** One selected assembly per episode creates 532 episode-level genome profiles. The temporal layer is 371 adjacent transitions from 139 residents.

**Emphasise:** Name the assembly policy and analysis unit before discussing results.

**Boundary:** Do not mix source/QC genomes into the analytical cohort.

### Slide 03 — Selected cohort

**Say:** The release cohort is 532 episodes, 161 residents, 16 operational UTI, and 516 operational Not_UTI.

**Emphasise:** There are 227 binary VF features in the current release.

**Boundary:** Not_UTI is a mixed comparator, not one biological state.

### Slide 04 — Numbered pipeline

**Say:** Walk through the validated clinical-to-VF handoffs and point to the RQ01–RQ10 release layer.

**Emphasise:** Each stage consumes the same selected cohort.

**Boundary:** Do not introduce an alternate assembly path.

### Slide 05 — VF representation

**Say:** Each selected episode is represented by 227 binary VF features, then summarised through curated modules and longitudinal comparisons.

**Emphasise:** Genes remain the underlying evidence even when modules are shown.

**Boundary:** Presence/absence is not expression, activity, or causality.

### Slide 06 — Gene repertoire

**Say:** Read the ranked genes as a descriptive view of what is commonly detected in the selected cohort.

**Emphasise:** The denominator is the selected episode-level genome set.

**Boundary:** Do not interpret prevalence ranking as a UTI association.

### Slide 07 — Curated modules

**Say:** Explain that modules organise genes into interpretable systems for navigation and summary.

**Emphasise:** Curation improves readability without changing binary gene evidence.

**Boundary:** Modules are not validated disease-causality scores.

### Slide 08 — Pair layers

**Say:** Separate 893 direct within-resident pairs from 371 adjacent transitions. The first is broad; the second is temporally ordered.

**Emphasise:** 140 adjacent transitions are at or below the operational 25-SNP threshold.

**Boundary:** Always name the pair unit beside the count.

### Slide 09 — Gain and loss

**Say:** Use the 371 adjacent transitions to identify candidate VF gains and losses for follow-up.

**Emphasise:** Interpret change with lineage, genome distance, and QC context.

**Boundary:** A detected change may reflect replacement or calling differences as well as true gene-content change.

### Slide 10 — MLST context

**Say:** Preferred MLST labels are available for 493 selected episodes across 76 distinct labels.

**Emphasise:** Use lineage as a diagnostic layer before status interpretation.

**Boundary:** ST agreement does not prove the same strain.

### Slide 11 — Exploratory clinical overlay

**Say:** The selected phenotype split is 16 operational UTI versus 516 operational Not_UTI episodes.

**Emphasise:** The clinical model is hypothesis-generating and participant-aware.

**Boundary:** No causal or confirmatory VF status claim is supported by the sparse phenotype denominator.

### Slide 12 — Aggregate focused transitions

**Say:** There are 9 focused Not_UTI-to-UTI adjacent transitions; 5 are at or below 25 SNPs. The casebook is 9/9/0 cases/linked/missing.

**Emphasise:** Report the aggregate casebook, not a single anecdote.

**Boundary:** Low genome distance cannot identify host-state change or bacterial regulation by itself.

### Slide 13 — Evidence, boundary, next steps

**Say:** Land on three messages: the selected longitudinal evidence is coherent; the clinical phenotype remains sparse; and follow-up must be lineage-aware.

**Emphasise:** Ask what extra evidence would distinguish bacterial change from host-state and sampling effects.

**Boundary:** Keep causal interpretation outside the supported claim set.

### Slide 14 — Operational handover

**Say:** Point to the selected manifest, clinical key, VF matrix, pair tables, transition casebook, and RQ01–RQ10 runner.

**Emphasise:** The release claim registry is the provenance anchor.

**Boundary:** Do not use an unregistered intermediate as a release source.

### Slide 15 — Denominator audit

**Say:** The source clinical context is 583 episodes and 166 residents (18 UTI, 565 Not_UTI). It is retained only for attrition/QC context.

**Emphasise:** Analytical claims remain 532/161/16/516.

**Boundary:** Never substitute source/QC counts into analytical claims.

### Slide 16 — Module prevalence annotation

**Say:** Use module prevalence by status as a descriptive clinical annotation view.

**Emphasise:** Keep the operational phenotype denominator visible.

**Boundary:** Repeated measures and lineage are not fully resolved by a simple prevalence plot.

### Slide 17 — Complete mechanism casebook

**Say:** The focused mechanism casebook contains 9 cases, all 9 linked, with 0 missing endpoints.

**Emphasise:** Use the evidence matrix to show heterogeneity across cases.

**Boundary:** Buckets organise evidence; they do not prove mechanism.

### Slide 18 — Stable strain and changing clinical state

**Say:** Among the 9 focused transitions, 5 are at or below 25 SNPs.

**Emphasise:** Stable genomic context can motivate host-state or regulation hypotheses.

**Boundary:** It cannot establish either hypothesis without additional evidence.

### Slide 19 — Robustness boundary

**Say:** The operational phenotype contains only 16 UTI episodes, and the 17-row near-miss audit is separate.

**Emphasise:** Use robustness diagnostics to define uncertainty.

**Boundary:** Near-miss rows are not operational UTI cases.

### Slide 20 — Lineage diagnostic

**Say:** Use the VF-profile PCoA by preferred sequence type to assess lineage structure.

**Emphasise:** Lineage structure should be considered before interpreting status patterns.

**Boundary:** Clustering is diagnostic context, not causal evidence.

### Slide 21 — Accessory context

**Say:** Keep accessory and mobile-element information as transition-level context when it helps interpret a specific pattern.

**Emphasise:** The VF review remains the scientific centre.

**Boundary:** No dedicated AMR association or causal accessory claim is made.

### Slide 22 — Release registry and sources

**Say:** Close the appendix by showing the registry path, SHA-256, selected denominators, pair layers, focused casebook, near-miss audit, and RQ01–RQ10 marker.

**Emphasise:** Use the registry for provenance questions.

**Boundary:** Do not answer release-count questions from memory or an earlier deck.
