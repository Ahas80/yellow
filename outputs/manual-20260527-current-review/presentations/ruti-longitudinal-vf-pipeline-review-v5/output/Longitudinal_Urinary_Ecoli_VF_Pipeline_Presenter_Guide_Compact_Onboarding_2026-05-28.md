# Longitudinal urinary E. coli virulence-factor review

Canonical Longcycler-only presenter guide for `v5-compact`.

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

### Slide 02 — Plain-English loop

**Say:** Walk from clinical episode to selected genome, binary VF row, longitudinal comparison, and finally the clinical overlay.

**Emphasise:** The genome source is selected QC-passing Longcycler only.

**Boundary:** Do not collapse episodes, residents, pairs, and transitions into one unit.

### Slide 03 — Denominator units

**Say:** Explain why the counts change: 583 source episodes retained only for attrition/QC context, 532 selected analytical episodes, 893 direct pairs, and 371 adjacent transitions.

**Emphasise:** Name the unit whenever a number appears.

**Boundary:** Source/QC counts do not replace analytical counts.

### Slide 04 — Operational phenotype

**Say:** Explain the project annotation that yields 16 operational UTI and 516 operational Not_UTI selected episodes.

**Emphasise:** Not_UTI is a mixed comparator.

**Boundary:** The operational label is not a universal clinical diagnosis.

### Slide 05 — Selected cohort

**Say:** The release cohort is 532 episodes, 161 residents, 16 operational UTI, and 516 operational Not_UTI.

**Emphasise:** There are 227 binary VF features in the current release.

**Boundary:** Not_UTI is a mixed comparator, not one biological state.

### Slide 06 — VF representation

**Say:** Each selected episode is represented by 227 binary VF features, then summarised through curated modules and longitudinal comparisons.

**Emphasise:** Genes remain the underlying evidence even when modules are shown.

**Boundary:** Presence/absence is not expression, activity, or causality.

### Slide 07 — Paired repertoire views

**Say:** Show gene-level prevalence and curated module structure together as two descriptive views of the same selected cohort.

**Emphasise:** The current release contains 227 binary VF features.

**Boundary:** Neither plot is a causal clinical model.

### Slide 08 — Pair-layer stability

**Say:** Describe VF similarity using 893 direct pairs and 371 adjacent transitions from 139 residents.

**Emphasise:** Use adjacent transitions for temporal questions.

**Boundary:** Do not reuse an earlier plot denominator.

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

### Slide 15 — Numbered pipeline

**Say:** Walk through the validated clinical-to-VF handoffs and point to the RQ01–RQ10 release layer.

**Emphasise:** Each stage consumes the same selected cohort.

**Boundary:** Do not introduce an alternate assembly path.

### Slide 16 — Complete mechanism casebook

**Say:** The focused mechanism casebook contains 9 cases, all 9 linked, with 0 missing endpoints.

**Emphasise:** Use the evidence matrix to show heterogeneity across cases.

**Boundary:** Buckets organise evidence; they do not prove mechanism.

### Slide 17 — Robustness and lineage diagnostics

**Say:** Keep the robustness plot beside the lineage PCoA. The 17-row near-miss audit remains separate from operational UTI cases.

**Emphasise:** These are backup interpretation aids.

**Boundary:** Neither diagnostic converts the sparse phenotype analysis into a confirmatory claim.
