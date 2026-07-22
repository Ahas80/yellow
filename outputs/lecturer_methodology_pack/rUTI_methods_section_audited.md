# Audited methods section

Clinical episodes were classified before genomic analysis using the current versioned operational UTI phenotype, which combines implemented culture-support and catheter-aware symptom rules. This operational phenotype is not presented as a reconstruction of the full published protocol, and the Not_UTI category is not interpreted as a healthy control.

Genomic analyses were restricted to an exact cohort of 532 genome-linked episodes from 161 residents. Each analytical episode was represented by one selected canonical Longcycler assembly that passed the implemented assembly-QC screen. The selected manifest retained the input path and FASTA SHA-256 digest, and the episode keys matched the analytical clinical cohort exactly. The cohort contained 16 operational UTI and 516 operational Not_UTI episodes. The full clinical source (583 episodes from 166 residents) was retained only for labelled attrition/QC context.

Virulence-factor screening produced 227 binary VFDB features on the selected cohort. MLST calls were used as lineage context; provider calls were accepted when key-linked to selected inputs, with labelled local fallback where required. Direct pairwise genomic comparisons covered all 893 unordered within-resident selected-episode pairs. Pair-specific direct distance evidence was treated as primary, so graph transitivity or ST agreement could not overrule a direct distance above the operational ≤25-SNP boundary.

Longitudinal comparisons were rebuilt from the selected cohort in clinical time order. This produced 371 adjacent transitions from 139 residents, of which 140 were at or below 25 SNPs. There were 9 Not_UTI-to-UTI transitions; 5 met the operational SNP boundary. The focused mechanism casebook contained 9 cases, all 9 linked and 0 missing. The 17 near-miss rows were retained as a separate audit population and were not counted as operational UTI cases.

All analyses were interpreted as exploratory and observational. RQ01–RQ10 formed the numbered release layer. No causal UTI mechanism, antibiotic effect, demographic or host-factor effect, named mutation conclusion, or cross-method distance equivalence was inferred.
