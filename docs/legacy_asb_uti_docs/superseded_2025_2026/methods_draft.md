# Methods Draft

## Longitudinal Analysis & Strain Tracking
To reconstruct patient infection timelines, we analyzed 276 clinical episodes from 87 participants. We defined strong "Same Strain" persistence using the prior YELLOW study SNP threshold: 0-25 core-genome SNPs. Pairs with >25 SNPs were treated as above the same-strain SNP threshold, and pairs with missing SNPs were not promoted to same-strain status by ST alone. ST was reported separately as secondary lineage context (`Same ST`, `Different ST`, or `Missing ST evidence`), while pairwise `Different` calls and different-ST pairs were used to flag likely replacement when SNPs did not support same-strain persistence. Pairwise comparisons were performed using `MASH` for rapid screening and `nucmer` for precise SNP counting. Isolates were assigned to global "Strain IDs" using a graph-based clustering approach (igraph), where nodes represented isolates and edges represented "Same Strain" relationships.

## Phenotype Switch Characterization
We identified "Phenotype Switch" events where a patient transitioned between Asymptomatic Bacteriuria (ASB) and Symptomatic UTI (or vice versa) while carrying the *same* bacterial strain. For these candidates, we performed deep genomic characterization:
1.  **Variant Calling**: We aligned the "To" isolate assembly against the "From" isolate assembly using `nucmer` (MUMmer4) to identify Single Nucleotide Polymorphisms (SNPs) and Indels.
2.  **Annotation**: Variants were annotated using the Prokka-generated GFF files of the reference isolate. We distinguished between Intergenic variants (potential promoter mutations) and Coding Sequence (CDS) variants.
3.  **Gene Content**: We compared the presence/absence of Virulence Factors (VFDB) and Plasmid Replicons (PlasmidFinder) to determine if the phenotype switch was driven by gene gain/loss.

## Virulence Factor Marker And Module Framework
Because no validated UTI-specific VF score exists for this cohort, composite VF measures were treated as supplementary rather than primary endpoints. The primary biological interpretation used prespecified marker groups, module-level VF repertoire patterns, lineage context, and within-resident longitudinal VF stability.

We replaced the previous UPEC-candidate gene-count endpoint with two explicitly labelled supplementary measures. First, an ExPEC-like marker classifier used the conventional "at least two of five marker groups" rule: P fimbriae (`papA` or `papC`), S/F1C fimbriae (`sfaD`, `sfaE`, `focD`, or `focE`), Afa/Dr adhesins (`afaB*`, `afaC*`, `draB`, or `draC`), aerobactin receptor (`iutA`), and capsule marker (`kpsM`). The available `kpsM` call was interpreted as a `kpsM` proxy and was not described as confirmed `kpsM II` unless subtype annotation supported that claim. Second, a UPEC system count gave one point per prespecified UPEC-associated biological module present, rather than one point per gene. Total and curated VF gene counts were retained only as descriptive burden measures.

## Lineage Analysis
We determined the Multi-Locus Sequence Type (ST) of all isolates using `mlst`. We calculated the "UTI Risk" for each ST as the proportion of episodes manifesting as symptomatic UTI. Fisher's Exact Test (with FDR correction) was used to assess if any specific lineage was significantly associated with symptomatic infection compared to the rest of the cohort.
