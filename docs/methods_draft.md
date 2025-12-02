# Methods Draft

## Longitudinal Analysis & Strain Tracking
To reconstruct patient infection timelines, we analyzed 276 clinical episodes from 87 participants. We defined "Same Strain" persistence based on a strict genomic threshold: >99.9% Average Nucleotide Identity (ANI) and <50 SNPs difference in the core genome. Pairwise comparisons were performed using `MASH` for rapid screening and `nucmer` for precise SNP counting. Isolates were assigned to global "Strain IDs" using a graph-based clustering approach (igraph), where nodes represented isolates and edges represented "Same Strain" relationships.

## Phenotype Switch Characterization
We identified "Phenotype Switch" events where a patient transitioned between Asymptomatic Bacteriuria (ASB) and Symptomatic UTI (or vice versa) while carrying the *same* bacterial strain. For these candidates, we performed deep genomic characterization:
1.  **Variant Calling**: We aligned the "To" isolate assembly against the "From" isolate assembly using `nucmer` (MUMmer4) to identify Single Nucleotide Polymorphisms (SNPs) and Indels.
2.  **Annotation**: Variants were annotated using the Prokka-generated GFF files of the reference isolate. We distinguished between Intergenic variants (potential promoter mutations) and Coding Sequence (CDS) variants.
3.  **Gene Content**: We compared the presence/absence of Virulence Factors (VFDB) and Plasmid Replicons (PlasmidFinder) to determine if the phenotype switch was driven by gene gain/loss.

## Lineage Analysis
We determined the Multi-Locus Sequence Type (ST) of all isolates using `mlst`. We calculated the "UTI Risk" for each ST as the proportion of episodes manifesting as symptomatic UTI. Fisher's Exact Test (with FDR correction) was used to assess if any specific lineage was significantly associated with symptomatic infection compared to the rest of the cohort.
