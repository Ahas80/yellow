# rUTIs report catalog

## Current source of truth

The active release is the registry-bound July 2026 Longcycler-only analysis:

- 532 analytical episodes from 161 residents.
- 16 operational UTI and 516 operational Not_UTI episodes.
- 83/83 pipeline verification checks passed before cleanup.
- 580/580 communication-deliverable checks passed before cleanup.

Start with these files:

1. `FOLDER_MAP.md` — active paths, release contract, denominators, and interpretation guardrails.
2. `docs/LECTURER_README.md` — concise current explanation for presenting the work.
3. `results/summary/final_key_results_summary.md` — current analytical results summary.
4. `results/pipeline/longcycler_release_claim_registry.json` — machine-readable release claims.
5. `results/qc/longcycler_only_pipeline_verification.txt` and `results/qc/longcycler_release_deliverables_verification.txt` — release gates.

For the supporting genomic-AMR layer, use
`results/summary/table_13_genomic_amr_summary.csv` for numeric results,
`results/amr/interpretation_report.md` for the plain-language interpretation,
and the three registered supplementary plots under `plots/amr/`. These are
genomic determinant and predicted-phenotype results, not phenotypic AST.

The older date-like names under `outputs/manual-20260526-*` and
`outputs/manual-20260527-*` do not mean those decks are stale. Their contents
were regenerated and registry-verified on 15 July 2026; keep them at their
current paths.

## Categories

| Category | Entries | Size | How to use |
|---|---:|---:|---|
| Current authoritative | 25 | 0.21 MiB | Use for current findings, denominators, methods, and release status. |
| Current deliverables | 30 | 14.31 MiB | Current decks, presenter guides, codebook, and communication companions. |
| Useful technical reference | 29,097 | 58.08 MiB | Retain for technical interpretation, QA, and reproducibility. Most are low-level tool reports rather than human-facing summaries. |
| Historical but useful | 832 | 197.20 MiB | Retain for development, cleanup, QA-render, and provenance history; do not treat as the current results summary. |
| Superseded — do not use for current findings | 120 | 4.28 MiB | Preserved only for traceability because the denominator, terminology, workflow state, or interpretation is obsolete. |
| Temporary/valueless | 186 | 1.52 MiB | Deleted system metadata, caches, temporary inspection ledgers, broken symlinks, and their empty directories. |

The exhaustive 30,290-row file-level catalog is
`archive/cleanup_2026-07-15/manifests/report_inventory.csv`. It records each
path, format, size, modification time, category, usefulness, reason, canonical
replacement, and cleanup action.

## Current deliverables

The verified communication package remains in place and includes:

- 11 canonical PowerPoint decks across the lecturer, methods-summary, onboarding, scientific-review, VF-review, and longitudinal-review packages.
- Presenter guides and current Markdown companions for the detailed, onboarding, and compact variants.
- `outputs/codebooks/vf_analysis_ready_lay_codebook.docx` and its delivered PDF.
- The current methodology register, audit findings, provenance tables, handout, talking points, and flowchart counts.

Use `results/qc/longcycler_release_deliverables_verification.csv` for the exact
verified inventory and `results/qc/longcycler_release_visual_qa.csv` for the
visual-QA ledger.

## Superseded reports

Superseded scientific documents now live under
`docs/legacy_asb_uti_docs/superseded_2025_2026/`. Historical figure guidance is
under `docs/legacy_asb_uti_docs/figures/`. These include the older 276/361 and
579/550 cohort narratives, legacy ASB-vs-UTI interpretation, old WGS restart
notes, and pre-release methods or figure guidance.

Obsolete generated summaries and retired participant TeX reports now live under
`archive/cleanup_2026-07-15/superseded_generated_reports/`. Their canonical
replacement is the current key-results summary, folder map, and lecturer guide
listed above.

## Cleanup records

The full movement and deletion report is
`archive/cleanup_2026-07-15/CLEANUP_REPORT.md`. Machine-readable records are in
`archive/cleanup_2026-07-15/manifests/`, including verified move hashes,
deletion reasons, report classifications, filesystem inventories, Git-status
snapshots, and post-cleanup verification outputs.
