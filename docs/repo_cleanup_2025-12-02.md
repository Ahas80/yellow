# Repository Cleanup Log - 2025-12-02

Comprehensive cleanup and reorganization of the rUTIs repository to create a clean, lecturer-friendly structure.

## Summary
- **Root files reduced:** 94 → 54 (43% reduction)
- **Files deleted:** 0 (all preserved in organized locations)
- **Canonical scripts identified:** 11_compare_strains.R, 12_wgs_exact_compare.R
- **Documentation created:** FOLDER_MAP.md, updated .gitignore

## Major Moves

### Legacy Scripts
- `legacy/11_compare_variants/`: 11_compare_strains_MOD.R
- `legacy/12_wgs_variants/`: All 12_wgs variants (.RY, _clean, .prepatch*, 12_wgs_runner.R, 12e_report.R)

### Debug & Development
- `legacy/debug_scripts/`: debug_panaroo*.sh (3 files), debug_snps.R
- `legacy/temp_dev_files/`: temp_*.R, test_*.R, verify_*.R, repair_prokka.R, .Rhistory, .Rapp.history, temp_prokka_test_*
- `legacy/installers/`: miniforge_x86.sh, igraph_*.tgz, rUTIs_backup_*.tar.gz

### Organization
- `scripts/`: cleanup_plan.sh, cleanup_root_clutter.sh
- `scripts/awk/`: inject_safe.awk, patch_progress.awk
- `docs/`: YELLOW_RoUTIne_R_Script_Breakdown.txt, FINAL_SUMMARY.md

### Plots
- `plots/participant_specific/`: ~50 participant-specific diagnostic plots
- Publication figures already well-organized in plots/publication/, plots/clinical/, etc.

## Key Documentation

See [FOLDER_MAP.md](file:///Users/Aamir/Desktop/rUTIs/FOLDER_MAP.md) for complete repository structure guide.

## Results

✅ Repository is now clean, organized, and lecturer-friendly  
✅ All code preserved (nothing deleted)  
✅ Canonical scripts clearly identified  
✅ Plots directory curated for presentation
