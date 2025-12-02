# Timepoint vs Isolate Terminology Clarification

## The Issue
Some participants appeared to have 6-10 "timepoints" in Panaroo selection plots, which was confusing since the study design only has T0-T4 + Uricult (max 6 biological timepoints).

## The Explanation
**Multiple assemblers per sample create multiple isolates:**

- Each biological sample is assembled using **2 different assemblers** (e.g., Flye and Unicycler)
- Each assembly creates a separate GFF file for Panaroo
- Example: Participant 20031 has:
  - 4 biological timepoints: T0, T1, T2, Uricult
  - 8 assembly files (isolates): 4 samples × 2 assemblers = 8 entries

## The Fix
Updated `13_visualise_panaroo_selection.R` to:
- Changed plot labels from "timepoints" → "isolates"
- Added subtitle: "Isolates = assembly files passing QC (note: multiple assemblers per sample)"
- Added comments clarifying that counts represent assembly files, not unique biological samples

## Scripts Reviewed
✅ **Script 13** - Updated (Panaroo counts isolates)  
✅ **Script 15** - No change needed (uses biological timepoints from status_map)  
✅ **Script 21** - No change needed (uses biological timepoints from status_map)  

## Key Distinction
- **Isolates**: Assembly files (can be 2+ per biological sample due to multiple assemblers)
- **Timepoints**: Biological sampling occasions (T0, T1, T2, etc.)

In Panaroo analysis specifically, we're counting **isolates** (assembly files), not **timepoints** (sampling occasions).
