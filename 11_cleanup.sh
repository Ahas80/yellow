#!/usr/bin/env bash
# ==============================================================================
# 11_cleanup.sh
# ------------------------------------------------------------------------------
# Purpose: Organize output files into docs/ and results/legacy/ and clean up.
# Inputs: Various files in root.
# Outputs: Organized directories.
# ==============================================================================

set -euo pipefail
shopt -s nullglob

# 1. Setup Directories
mkdir -p docs/slides
mkdir -p docs/papers
mkdir -p results/legacy

# 2. Move Files (Safe Moves)
move_if_exists() {
    if [ -e "$1" ]; then
        echo "Moving $1 to $2"
        mv -i "$1" "$2"
    fi
}

# Slides
move_if_exists "Virulence factors.pptx" "docs/slides/"
move_if_exists "Poster_upec.pptx" "docs/slides/"
move_if_exists "Longread_UTI_interim.pptx" "docs/slides/"
move_if_exists "YELLOW_RoUTIne_Key_Results_Slide copy.pptx" "docs/slides/"
move_if_exists '~$Virulence factors.pptx' "docs/slides/"

# Papers & Docs
move_if_exists "Recurrent E. coli Urinary Tract Infections in Nursing Homes- Insight in Sequence Types and Antibiotic Resistance Patterns.pdf" "docs/papers/"
move_if_exists "Bacterial virulence phenotypes of Escherichia coli and host susceptibility determines risk for urinary tract infections.pdf" "docs/papers/"
move_if_exists "YELLOW RoUTIne.pdf" "docs/papers/"
move_if_exists "Virulence Mechanisms of Common Uropathogens and Their Intracellular Localisation within Urothelial Cells.pdf" "docs/papers/"
move_if_exists "Whole-genome characterisation of Escherichia coli isolates from older women with urinary tract infection or asymptomatic bacteriuria.docx" "docs/papers/"
move_if_exists "e.coli_manuscript_vs7_MB_CS_MvdB.docx" "docs/papers/"
move_if_exists "dkac307.pdf" "docs/papers/"
move_if_exists "Analysis of Enterobacter hormaechei.docx" "docs/papers/"
move_if_exists '~$ole-genome characterisation of Escherichia coli isolates from older women with urinary tract infection or asymptomatic bacteriuria.docx' "docs/papers/"
move_if_exists "Rplots.pdf" "docs/papers/"

# Move Plots to Papers (if they exist in plots/)
move_if_exists "plots/vf/variable_gene_heatmap.pdf" "docs/papers/"
move_if_exists "plots/clinical/overview_within_person_status_all_samples.pdf" "docs/papers/"
move_if_exists "plots/clinical/overview_plots.pdf" "docs/papers/"
move_if_exists "plots/clinical/within_person_status_all_samples_single_panel.pdf" "docs/papers/"
move_if_exists "plots/clinical/overview_status_story.pdf" "docs/papers/"
move_if_exists "plots/clinical/within_person_status_pid_by_tp_subdivided.pdf" "docs/papers/"
move_if_exists "plots/clinical/within_person_status_all_individuals_bubbles.pdf" "docs/papers/"
move_if_exists "plots/mlst/top20_STs.pdf" "docs/papers/"
move_if_exists "plots/plasmids/ST_vs_replicon_network.pdf" "docs/papers/"
move_if_exists "plots/plasmids/replicon_cooccurrence.pdf" "docs/papers/"
move_if_exists "plots/plasmids/replicon_heatmap.pdf" "docs/papers/"

# Legacy Data
move_if_exists "vfdb_hits_all_assemblies.csv" "results/legacy/"
move_if_exists "vf_hits.rds" "results/legacy/"
move_if_exists "vfdb_presence_absence_matrix.csv" "results/legacy/"
move_if_exists "vfdb_gain_loss_by_participant.csv" "results/legacy/"

# 3. Delete Temporary Files
rm -f triangle.ids triangle.sorted

echo "✓ Cleanup complete."
