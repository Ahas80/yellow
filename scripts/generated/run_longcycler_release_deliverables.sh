#!/usr/bin/env bash
# Regenerate the canonical Longcycler-only presentation and handout release.
#
# This script is intentionally separate from RUN_COMPLETE_ANALYSIS.sh. It may
# run only after the final analysis marker and audited claim registry exist.
# Every presentation workspace is external to the repository; only named final
# artifacts and their generator-owned support files are written under outputs/.

set -Eeuo pipefail

SCRIPT_PATH="${BASH_SOURCE[0]}"
PROJECT_ROOT="$(cd "$(dirname "${SCRIPT_PATH}")/../.." && pwd -P)"
cd "${PROJECT_ROOT}"

NODE_BIN="${NODE_BIN:-/Users/Aamir/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node}"
PYTHON_BIN="${PYTHON_BIN:-/Users/Aamir/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3}"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
DOC_SKILL_DIR="${DOC_SKILL_DIR:-/Users/Aamir/.codex/plugins/cache/openai-primary-runtime/documents/26.709.11516/skills/documents}"
PDFTOTEXT_BIN="${PDFTOTEXT_BIN:-/Users/Aamir/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/poppler/bin/pdftotext}"

REGISTRY="${PROJECT_ROOT}/results/pipeline/longcycler_release_claim_registry.json"
FINAL_MARKER="${PROJECT_ROOT}/results/pipeline/RUN_COMPLETE.txt"
FAILED_MARKER="${PROJECT_ROOT}/results/pipeline/RUN_FAILED.txt"
RQ_MARKER="${PROJECT_ROOT}/results/research_questions/RUN_COMPLETE.txt"
DELIVERABLE_VERIFIER="${PROJECT_ROOT}/scripts/verify_longcycler_release_deliverables.R"
RELEASE_DOC_WRITER="${PROJECT_ROOT}/scripts/generated/write_longcycler_release_docs.R"
DELIVERABLE_COMPLETE_MARKER="${PROJECT_ROOT}/results/pipeline/DELIVERABLES_COMPLETE.txt"
DELIVERABLE_FAILED_MARKER="${PROJECT_ROOT}/results/pipeline/DELIVERABLES_FAILED.txt"
DELIVERABLES_RESUME_FROM_ORDER="${DELIVERABLES_RESUME_FROM_ORDER:-1}"

V3_ROOT="outputs/manual-20260527-current-review/presentations/ruti-longitudinal-vf-pipeline-review-v3"
V4_ROOT="outputs/manual-20260527-current-review/presentations/ruti-longitudinal-vf-pipeline-review-v4"
V5_ROOT="outputs/manual-20260527-current-review/presentations/ruti-longitudinal-vf-pipeline-review-v5"

print_manifest() {
  cat <<'MANIFEST'
ORDER  KIND              GENERATOR / WRAPPER                                      CANONICAL OUTPUT
01     deck              scripts/generated/build_lecturer_methodology_deck.mjs    outputs/lecturer_methodology_pack/rUTI_complete_methodology_for_lecturer.pptx
02     deck (group B)    scripts/generated/build_longcycler_methods_summary_decks.mjs  outputs/longcycler_only_methods_summary/Longcycler_only_methods_summary.pptx
03     deck (group B)    scripts/generated/build_longcycler_methods_summary_decks.mjs  outputs/longcycler_only_methods_summary/Longcycler_only_analysis_flowchart.pptx
04     deck (group B)    scripts/generated/build_longcycler_methods_summary_decks.mjs  outputs/longcycler_only_methods_summary/Longcycler_only_methods_summary_with_flowchart.pptx
05     deck (group C)    scripts/generated/build_longcycler_core_review_decks.mjs outputs/manual-20260526-ruti-onboarding/presentations/ruti-clinical-genomic-onboarding/output/ruti-clinical-genomic-onboarding-review.pptx
06     deck (group C)    scripts/generated/build_longcycler_core_review_decks.mjs outputs/manual-20260527-current-review/presentations/ruti-current-results-scientific-review/output/rUTI_Current_Results_Scientific_Review_2026-05-27.pptx
07     deck (group C)    scripts/generated/build_longcycler_core_review_decks.mjs outputs/manual-20260527-current-review/presentations/ruti-current-results-scientific-review-v2/output/rUTI_Current_Results_VF_Focused_Review_2026-05-27.pptx
08     deck              v3/scripts/build_vf_pipeline_deck.mjs                    v3/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-27.pptx
09     deck              v4/scripts/build_onboarding_safe_deck.mjs                v4/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_With_Onboarding_2026-05-28.pptx
10     deck              v5/scripts/build_full_review_deck.mjs                    v5/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-28.pptx
11     deck              v5/scripts/build_editable_compact_deck.mjs               v5/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_Compact_Onboarding_2026-05-28.pptx
12     companion pack    lecturer generator                                       lecturer MD/CSV provenance and methods companions
13     companion pack    methods-summary generator                                handout, lay explanation, talking points, and count CSVs
14     presenter guide   v5/scripts/build_longcycler_presenter_guide.py --variant v3         v3 detailed guide DOCX + MD
15     presenter guide   v5/scripts/build_longcycler_presenter_guide.py --variant v4         v4 onboarding guide DOCX + MD
16     presenter guide   v5/scripts/build_longcycler_presenter_guide.py --variant v5-compact v5 compact guide DOCX + MD
17     guide render QA   documents/render_docx.py                                              external full-page PNGs for all three presenter guides
18     codebook          scripts/generated/build_vf_analysis_ready_codebook_docx.py          outputs/codebooks/vf_analysis_ready_lay_codebook.docx
19     codebook render   documents/render_docx.py --emit_pdf                                  outputs/codebooks/vf_analysis_ready_lay_codebook.pdf
20     release docs      scripts/generated/write_longcycler_release_docs.R                    11 registry-bound current Markdown documents
21     cleanup           exact obsolete-path manifest                                         repaired/backup/image-only decks and superseded presenter scripts
22     final gate        scripts/verify_longcycler_release_deliverables.R                     deliverable verification report
MANIFEST
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

should_run_order() {
  local order="$1"
  (( order >= DELIVERABLES_RESUME_FROM_ORDER ))
}

require_nonempty() {
  [[ -s "$1" ]] || die "Required non-empty file is missing: $1"
}

require_command() {
  command -v "$1" >/dev/null 2>&1 || die "Required command is unavailable: $1"
}

install_atomic() {
  local source="$1"
  local target="$2"
  local temp
  require_nonempty "${source}"
  mkdir -p "$(dirname "${target}")"
  temp="$(dirname "${target}")/.${target##*/}.$$.tmp"
  cp -- "${source}" "${temp}"
  mv -f -- "${temp}" "${target}"
}

sha256_file() {
  "${NODE_BIN}" --input-type=module -e \
    'import fs from "node:fs"; import crypto from "node:crypto"; process.stdout.write(crypto.createHash("sha256").update(fs.readFileSync(process.argv[1])).digest("hex"));' \
    "$1"
}

is_canonical_output() {
  case "$1" in
    "outputs/lecturer_methodology_pack/rUTI_complete_methodology_for_lecturer.pptx"|\
    "outputs/longcycler_only_methods_summary/Longcycler_only_methods_summary.pptx"|\
    "outputs/longcycler_only_methods_summary/Longcycler_only_analysis_flowchart.pptx"|\
    "outputs/longcycler_only_methods_summary/Longcycler_only_methods_summary_with_flowchart.pptx"|\
    "outputs/manual-20260526-ruti-onboarding/presentations/ruti-clinical-genomic-onboarding/output/ruti-clinical-genomic-onboarding-review.pptx"|\
    "outputs/manual-20260527-current-review/presentations/ruti-current-results-scientific-review/output/rUTI_Current_Results_Scientific_Review_2026-05-27.pptx"|\
    "outputs/manual-20260527-current-review/presentations/ruti-current-results-scientific-review-v2/output/rUTI_Current_Results_VF_Focused_Review_2026-05-27.pptx"|\
    "${V3_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-27.pptx"|\
    "${V4_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_With_Onboarding_2026-05-28.pptx"|\
    "${V5_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-28.pptx"|\
    "${V5_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_Compact_Onboarding_2026-05-28.pptx") return 0 ;;
    *) return 1 ;;
  esac
}

safe_remove_obsolete() {
  local relative="$1"
  local absolute
  [[ "${relative}" != /* ]] || die "Cleanup path must be project-relative: ${relative}"
  [[ "${relative}" != *".."* ]] || die "Cleanup path may not contain '..': ${relative}"
  is_canonical_output "${relative}" && die "Refusing to remove canonical output: ${relative}"
  absolute="${PROJECT_ROOT}/${relative}"
  if [[ -e "${absolute}" || -L "${absolute}" ]]; then
    rm -f -- "${absolute}"
    echo "Removed obsolete file: ${relative}"
  fi
}

safe_remove_obsolete_tree() {
  local relative="$1"
  local absolute
  case "${relative}" in
    "${V3_ROOT}/template-starter-preview"|"${V3_ROOT}/layout"|"${V3_ROOT}/template-inspect"|\
    "${V3_ROOT}/qa"|"${V3_ROOT}/slides"|"${V3_ROOT}/preview"|"${V3_ROOT}/template-starter-layout"|\
    "${V4_ROOT}/template-starter-preview"|"${V4_ROOT}/layout"|"${V4_ROOT}/template-inspect"|\
    "${V4_ROOT}/qa"|"${V4_ROOT}/slides"|"${V4_ROOT}/preview"|"${V4_ROOT}/template-starter-layout"|\
    "${V5_ROOT}/layout"|"${V5_ROOT}/qa"|"${V5_ROOT}/preview") ;;
    *) die "Refusing to remove unapproved generated tree: ${relative}" ;;
  esac
  absolute="${PROJECT_ROOT}/${relative}"
  if [[ -d "${absolute}" ]]; then
    rm -rf -- "${absolute}"
    echo "Removed obsolete generated tree: ${relative}"
  fi
}

# Invalidate any prior deliverable PASS before the first prerequisite check.
# Otherwise an early failure (for example a missing registry or runtime) could
# leave a stale completion marker that appears to bless the failed rerun.
mkdir -p "$(dirname "${DELIVERABLE_COMPLETE_MARKER}")"
rm -f -- "${DELIVERABLE_COMPLETE_MARKER}" "${DELIVERABLE_FAILED_MARKER}"
WORKSPACE_BASE=""
on_exit() {
  local status=$?
  if [[ ${status} -ne 0 ]]; then
    {
      echo "Longcycler-only deliverable regeneration: FAILED"
      echo "Ended: $(date)"
      echo "Scratch workspace: ${WORKSPACE_BASE:-not-created}"
    } > "${DELIVERABLE_FAILED_MARKER}"
    echo "Deliverable regeneration failed; scratch: ${WORKSPACE_BASE:-not-created}" >&2
  fi
}
trap on_exit EXIT

require_command "${NODE_BIN}"
require_command "${PYTHON_BIN}"
require_command "${RSCRIPT_BIN}"
require_command "${PDFTOTEXT_BIN}"
export PDFTOTEXT_BIN
[[ "${DELIVERABLES_RESUME_FROM_ORDER}" =~ ^[0-9]+$ ]] || \
  die "DELIVERABLES_RESUME_FROM_ORDER must be an integer from 1 to 22"
(( DELIVERABLES_RESUME_FROM_ORDER >= 1 && DELIVERABLES_RESUME_FROM_ORDER <= 22 )) || \
  die "DELIVERABLES_RESUME_FROM_ORDER must be between 1 and 22"
require_nonempty "${REGISTRY}"
require_nonempty "${FINAL_MARKER}"
require_nonempty "${RQ_MARKER}"
require_nonempty "${DELIVERABLE_VERIFIER}"
require_nonempty "${RELEASE_DOC_WRITER}"
[[ ! -e "${FAILED_MARKER}" ]] || die "Analysis failure marker still exists: ${FAILED_MARKER}"
grep -Fqx "Complete analysis: PASS" "${FINAL_MARKER}" || die "Final marker does not record PASS"
[[ ! "${REGISTRY}" -nt "${FINAL_MARKER}" ]] || die "Claim registry is newer than the final analysis marker"
REGISTRY_SHA256_START="$(sha256_file "${REGISTRY}")"
[[ "${REGISTRY_SHA256_START}" =~ ^[0-9a-f]{64}$ ]] || die "Could not hash the final claim registry"

SCRATCH_ROOT="${SCRATCH_ROOT:-$("${NODE_BIN}" -p "require('node:os').tmpdir()") }"
SCRATCH_ROOT="${SCRATCH_ROOT% }"
mkdir -p "${SCRATCH_ROOT}"
SCRATCH_ROOT="$(cd "${SCRATCH_ROOT}" && pwd -P)"
case "${SCRATCH_ROOT}/" in
  "${PROJECT_ROOT}/"*) die "Scratch root must be outside the repository: ${SCRATCH_ROOT}" ;;
esac
WORKSPACE_BASE="$(mktemp -d "${SCRATCH_ROOT}/ruti-longcycler-deliverables.XXXXXX")"
export PYTHONPYCACHEPREFIX="${WORKSPACE_BASE}/pycache"
export TMPDIR="${TMPDIR:-/private/tmp}"

echo "Ordered regeneration manifest:"
print_manifest
print_manifest > "${WORKSPACE_BASE}/ordered-command-manifest.txt"
echo "External scratch workspace: ${WORKSPACE_BASE}"
echo "Resume from manifest order: ${DELIVERABLES_RESUME_FROM_ORDER}"

# Static syntax and implementation-policy gates run before any final artifact is touched.
bash -n "${SCRIPT_PATH}"
node_sources=(
  scripts/generated/longcycler_release_presentation_common.mjs
  scripts/generated/build_lecturer_methodology_deck.mjs
  scripts/generated/build_longcycler_methods_summary_decks.mjs
  scripts/generated/build_longcycler_core_review_decks.mjs
  scripts/generated/create_lecturer_template_map.mjs
  scripts/generated/create_core_review_template_map.mjs
  "${V3_ROOT}/scripts/build_vf_pipeline_deck.mjs"
  "${V4_ROOT}/scripts/build_onboarding_safe_deck.mjs"
  "${V5_ROOT}/scripts/longcycler_review_deck_engine.mjs"
  "${V5_ROOT}/scripts/build_full_review_deck.mjs"
  "${V5_ROOT}/scripts/build_editable_compact_deck.mjs"
)
for source in "${node_sources[@]}"; do
  require_nonempty "${PROJECT_ROOT}/${source}"
  "${NODE_BIN}" --check "${PROJECT_ROOT}/${source}"
done

python_sources=(
  "${V5_ROOT}/scripts/build_longcycler_presenter_guide.py"
  scripts/generated/build_vf_analysis_ready_codebook_docx.py
)
for source in "${python_sources[@]}"; do
  require_nonempty "${PROJECT_ROOT}/${source}"
  "${PYTHON_BIN}" -m py_compile "${PROJECT_ROOT}/${source}"
done
for source in "${RELEASE_DOC_WRITER}" "${DELIVERABLE_VERIFIER}"; do
  "${RSCRIPT_BIN}" -e 'invisible(parse(file = commandArgs(trailingOnly = TRUE)[[1L]]))' "${source}"
done

if grep -Eiq 'pptxgenjs|python-pptx|from[[:space:]]+pptx|import[[:space:]]+pptx' "${node_sources[@]}" "${python_sources[@]}"; then
  die "A presentation generator references a forbidden PPTX implementation"
fi

# Reconfirm the final analytical release immediately before presentation work.
"${RSCRIPT_BIN}" scripts/verify_longcycler_only_pipeline.R --stage final

# Commands A-C generate decks 01-07 and their generator-owned companion files.
if should_run_order 1; then
  "${NODE_BIN}" scripts/generated/build_lecturer_methodology_deck.mjs \
    --project-root "${PROJECT_ROOT}" --registry "${REGISTRY}" \
    --workspace "${WORKSPACE_BASE}/01-lecturer"
fi

if should_run_order 2; then
  "${NODE_BIN}" scripts/generated/build_longcycler_methods_summary_decks.mjs \
    --project-root "${PROJECT_ROOT}" --registry "${REGISTRY}" \
    --workspace "${WORKSPACE_BASE}/02-methods-summary"
fi

if should_run_order 5; then
  "${NODE_BIN}" scripts/generated/build_longcycler_core_review_decks.mjs \
    --project-root "${PROJECT_ROOT}" --registry "${REGISTRY}" \
    --workspace "${WORKSPACE_BASE}/03-core-reviews"
fi

# Commands D-G generate decks 08-11 at exact canonical destinations.
if should_run_order 8; then
  "${NODE_BIN}" "${V3_ROOT}/scripts/build_vf_pipeline_deck.mjs" \
    --project-root "${PROJECT_ROOT}" --registry "${REGISTRY}" \
    --workspace "${WORKSPACE_BASE}/04-v3-review" \
    --output "${PROJECT_ROOT}/${V3_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-27.pptx"
fi

if should_run_order 9; then
  "${NODE_BIN}" "${V4_ROOT}/scripts/build_onboarding_safe_deck.mjs" \
    --project-root "${PROJECT_ROOT}" --registry "${REGISTRY}" \
    --workspace "${WORKSPACE_BASE}/05-v4-onboarding" \
    --output "${PROJECT_ROOT}/${V4_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_With_Onboarding_2026-05-28.pptx"
fi

if should_run_order 10; then
  "${NODE_BIN}" "${V5_ROOT}/scripts/build_full_review_deck.mjs" \
    --project-root "${PROJECT_ROOT}" --registry "${REGISTRY}" \
    --workspace "${WORKSPACE_BASE}/06-v5-full" \
    --output "${PROJECT_ROOT}/${V5_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-28.pptx"
fi

if should_run_order 11; then
  "${NODE_BIN}" "${V5_ROOT}/scripts/build_editable_compact_deck.mjs" \
    --project-root "${PROJECT_ROOT}" --registry "${REGISTRY}" \
    --workspace "${WORKSPACE_BASE}/07-v5-compact" \
    --output "${PROJECT_ROOT}/${V5_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_Compact_Onboarding_2026-05-28.pptx"
fi

# Companion presenter guides use the same registry and canonical deck variants.
GUIDE_GENERATOR="${V5_ROOT}/scripts/build_longcycler_presenter_guide.py"
V3_GUIDE_DOCX="${V3_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Detailed_Presenter_Script_2026-05-28.docx"
V4_GUIDE_DOCX="${V4_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Presenter_Guide_With_Onboarding_2026-05-28.docx"
V5_GUIDE_DOCX="${V5_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Presenter_Guide_Compact_Onboarding_2026-05-28.docx"
if should_run_order 14; then
  "${PYTHON_BIN}" "${GUIDE_GENERATOR}" --variant v3 --registry "${REGISTRY}" \
    --output-docx "${V3_GUIDE_DOCX}" \
    --output-md "${V3_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Detailed_Presenter_Script_2026-05-28.md"
  "${PYTHON_BIN}" "${GUIDE_GENERATOR}" --variant v4 --registry "${REGISTRY}" \
    --output-docx "${V4_GUIDE_DOCX}" \
    --output-md "${V4_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Presenter_Guide_With_Onboarding_2026-05-28.md"
  "${PYTHON_BIN}" "${GUIDE_GENERATOR}" --variant v5-compact --registry "${REGISTRY}" \
    --output-docx "${V5_GUIDE_DOCX}" \
    --output-md "${V5_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Presenter_Guide_Compact_Onboarding_2026-05-28.md"
fi

# Render every presenter guide into retained external scratch. A successful
# render is a structural gate; full-size page inspection remains mandatory.
if should_run_order 17; then
  guide_render_index=0
  for guide_docx in "${V3_GUIDE_DOCX}" "${V4_GUIDE_DOCX}" "${V5_GUIDE_DOCX}"; do
    guide_render_index=$((guide_render_index + 1))
    guide_render_dir="${WORKSPACE_BASE}/08-guide-${guide_render_index}"
    "${PYTHON_BIN}" "${DOC_SKILL_DIR}/render_docx.py" \
      "${PROJECT_ROOT}/${guide_docx}" --output_dir "${guide_render_dir}"
    require_nonempty "${guide_render_dir}/page-1.png"
  done
fi

# Build and externally render the codebook; install only its canonical root PDF.
CODEBOOK_RENDER_DIR="${WORKSPACE_BASE}/09-codebook-render"
if should_run_order 18; then
  "${PYTHON_BIN}" scripts/generated/build_vf_analysis_ready_codebook_docx.py
fi
if should_run_order 19; then
  "${PYTHON_BIN}" "${DOC_SKILL_DIR}/render_docx.py" \
    "${PROJECT_ROOT}/outputs/codebooks/vf_analysis_ready_lay_codebook.docx" \
    --output_dir "${CODEBOOK_RENDER_DIR}" --emit_pdf
  install_atomic "${CODEBOOK_RENDER_DIR}/vf_analysis_ready_lay_codebook.pdf" \
    "${PROJECT_ROOT}/outputs/codebooks/vf_analysis_ready_lay_codebook.pdf"
fi

# Rewrite all current prose documentation from the same audited registry.
if should_run_order 20; then
  "${RSCRIPT_BIN}" "${RELEASE_DOC_WRITER}"
fi

required_outputs=(
  outputs/lecturer_methodology_pack/rUTI_complete_methodology_for_lecturer.pptx
  outputs/lecturer_methodology_pack/layperson_to_technical_talking_points.md
  outputs/lecturer_methodology_pack/methodology_audit_findings.md
  outputs/lecturer_methodology_pack/rUTI_methods_section_audited.md
  outputs/lecturer_methodology_pack/numbered_R_script_methods_register.csv
  outputs/lecturer_methodology_pack/numbered_R_script_methods_register.md
  outputs/lecturer_methodology_pack/presentation_number_provenance.csv
  outputs/longcycler_only_methods_summary/Longcycler_only_methods_summary.pptx
  outputs/longcycler_only_methods_summary/Longcycler_only_analysis_flowchart.pptx
  outputs/longcycler_only_methods_summary/Longcycler_only_methods_summary_with_flowchart.pptx
  outputs/longcycler_only_methods_summary/Longcycler_only_methods_handout.md
  outputs/longcycler_only_methods_summary/Longcycler_only_methods_layperson_explanation.md
  outputs/longcycler_only_methods_summary/Longcycler_only_methods_talking_points.md
  outputs/longcycler_only_methods_summary/longcycler_only_methods_counts.csv
  outputs/longcycler_only_methods_summary/flowchart_counts.csv
  outputs/manual-20260526-ruti-onboarding/presentations/ruti-clinical-genomic-onboarding/output/ruti-clinical-genomic-onboarding-review.pptx
  outputs/manual-20260527-current-review/presentations/ruti-current-results-scientific-review/output/rUTI_Current_Results_Scientific_Review_2026-05-27.pptx
  outputs/manual-20260527-current-review/presentations/ruti-current-results-scientific-review-v2/output/rUTI_Current_Results_VF_Focused_Review_2026-05-27.pptx
  "${V3_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-27.pptx"
  "${V3_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Detailed_Presenter_Script_2026-05-28.docx"
  "${V3_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Detailed_Presenter_Script_2026-05-28.md"
  "${V4_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_With_Onboarding_2026-05-28.pptx"
  "${V4_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Presenter_Guide_With_Onboarding_2026-05-28.docx"
  "${V4_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Presenter_Guide_With_Onboarding_2026-05-28.md"
  "${V5_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-28.pptx"
  "${V5_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_Compact_Onboarding_2026-05-28.pptx"
  "${V5_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Presenter_Guide_Compact_Onboarding_2026-05-28.docx"
  "${V5_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Presenter_Guide_Compact_Onboarding_2026-05-28.md"
  outputs/codebooks/vf_analysis_ready_lay_codebook.docx
  outputs/codebooks/vf_analysis_ready_lay_codebook.pdf
  FOLDER_MAP.md
  CODE_REVIEW_RECONCILIATION_README.md
  docs/LECTURER_README.md
  docs/PIPELINE_FAILURE_LOG.md
  docs/VF_abstract_draft.md
  docs/VF_merge_diagnostics.md
  docs/VF_verification_report.md
  docs/figures/timepoint_vs_isolate_clarification.md
  docs/legacy_asb_uti_docs/VF_provenance_map.md
  docs/workflow_case_count_flowchart.md
  docs/workflow_flowchart.md
)
for output in "${required_outputs[@]}"; do
  require_nonempty "${PROJECT_ROOT}/${output}"
done

# Exact-only cleanup: no globs, directories, or canonical outputs are accepted.
obsolete_paths=(
  "${V3_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-27_EDITABLE_GENERATED_REPAIR_WARNING.pptx"
  "${V3_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-27_POWERPOINT_SAFE_IMAGE_DECK.pptx"
  "${V3_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-27_REPAIRED.pptx"
  "${V3_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-27_REPAIRED_IDS.pptx"
  "${V3_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-27_pre_id_repair_backup.pptx"
  "${V3_ROOT}/qa/roundtrip/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_roundtrip.pptx"
  "${V3_ROOT}/template-starter.pptx"
  "${V4_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-27.pptx"
  "${V4_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-27_EDITABLE_GENERATED_REPAIR_WARNING.pptx"
  "${V4_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-27_POWERPOINT_SAFE_IMAGE_DECK.pptx"
  "${V4_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-27_REPAIRED.pptx"
  "${V4_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-27_REPAIRED_IDS.pptx"
  "${V4_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_2026-05-27_pre_id_repair_backup.pptx"
  "${V4_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_With_Onboarding_BASE_SHIFTED_EDITABLE_GENERATED_REPAIR_WARNING.pptx"
  "${V4_ROOT}/qa/roundtrip/Longitudinal_Urinary_Ecoli_VF_Pipeline_Review_roundtrip.pptx"
  "${V4_ROOT}/template-starter.pptx"
  "${V3_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Presenter_Script_2026-05-27.docx"
  "${V3_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Presenter_Script_2026-05-27_REPAIRED.docx"
  "${V4_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Detailed_Presenter_Script_2026-05-28.docx"
  "${V4_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Detailed_Presenter_Script_2026-05-28.md"
  "${V4_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Detailed_Presenter_Script_With_Onboarding_2026-05-28.docx"
  "${V4_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Detailed_Presenter_Script_With_Onboarding_2026-05-28.md"
  "${V4_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Presenter_Script_2026-05-27.docx"
  "${V4_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Presenter_Script_2026-05-27_BEFORE_ONBOARDING_BACKUP.docx"
  "${V4_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Presenter_Script_2026-05-27_REPAIRED.docx"
  "${V4_ROOT}/output/Longitudinal_Urinary_Ecoli_VF_Pipeline_Presenter_Script_2026-05-28.docx"
  "${V4_ROOT}/output/~\$ngitudinal_Urinary_Ecoli_VF_Pipeline_Presenter_Script_2026-05-27.docx"
  "${V3_ROOT}/scripts/build_presenter_script.py"
  "${V3_ROOT}/scripts/build_detailed_presenter_script.py"
  "${V4_ROOT}/scripts/build_presenter_script.py"
  "${V4_ROOT}/scripts/build_detailed_presenter_script.py"
  "${V5_ROOT}/scripts/build_compact_presenter_guide.py"
)
if should_run_order 21; then
  for obsolete in "${obsolete_paths[@]}"; do
    safe_remove_obsolete "${obsolete}"
  done
fi

# These historical in-repository render/layout workspaces contain stale release
# claims. Current QA is regenerated in the external scratch workspace above;
# remove the old trees so the active outputs hierarchy is Longcycler-only.
obsolete_generated_trees=(
  "${V3_ROOT}/template-starter-preview"
  "${V3_ROOT}/layout"
  "${V3_ROOT}/template-inspect"
  "${V3_ROOT}/qa"
  "${V3_ROOT}/slides"
  "${V3_ROOT}/preview"
  "${V3_ROOT}/template-starter-layout"
  "${V4_ROOT}/template-starter-preview"
  "${V4_ROOT}/layout"
  "${V4_ROOT}/template-inspect"
  "${V4_ROOT}/qa"
  "${V4_ROOT}/slides"
  "${V4_ROOT}/preview"
  "${V4_ROOT}/template-starter-layout"
  "${V5_ROOT}/layout"
  "${V5_ROOT}/qa"
  "${V5_ROOT}/preview"
)
if should_run_order 21; then
  for obsolete_tree in "${obsolete_generated_trees[@]}"; do
    safe_remove_obsolete_tree "${obsolete_tree}"
  done
fi

# The deliverable verifier is deliberately the final command and must fail closed.
REGISTRY_SHA256_CURRENT="$(sha256_file "${REGISTRY}")"
[[ "${REGISTRY_SHA256_CURRENT}" == "${REGISTRY_SHA256_START}" ]] || \
  die "Claim-registry bytes changed during deliverable regeneration"
"${RSCRIPT_BIN}" "${DELIVERABLE_VERIFIER}"
REGISTRY_SHA256_END="$(sha256_file "${REGISTRY}")"
[[ "${REGISTRY_SHA256_END}" == "${REGISTRY_SHA256_START}" ]] || \
  die "Claim-registry bytes changed during final deliverable verification"
DELIVERABLE_COMPLETE_TMP="${DELIVERABLE_COMPLETE_MARKER}.$$.tmp"
{
  echo "Longcycler-only deliverable regeneration: PASS"
  echo "Ended: $(date)"
  echo "Claim registry SHA-256: ${REGISTRY_SHA256_START}"
  echo "Verification: results/qc/longcycler_release_deliverables_verification.csv"
  echo "Retained QA scratch: ${WORKSPACE_BASE}"
} > "${DELIVERABLE_COMPLETE_TMP}"
rm -f -- "${DELIVERABLE_FAILED_MARKER}"
mv -f -- "${DELIVERABLE_COMPLETE_TMP}" "${DELIVERABLE_COMPLETE_MARKER}"
trap - EXIT
echo "Longcycler-only deliverable regeneration: PASS"
echo "Retained scratch workspace: ${WORKSPACE_BASE}"
