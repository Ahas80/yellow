#!/usr/bin/env bash
set -euo pipefail
b="results/wgs"
q="results/qc"

echo "== Canonical WGS artifact spot-checks =="
check() { compgen -G "$1" >/dev/null && echo "[OK] $2" || echo "[MISS] $2"; }
check "$b/qc_summary.csv"                         "12a: WGS QC summary"
check "$q/canonical_assembly_selection.csv"        "12a: canonical assembly selection"
check "$b/core/parsnp_out/parsnp.xmfa"             "12b: Parsnp XMFA"
check "$b/core/parsnp_out/parsnp.fasta"            "12b: Parsnp FASTA alignment"
check "$b/pan/gene_presence_absence.csv"           "12c: Panaroo gene_presence_absence"
check "$b/pan/panaroo_input_manifest.csv"          "12c: Panaroo input manifest"
check "$b/pan/panaroo_staleness_report.txt"        "12c: Panaroo staleness report"
check "results/vf/vf_pa_all.csv"                   "02: VF presence/absence matrix"
check "results/mlst/mlst_with_meta.csv"            "06: local MLST provenance table"
check "results/mlst/mlst_provider_preferred.csv"   "MLST: active provider-preferred episode table"
check "results/mlst/mlst_provider_source_audit.csv" "MLST: provider source audit"

echo
echo "== Canonical selection summary =="
if [[ -f "$q/canonical_assembly_selection.csv" ]]; then
  awk -F, '
    NR==1 {
      for (i=1; i<=NF; i++) h[$i]=i
      next
    }
    {
      sel=tolower($(h["selected_canonical"]))
      qc=tolower($(h["QC_PASS"]))
      if (sel=="true" || sel=="1") selected++
      if (qc=="true" || qc=="1") qcp++
      total++
    }
    END { printf "rows=%d QC_PASS=%d selected_canonical=%d\n", total, qcp, selected }
  ' "$q/canonical_assembly_selection.csv"
else
  echo "No canonical assembly selection file."
fi

echo
echo "== Panaroo status =="
if [[ -f "$b/pan/panaroo_staleness_report.txt" ]]; then
  grep '^Status:' "$b/pan/panaroo_staleness_report.txt" || true
else
  echo "No Panaroo staleness report."
fi

echo
echo "== Last 3 WGS logs (tail 20) =="
recent_logs="$(find "$b" -path '*/logs/*.log' -type f 2>/dev/null | sort -r | sed -n '1,3p' || true)"
printf '%s\n' "$recent_logs" | while read -r f; do
  [[ -n "$f" ]] || continue
  echo "--- $f ---"; tail -n 20 "$f"; echo
done
