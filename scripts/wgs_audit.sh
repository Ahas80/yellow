#!/usr/bin/env bash
set -euo pipefail
b="results/wgs"; p="$b/reports/progress.csv"

echo "== Progress heartbeat summary (stage,status,count) =="
if [[ -f "$p" ]]; then
  awk -F, 'NR>1{g[$1","$3]++} END{for(k in g) print k","g[k]}' "$p" | sort -t, -k1,1 -k2,2
else
  echo "No progress.csv yet ($p)."
fi

echo; echo "== Quick artifact spot-checks =="
check() { compgen -G "$1" >/dev/null && echo "[OK] $2" || echo "[MISS] $2"; }
check "$b/core/joint_*.vcf.gz"      "A: joint VCF"
check "$b/kmer/mash_all_vs_all.tab" "B: Mash all-vs-all"
check "$b/pangenome/pan_jaccard_pairs.csv" "C: Pan-Jaccard pairs"
check "$b/sv/sv_pairs.csv"          "D: SV pairs"
check "$b/vaf/het_snp_summary.csv"  "E: VAF het summary"
check "$b/plasmids/plasmid_clusters.csv" "F: Plasmid clusters"
echo; echo "== Last 3 logs (tail 20) =="
ls -1t $b/logs/*.log 2>/dev/null | head -3 | while read -r f; do
  echo "--- $f ---"; tail -n 20 "$f"; echo
done
