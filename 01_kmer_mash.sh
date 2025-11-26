#!/usr/bin/env bash
# ==============================================================================
# 01_kmer_mash.sh
# ------------------------------------------------------------------------------
# Purpose: Performs k-mer based comparison of assemblies using Mash.
#          Generates distance matrices and identifies near-identical clusters.
# Inputs:  
#   - assemblies.list (default) or specified file
# Outputs: 
#   - results/wgs/kmer/ (asm.list, mash_all_vs_all.tab, mash_matrix.tsv, etc.)
# ==============================================================================

set -euo pipefail

# 1. Configuration
# ------------------------------------------------------------------------------
DIR_ROOT="${PWD}"
DIR_RESULTS="${DIR_ROOT}/results"
DIR_WGS="${DIR_RESULTS}/wgs"
DIR_KMER="${DIR_WGS}/kmer"
DEFAULT_ASM_LIST="${DIR_ROOT}/assemblies.list"

THREADS=8
ASM_LIST="$DEFAULT_ASM_LIST"

# 2. Parse Arguments
# ------------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --asm-list) ASM_LIST="$2"; shift 2;;
    --threads|-p) THREADS="$2"; shift 2;;
    *) echo "Usage: $0 [--asm-list FILE] [--threads N]" >&2; exit 2;;
  esac
done

# 3. Setup
# ------------------------------------------------------------------------------
mkdir -p "$DIR_KMER"

# 3. Auto-Generate Assembly List if Missing
# ------------------------------------------------------------------------------
mkdir -p "$DIR_KMER"

if [[ ! -s "$ASM_LIST" ]]; then
  echo "[info] Assembly list missing or empty: $ASM_LIST"
  echo "[info] Attempting to generate from metadata or directory..."

  METADATA_FILE="${DIR_ROOT}/assembly_metadata.csv"
  FASTAS_DIR="${DIR_ROOT}/ont-yellow-routine-fastas"

  if [[ -f "$METADATA_FILE" ]]; then
    echo "[info] Generating from $METADATA_FILE..."
    # Extract 'full_path' column (assuming it's the last column or named full_path)
    # We'll try to find the column index or just grep for fasta paths
    # Simple robust approach: look for .fasta in the csv
    grep -oE "[^,]+\.(fasta|fa|fna)(\.gz)?" "$METADATA_FILE" | while read -r fn; do
        # Check if it's a full path or relative
        if [[ -f "$fn" ]]; then
            echo "$fn"
        elif [[ -f "$DIR_ROOT/$fn" ]]; then
            echo "$DIR_ROOT/$fn"
        elif [[ -f "$FASTAS_DIR/$(basename "$fn")" ]]; then
            echo "$FASTAS_DIR/$(basename "$fn")"
        fi
    done | sort -u > "$ASM_LIST"
  
  elif [[ -d "$FASTAS_DIR" ]]; then
    echo "[info] Generating from directory $FASTAS_DIR..."
    find "$FASTAS_DIR" -type f -name "*.fasta" -o -name "*.fa" -o -name "*.fna" | sort > "$ASM_LIST"
  fi

  # Re-check
  if [[ ! -s "$ASM_LIST" ]]; then
     echo "[error] Could not generate assembly list. Please create $ASM_LIST manually." >&2
     exit 1
  fi
  echo "[info] Generated $ASM_LIST with $(wc -l < "$ASM_LIST" | xargs) entries."
fi

# 4. Build Absolute FASTA List
# ------------------------------------------------------------------------------
echo "[1/6] Building absolute FASTA list from: $ASM_LIST"
rm -f "$DIR_KMER/asm.list"
tmp_list="$(mktemp)"

# Trap to clean up temp file on exit
trap 'rm -f "$tmp_list"' EXIT

while IFS= read -r P || [[ -n "$P" ]]; do
  [[ -z "$P" || "$P" =~ ^# ]] && continue
  
  if [[ -d "$P" ]]; then
    # Find FASTAs in directory
    find "$P" -type f -size +0 \( -iname "*.fa" -o -iname "*.fna" -o -iname "*.fasta" -o -iname "*.fa.gz" -o -iname "*.fna.gz" -o -iname "*.fasta.gz" \) >> "$tmp_list"
  elif [[ -f "$P" ]]; then
    # Check extension if it's a file
    if echo "$P" | grep -Eiq '\.(fa|fna|fasta)(\.gz)?$'; then
      echo "$P" >> "$tmp_list"
    else
      echo "[warn] Skipping non-FASTA file: $P" >&2
    fi
  else
    echo "[warn] Path not found: $P" >&2
  fi
done < "$ASM_LIST"

# Canonicalize paths and unique sort
awk 'NF' "$tmp_list" | while read -r f; do
  d=$(cd "$(dirname "$f")" 2>/dev/null && pwd)
  [[ -n "$d" ]] && echo "$d/$(basename "$f")"
done | sort -u > "$DIR_KMER/asm.list"

# Validate count
N=$(wc -l < "$DIR_KMER/asm.list" | xargs || echo 0)
echo "[info] Found $N unique FASTA files"
if (( N < 2 )); then
  echo "[error] Need at least 2 FASTAs to perform comparison."
  exit 1
fi

# 5. Check Prerequisites
# ------------------------------------------------------------------------------
echo "[2/6] Checking for Mash..."
if ! command -v mash >/dev/null 2>&1; then
  echo "[info] Mash not found. Attempting installation via conda..."
  conda install -y -c bioconda -c conda-forge mash || { echo "[error] Failed to install mash."; exit 1; }
fi

# 6. Run Mash
# ------------------------------------------------------------------------------
echo "[3/6] Running Mash sketch..."
# Clean old outputs
rm -f "$DIR_KMER/asm.msh" "$DIR_KMER/mash_all_vs_all.tab" "$DIR_KMER/mash_matrix.tsv" \
      "$DIR_KMER/labels.tsv" "$DIR_KMER/duplicates.tsv" "$DIR_KMER/near_identical_clusters.tsv"

mash sketch -p "$THREADS" -o "$DIR_KMER/asm" -l "$DIR_KMER/asm.list"

echo "[4/6] Running Mash dist (all-vs-all)..."
mash dist -p "$THREADS" "$DIR_KMER/asm.msh" "$DIR_KMER/asm.msh" > "$DIR_KMER/mash_all_vs_all.tab"

# 7. Post-Processing (Python)
# ------------------------------------------------------------------------------
echo "[5/6] Generating matrix and reports..."
python3 - <<PY
import os, sys, collections

# Config from bash
base_dir = "$DIR_KMER"
asm_list_path = os.path.join(base_dir, "asm.list")
pairs_tab_path = os.path.join(base_dir, "mash_all_vs_all.tab")

# 5a) Read paths and detect duplicates
try:
    with open(asm_list_path) as f:
        paths = [ln.strip() for ln in f if ln.strip()]
except FileNotFoundError:
    print(f"[error] Could not read {asm_list_path}")
    sys.exit(1)

basenames = [os.path.basename(p) for p in paths]
parents   = [os.path.basename(os.path.dirname(p)) for p in paths]
cnt = collections.Counter(basenames)

# Write duplicates report
with open(os.path.join(base_dir, "duplicates.tsv"), "w") as fo:
    fo.write("basename\tcount\tpaths\n")
    for b, c in cnt.items():
        if c > 1:
            ps = [p for p in paths if os.path.basename(p) == b]
            fo.write(f"{b}\t{c}\t{';'.join(ps)}\n")

# 5b) Create unique labels
used = set()
labels = {}
for p, b, pa in zip(paths, basenames, parents):
    lab = b if cnt[b] == 1 else f"{pa}__{b}"
    # De-dup if still clashing
    orig = lab
    i = 2
    while lab in used:
        lab = f"{orig}~{i}"
        i += 1
    labels[p] = lab
    used.add(lab)

with open(os.path.join(base_dir, "labels.tsv"), "w") as fo:
    fo.write("label\tpath\n")
    for p in paths:
        fo.write(f"{labels[p]}\t{p}\n")

# 5c) Build Matrix
pairs = {}
names = set()

try:
    with open(pairs_tab_path) as fh:
        for ln in fh:
            s = ln.strip().split("\t")
            if len(s) < 3: continue
            q, r, d = s[0], s[1], float(s[2])
            if q not in labels or r not in labels:
                continue
            ql, rl = labels[q], labels[r]
            names.update([ql, rl])
            pairs[(ql, rl)] = d
            pairs[(rl, ql)] = d
except FileNotFoundError:
    print(f"[error] Could not read {pairs_tab_path}")
    sys.exit(1)

names = sorted(names)
mat_path = os.path.join(base_dir, "mash_matrix.tsv")
with open(mat_path, "w") as fo:
    fo.write("\t" + "\t".join(names) + "\n")
    for a in names:
        row = []
        for b in names:
            if a == b:
                row.append("0")
            else:
                row.append(str(pairs.get((a, b), "NA")))
        fo.write(a + "\t" + "\t".join(row) + "\n")

# 5d) Near-identical clusters (d <= 0.001)
th = 0.001
from collections import defaultdict, deque
G = defaultdict(set)
for (a, b), d in pairs.items():
    if a != b and d != "NA" and float(d) <= th:
        G[a].add(b)
        G[b].add(a)

seen = set()
clusters = []
for n in sorted(G):
    if n in seen: continue
    q = deque([n])
    seen.add(n)
    comp = []
    while q:
        x = q.popleft()
        comp.append(x)
        for y in G[x]:
            if y not in seen:
                seen.add(y)
                q.append(y)
    clusters.append(sorted(comp))

with open(os.path.join(base_dir, "near_identical_clusters.tsv"), "w") as fo:
    fo.write("# cluster_id\tcount\tlabels\n")
    for i, c in enumerate(sorted(clusters, key=len, reverse=True), 1):
        fo.write(f"group{i}\t{len(c)}\t{','.join(c)}\n")

print(f"[ok] Generated matrix: {mat_path}")
PY

echo "[6/6] Done. Results in $DIR_KMER"
