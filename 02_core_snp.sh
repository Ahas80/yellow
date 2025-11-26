#!/usr/bin/env bash
# ==============================================================================
# 02_core_snp.sh
# ------------------------------------------------------------------------------
# Purpose: Performs core genome SNP analysis using Parsnp.
#          Generates core alignment, SNP distance matrix, and clusters.
# Inputs:
#   - results/wgs/kmer/asm.list (from step 01)
# Outputs:
#   - results/wgs/parsnp/ (core.aln.fasta, snp-dists.tsv, etc.)
# ==============================================================================

set -euo pipefail

# 1. Configuration
# ------------------------------------------------------------------------------
DIR_ROOT="${PWD}"
DIR_RESULTS="${DIR_ROOT}/results"
DIR_WGS="${DIR_RESULTS}/wgs"
DIR_PARSNP="${DIR_WGS}/parsnp"
DIR_INPUT="${DIR_PARSNP}/input"
DEFAULT_ASM_LIST="${DIR_WGS}/kmer/asm.list"

THREADS=8
REF="!"     # "!" means auto-select reference
FORCE=0
MAX_N=0     # 0 = all
SNP_TH=10   # Threshold for near-identical clusters

ASM_LIST="$DEFAULT_ASM_LIST"

# 2. Parse Arguments
# ------------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --asm-list) ASM_LIST="$2"; shift 2;;
    --threads|-p) THREADS="$2"; shift 2;;
    --ref) REF="$2"; shift 2;;
    --force) FORCE=1; shift;;
    --max-n) MAX_N="$2"; shift 2;;
    --snp-threshold) SNP_TH="$2"; shift 2;;
    *) echo "Usage: $0 [--asm-list FILE] [--threads N] [--ref REF.fa] [--force] [--max-n N] [--snp-threshold N]" >&2; exit 2;;
  esac
done

# 3. Setup Output Directories
# ------------------------------------------------------------------------------
mkdir -p "$DIR_PARSNP"
mkdir -p "$DIR_INPUT"

if [[ ! -s "$ASM_LIST" ]]; then
  echo "[error] Assembly list missing or empty: $ASM_LIST" >&2
  exit 1
fi

if [[ $FORCE -eq 1 ]]; then
  echo "[info] Force mode: cleaning previous results..."
  rm -rf "$DIR_INPUT" "$DIR_PARSNP"/parsnp.* "$DIR_PARSNP"/core.aln.fasta "$DIR_PARSNP"/snp-dists.tsv "$DIR_PARSNP"/snp_clusters.tsv
  mkdir -p "$DIR_INPUT"
fi

# 4. Prepare Input Directory (Symlinks)
# ------------------------------------------------------------------------------
echo "[1/7] Linking inputs to: $DIR_INPUT"
# Clean existing links to ensure fresh state
find "$DIR_INPUT" -type l -delete

count=0
while IFS= read -r f || [[ -n "$f" ]]; do
  [[ -z "$f" ]] && continue
  if [[ ! -f "$f" ]]; then
    echo "[warn] Missing file: $f" >&2
    continue
  fi
  
  ((count++))
  # Use 5-digit prefix to preserve order if needed, though parsnp reorders
  ln -sf "$f" "$DIR_INPUT/$(printf "%05d" "$count")__$(basename "$f")"
  
  if [[ "$MAX_N" -gt 0 && "$count" -ge "$MAX_N" ]]; then
    echo "[info] Reached limit --max-n=$MAX_N"
    break
  fi
done < "$ASM_LIST"

echo "[info] Linked $count FASTAs"

# 5. Rosetta / Conda Environment Setup (Apple Silicon Support)
# ------------------------------------------------------------------------------
echo "[2/7] Checking Rosetta..."
if [[ "$(uname -m)" == "arm64" ]]; then
    softwareupdate --install-rosetta --agree-to-license >/dev/null 2>&1 || true
fi

echo "[3/7] Running Parsnp in x86 environment..."
# We execute the heavy lifting inside a subshell with the correct architecture
arch -x86_64 /bin/bash -lc '
  set -e
  
  # Source conda
  # Try standard locations
  if [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
  elif [[ -f "$HOME/opt/miniconda3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/opt/miniconda3/etc/profile.d/conda.sh"
  elif [[ -f "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh" ]]; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
  else
    echo "[error] Could not find conda.sh"
    exit 1
  fi

  # Create env if missing
  if ! conda env list | grep -q "^asm-snp-x86 "; then
    echo "[info] Creating conda env: asm-snp-x86"
    CONDA_SUBDIR=osx-64 conda create -y -n asm-snp-x86 -c bioconda -c conda-forge parsnp harvesttools snp-dists
  fi
  
  conda activate asm-snp-x86
  conda config --env --set subdir osx-64

  PAR_OUT="'"$DIR_PARSNP"'"
  INP="'"$DIR_INPUT"'"
  THREADS="'"$THREADS"'"
  REF="'"$REF"'"

  echo "[4/7] Executing Parsnp..."
  # Note: parsnp outputs to current dir if not careful, but -o handles it.
  # We explicitly cd to output dir to capture any stray files, then move back?
  # Better: trust -o but monitor.
  
  if [[ "$REF" == "!" ]]; then
    parsnp -r "!" -d "$INP" -o "$PAR_OUT" -p "$THREADS"
  else
    parsnp -r "$REF" -d "$INP" -o "$PAR_OUT" -p "$THREADS"
  fi

  echo "[5/7] Extracting core alignment..."
  if [[ -f "$PAR_OUT/parsnp.xmfa" ]]; then
      # Try harvesttools first
      set +e
      harvesttools -x -i "$PAR_OUT/parsnp.xmfa" -M "$PAR_OUT/core.aln.fasta"
      rc=$?
      set -e
      
      if [[ $rc -ne 0 || ! -s "$PAR_OUT/core.aln.fasta" ]]; then
          echo "[warn] harvesttools failed (exit code $rc). Attempting Python fallback..."
          
          python3 - <<PY
import sys, re, os

# Bash variables injected via environment or just hardcoded if we use HEREDOC with variable expansion
# But we want to avoid expansion of python variables.
# So we use quoted HEREDOC 'PY' and pass variables via environment or arguments.
# Or we use unquoted HEREDOC and escape python variables.
# Simplest: Pass variables as arguments.

xmfa_path = "$PAR_OUT/parsnp.xmfa"
out_path  = "$PAR_OUT/core.aln.fasta"

def parse_xmfa(path, out):
    seqs = {}
    current_header = None
    
    try:
        with open(path) as f:
            blocks = {} # idx -> list of seq chunks
            
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    # > 1:243-599 + ...
                    m = re.match(r"> ?(\d+):", line)
                    if m:
                        idx = int(m.group(1))
                        current_header = idx
                        if idx not in blocks: blocks[idx] = []
                    else:
                        current_header = None
                elif line.startswith("="):
                    current_header = None
                else:
                    if current_header is not None:
                        blocks[current_header].append(line)
                        
        # Concatenate
        final_seqs = {}
        for idx, chunks in blocks.items():
            final_seqs[idx] = "".join(chunks)
            
        with open(out, "w") as fo:
            for idx in sorted(final_seqs.keys()):
                fo.write(f">genome_{idx}\n{final_seqs[idx]}\n")
                
        print(f"[ok] Extracted {len(final_seqs)} sequences via Python.")
        
    except Exception as e:
        print(f"[error] Python fallback failed: {e}")
        sys.exit(1)

parse_xmfa(xmfa_path, out_path)
PY
      fi
  else
      echo "[error] parsnp.xmfa not found. Parsnp failed?"
      exit 1
  fi

  echo "[6/7] Calculating SNP distance matrix..."
  snp-dists "$PAR_OUT/core.aln.fasta" > "$PAR_OUT/snp-dists.tsv"
'

# 6. Post-Processing: Clustering
# ------------------------------------------------------------------------------
echo "[7/7] Clustering near-identical isolates (threshold=$SNP_TH SNPs)..."
python3 - <<PY
import os, sys, collections
from collections import defaultdict, deque

par_dir = "$DIR_PARSNP"
tab_path = os.path.join(par_dir, "snp-dists.tsv")
out_path = os.path.join(par_dir, "snp_clusters.tsv")
threshold = int("$SNP_TH")

if not os.path.exists(tab_path):
    print(f"[error] {tab_path} not found")
    sys.exit(1)

# Read Header
with open(tab_path) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    labels = hdr[1:]

# Parse Matrix
D = defaultdict(dict)
with open(tab_path) as fh:
    next(fh) # skip header
    for ln in fh:
        s = ln.rstrip("\n").split("\t")
        a = s[0]
        for b, val in zip(labels, s[1:]):
            if a == b: continue
            if val.isdigit():
                D[a][b] = int(val)

# Cluster
G = defaultdict(set)
for a in labels:
    for b, v in D[a].items():
        if v <= threshold:
            G[a].add(b)
            G[b].add(a)

seen = set()
clusters = []
for n in labels:
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
    if len(comp) > 1:
        clusters.append(sorted(comp))

# Write Output
with open(out_path, "w") as fo:
    fo.write("# cluster_id\tcount\tlabels\n")
    for i, c in enumerate(sorted(clusters, key=len, reverse=True), 1):
        fo.write(f"group{i}\t{len(c)}\t{','.join(c)}\n")

print(f"[ok] Wrote clusters to {out_path}")
PY

echo "Done. Outputs in $DIR_PARSNP"
