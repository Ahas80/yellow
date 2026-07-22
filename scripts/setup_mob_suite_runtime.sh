#!/usr/bin/env bash
set -euo pipefail

# Reproducible local runtime for 09b_mob_plasmid_reconstruction.R. The runtime
# and downloaded databases live under ignored data/ and are never committed.

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUNTIME_ROOT="${MOB_RUNTIME_ROOT:-${ROOT}/data/mob_suite_runtime}"
PREFIX="${MOB_RUNTIME_PREFIX:-${RUNTIME_ROOT}/env}"
MAMBA_BIN="${MAMBA_BIN:-/Users/Aamir/miniforge_x86/bin/mamba}"
MOB_SUITE_VERSION="${MOB_SUITE_VERSION:-3.1.9}"
BLAST_VERSION="${MOB_BLAST_VERSION:-2.15.0}"
MASH_VERSION="${MOB_MASH_VERSION:-2.3}"
ZENODO_RECORD="${MOB_DATABASE_ZENODO_RECORD:-10304948}"

EXPECTED_NAMES=(
  "ncbi_plasmid_full_seqs.fas"
  "clusters.txt"
  "rep.dna.fas"
  "mob.proteins.faa"
  "mpf.proteins.faa"
  "orit.fas"
  "repetitive.dna.fas"
)

expected_hash() {
  case "$1" in
    "ncbi_plasmid_full_seqs.fas")
      echo "939ca451e97f9a5e6ea3e38286d5d11c620aafdefe693705c2611f4ed98355d5"
      ;;
    "clusters.txt")
      echo "2811a3f1af5a632d985f6e64b6c3be46ec249acc4bf8f76b0bc3bdf569937e11"
      ;;
    "rep.dna.fas")
      echo "fdc10d866d8fdc3c75db0f945c8b52797392abb7b7bf389c5de36a2429500eb0"
      ;;
    "mob.proteins.faa")
      echo "e5db0e07f9e94f7252f8adceaba4355851331b1ea4e60f51ce90af97490c47a9"
      ;;
    "mpf.proteins.faa")
      echo "fc6a9f78465271826120659e46c6876aa1b0f17baa0b806a525104b2e41f12fa"
      ;;
    "orit.fas")
      echo "821b2d39ff960f9c05dd467a3f9694d108b2e45c368923d21dfbf9a7383b921d"
      ;;
    "repetitive.dna.fas")
      echo "bf039b8cf672ffed353ea35d8088e9567ea9d16191003bb10a5ec5f2157d3bf7"
      ;;
    *)
      echo "No expected MOB database hash is registered for $1" >&2
      return 1
      ;;
  esac
}

if [[ ! -x "${MAMBA_BIN}" ]]; then
  echo "mamba was not found at ${MAMBA_BIN}" >&2
  exit 1
fi

mkdir -p "${RUNTIME_ROOT}"
mkdir -p "${RUNTIME_ROOT}/pkgs" "${RUNTIME_ROOT}/cache"
export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-${RUNTIME_ROOT}/pkgs}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${RUNTIME_ROOT}/cache}"

runtime_ready=false
if [[ -x "${PREFIX}/bin/mob_recon" &&
      -x "${PREFIX}/bin/blastn" &&
      -x "${PREFIX}/bin/mash" ]]; then
  mob_version="$("${PREFIX}/bin/mob_recon" --version 2>&1 || true)"
  blast_version="$("${PREFIX}/bin/blastn" -version 2>&1 | head -n 1 || true)"
  mash_version="$("${PREFIX}/bin/mash" --version 2>&1 || true)"
  if [[ "${mob_version}" == *"${MOB_SUITE_VERSION}"* &&
        "${blast_version}" == *"${BLAST_VERSION}"* &&
        "${mash_version}" == "${MASH_VERSION}" ]]; then
    runtime_ready=true
  fi
fi

if [[ "${runtime_ready}" != true ]]; then
  "${MAMBA_BIN}" create -y -p "${PREFIX}" \
    --override-channels -c conda-forge -c bioconda --strict-channel-priority \
    "mob_suite=${MOB_SUITE_VERSION}" \
    "blast=${BLAST_VERSION}" \
    "mash=${MASH_VERSION}"
fi

db_dir="$("${PREFIX}/bin/python" -c \
  'from pathlib import Path; import mob_suite; print(Path(mob_suite.__file__).resolve().parent / "databases")')"
if [[ ! -d "${db_dir}" ]]; then
  echo "MOB-suite database directory is missing: ${db_dir}" >&2
  exit 1
fi

hash_file() {
  shasum -a 256 "$1" | awk '{print $1}'
}

for name in "${EXPECTED_NAMES[@]}"; do
  path="${db_dir}/${name}"
  if [[ ! -f "${path}" ]]; then
    echo "Required MOB-suite database file is missing: ${path}" >&2
    exit 1
  fi
  observed="$(hash_file "${path}")"
  expected="$(expected_hash "${name}")"
  if [[ "${observed}" != "${expected}" ]]; then
    echo "MOB-suite database hash mismatch for ${name}" >&2
    echo "observed=${observed}" >&2
    echo "expected=${expected}" >&2
    exit 1
  fi
done

{
  echo "MOB-suite runtime installed"
  echo "prefix=${PREFIX}"
  echo "mob_suite_version=$("${PREFIX}/bin/mob_recon" --version 2>&1)"
  echo "blast_version=$("${PREFIX}/bin/blastn" -version 2>&1 | head -n 1)"
  echo "mash_version=$("${PREFIX}/bin/mash" --version 2>&1)"
  echo "database_zenodo_record=${ZENODO_RECORD}"
  echo "database_dir=${db_dir}"
  for name in $(printf '%s\n' "${EXPECTED_NAMES[@]}" | sort); do
    echo "sha256_${name}=$(expected_hash "${name}")"
  done
} > "${RUNTIME_ROOT}/INSTALL_COMPLETE.txt"

echo "Pinned MOB-suite runtime is ready at ${PREFIX}"
