#!/usr/bin/env bash
set -euo pipefail

# Reproducible local runtime for numbered script 29. Database repositories are
# cloned once and then left at their recorded commits; this script never pulls
# an existing database checkout.

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUNTIME_ROOT="${AMR_RUNTIME_ROOT:-${ROOT}/data/amr_runtime}"
PREFIX="${AMR_RUNTIME_PREFIX:-${RUNTIME_ROOT}/env}"
DB_ROOT="${AMR_DATABASE_ROOT:-${RUNTIME_ROOT}/databases}"
MAMBA_BIN="${MAMBA_BIN:-/Users/Aamir/miniforge_x86/bin/mamba}"
AMRFINDER_VERSION="${AMRFINDER_VERSION:-4.2.7}"
AMRFINDER_DB_VERSION="${AMRFINDER_DB_VERSION:-2026-05-15.1}"
RESFINDER_VERSION="${RESFINDER_VERSION:-4.7.2}"
RESFINDER_DB_COMMIT="${RESFINDER_DB_COMMIT:-eecf0aa207594fe6d51badf808473de62b28cb06}"
POINTFINDER_DB_COMMIT="${POINTFINDER_DB_COMMIT:-44ce624a806c6d2b70f7e39841a5f9cb4d9010aa}"
AMRFINDER_DB_ROOT="${AMRFINDER_DB_ROOT:-${DB_ROOT}/amrfinderplus}"
AMRFINDER_DB="${AMRFINDER_DB:-${AMRFINDER_DB_ROOT}/${AMRFINDER_DB_VERSION}}"

if [[ ! -x "${MAMBA_BIN}" ]]; then
  echo "mamba was not found at ${MAMBA_BIN}" >&2
  exit 1
fi

mkdir -p "${RUNTIME_ROOT}" "${DB_ROOT}"

runtime_ready=false
if [[ -x "${PREFIX}/bin/amrfinder" && -x "${PREFIX}/bin/python" && -x "${PREFIX}/bin/kma" ]]; then
  amrfinder_version="$("${PREFIX}/bin/amrfinder" --version 2>/dev/null || true)"
  resfinder_version="$("${PREFIX}/bin/python" -c "import importlib.metadata as m; print(m.version('resfinder'))" 2>/dev/null || true)"
  if [[ "${amrfinder_version}" == "${AMRFINDER_VERSION}" &&
        "${resfinder_version}" == "${RESFINDER_VERSION}" ]]; then
    runtime_ready=true
  fi
fi

if [[ "${runtime_ready}" != true ]]; then
  "${MAMBA_BIN}" create -y -p "${PREFIX}" \
    -c conda-forge -c bioconda \
    "ncbi-amrfinderplus=${AMRFINDER_VERSION}" \
    "resfinder=${RESFINDER_VERSION}" \
    kma \
    blast
fi

if [[ ! -d "${DB_ROOT}/resfinder_db/.git" ]]; then
  git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git \
    "${DB_ROOT}/resfinder_db"
fi

if [[ ! -d "${DB_ROOT}/pointfinder_db/.git" ]]; then
  git clone https://bitbucket.org/genomicepidemiology/pointfinder_db.git \
    "${DB_ROOT}/pointfinder_db"
fi

for db in resfinder_db pointfinder_db; do
  if [[ -n "$(git -C "${DB_ROOT}/${db}" status --porcelain)" ]]; then
    echo "${db} contains local changes; refusing to alter the pinned database checkout." >&2
    exit 1
  fi
done

git -C "${DB_ROOT}/resfinder_db" checkout --detach "${RESFINDER_DB_COMMIT}"
git -C "${DB_ROOT}/pointfinder_db" checkout --detach "${POINTFINDER_DB_COMMIT}"

if [[ -f "${DB_ROOT}/resfinder_db/INSTALL.py" &&
      ! -f "${DB_ROOT}/resfinder_db/all.comp.b" ]]; then
  (
    cd "${DB_ROOT}/resfinder_db"
    "${PREFIX}/bin/python" INSTALL.py "${PREFIX}/bin/kma" non_interactive
  )
fi

if [[ -f "${DB_ROOT}/pointfinder_db/INSTALL.py" &&
      ! -f "${DB_ROOT}/pointfinder_db/escherichia_coli/escherichia_coli.comp.b" ]]; then
  (
    cd "${DB_ROOT}/pointfinder_db"
    "${PREFIX}/bin/python" INSTALL.py "${PREFIX}/bin/kma" non_interactive
  )
fi

if [[ ! -f "${AMRFINDER_DB}/database_format_version.txt" ]]; then
  mkdir -p "${AMRFINDER_DB_ROOT}"
  source_root="${AMRFINDER_DB_SOURCE_ROOT:-}"
  if [[ -z "${source_root}" ]]; then
    for candidate in \
      "${PREFIX}/share/amrfinderplus/data" \
      "${PREFIX}/bin/data" \
      "/Users/Aamir/miniforge_x86/share/amrfinderplus/data"; do
      if [[ -f "${candidate}/${AMRFINDER_DB_VERSION}/database_format_version.txt" ]]; then
        source_root="${candidate}"
        break
      fi
    done
  fi
  if [[ -n "${source_root}" &&
        -f "${source_root}/${AMRFINDER_DB_VERSION}/database_format_version.txt" ]]; then
    cp -R "${source_root}/${AMRFINDER_DB_VERSION}" "${AMRFINDER_DB_ROOT}/"
  else
    "${PREFIX}/bin/amrfinder_update" -d "${AMRFINDER_DB_ROOT}"
  fi
fi

if [[ ! -f "${AMRFINDER_DB}/database_format_version.txt" ]]; then
  echo "Pinned AMRFinderPlus database ${AMRFINDER_DB_VERSION} is unavailable." >&2
  echo "Set AMRFINDER_DB_SOURCE_ROOT to a directory containing that version." >&2
  exit 1
fi

{
  echo "AMR runtime installed"
  echo "prefix=${PREFIX}"
  echo "amrfinder_version=$("${PREFIX}/bin/amrfinder" --version 2>&1)"
  echo "amrfinder_db=${AMRFINDER_DB}"
  echo "amrfinder_db_version=${AMRFINDER_DB_VERSION}"
  echo "amrfinder_db_report=$("${PREFIX}/bin/amrfinder" -d "${AMRFINDER_DB}" -V 2>&1 | tr '\n' ' ')"
  echo "resfinder_version=$("${PREFIX}/bin/python" -c "import importlib.metadata as m; print(m.version('resfinder'))")"
  echo "resfinder_db_commit=$(git -C "${DB_ROOT}/resfinder_db" rev-parse HEAD)"
  echo "pointfinder_db_commit=$(git -C "${DB_ROOT}/pointfinder_db" rev-parse HEAD)"
} > "${RUNTIME_ROOT}/INSTALL_COMPLETE.txt"

echo "Pinned AMR runtime is ready at ${PREFIX}"
