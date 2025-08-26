#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   workflow/01_fetch_and_subsample.sh <SRR_ACCESSION> [N_READS] [SEED]
# Requires: curl, awk, seqtk, pigz   (aria2c optional but faster)

SRR="${1:?Provide SRR accession, e.g., SRR32733476}"
NREAD="${2:-1500000}"
SEED="${3:-42}"

# repo-relative paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="$(cd "${SCRIPT_DIR}" && pwd)"
RAW="${WORKDIR}/data/raw"
PROC="${WORKDIR}/data/processed"

mkdir -p "${RAW}" "${PROC}"

have() { command -v "$1" >/dev/null 2>&1; }
need() { have "$1" || { echo "[ERR] Missing dependency: $1" >&2; exit 2; }; }

need curl; need awk; need seqtk; need pigz

echo "[INFO] Resolving ENA URL for ${SRR}…"
report="$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=fastq_ftp,library_layout&format=tsv" | awk 'NR==2{print}')"
fastq_ftp="$(echo "$report" | cut -f1)"
layout="$(echo "$report" | cut -f2)"
urls="$(echo "$fastq_ftp" | sed 's/^ftp/https/; s/; /\n/g; s/;/\n/g')"

[[ -n "${urls// }" ]] || { echo "[ERR] Could not resolve FASTQ for ${SRR}" >&2; exit 2; }

out="${RAW}/${SRR}.fastq.gz"
dl_url="$(echo "$urls" | head -n1)"   # if paired, we take R1 only (demo keeps size small)
echo "[DL] ${dl_url}"
if have aria2c; then
  aria2c -x8 -s8 -o "$(basename "$out")" -d "$RAW" "$dl_url"
else
  curl -L -o "$out" "$dl_url"
fi

echo "[SUBSAMPLE] ${out} → ${NREAD} reads (SEED=${SEED})"
seqtk sample -s"${SEED}" "$out" "${NREAD}" | pigz -c > "${PROC}/${SRR}.sub.fastq.gz"

echo "[DONE] ${PROC}/${SRR}.sub.fastq.gz"
ls -lh "${PROC}/${SRR}.sub.fastq.gz"
