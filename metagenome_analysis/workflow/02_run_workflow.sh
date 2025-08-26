#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   workflow/02_run_workflow.sh <SRR_ACCESSION> [--threads N] [--trimmomatic-dir DIR] \
#                               [--host-db-prefix PATH] [--mp-db-dir DIR] [--install-metaphlan-db]
#
# Examples:
#   workflow/02_run_workflow.sh SRR32733476
#   workflow/02_run_workflow.sh SRR32733476 --threads 8 --trimmomatic-dir /usr/share/trimmomatic
#
# Expects: kneaddata (>=0.12), bowtie2, trimmomatic (jar), metaphlan (>=4.1), python3+pandas
# Places everything under workflow/: db/, qc/, results/, data/

SRR="${1:?Provide SRR accession, e.g., SRR32733476}"; shift || true

THREADS=4
TRIMDIR=""              # dir containing trimmomatic.jar (auto-detect if empty)
HOST_DB_PREFIX=""       # e.g., workflow/db/host/hg37dec_v0.1  (no .1.bt2)
MP_DB_DIR=""            # e.g., workflow/db/metaphlan
INSTALL_MPDB=0          # pass --install-metaphlan-db to download (LARGE)

while [[ $# -gt 0 ]]; do
  case "$1" in
    --threads) THREADS="${2:?}"; shift 2;;
    --trimmomatic-dir) TRIMDIR="${2:?}"; shift 2;;
    --host-db-prefix) HOST_DB_PREFIX="${2:?}"; shift 2;;
    --mp-db-dir) MP_DB_DIR="${2:?}"; shift 2;;
    --install-metaphlan-db) INSTALL_MPDB=1; shift;;
    *) echo "[ERR] Unknown option: $1" >&2; exit 2;;
  esac
done

have() { command -v "$1" >/dev/null 2>&1; }
need() { have "$1" || { echo "[ERR] Missing dependency: $1" >&2; exit 2; }; }

need kneaddata; need bowtie2; need metaphlan; need python3

# ensure pandas for the summary step
python3 - <<'PY' 2>/dev/null || python3 -m pip -q install pandas
import pandas as pd
print("ok")
PY

# repo-relative paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="$(cd "${SCRIPT_DIR}" && pwd)"
PROC="${WORKDIR}/data/processed"
DBHOST_DIR="${WORKDIR}/db/host"
MPDB="${MP_DB_DIR:-${WORKDIR}/db/metaphlan}"
KNEAD="${WORKDIR}/qc/kneaddata"
RES="${WORKDIR}/results"

mkdir -p "$KNEAD" "$RES" "$MPDB"

IN_SE="${PROC}/${SRR}.sub.fastq.gz"
[[ -s "$IN_SE" ]] || { echo "[ERR] Missing ${IN_SE}. Run workflow/01_fetch_and_subsample.sh first." >&2; exit 2; }

# Trimmomatic directory
if [[ -z "${TRIMDIR}" ]]; then
  if [[ -f /usr/share/java/trimmomatic.jar ]]; then
    TRIMDIR="/usr/share/java"
  elif [[ -d /usr/share/trimmomatic ]]; then
    TRIMDIR="/usr/share/trimmomatic"
  else
    echo "[ERR] Provide --trimmomatic-dir pointing to directory containing trimmomatic.jar" >&2
    exit 2
  fi
fi
[[ -f "${TRIMDIR}/trimmomatic.jar" ]] || { echo "[ERR] trimmomatic.jar not found in ${TRIMDIR}" >&2; exit 2; }

# Host Bowtie2 DB prefix
if [[ -z "${HOST_DB_PREFIX}" ]]; then
  p="$(ls "${DBHOST_DIR}"/*.1.bt2 2>/dev/null | head -n1 || true)"
  [[ -n "$p" ]] || { echo "[ERR] No Bowtie2 index found under ${DBHOST_DIR}. Place hg37dec_v0.1.*.bt2 (or equivalent) there." >&2; exit 2; }
  HOST_DB_PREFIX="${p%.1.bt2}"
fi

echo "[KNEADDATA] input=${IN_SE}"
echo "[KNEADDATA] host_db_prefix=${HOST_DB_PREFIX}"
echo "[KNEADDATA] trimmomatic_dir=${TRIMDIR}"
kneaddata \
  --unpaired "${IN_SE}" \
  -db "${HOST_DB_PREFIX}" \
  --trimmomatic "${TRIMDIR}" \
  -t "${THREADS}" \
  --bypass-trf \
  --remove-intermediate-output \
  --output-prefix "${SRR}.sub" \
  -o "${KNEAD}"

# Ensure gz output exists (kneaddata may output .fastq)
CLEAN="${KNEAD}/${SRR}.sub.fastq.gz"
if [[ ! -s "$CLEAN" && -s "${KNEAD}/${SRR}.sub.fastq" ]]; then
  gzip -c "${KNEAD}/${SRR}.sub.fastq" > "$CLEAN"
fi
[[ -s "$CLEAN" ]] || { echo "[ERR] Cleaned FASTQ not found in ${KNEAD}" >&2; exit 3; }

# MetaPhlAn DB
IDX=""
if [[ -f "${MPDB}/mpa_latest" ]]; then
  IDX="$(cat "${MPDB}/mpa_latest")"
elif [[ "${INSTALL_MPDB}" -eq 1 ]]; then
  echo "[MetaPhlAn] Installing DB into ${MPDB} (LARGE download)…"
  metaphlan --install --bowtie2db "${MPDB}" --nproc "${THREADS}"
  [[ -f "${MPDB}/mpa_latest" ]] && IDX="$(cat "${MPDB}/mpa_latest")"
else
  echo "[ERR] MetaPhlAn DB not found in ${MPDB}. Pre-populate it or rerun with --install-metaphlan-db" >&2
  exit 4
fi

OUT_TSV="${RES}/${SRR}_metaphlan.tsv"
echo "[MetaPhlAn] profiling → ${OUT_TSV}"
metaphlan "${CLEAN}" \
  --input_type fastq \
  --bowtie2db "${MPDB}" \
  -x "${IDX}" \
  --offline \
  --nproc "${THREADS}" \
  --add_viruses \
  -o "${OUT_TSV}"

# ---- Per-rank top20 CSVs (phylum/genus/species) ----
export MPA_TSV="${OUT_TSV}"
python3 - <<'PY'
import pandas as pd, os
tsv=os.environ["MPA_TSV"]
outdir=os.path.dirname(tsv)
srr=os.path.basename(tsv).split("_metaphlan.tsv")[0]

rows=[]; header=None
with open(tsv) as f:
    for line in f:
        if line.startswith("#clade_name"):
            header=line[1:].strip().split("\t"); continue
        if line.startswith("#") or not line.strip(): continue
        rows.append(line.rstrip("\n").split("\t"))
if header is None:
    header=["clade_name","NCBI_tax_id","relative_abundance","additional_species"]
df=pd.DataFrame(rows, columns=header)
df["relative_abundance"]=pd.to_numeric(df["relative_abundance"], errors="coerce").fillna(0.0)

def ranks(s):
    d={}
    for p in s.split("|"):
        if "__" in p:
            r,n=p.split("__",1); d[r]=n
    return d

rec=[]
for _,r in df.iterrows():
    d=ranks(r["clade_name"])
    rec.append({"kingdom":d.get("k",""),"phylum":d.get("p",""),"class":d.get("c",""),
                "order":d.get("o",""),"family":d.get("f",""),"genus":d.get("g",""),
                "species":d.get("s",""),"relative_abundance":r["relative_abundance"]})
tdf=pd.DataFrame(rec)

def collapse(level):
    return (tdf[tdf[level]!=""]
            .groupby(level,as_index=False)["relative_abundance"].sum()
            .sort_values("relative_abundance",ascending=False)
            .head(20))

collapse("phylum").to_csv(os.path.join(outdir, f"{srr}_metaphlan_top20_phylum.csv"), index=False)
collapse("genus").to_csv(os.path.join(outdir, f"{srr}_metaphlan_top20_genus.csv"), index=False)
collapse("species").to_csv(os.path.join(outdir, f"{srr}_metaphlan_top20_species.csv"), index=False)
print("[OK] Wrote top-20 CSVs in", outdir)
PY

echo "[DONE]"
ls -lh "${RES}" | sed -n '1,200p'
