# colorectal_biomarker_metagenomics — minimal metagenomic profiling

## 1) Data sources
This repo reproduces a small, offline-friendly metagenomics workflow using subsampled reads (<1 GB each) from public stool WGS datasets (2024–2025).

**Suggested sources (choose 1–2 samples to demo):**
ONCOBIOME cohorts (ENA): PRJEB72524, PRJEB72525, PRJEB72526, PRJEB72523 ;
NHSII (SRA): PRJNA1237248 ;
Cohort 6 (SRA): PRJNA1167935 

(Profiles/metadata from the study are on Zenodo: DOI 10.5281/zenodo.15069069)

Example used in this repo: SRR32733476 (stool WGS; single-end). Subsampled FASTQs are placed in workflow/data/processed/ and used as inputs. 
After DBs are cached locally under workflow/db/, the workflow can run offline.

**Repository layout: **
```bash
colorectal_biomarker_metagenomics/
├── README.md
├── metadata.yaml          
├── questions.yaml          
├── answers.yaml            
└── workflow/
    ├── data/
    │   ├── raw/            # full FASTQs (downloaded)
    │   └── processed/      # subsampled FASTQs (<1 GB each)
    ├── db/
    │   ├── host/           # Bowtie2 human index (e.g., hg37dec_v0.1.*.bt2)
    │   └── metaphlan/      # MetaPhlAn 4 DB cache (large; optional auto-install)
    ├── qc/
    │   └── kneaddata/      # cleaned FASTQs + KneadData logs
    ├── results/            # MetaPhlAn profiles + top-20 CSVs
    ├── 01_fetch_and_subsample.sh
    └── 02_run_workflow.sh
```
## 2) Installation
**Requirements (CLI):**
curl, awk, seqtk, pigz (for fetch + subsample) ;
kneaddata ≥ 0.12.0, bowtie2, Trimmomatic (the jar), metaphlan ≥ 4.1 ;
python3 with pandas

```bash
conda create -n meta-mini -y python=3.11
conda activate meta-mini
conda install -y -c bioconda kneaddata bowtie2 trimmomatic metaphlan seqtk pigz
python -m pip install pandas
```
Apt/pip (Debian/Ubuntu/Colab):
```bash
sudo apt-get update -y
sudo apt-get install -y bowtie2 trimmomatic seqtk pigz
python3 -m pip install "kneaddata==0.12.0" "metaphlan==4.1.1" pandas
```
## 3) Pre-processing / subsampling
Subsampled inputs look like:
```bash
workflow/data/processed/
└─ SRR32733476.sub.fastq.gz
```
To (re)download and subsample from ENA/SRA:
```bash
# from repo root
chmod +x workflow/01_fetch_and_subsample.sh
workflow/01_fetch_and_subsample.sh SRR32733476 1500000 42
```
By default this keeps ~1.5M reads (adjust second argument to control size). 

If the run is paired-end, the script intentionally uses R1 only to keep the demo tiny.

## 4) Databases (host + MetaPhlAn)
**Host (human) Bowtie2 index — place here:**
```bash
workflow/db/host/
  hg37dec_v0.1.1.bt2
  hg37dec_v0.1.2.bt2
  hg37dec_v0.1.3.bt2
  hg37dec_v0.1.4.bt2
  hg37dec_v0.1.rev.1.bt2
  hg37dec_v0.1.rev.2.bt2
```
Tip (one-liner install via KneadData helper):
```bash
kneaddata_database --download human_genome bowtie2 workflow/db/host
```
MetaPhlAn DB — either copy a prepared cache here:
```bash
workflow/db/metaphlan/
  mpa_latest
  bowtie2_indexes/...
  *.pkl, *.bt2, *.bz2, ...
```
or let the workflow script install it once with --install-metaphlan-db (large download).
## 5) Workflow
Step 1 — KneadData (trim + human decontam) → cleaned FASTQ
Step 2 — MetaPhlAn 4 → taxonomic profile + top-20 summaries

Run:
```bash
# from repo root
chmod +x workflow/02_run_workflow.sh

# Option A: Use existing MetaPhlAn DB in workflow/db/metaphlan/
workflow/02_run_workflow.sh SRR32733476 --threads 4 --trimmomatic-dir /usr/share/trimmomatic

# Option B: Install MetaPhlAn DB into workflow/db/metaphlan/ (one-time, large)
workflow/02_run_workflow.sh SRR32733476 --threads 4 --trimmomatic-dir /usr/share/trimmomatic --install-metaphlan-db
```
The script auto-detects the Bowtie2 host index prefix from workflow/db/host/*.1.bt2.
If Trimmomatic is not in /usr/share/trimmomatic or /usr/share/java, pass --trimmomatic-dir /path/to/dir-with-trimmomatic.jar.

## 6) Primary outputs
```bash
workflow/qc/kneaddata/
  SRR32733476.sub.fastq.gz         # cleaned, decontaminated reads
  SRR32733476.sub.log              # trimming/decontam log

workflow/results/
  SRR32733476_metaphlan.tsv        # full clade profile (MetaPhlAn v4)
  SRR32733476_metaphlan_top20_phylum.csv
  SRR32733476_metaphlan_top20_genus.csv
  SRR32733476_metaphlan_top20_species.csv
```
Quick checks:
```bash
zcat workflow/qc/kneaddata/SRR32733476.sub.fastq.gz | head
head -n 30 workflow/results/SRR32733476_metaphlan.tsv
```
## 7) Notes
- This minimal pipeline mirrors the paper’s per-sample steps (host decontam + taxonomic profiling).
- Cohort-level analyses (diversity, DA tests, ordinations) require multiple samples and study metadata.
- Files are kept <1 GB by subsampling so they’re Taiga/SEPAL-friendly; adjust N_READS in script 1 as needed.
- After first DB installs, re-runs are offline and quick; DB folders live under workflow/db/.


