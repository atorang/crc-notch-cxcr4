"""00_download_data.py

Pelka et al. scRNA-seq data acquisition. The Broad Single Cell Portal
(SCP1162) requires an authenticated, manual download, see
docs/DATA_ACCESS.md for the exact steps and files. This script's automated
default path is the GEO mirror (GSE178341).

Fails loudly if neither source is available.
"""
import sys
import tarfile
import urllib.request
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
cfg = yaml.safe_load((REPO_ROOT / "config" / "config.yaml").read_text())

raw_dir = REPO_ROOT / cfg["paths"]["data_raw"] / "scrnaseq"
raw_dir.mkdir(parents=True, exist_ok=True)

SCP_EXPECTED_FILES = [
    "scRNAseq_Pelka_Epithelial.h5ad",  # or the raw metatable + matrix.mtx.gz triple, see DATA_ACCESS.md
]
GEO_ACCESSION = cfg["data_sources"]["pelka"]["geo_accession"]
GEO_SUPP_BASE = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{GEO_ACCESSION[:-3]}nnn/{GEO_ACCESSION}/suppl"
# Confirmed live (2026-08-20) at the URL above:
#   GSE178341_crc10x_full_c295v4_submit.h5              (raw UMI counts — confirmed integer values)
#   GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz
#   GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz
GEO_FILES = [
    "GSE178341_crc10x_full_c295v4_submit.h5",
    "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz",
    "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz",
]


def scp_files_present():
    return all((raw_dir / f).exists() for f in SCP_EXPECTED_FILES)


def download_from_geo():
    print(f"SCP files not found manually; downloading from GEO {GEO_ACCESSION}...")
    for fname in GEO_FILES:
        dest = raw_dir / fname
        if dest.exists():
            print("already have", fname)
            continue
        url = f"{GEO_SUPP_BASE}/{fname}"
        print("Downloading", url, "...")
        try:
            urllib.request.urlretrieve(url, dest)
        except Exception as e:
            raise SystemExit(
                f"Failed to download {fname} from GEO ({e}). Either fix network access, or "
                f"follow the manual SCP1162 download in docs/DATA_ACCESS.md and place the "
                f"files in {raw_dir}."
            )
    print(
        "Downloaded raw GEO matrices (GSE178341_crc10x_full_c295v4_submit.h5 carries true raw "
        "UMI counts, confirmed by direct inspection, better than the SCP mirror, which is "
        "normalized-only). NOTE: this GEO bundle's cell/cluster metadata "
        "(*_cluster.csv.gz / *_metatables.csv.gz) has NOT been diffed here against the "
        "SCP-derived metatable_v3_fix_v3.tsv / crc10x_tSNE_cl_Epi.tsv the existing pipeline "
        "scripts actually read, that cross-check (matching ClusterMidway/donor_id/MMR-status "
        "labels cell-for-cell) still needs to be done and reported in docs/DATA_ACCESS.md "
        "before this GEO path can be trusted as a drop-in replacement for the manual SCP export."
    )


if __name__ == "__main__":
    if scp_files_present():
        print("Manually-downloaded SCP files found: nothing to do.")
        sys.exit(0)
    download_from_geo()
