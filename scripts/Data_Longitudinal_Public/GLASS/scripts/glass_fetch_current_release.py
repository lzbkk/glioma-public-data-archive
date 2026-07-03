#!/usr/bin/env python3

"""Fetch the minimum GLASS current-release inputs after access is approved."""

from __future__ import annotations

import os
from pathlib import Path

import synapseclient


PROJECT_ROOT = Path("/home/lzb/glioma")
OUT_DIR = PROJECT_ROOT / "Data_Longitudinal_Public" / "GLASS" / "data"
TABLE_DIR = OUT_DIR / "tables"
FILE_DIR = OUT_DIR / "expression"

TABLES = {
    "clinical_cases": "syn69931132",
    "clinical_surgeries": "syn69931127",
    "analysis_tumor_pairs": "syn69931139",
    "analysis_rna_silver_set": "syn69930927",
}

FILES = {
    "gene_tpm_matrix_all_samples.tsv": "syn69961520",
}


def main() -> None:
    token = os.environ.get("SYNAPSE_AUTH_TOKEN")
    if not token:
        raise SystemExit(
            "SYNAPSE_AUTH_TOKEN is unset. Accept the GLASS access terms in "
            "Synapse, create a personal access token with view/download scopes, "
            "and export it in the shell before running this script."
        )

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FILE_DIR.mkdir(parents=True, exist_ok=True)

    syn = synapseclient.Synapse()
    syn.login(authToken=token, silent=True)

    for name, syn_id in TABLES.items():
        result = syn.tableQuery(f"SELECT * FROM {syn_id}")
        result.asDataFrame().to_csv(TABLE_DIR / f"{name}.csv", index=False)

    for name, syn_id in FILES.items():
        syn.get(syn_id, downloadLocation=str(FILE_DIR), ifcollision="overwrite.local")

    print(f"GLASS current-release inputs written to: {OUT_DIR}")


if __name__ == "__main__":
    main()
