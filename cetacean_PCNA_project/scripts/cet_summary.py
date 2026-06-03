#!/usr/bin/env python3
"""
cet_summary.py
Purpose: to read all species' raw genome and rnaseq files and builds one summary CSV.

Usage (run from project root):
    python3 scripts/cet_summary.py
"""

import os
import csv
import glob
import pandas as pd
from pathlib import Path

TABLES_DIR  = Path("tables")
RNASEQ_DIR  = Path("rna_seq")
SUMMARY_OUT = TABLES_DIR / "cetacean_summary.csv"

rows = []

# Find all genome files
for genome_file in sorted(TABLES_DIR.glob("*_genomes.tsv")):
    common = genome_file.stem.replace("_genomes", "")
    rnaseq_file = RNASEQ_DIR / f"{common}_rnaseq_runinfo.csv"

    # Genome stats====
    gdf = pd.read_csv(genome_file, sep="\t")
    n_assemblies = len(gdf)

    # Prefer RefSeq (GCF_) as best assembly
    refseq = gdf[gdf["Assembly Accession"].str.startswith("GCF_", na=False)]
    best = refseq.iloc[0] if len(refseq) > 0 else gdf.iloc[0]

    best_acc   = best.get("Assembly Accession", "")
    best_level = best.get("Assembly Level", "")
    contig_n50 = best.get("Assembly Stats Contig N50", "")
    scaff_n50  = best.get("Assembly Stats Scaffold N50", "")
    total_len  = best.get("Assembly Stats Total Sequence Length", "")
    annotated  = "yes" if str(best.get("Annotation Name", "")).strip() not in ("", "nan") else "no"

    total_gbp = ""
    try:
        total_gbp = round(float(total_len) / 1e9, 3)
    except (ValueError, TypeError):
        pass

    # RNA-seq stats ======
    n_runs = bioprojects = platforms = date_range = ""

    if rnaseq_file.exists() and rnaseq_file.stat().st_size > 100:
        rdf = pd.read_csv(rnaseq_file, low_memory=False)

        n_runs      = len(rdf)
        bioprojects = rdf["BioProject"].nunique() if "BioProject" in rdf.columns else ""
        platforms   = ", ".join(rdf["Platform"].dropna().unique()) if "Platform" in rdf.columns else ""

        if "ReleaseDate" in rdf.columns:
            dates = pd.to_datetime(rdf["ReleaseDate"], errors="coerce").dropna()
            if not dates.empty:
                date_range = f"{dates.min().date()} – {dates.max().date()}"

    rows.append({
        "species":         common,
        "n_assemblies":    n_assemblies,
        "best_accession":  best_acc,
        "best_level":      best_level,
        "contig_n50_bp":   contig_n50,
        "scaffold_n50_bp": scaff_n50,
        "genome_size_gbp": total_gbp,
        "annotated":       annotated,
        "n_rnaseq_runs":   n_runs,
        "n_bioprojects":   bioprojects,
        "platforms":       platforms,
        "rnaseq_dates":    date_range,
    })

    print(f"  {common}: {n_assemblies} assemblies, {n_runs} RNA-seq runs")

# Write summary 
summary_df = pd.DataFrame(rows)
summary_df.to_csv(SUMMARY_OUT, index=False)
print(f"\nSummary written to {SUMMARY_OUT}")