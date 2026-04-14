#!/usr/bin/env python3
"""
Convert *_genome.tsv files to *_transcriptome.tsv using Ensembl 104 GTF.

Maps genomic coordinates to transcript-relative coordinates via exon overlap.
One genomic site may map to multiple transcripts (all mappings kept).

Usage:
    docker build -f Dockerfile.analysis -t baleen-analysis .
    docker run --rm -v "$(pwd):/data" baleen-analysis python /data/genome_to_transcriptome.py
    docker run --rm -v "$(pwd):/data" baleen-analysis python /data/genome_to_transcriptome.py \
        --dirs hek293t/m5c hela/m5c hek293t/pseudo-u hela/pseudo-u a549/pseudo-u
"""

import argparse
import pandas as pd
import pyranges as pr
import numpy as np
import urllib.request
import gzip
import shutil
from pathlib import Path


GTF_URL = "https://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz"
GTF_FILE = "Homo_sapiens.GRCh38.104.gtf"
DATA_DIR = Path("/data")


def download_gtf():
    """Download and decompress Ensembl 104 GTF if not present."""
    gtf_path = DATA_DIR / GTF_FILE
    if gtf_path.exists():
        print(f"GTF exists: {gtf_path}")
        return gtf_path

    gz_path = DATA_DIR / f"{GTF_FILE}.gz"
    print(f"Downloading {GTF_URL} ...")
    urllib.request.urlretrieve(GTF_URL, gz_path)
    print("Decompressing ...")
    with gzip.open(gz_path, "rb") as f_in, open(gtf_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    gz_path.unlink()
    print(f"GTF ready: {gtf_path}")
    return gtf_path


def build_exon_models(gtf_path):
    """Parse GTF exons and compute cumulative transcript offsets.

    Returns a DataFrame with columns:
        Chromosome, Start, End, Strand, transcript_id, gene_id,
        exon_number, exon_length, cum_offset
    Coordinates are 0-based half-open (pyranges convention).
    """
    print("Parsing GTF ...")
    gtf = pr.read_gtf(str(gtf_path))
    exons = gtf[gtf.Feature == "exon"].df.copy()

    exons = exons[
        ["Chromosome", "Start", "End", "Strand",
         "transcript_id", "gene_id", "exon_number"]
    ].copy()
    # Convert Categorical columns from GTF parsing to plain strings
    exons["Strand"] = exons["Strand"].astype(str)
    exons["Chromosome"] = exons["Chromosome"].astype(str)
    exons["exon_number"] = exons["exon_number"].astype(int)
    exons["exon_length"] = exons["End"] - exons["Start"]

    # Sort by transcript + exon_number (= 5'→3' transcript order)
    exons.sort_values(["transcript_id", "exon_number"], inplace=True)
    exons.reset_index(drop=True, inplace=True)

    # Cumulative offset = transcript position at the start of each exon
    exons["cum_offset"] = (
        exons.groupby("transcript_id")["exon_length"].cumsum()
        - exons["exon_length"]
    )

    n_tx = exons["transcript_id"].nunique()
    print(f"  {len(exons):,} exon records from {n_tx:,} transcripts")
    return exons


def convert_file(genome_tsv, exons_df):
    """Convert one genome TSV → transcriptome TSV via pyranges overlap join."""
    print(f"\n{'='*60}")
    print(f"Converting {genome_tsv.name}")
    sites = pd.read_csv(genome_tsv, sep="\t", dtype={"chr": str})
    n_input = len(sites)
    print(f"  Input: {n_input:,} sites")

    # Drop rows with invalid strand (must be + or -)
    sites = sites[sites["strand"].isin(["+", "-"])].copy()
    if len(sites) < n_input:
        print(f"  Dropped {n_input - len(sites)} sites with invalid strand")
        n_input = len(sites)

    # ---- prepare sites PyRanges ----
    extra_cols = [c for c in sites.columns if c not in ("chr", "start", "end", "strand")]

    sp = pd.DataFrame({
        "Chromosome": sites["chr"].astype(str),
        "Start": sites["start"].astype(int),
        "End": sites["end"].astype(int),
        "Strand": sites["strand"].astype(str),
    })
    for col in extra_cols:
        sp[f"s_{col}"] = sites[col].values
    sp["s_idx"] = np.arange(n_input)

    sites_pr = pr.PyRanges(sp)

    # ---- prepare exons PyRanges ----
    ep = exons_df[
        ["Chromosome", "Start", "End", "Strand",
         "transcript_id", "gene_id", "cum_offset"]
    ].copy()
    exons_pr = pr.PyRanges(ep)

    # ---- strand-aware overlap join ----
    joined = sites_pr.join(exons_pr, strandedness="same").df

    if joined.empty:
        print("  WARNING: no overlaps found — skipping")
        return None

    # After join, site coords are Start/End, exon coords are Start_b/End_b
    # Compute transcript-relative coordinates
    plus = joined["Strand"] == "+"

    tx_start = pd.Series(0, index=joined.index, dtype=int)
    # + strand: tx_pos = cum_offset + (site_start − exon_start)
    tx_start[plus] = (
        joined.loc[plus, "cum_offset"]
        + joined.loc[plus, "Start"]
        - joined.loc[plus, "Start_b"]
    ).astype(int)
    # − strand: tx_pos = cum_offset + (exon_end − site_end)
    tx_start[~plus] = (
        joined.loc[~plus, "cum_offset"]
        + joined.loc[~plus, "End_b"]
        - joined.loc[~plus, "End"]
    ).astype(int)

    tx_end = tx_start + (joined["End"] - joined["Start"]).astype(int)

    # ---- build output ----
    out = pd.DataFrame()
    out["chr"] = joined["transcript_id"].values
    out["start"] = tx_start.values
    out["end"] = tx_end.values
    out["strand"] = "."

    for col in extra_cols:
        out[col] = joined[f"s_{col}"].values

    out["gene_id"] = joined["gene_id"].values
    # Build genomic_position string (handle Categorical dtypes from pyranges)
    chrom_str = joined["Chromosome"].astype(object).astype(str)
    strand_str = joined["Strand"].astype(object).astype(str)
    out["genomic_position"] = (
        chrom_str + ":"
        + joined["Start"].astype(str) + ":"
        + joined["End"].astype(str) + ":"
        + strand_str
    )

    out.sort_values(["chr", "start"], inplace=True)
    out.reset_index(drop=True, inplace=True)

    # ---- write ----
    out_name = genome_tsv.name.replace("_genome.tsv", "_transcriptome.tsv")
    out_path = genome_tsv.parent / out_name
    out.to_csv(out_path, sep="\t", index=False)

    n_output = len(out)
    n_tx = out["chr"].nunique()
    n_mapped = joined["s_idx"].nunique()
    print(f"  Mapped sites: {n_mapped:,}/{n_input:,} ({100*n_mapped/n_input:.1f}%)")
    print(f"  Output: {n_output:,} rows across {n_tx:,} transcripts → {out_name}")
    print(f"  Avg transcripts/site: {n_output/n_mapped:.1f}")

    return out_path


def parse_args():
    parser = argparse.ArgumentParser(description="Convert genome TSVs to transcriptome coordinates")
    parser.add_argument(
        "--dirs", nargs="+", default=["hek293t/m6a"],
        help="Directories to process (relative to /data). Default: hek293t/m6a",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    gtf_path = download_gtf()
    exons_df = build_exon_models(gtf_path)

    for d in args.dirs:
        target_dir = DATA_DIR / d
        genome_files = sorted(target_dir.glob("*_genome.tsv"))
        print(f"\n{'='*60}")
        print(f"Directory: {d} — {len(genome_files)} genome files")
        for f in genome_files:
            print(f"  {f.name}")

        for gf in genome_files:
            convert_file(gf, exons_df)

    print(f"\n{'='*60}")
    print("All conversions complete.")


if __name__ == "__main__":
    main()
