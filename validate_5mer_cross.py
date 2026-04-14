#!/usr/bin/env python3
"""
Cross-validate 5-mers between genome FASTA and transcriptome FASTA.

For each *_transcriptome.tsv, extract the 5-mer from BOTH:
  1. Genome FASTA — using genomic_position (chr:start:end:strand), with RC for minus strand
  2. Transcriptome FASTA — using transcript_id (chr column) + start, direct extraction (no RC)

Then compare the two extractions to verify coordinate conversion correctness.

Requires gffread to build the transcriptome FASTA from genome + GTF.

Usage:
    docker run --rm -v "$(pwd):/data" baleen-analysis python /data/validate_5mer_cross.py
    docker run --rm -v "$(pwd):/data" baleen-analysis python /data/validate_5mer_cross.py \
        --dirs hek293t/m5c hela/m5c hek293t/pseudo-u hela/pseudo-u a549/pseudo-u
"""

import argparse
import subprocess
import sys
import pandas as pd
import pysam
import urllib.request
import gzip
import shutil
from pathlib import Path

GENOME_URL = (
    "https://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/"
    "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
)
GENOME_FILE = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_URL = (
    "https://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/"
    "Homo_sapiens.GRCh38.104.gtf.gz"
)
GTF_FILE = "Homo_sapiens.GRCh38.104.gtf"
TX_FASTA = "transcriptome.fa"

DATA_DIR = Path("/data")

COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq):
    return seq.translate(COMP)[::-1]


def download_and_decompress(url, filename):
    """Download a .gz file and decompress it if the target doesn't exist."""
    path = DATA_DIR / filename
    if path.exists():
        print(f"  Exists: {path}")
        return path

    gz_path = DATA_DIR / f"{filename}.gz"
    if not gz_path.exists():
        print(f"  Downloading {url} ...")
        urllib.request.urlretrieve(url, gz_path)

    print(f"  Decompressing {gz_path} ...")
    with gzip.open(gz_path, "rb") as f_in, open(path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    gz_path.unlink()
    print(f"  Ready: {path}")
    return path


def build_transcriptome_fasta(genome_fa, gtf_path):
    """Run gffread to generate transcriptome FASTA from genome + GTF."""
    tx_fa = DATA_DIR / TX_FASTA
    if tx_fa.exists():
        print(f"  Transcriptome FASTA exists: {tx_fa}")
        return tx_fa

    print("  Running gffread ...")
    cmd = [
        "gffread", str(gtf_path),
        "-g", str(genome_fa),
        "-w", str(tx_fa),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  gffread stderr: {result.stderr}", file=sys.stderr)
        sys.exit(1)
    print(f"  Transcriptome FASTA ready: {tx_fa}")
    return tx_fa


def open_fasta(fa_path):
    """Index (if needed) and open a FASTA for random access."""
    fai_path = Path(str(fa_path) + ".fai")
    if not fai_path.exists():
        print(f"  Indexing {fa_path.name} ...")
        pysam.faidx(str(fa_path))
    return pysam.FastaFile(str(fa_path))


def extract_5mer_genome(fa, chrom, start, strand):
    """Extract 5-mer from genome FASTA centered on the site. RC for minus strand."""
    left = start - 2
    right = start + 3
    if left < 0:
        return None
    try:
        seq = fa.fetch(chrom, left, right).upper()
    except (KeyError, ValueError):
        return None
    if len(seq) != 5:
        return None
    if strand == "-":
        seq = reverse_complement(seq)
    return seq


def extract_5mer_transcriptome(fa, transcript_id, start):
    """Extract 5-mer from transcriptome FASTA. Coords are already 5'->3', no RC."""
    left = start - 2
    right = start + 3
    if left < 0:
        return None
    try:
        seq = fa.fetch(transcript_id, left, right).upper()
    except (KeyError, ValueError):
        return None
    if len(seq) != 5:
        return None
    return seq


def validate_file(genome_fa, tx_fa, tx_tsv):
    """Cross-validate 5-mers for one transcriptome TSV file."""
    print(f"\n--- {tx_tsv.name} ---")
    df = pd.read_csv(tx_tsv, sep="\t", dtype={"chr": str})

    if "genomic_position" not in df.columns:
        print("  No genomic_position column — skipping")
        return 0, 0, 0

    sample_gp = str(df["genomic_position"].iloc[0])
    if sample_gp.count(":") != 3:
        print("  genomic_position format not chr:start:end:strand — skipping")
        return 0, 0, 0

    match = 0
    mismatch = 0
    skipped = 0
    mismatch_examples = []

    for _, row in df.iterrows():
        parts = str(row["genomic_position"]).split(":")
        g_chrom, g_start, g_strand = parts[0], int(parts[1]), parts[3]

        tx_id = str(row["chr"])
        tx_start = int(row["start"])

        genome_5mer = extract_5mer_genome(genome_fa, g_chrom, g_start, g_strand)
        tx_5mer = extract_5mer_transcriptome(tx_fa, tx_id, tx_start)

        if genome_5mer is None or tx_5mer is None:
            skipped += 1
            continue

        if genome_5mer == tx_5mer:
            match += 1
        else:
            mismatch += 1
            if len(mismatch_examples) < 5:
                mismatch_examples.append(
                    f"    {tx_id}:{tx_start}  genomic={row['genomic_position']}  "
                    f"genome_5mer={genome_5mer}  tx_5mer={tx_5mer}"
                )

    total = match + mismatch
    print(f"  Compared: {total:,}")
    if total:
        print(f"  Match: {match:,} ({100 * match / total:.1f}%)")
    if mismatch:
        print(f"  Mismatch: {mismatch:,}")
        for ex in mismatch_examples:
            print(ex)
    if skipped:
        print(f"  Skipped (could not extract): {skipped:,}")

    return match, mismatch, skipped


def parse_args():
    parser = argparse.ArgumentParser(description="Cross-validate 5-mers between genome and transcriptome FASTA")
    parser.add_argument(
        "--dirs", nargs="+", default=["hek293t/m6a"],
        help="Directories to process (relative to /data). Default: hek293t/m6a",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    print("=" * 60)
    print("CROSS-VALIDATION: genome FASTA vs transcriptome FASTA 5-mers")
    print("=" * 60)

    print("\n[1/3] Downloading reference files ...")
    genome_path = download_and_decompress(GENOME_URL, GENOME_FILE)
    gtf_path = download_and_decompress(GTF_URL, GTF_FILE)

    print("\n[2/3] Building transcriptome FASTA ...")
    tx_fasta_path = build_transcriptome_fasta(genome_path, gtf_path)

    print("\n[3/3] Cross-validating 5-mers ...")
    genome_fa = open_fasta(genome_path)
    tx_fa = open_fasta(tx_fasta_path)
    print(f"  Genome: {genome_fa.nreferences} sequences")
    print(f"  Transcriptome: {tx_fa.nreferences} sequences")

    total_match = 0
    total_mismatch = 0
    total_skipped = 0

    for d in args.dirs:
        target_dir = DATA_DIR / d
        tx_files = sorted(target_dir.glob("*_transcriptome.tsv"))
        if tx_files:
            print(f"\n--- {d} ---")
        for tsv in tx_files:
            m, mm, s = validate_file(genome_fa, tx_fa, tsv)
            total_match += m
            total_mismatch += mm
            total_skipped += s

    genome_fa.close()
    tx_fa.close()

    print("\n" + "=" * 60)
    total_compared = total_match + total_mismatch
    print(f"OVERALL: {total_compared:,} compared, "
          f"{total_match:,} match, {total_mismatch:,} mismatch, "
          f"{total_skipped:,} skipped")
    if total_compared:
        print(f"Match rate: {100 * total_match / total_compared:.1f}%")
    print("=" * 60)

    if total_compared > 0:
        mismatch_rate = total_mismatch / total_compared
        if mismatch_rate > 0.01:
            print(f"\nFAILED: mismatch rate {100*mismatch_rate:.2f}% exceeds 1% threshold!")
            sys.exit(1)
        elif total_mismatch > 0:
            print(f"\nPASSED: {total_mismatch:,} mismatches ({100*mismatch_rate:.3f}%) "
                  f"within 1% tolerance (likely exon boundary edge cases).")
        else:
            print("\nPASSED: all 5-mers match between genome and transcriptome FASTA.")
    else:
        print("\nNo comparisons made.")


if __name__ == "__main__":
    main()
