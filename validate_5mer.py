#!/usr/bin/env python3
"""
Validate 5-mer motifs by extracting sequences from Ensembl 104 genome.

For genome files with a Motif column: extract 5-mer at each site and compare.
For transcriptome files: use genomic_position to extract and verify.

Usage:
    docker run --rm -v "$(pwd):/data" baleen-analysis python /data/validate_5mer.py
"""

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
DATA_DIR = Path("/data")
M6A_DIR = DATA_DIR / "hek293t" / "m6a"

COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq):
    return seq.translate(COMP)[::-1]


def download_genome():
    """Download and decompress Ensembl 104 primary assembly if not present."""
    fa_path = DATA_DIR / GENOME_FILE
    if fa_path.exists():
        print(f"Genome exists: {fa_path}")
        return fa_path

    gz_path = DATA_DIR / f"{GENOME_FILE}.gz"
    if not gz_path.exists():
        print(f"Downloading genome ({GENOME_URL}) ...")
        urllib.request.urlretrieve(GENOME_URL, gz_path)

    print("Decompressing (this takes a while) ...")
    with gzip.open(gz_path, "rb") as f_in, open(fa_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    gz_path.unlink()
    print(f"Genome ready: {fa_path}")
    return fa_path


def index_genome(fa_path):
    """Create .fai index if missing, then open for random access."""
    fai_path = Path(str(fa_path) + ".fai")
    if not fai_path.exists():
        print("Indexing genome ...")
        pysam.faidx(str(fa_path))
    return pysam.FastaFile(str(fa_path))


def extract_5mer(fa, chrom, start, end, strand):
    """Extract the 5-mer centered on the site (single-base, 0-based).

    The m6A site is at position 2 (0-indexed) of the 5-mer (= the 'A' in DRACH).
    So 5-mer spans [start-2, start+3) on the + strand sense.
    """
    pos = start  # 0-based position of the modified base
    left = pos - 2
    right = pos + 3

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


def validate_genome_file(fa, genome_tsv):
    """Validate 5-mers for a genome file that has a Motif column."""
    print(f"\n--- {genome_tsv.name} ---")
    df = pd.read_csv(genome_tsv, sep="\t", dtype={"chr": str})

    if "Motif" not in df.columns:
        print("  No Motif column — skipping")
        return

    total = 0
    match = 0
    mismatch = 0
    skipped = 0
    mismatch_examples = []

    for _, row in df.iterrows():
        motif = str(row["Motif"]).strip()
        if not motif or motif == "nan":
            skipped += 1
            continue

        total += 1
        extracted = extract_5mer(fa, str(row["chr"]), int(row["start"]),
                                 int(row["end"]), str(row["strand"]))
        if extracted is None:
            skipped += 1
            continue

        if extracted == motif.upper():
            match += 1
        else:
            mismatch += 1
            if len(mismatch_examples) < 5:
                mismatch_examples.append(
                    f"    {row['chr']}:{row['start']}:{row['strand']}  "
                    f"expected={motif}  got={extracted}"
                )

    print(f"  Total with Motif: {total:,}")
    print(f"  Match: {match:,} ({100*match/total:.1f}%)" if total else "  No motifs")
    if mismatch:
        print(f"  Mismatch: {mismatch:,}")
        for ex in mismatch_examples:
            print(ex)
    if skipped:
        print(f"  Skipped (no motif/out of range): {skipped:,}")


def validate_transcriptome_file(fa, tx_tsv):
    """Validate 5-mers for a transcriptome file using genomic_position."""
    print(f"\n--- {tx_tsv.name} ---")
    df = pd.read_csv(tx_tsv, sep="\t")

    if "genomic_position" not in df.columns:
        print("  No genomic_position column — skipping (original m6anet?)")
        return

    # Check if genomic_position is in our chr:start:end:strand format
    sample_gp = str(df["genomic_position"].iloc[0])
    if sample_gp.count(":") != 3:
        print(f"  genomic_position format not chr:start:end:strand — skipping")
        return

    if "Motif" not in df.columns:
        print("  No Motif column — extracting 5-mers from genomic_position only")
        # Just show we CAN extract 5-mers; sample 10 rows
        sample = df.head(10)
        for _, row in sample.iterrows():
            parts = str(row["genomic_position"]).split(":")
            chrom, start, end, strand = parts[0], int(parts[1]), int(parts[2]), parts[3]
            mer = extract_5mer(fa, chrom, start, end, strand)
            print(f"  {row['chr']}:{row['start']}  genomic={row['genomic_position']}  "
                  f"5-mer={mer}")
        print(f"  ... ({len(df):,} total rows)")
        return

    # Has both Motif and genomic_position: validate
    total = 0
    match = 0
    mismatch = 0
    skipped = 0
    mismatch_examples = []

    for _, row in df.iterrows():
        motif = str(row["Motif"]).strip()
        if not motif or motif == "nan":
            skipped += 1
            continue

        total += 1
        parts = str(row["genomic_position"]).split(":")
        chrom, start, end, strand = parts[0], int(parts[1]), int(parts[2]), parts[3]
        extracted = extract_5mer(fa, chrom, start, end, strand)

        if extracted is None:
            skipped += 1
            continue

        if extracted == motif.upper():
            match += 1
        else:
            mismatch += 1
            if len(mismatch_examples) < 5:
                mismatch_examples.append(
                    f"    {row['chr']}:{row['start']}  "
                    f"genomic={row['genomic_position']}  "
                    f"expected={motif}  got={extracted}"
                )

    print(f"  Total with Motif: {total:,}")
    print(f"  Match: {match:,} ({100*match/total:.1f}%)" if total else "  No motifs")
    if mismatch:
        print(f"  Mismatch: {mismatch:,}")
        for ex in mismatch_examples:
            print(ex)
    if skipped:
        print(f"  Skipped: {skipped:,}")


def main():
    fa_path = download_genome()
    fa = index_genome(fa_path)
    print(f"Genome loaded: {fa.nreferences} chromosomes")

    # --- Validate genome files ---
    print("\n" + "=" * 60)
    print("GENOME FILE VALIDATION (Motif vs extracted 5-mer)")
    print("=" * 60)
    for tsv in sorted(M6A_DIR.glob("*_genome.tsv")):
        validate_genome_file(fa, tsv)

    # --- Validate transcriptome files ---
    print("\n" + "=" * 60)
    print("TRANSCRIPTOME FILE VALIDATION (genomic_position → 5-mer)")
    print("=" * 60)
    for tsv in sorted(M6A_DIR.glob("*_transcriptome.tsv")):
        validate_transcriptome_file(fa, tsv)

    fa.close()
    print("\nDone!")


if __name__ == "__main__":
    main()
