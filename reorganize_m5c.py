#!/usr/bin/env python3
"""
Reorganize m5C raw files into {tool}_genome.tsv format.

Standardized first 5 columns: chr, start, end, strand, label
Followed by source-specific columns.

Directories processed:
  - hek293t/m5c/
  - hela/m5c/
"""

import pandas as pd
import numpy as np
import os

DATA = '/data'


def write_tsv(df, outdir, name, expected_rows=None):
    os.makedirs(outdir, exist_ok=True)
    path = f'{outdir}/{name}'
    df.to_csv(path, sep='\t', index=False)
    n = len(df)
    check = f" (expected {expected_rows})" if expected_rows and n != expected_rows else ""
    print(f"  {name}: {n:,} rows{check}")
    assert list(df.columns[:5]) == ['chr', 'start', 'end', 'strand', 'label'], \
        f"Bad columns in {name}: {list(df.columns[:5])}"
    if expected_rows:
        assert n == expected_rows, f"Row count mismatch in {name}: {n} != {expected_rows}"


def process_gse122254():
    """hek293t/m5c/gse122254.tsv → gse122254_genome.tsv

    260 rows, 1-based no chr prefix. All rows are m5C sites (if m5C site = TRUE).
    """
    print("Loading hek293t/m5c/gse122254.tsv ...")
    indir = f'{DATA}/hek293t/m5c'
    df = pd.read_csv(f'{indir}/gse122254.tsv', sep='\t')

    out = pd.DataFrame()
    out['chr'] = df['Chromosome'].astype(str)
    out['start'] = df['Position'].astype(int) - 1  # 1-based → 0-based
    out['end'] = df['Position'].astype(int)
    out['strand'] = df['Strand']
    out['label'] = df['if m5C site'].map({True: 1, 'TRUE': 1, False: 0, 'FALSE': 0}).fillna(1).astype(int)
    for col in df.columns:
        if col not in ('Chromosome', 'Position', 'Strand', 'if m5C site'):
            out[col] = df[col].values

    write_tsv(out, f'{DATA}/hek293t/m5c', 'gse122254_genome.tsv', expected_rows=260)


def process_gse225614(cell_line):
    """GSE225614 WT sites → gse225614_genome.tsv

    Drop rows where chromosome=".". 1-based positions. Positive-only (label=NA).
    """
    if cell_line == 'hek293t':
        filename = 'GSE225614_HEK293T-WT_sites.tsv'
    else:
        filename = 'GSE225614_HeLa-WT_sites.tsv'

    indir = f'{DATA}/{cell_line}/m5c'
    print(f"Loading {cell_line}/m5c/{filename} ...")
    df = pd.read_csv(f'{indir}/{filename}', sep='\t')

    # Drop rRNA rows (chromosome=".")
    n_before = len(df)
    df = df[df['chromosome'] != '.'].copy()
    print(f"  Dropped {n_before - len(df)} rRNA rows (chromosome='.')")

    out = pd.DataFrame()
    out['chr'] = df['chromosome'].astype(str)
    out['start'] = df['position'].astype(int) - 1  # 1-based → 0-based
    out['end'] = df['position'].astype(int)
    out['strand'] = df['strand']
    out['label'] = pd.NA  # positive-only dataset
    for col in df.columns:
        if col not in ('chromosome', 'position', 'strand'):
            out[col] = df[col].values

    write_tsv(out, indir, 'gse225614_genome.tsv')


def process_gse140995():
    """hela/m5c/GSE140995_transcriptome-wide_sites.xlsx → gse140995_genome.tsv

    Already 0-based BED with chr prefix. Strip chr prefix. label=NA.
    """
    print("Loading hela/m5c/GSE140995_transcriptome-wide_sites.xlsx ...")
    indir = f'{DATA}/hela/m5c'
    df = pd.read_excel(f'{indir}/GSE140995_transcriptome-wide_sites.xlsx')

    # Drop rows with missing coordinates
    df = df.dropna(subset=['Chrom', 'Start', 'End', 'Strand']).copy()

    out = pd.DataFrame()
    out['chr'] = df['Chrom'].astype(str).str.replace('chr', '', regex=False)
    out['start'] = df['Start'].astype(int)
    out['end'] = df['End'].astype(int)
    out['strand'] = df['Strand']
    out['label'] = pd.NA
    for col in df.columns:
        if col not in ('Chrom', 'Start', 'End', 'Strand'):
            out[col] = df[col].values

    write_tsv(out, indir, 'gse140995_genome.tsv')


def process_gse93749():
    """hela/m5c/GSE93749_hg19_Human_m5C_sites_information.xls → gse93749_genome.tsv

    hg19 coordinates → liftover to hg38. Strip chr prefix. 1-based → 0-based.
    Label from siCTRL Status_rep1: "1"→1, "/"→0.
    """
    from pyliftover import LiftOver

    print("Loading hela/m5c/GSE93749_hg19_Human_m5C_sites_information.xls ...")
    indir = f'{DATA}/hela/m5c'
    df = pd.read_excel(
        f'{indir}/GSE93749_hg19_Human_m5C_sites_information.xls',
        header=1,
    )

    print(f"  Raw rows: {len(df):,}")

    # Label from siCTRL Status_rep1: "1"→1, "/"→0
    df['label_raw'] = df['Status_rep1'].astype(str).map({'1': 1, '/': 0})

    # Liftover hg19 → hg38
    lo = LiftOver('hg19', 'hg38')
    chrom_hg38 = []
    pos_hg38 = []
    strand_hg38 = []
    failed = 0

    for _, row in df.iterrows():
        chrom = str(row['Chromosome'])
        if not chrom.startswith('chr'):
            chrom = 'chr' + chrom
        pos = int(row['Position']) - 1  # 1-based → 0-based for liftover
        strand = str(row['Strand'])

        result = lo.convert_coordinate(chrom, pos, strand)
        if result and len(result) > 0:
            new_chrom, new_pos, new_strand, _ = result[0]
            chrom_hg38.append(new_chrom.replace('chr', ''))
            pos_hg38.append(new_pos)
            strand_hg38.append(new_strand)
        else:
            chrom_hg38.append(None)
            pos_hg38.append(None)
            strand_hg38.append(None)
            failed += 1

    df['chr_hg38'] = chrom_hg38
    df['pos_hg38'] = pos_hg38
    df['strand_hg38'] = strand_hg38

    # Drop failed liftover rows
    n_before = len(df)
    df = df.dropna(subset=['chr_hg38']).copy()
    print(f"  Liftover: {failed} failed, {len(df):,} succeeded")

    out = pd.DataFrame()
    out['chr'] = df['chr_hg38'].astype(str)
    out['start'] = df['pos_hg38'].astype(int)
    out['end'] = df['pos_hg38'].astype(int) + 1
    out['strand'] = df['strand_hg38'].astype(str)
    out['label'] = df['label_raw'].values
    for col in df.columns:
        if col not in ('Chromosome', 'Position', 'Strand', 'label_raw',
                        'chr_hg38', 'pos_hg38', 'strand_hg38'):
            out[col] = df[col].values

    write_tsv(out, indir, 'gse93749_genome.tsv')


def main():
    print("=" * 60)
    print("Reorganizing m5C datasets")
    print("=" * 60)

    # hek293t/m5c
    print("\n--- hek293t/m5c ---")
    process_gse122254()
    process_gse225614('hek293t')

    # hela/m5c
    print("\n--- hela/m5c ---")
    process_gse225614('hela')
    process_gse140995()
    process_gse93749()

    print("\nDone!")


if __name__ == '__main__':
    main()
