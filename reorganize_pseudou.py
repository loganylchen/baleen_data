#!/usr/bin/env python3
"""
Reorganize pseudo-U (BID-seq) raw files into bidseq_genome.tsv format.

All 3 files share identical GSE179798 BID-seq schema (header at Excel row 4).

Standardized first 5 columns: chr, start, end, strand, label
Followed by source-specific columns.

Directories processed:
  - hek293t/pseudo-u/
  - hela/pseudo-u/
  - a549/pseudo-u/
"""

import pandas as pd
import numpy as np
import os

DATA = '/data'

DATASETS = [
    {
        'cell_line': 'hek293t',
        'filename': 'GSE179798_HEK293T_mRNA_WT_BID-seq.xlsx',
        'expected_rows': 543,
    },
    {
        'cell_line': 'hela',
        'filename': 'GSE179798_HeLa_mRNA_WT_BID-seq.xlsx',
        'expected_rows': 575,
    },
    {
        'cell_line': 'a549',
        'filename': 'GSE179798_A549_mRNA_WT_BID-seq.xlsx',
        'expected_rows': 922,
    },
]


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


def process_bidseq(cell_line, filename, expected_rows):
    """GSE179798 BID-seq xlsx → bidseq_genome.tsv

    Strip chr prefix. pos (1-based) → start=pos-1, end=pos. label=NA (positive-only).
    """
    indir = f'{DATA}/{cell_line}/pseudo-u'
    print(f"Loading {cell_line}/pseudo-u/{filename} ...")
    df = pd.read_excel(f'{indir}/{filename}', header=3)

    out = pd.DataFrame()
    out['chr'] = df['chr'].astype(str).str.replace('chr', '', regex=False)
    out['start'] = df['pos'].astype(int) - 1  # 1-based → 0-based
    out['end'] = df['pos'].astype(int)
    out['strand'] = df['strand'].str.strip()
    out['label'] = pd.NA  # positive-only dataset
    for col in df.columns:
        if col not in ('chr', 'pos', 'strand'):
            out[col] = df[col].values

    write_tsv(out, indir, 'bidseq_genome.tsv', expected_rows=expected_rows)


def main():
    print("=" * 60)
    print("Reorganizing pseudo-U (BID-seq) datasets")
    print("=" * 60)

    for ds in DATASETS:
        print(f"\n--- {ds['cell_line']}/pseudo-u ---")
        process_bidseq(ds['cell_line'], ds['filename'], ds['expected_rows'])

    print("\nDone!")


if __name__ == '__main__':
    main()
