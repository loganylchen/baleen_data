#!/usr/bin/env python3
"""
Reorganize FTO-m m6A SAC-seq raw files (GSE211303) into {tool}_genome.tsv format.

Two replicates use different controls:
  - ftom.ftop : FTO-minus (test) vs FTO-plus (demethylated control)
  - ftom.ivt  : FTO-minus (test) vs IVT (unmodified control)

The .hits.txt files contain only called positive sites; label=NA (positive-only).

Standardized first 5 columns: chr, start, end, strand, label
Followed by source-specific columns from the raw file.

Directories processed:
  - hela/m6a/
"""

import pandas as pd
import os

DATA = '/data'

DATASETS = [
    {
        'cell_line': 'hela',
        'filename': 'GSE211303_hela.polya.wt.ftom.ftop.rep1.deep.hits.txt',
        'outname': 'ftom_ftop_genome.tsv',
        'expected_rows': 40096,
    },
    {
        'cell_line': 'hela',
        'filename': 'GSE211303_hela.polya.wt.ftom.ivt.rep1.deep.hits.txt',
        'outname': 'ftom_ivt_genome.tsv',
        'expected_rows': 69834,
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


def process_ftom(cell_line, filename, outname, expected_rows):
    """GSE211303 .hits.txt → {outname}

    `pos` column format: 'chr1_632652_-' (chrom_pos_strand, 1-based).
    Convert to 0-based BED: start=pos-1, end=pos. Strip 'chr' prefix.
    """
    indir = f'{DATA}/{cell_line}/m6a'
    print(f"Loading {cell_line}/m6a/{filename} ...")
    df = pd.read_csv(f'{indir}/{filename}', sep='\t')

    # Parse 'chrN_POS_STRAND' from `pos` column
    parts = df['pos'].astype(str).str.rsplit('_', n=2, expand=True)
    parts.columns = ['_chr', '_pos', '_strand']

    out = pd.DataFrame()
    out['chr'] = parts['_chr'].str.replace('chr', '', regex=False)
    out['start'] = parts['_pos'].astype(int) - 1  # 1-based → 0-based
    out['end'] = parts['_pos'].astype(int)
    out['strand'] = parts['_strand']
    out['label'] = pd.NA  # positive-only dataset (called hits)
    for col in df.columns:
        if col != 'pos':
            out[col] = df[col].values

    write_tsv(out, indir, outname, expected_rows=expected_rows)


def main():
    print("=" * 60)
    print("Reorganizing FTO-m m6A (GSE211303) datasets")
    print("=" * 60)

    for ds in DATASETS:
        print(f"\n--- {ds['cell_line']}/m6a ---")
        process_ftom(ds['cell_line'], ds['filename'], ds['outname'], ds['expected_rows'])

    print("\nDone!")


if __name__ == '__main__':
    main()
