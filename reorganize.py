#!/usr/bin/env python3
"""
Reorganize baleen_data into hek293t/m6a/{tool}_genome.tsv format.

Standardized first 5 columns: chr, start, end, strand, label
Followed by source-specific columns.
"""

import pandas as pd
import numpy as np
import os

DATA = '/data'
OUT = f'{DATA}/hek293t/m6a'
os.makedirs(OUT, exist_ok=True)


def write_tsv(df, name, expected_rows=None):
    path = f'{OUT}/{name}'
    df.to_csv(path, sep='\t', index=False)
    n = len(df)
    check = f" (expected {expected_rows})" if expected_rows and n != expected_rows else ""
    print(f"  {name}: {n:,} rows{check}")
    assert list(df.columns[:5]) == ['chr', 'start', 'end', 'strand', 'label'], \
        f"Bad columns in {name}: {list(df.columns[:5])}"
    if expected_rows:
        assert n == expected_rows, f"Row count mismatch in {name}: {n} != {expected_rows}"


def process_author():
    """m6ACE-seq raw from author → m6aceseq_genome.tsv"""
    print("Loading Author (m6ACE-seq)...")
    df = pd.read_excel(f'{DATA}/hek293t_m6ace/m6ace_hek293t_fromAuthor.xlsx')
    df['Chr'] = df['Chr'].str.replace('chr', '', regex=False)
    df['WT RML average'] = df[['WT rep1 RML', 'WT rep2 RML', 'WT rep3 RML']].mean(axis=1)
    df['Mettl3-KO RML average'] = df[['Mettl3-KO rep1 RML', 'Mettl3-KO rep2 RML', 'Mettl3-KO rep3 RML']].mean(axis=1)
    df['WT-KO'] = df['WT RML average'] / df['Mettl3-KO RML average']
    df['mettl3-m6a'] = (df['WT-KO'] >= 4).astype(int)

    # Standardize first 5 columns
    out = pd.DataFrame()
    out['chr'] = df['Chr']
    out['start'] = df['Start']
    out['end'] = df['End']
    out['strand'] = df['Strand']
    out['label'] = df['mettl3-m6a']
    # Append original columns (skip Chr/Start/End/Strand which are now standardized)
    for col in df.columns:
        if col not in ('Chr', 'Start', 'End', 'Strand', 'mettl3-m6a'):
            out[col] = df[col].values

    write_tsv(out, 'm6aceseq_genome.tsv', expected_rows=48027)


def process_xpore():
    """Paper/xPore → xpore_genome.tsv"""
    print("Loading Paper/xPore...")
    df = pd.read_excel(f'{DATA}/hek293t_m6ace/m6ace_hek293t_fromPaper.xlsx')
    df['Chr'] = df['Chr'].str.replace('chr', '', regex=False)
    df['WT RML average'] = df[['WT rep1 RML', 'WT rep2 RML', 'WT rep3 RML']].mean(axis=1)
    df['Mettl3-KO RML average'] = df[['Mettl3-KO rep1 RML', 'Mettl3-KO rep2 RML', 'Mettl3-KO rep3 RML']].mean(axis=1)
    df['WT-KO'] = df['WT RML average'] / df['Mettl3-KO RML average']
    df['mettl3-m6a'] = (df['WT-KO'] >= 4).astype(int)

    out = pd.DataFrame()
    out['chr'] = df['Chr']
    out['start'] = df['Start']
    out['end'] = df['End']
    out['strand'] = df['Strand']
    out['label'] = df['mettl3-m6a']
    for col in df.columns:
        if col not in ('Chr', 'Start', 'End', 'Strand', 'mettl3-m6a'):
            out[col] = df[col].values

    write_tsv(out, 'xpore_genome.tsv', expected_rows=15073)


def process_m6anet():
    """m6anet → m6anet_transcriptome.tsv (transcript coordinates)"""
    print("Loading m6anet...")
    df = pd.read_csv(f'{DATA}/hek293t_m6ace/m6ace_hek293t_fromm6anet.csv')
    # Columns: gene_id, genomic_position, transcript_id, transcript_position,
    #          HCT116_model, n_reads, y_true, HEK293T_model, Arabidopsis_model

    out = pd.DataFrame()
    out['transcript_id'] = df['transcript_id']
    out['start'] = df['transcript_position']
    out['end'] = df['transcript_position'] + 1
    out['strand'] = '.'  # unknown from transcript coords
    out['label'] = df['y_true'].astype(int)
    for col in df.columns:
        if col not in ('transcript_id', 'transcript_position', 'y_true'):
            out[col] = df[col].values

    # Rename first col to match standard (chr equivalent for transcriptome)
    out.rename(columns={'transcript_id': 'chr'}, inplace=True)

    write_tsv(out, 'm6anet_transcriptome.tsv', expected_rows=15790)


def process_glory():
    """GLORY → glory_genome.tsv"""
    print("Loading GLORY...")
    df = pd.read_excel(f'{DATA}/HEK293T/glory.xlsx')
    # Columns include: Chr, Strand, Sites, Gene, AGCov_rep1/2, m6A_level_rep1/2, Cluster_info

    out = pd.DataFrame()
    out['chr'] = df['Chr'].astype(str).str.replace('chr', '', regex=False)
    out['start'] = df['Sites'].astype(int) - 1
    out['end'] = df['Sites'].astype(int)
    out['strand'] = df['Strand']
    out['label'] = pd.NA  # GLORY has no m6A label
    for col in df.columns:
        if col not in ('Chr', 'Strand', 'Sites'):
            out[col] = df[col].values

    write_tsv(out, 'glory_genome.tsv', expected_rows=170240)


def process_gt(filename, outname, expected_rows):
    """GT files → gt_*.tsv"""
    print(f"Loading {filename}...")
    df = pd.read_csv(f'{DATA}/hek293t_m6ace/{filename}')

    out = pd.DataFrame()
    out['chr'] = df['Chr']
    out['start'] = df['Start']
    out['end'] = df['End']
    out['strand'] = df['Strand']
    out['label'] = df['mettl3-m6a']
    for col in df.columns:
        if col not in ('Chr', 'Start', 'End', 'Strand', 'mettl3-m6a'):
            out[col] = df[col].values

    write_tsv(out, outname, expected_rows=expected_rows)


def main():
    print(f"Output directory: {OUT}\n")

    process_author()
    process_xpore()
    process_m6anet()
    process_glory()

    process_gt('m6ace_hek293t_GT.csv', 'gt_genome.tsv', 44186)
    process_gt('m6ace_hek293t_GT_positive.csv', 'gt_positive_genome.tsv', 13768)
    process_gt('m6ace_hek293t_GT_negative.csv', 'gt_negative_genome.tsv', 30418)
    process_gt('m6ace_hek293t_GT_extend.csv', 'gt_extend_genome.tsv', 48802)

    print(f"\nDone! All files written to {OUT}")


if __name__ == '__main__':
    main()
