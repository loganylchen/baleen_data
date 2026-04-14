#!/usr/bin/env python3
"""
Filter RMBase v3 BED files by cell line and produce standardized *_genome.tsv.

RMBase col29 BED columns (0-indexed):
  0: chr, 1: start, 2: end, 3: site_id, 4: score, 5: strand,
  6: mod_type, 7: n_supporting, 8: GSE, 9: GSM, 10: PMID,
  11: cell_line (comma-separated), 12: method,
  13: gene_id, 14: transcript_id, 15: gene_name, 16: gene_type,
  17: region, 18: sequence_context, 19: score_value,
  20-28: additional (mostly na)

Output: chr/start/end/strand/label + source columns.
"""

import pandas as pd
import os

DATA = '/home/logan/Projects/baleen_data'
RMBASE = f'{DATA}/rmbase'

COLNAMES = [
    'chr', 'start', 'end', 'site_id', 'score', 'strand',
    'mod_type', 'n_supporting', 'GSE', 'GSM', 'PMID',
    'cell_line', 'method',
    'gene_id', 'transcript_id', 'gene_name', 'gene_type',
    'region', 'sequence_context', 'score_value',
    'col21', 'col22', 'col23', 'col24', 'col25', 'col26', 'col27', 'col28', 'col29',
]

# Map rmbase filenames → (modification_name, directory_suffix)
RMBASE_FILES = {
    'human.hg38.m6A.result.col29.bed': ('m6A', 'm6a'),
    'human.hg38.m5C.result.col29.bed': ('m5C', 'm5c'),
    'human.hg38.m1A.result.col29.bed': ('m1A', 'm1a'),
    'human.hg38.m7G.result.col29.bed': ('m7G', 'm7g'),
    'human.hg38.Nm.result.col29.bed':  ('Nm',  'nm'),
    'human.hg38.Pseudo.result.col29.bed': ('Pseudo-U', 'pseudo-u'),
}

# Cell line matching: rmbase_token → our directory name
CELL_LINES = {
    'HEK293T': 'hek293t',
    'HeLa': 'hela',
    'A549': 'a549',
}


def cell_line_matches(cell_line_field, target):
    """Check if target cell line appears in the comma-separated cell_line field."""
    if pd.isna(cell_line_field):
        return False
    return target in str(cell_line_field).split(',')


def process_file(filename, mod_name, dir_suffix):
    """Filter one rmbase BED file by cell line and write _genome.tsv files."""
    filepath = f'{RMBASE}/{filename}'
    print(f"\n{'='*60}")
    print(f"Loading {filename} ({mod_name}) ...")

    df = pd.read_csv(filepath, sep='\t', header=None, names=COLNAMES,
                     dtype={'chr': str}, low_memory=False)
    print(f"  Total rows: {len(df):,}")

    # Strip chr prefix from chromosome
    df['chr'] = df['chr'].str.replace('chr', '', regex=False)

    for cl_token, cl_dir in CELL_LINES.items():
        mask = df['cell_line'].apply(lambda x: cell_line_matches(x, cl_token))
        subset = df[mask].copy()

        if len(subset) == 0:
            print(f"  {cl_token}: 0 rows — skipping")
            continue

        outdir = f'{DATA}/{cl_dir}/{dir_suffix}'
        os.makedirs(outdir, exist_ok=True)

        # Build standardized output — keep it lean
        out = pd.DataFrame()
        out['chr'] = subset['chr'].values
        out['start'] = subset['start'].values
        out['end'] = subset['end'].values
        out['strand'] = subset['strand'].values
        out['label'] = pd.NA  # RMBase sites are positive-only
        out['n_supporting'] = subset['n_supporting'].values
        out['method'] = subset['method'].values
        out['gene_name'] = subset['gene_name'].values

        outpath = f'{outdir}/rmbase_genome.tsv'
        out.to_csv(outpath, sep='\t', index=False)
        print(f"  {cl_token}: {len(out):,} rows → {cl_dir}/{dir_suffix}/rmbase_genome.tsv")


def main():
    print("Reorganizing RMBase v3 data by cell line")

    for filename, (mod_name, dir_suffix) in RMBASE_FILES.items():
        filepath = f'{RMBASE}/{filename}'
        if not os.path.exists(filepath):
            print(f"\nSkipping {filename} — not found")
            continue
        process_file(filename, mod_name, dir_suffix)

    print(f"\n{'='*60}")
    print("Done!")


if __name__ == '__main__':
    main()
