#!/usr/bin/env python3
"""
Analyze HEK293T m6A benchmark ground truth datasets.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_memberships
from collections import Counter
from itertools import combinations
import warnings, os
warnings.filterwarnings('ignore')

OUTPUT_DIR = '/data/output'
os.makedirs(OUTPUT_DIR, exist_ok=True)


def make_site_key(df, chr_col='Chr', start_col='Start', end_col='End', strand_col='Strand'):
    return set(
        df[chr_col].astype(str) + ':' + df[start_col].astype(str) + ':' +
        df[end_col].astype(str) + ':' + df[strand_col].astype(str)
    )


def make_upset(site_sets, title, filename):
    """Create UpSet plot from dict of {name: set_of_site_keys}."""
    names = list(site_sets.keys())
    all_sites = set()
    for s in site_sets.values():
        all_sites |= s

    # Build boolean membership matrix
    rows = []
    for site in all_sites:
        rows.append({n: site in site_sets[n] for n in names})

    df = pd.DataFrame(rows)
    # Group by membership pattern and count
    counts = df.groupby(names).size()

    fig = plt.figure(figsize=(16, 8))
    upset = UpSet(counts, show_counts=True, sort_by='cardinality', show_percentages=True)
    upset.plot(fig=fig)
    plt.suptitle(title, fontsize=14, y=1.02)
    plt.savefig(f'{OUTPUT_DIR}/{filename}', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {OUTPUT_DIR}/{filename}")


def load_datasets():
    datasets = {}

    print("Loading Author (m6ACE-seq raw)...")
    author_df = pd.read_excel('/data/hek293t_m6ace/m6ace_hek293t_fromAuthor.xlsx')
    author_df['Chr'] = author_df['Chr'].str.replace('chr', '')
    author_df['WT RML average'] = author_df[['WT rep1 RML', 'WT rep2 RML', 'WT rep3 RML']].mean(axis=1)
    author_df['Mettl3-KO RML average'] = author_df[['Mettl3-KO rep1 RML', 'Mettl3-KO rep2 RML', 'Mettl3-KO rep3 RML']].mean(axis=1)
    author_df['WT-KO'] = author_df['WT RML average'] / author_df['Mettl3-KO RML average']
    author_df['mettl3-m6a'] = 0
    author_df.loc[author_df['WT-KO'] >= 4, 'mettl3-m6a'] = 1
    datasets['Author_m6ACE'] = author_df

    print("Loading Paper/xPore...")
    paper_df = pd.read_excel('/data/hek293t_m6ace/m6ace_hek293t_fromPaper.xlsx')
    paper_df['Chr'] = paper_df['Chr'].str.replace('chr', '')
    paper_df['WT RML average'] = paper_df[['WT rep1 RML', 'WT rep2 RML', 'WT rep3 RML']].mean(axis=1)
    paper_df['Mettl3-KO RML average'] = paper_df[['Mettl3-KO rep1 RML', 'Mettl3-KO rep2 RML', 'Mettl3-KO rep3 RML']].mean(axis=1)
    paper_df['WT-KO'] = paper_df['WT RML average'] / paper_df['Mettl3-KO RML average']
    paper_df['mettl3-m6a'] = 0
    paper_df.loc[paper_df['WT-KO'] >= 4, 'mettl3-m6a'] = 1
    datasets['Paper_xPore'] = paper_df

    print("Loading m6anet...")
    datasets['m6anet'] = pd.read_csv('/data/hek293t_m6ace/m6ace_hek293t_fromm6anet.csv')

    print("Loading GLORY...")
    datasets['GLORY'] = pd.read_excel('/data/HEK293T/glory.xlsx')

    print("Loading M6A.csv (Author+GLORY merged)...")
    datasets['M6A_merged'] = pd.read_csv('/data/HEK293T/M6A.csv')

    print("Loading derived GT files...")
    datasets['GT'] = pd.read_csv('/data/hek293t_m6ace/m6ace_hek293t_GT.csv')
    datasets['GT_positive'] = pd.read_csv('/data/hek293t_m6ace/m6ace_hek293t_GT_positive.csv')
    datasets['GT_negative'] = pd.read_csv('/data/hek293t_m6ace/m6ace_hek293t_GT_negative.csv')
    datasets['GT_extend'] = pd.read_csv('/data/hek293t_m6ace/m6ace_hek293t_GT_extend.csv')

    return datasets


def print_summary(datasets):
    print("\n" + "=" * 80)
    print("SUMMARY OF HEK293T m6A BENCHMARK GROUND TRUTH DATASETS")
    print("=" * 80)

    print("\n--- SOURCE DATASETS ---\n")

    df = datasets['Author_m6ACE']
    n_pos = (df['mettl3-m6a'] == 1).sum()
    n_neg = (df['mettl3-m6a'] == 0).sum()
    print(f"1. Author (m6ACE-seq raw from author)")
    print(f"   Total sites: {len(df):,}  |  Positive (WT-KO>=4): {n_pos:,}  |  Negative: {n_neg:,}")

    df = datasets['Paper_xPore']
    n_pos = (df['mettl3-m6a'] == 1).sum()
    print(f"\n2. Paper/xPore (supplementary data 4 = xPore zenodo)")
    print(f"   Total sites: {len(df):,}  |  ALL positive (WT-KO>=4): {n_pos:,}")
    print(f"   Note: Paper and xPore are IDENTICAL (15,073 sites)")

    df = datasets['m6anet']
    n_true = (df['y_true'] == True).sum()
    n_false = (df['y_true'] == False).sum()
    print(f"\n3. m6anet (ML-predicted m6A sites)")
    print(f"   Total sites: {len(df):,}  |  y_true=True: {n_true:,}  |  y_true=False: {n_false:,}")
    print(f"   Uses transcript coordinates (not genomic)")

    df = datasets['GLORY']
    print(f"\n4. GLORY")
    print(f"   Total sites: {len(df):,}")
    cluster_counts = df['Cluster_info'].value_counts()
    for k, v in cluster_counts.items():
        print(f"   Cluster_info={k}: {v:,}")

    df = datasets['M6A_merged']
    n_pos = (df['mettl3-m6a'] == 1).sum()
    n_neg = (df['mettl3-m6a'] == 0).sum()
    print(f"\n5. M6A.csv (Author + GLORY merged)")
    print(f"   Total sites: {len(df):,}  |  Positive: {n_pos:,}  |  Negative: {n_neg:,}")
    n_with_glory = df.dropna(subset=['Cluster_info']).shape[0]
    print(f"   Sites with GLORY annotation: {n_with_glory:,}")

    print("\n--- DERIVED GROUND TRUTH FILES ---\n")

    for name, label in [('GT', 'consensus'), ('GT_positive', 'positive only'),
                         ('GT_negative', 'negative only'), ('GT_extend', 'extended')]:
        df = datasets[name]
        if 'mettl3-m6a' in df.columns:
            n_pos = (df['mettl3-m6a'] == 1.0).sum()
            n_neg = (df['mettl3-m6a'] == 0.0).sum()
            print(f"  {name} ({label}): {len(df):,} total  |  pos={n_pos:,}  |  neg={n_neg:,}")
        else:
            print(f"  {name} ({label}): {len(df):,} total")


def compute_overlaps_and_plot(datasets):
    # ============================================================
    # Build all site key sets
    # ============================================================
    author_sites = make_site_key(datasets['Author_m6ACE'])
    paper_sites = make_site_key(datasets['Paper_xPore'])
    m6a_merged_sites = make_site_key(datasets['M6A_merged'])

    # GLORY: uses 'chr1' format + Sites (1-based position = End in other datasets)
    glory_df = datasets['GLORY'].copy()
    glory_df['Chr'] = glory_df['Chr'].astype(str).str.replace('chr', '', regex=False)
    glory_df['Start'] = glory_df['Sites'].astype(int) - 1
    glory_df['End'] = glory_df['Sites'].astype(int)

    # Debug: build one key manually for the known matching site
    m6a_with_sites = datasets['M6A_merged'].dropna(subset=['Sites'])
    if len(m6a_with_sites) > 0:
        row = m6a_with_sites.iloc[0]
        m6a_key = f"{row['Chr']}:{int(row['Start'])}:{int(row['End'])}:{row['Strand']}"
        glory_match = glory_df[glory_df['End'] == int(row['Sites'])]
        if len(glory_match) > 0:
            gr = glory_match.iloc[0]
            glory_key = f"{gr['Chr']}:{gr['Start']}:{gr['End']}:{gr['Strand']}"
            print(f"\nDebug key comparison:")
            print(f"  M6A_merged key: '{m6a_key}' (types: Chr={type(row['Chr'])}, Start={type(row['Start'])})")
            print(f"  GLORY key:      '{glory_key}' (types: Chr={type(gr['Chr'])}, Start={type(gr['Start'])})")
            print(f"  Match: {m6a_key == glory_key}")

    glory_sites = make_site_key(glory_df)
    print(f"GLORY unique site keys: {len(glory_sites):,}")
    print(f"GLORY ∩ M6A_merged: {len(glory_sites & m6a_merged_sites):,}")

    # ============================================================
    # PLOT 1: All source datasets (Author, Paper/xPore, GLORY, M6A_merged)
    # ============================================================
    print("\n" + "=" * 80)
    print("OVERLAP ANALYSIS")
    print("=" * 80)

    site_sets = {
        'Author_m6ACE': author_sites,
        'Paper_xPore': paper_sites,
        'M6A_merged': m6a_merged_sites,
        'GLORY': glory_sites,
    }

    print("\n--- Dataset Sizes ---")
    for name, sites in site_sets.items():
        print(f"  {name}: {len(sites):,}")

    print("\n--- Pairwise Overlaps ---")
    for (n1, s1), (n2, s2) in combinations(site_sets.items(), 2):
        overlap = len(s1 & s2)
        pct1 = overlap / len(s1) * 100 if len(s1) > 0 else 0
        pct2 = overlap / len(s2) * 100 if len(s2) > 0 else 0
        print(f"  {n1} ∩ {n2}: {overlap:,}  ({pct1:.1f}% / {pct2:.1f}%)")

    make_upset(site_sets,
               'HEK293T m6A: All Sites Overlap (Genomic Coordinates)',
               'upset_all_sites.png')

    # ============================================================
    # PLOT 2: Positive sites only (Author_pos, Paper_pos, M6A_merged_pos, GLORY)
    # ============================================================
    author_pos = make_site_key(datasets['Author_m6ACE'][datasets['Author_m6ACE']['mettl3-m6a'] == 1])
    paper_pos = make_site_key(datasets['Paper_xPore'][datasets['Paper_xPore']['mettl3-m6a'] == 1])
    m6a_merged_pos = make_site_key(datasets['M6A_merged'][datasets['M6A_merged']['mettl3-m6a'] == 1])

    pos_sets = {
        'Author_pos': author_pos,
        'Paper_pos': paper_pos,
        'M6A_merged_pos': m6a_merged_pos,
        'GLORY_all': glory_sites,
    }

    print("\n--- Positive Sites ---")
    for name, sites in pos_sets.items():
        print(f"  {name}: {len(sites):,}")

    make_upset(pos_sets,
               'HEK293T m6A: Positive Sites Overlap',
               'upset_positive_sites.png')

    # ============================================================
    # PLOT 3: Source positive vs derived GT
    # ============================================================
    gt_pos_sites = make_site_key(datasets['GT_positive'])
    gt_neg_sites = make_site_key(datasets['GT_negative'])

    compare_sets = {
        'Author_pos': author_pos,
        'Paper_pos': paper_pos,
        'GT_positive': gt_pos_sites,
        'GT_negative': gt_neg_sites,
    }

    print("\n--- Source vs Derived GT ---")
    for name, sites in compare_sets.items():
        print(f"  {name}: {len(sites):,}")

    make_upset(compare_sets,
               'HEK293T m6A: Source Positive vs Derived GT',
               'upset_source_vs_gt.png')

    # ============================================================
    # PLOT 4: Positive sites with GLORY cluster info
    # ============================================================
    m6a_with_glory = datasets['M6A_merged'].dropna(subset=['Cluster_info'])
    m6a_glory_pos = m6a_with_glory[m6a_with_glory['mettl3-m6a'] == 1]
    m6a_glory_neg = m6a_with_glory[m6a_with_glory['mettl3-m6a'] == 0]
    m6a_cluster = m6a_glory_pos[m6a_glory_pos['Cluster_info'] != 'Non-cluster']
    m6a_noncluster = m6a_glory_pos[m6a_glory_pos['Cluster_info'] == 'Non-cluster']

    cluster_sets = {}
    cluster_sets['Author_pos'] = author_pos
    cluster_sets['Paper_pos'] = paper_pos
    if len(m6a_cluster) > 0:
        cluster_sets['GLORY_Cluster'] = make_site_key(m6a_cluster)
    if len(m6a_noncluster) > 0:
        cluster_sets['GLORY_NonCluster'] = make_site_key(m6a_noncluster)

    print("\n--- Positive Sites + GLORY Cluster ---")
    for name, sites in cluster_sets.items():
        print(f"  {name}: {len(sites):,}")

    make_upset(cluster_sets,
               'HEK293T m6A: Positive Sites + GLORY Cluster Status',
               'upset_glory_cluster.png')

    # ============================================================
    # GT file relationships
    # ============================================================
    gt_sites = make_site_key(datasets['GT'])
    gt_ext_sites = make_site_key(datasets['GT_extend'])

    print("\n--- GT File Relationships ---")
    print(f"  GT: {len(gt_sites):,}")
    print(f"  GT_positive: {len(gt_pos_sites):,}")
    print(f"  GT_negative: {len(gt_neg_sites):,}")
    print(f"  GT_extend: {len(gt_ext_sites):,}")
    print(f"  GT_positive ∪ GT_negative = {len(gt_pos_sites | gt_neg_sites):,} (should = GT: {len(gt_sites):,})")
    print(f"  GT_positive ∩ GT_negative = {len(gt_pos_sites & gt_neg_sites):,} (should = 0)")
    print(f"  GT ⊂ GT_extend: {gt_sites.issubset(gt_ext_sites)}")
    print(f"  GT_extend \\ GT = {len(gt_ext_sites - gt_sites):,} extra sites")


def main():
    datasets = load_datasets()
    print_summary(datasets)
    compute_overlaps_and_plot(datasets)
    print("\n" + "=" * 80)
    print(f"DONE! All plots saved to {OUTPUT_DIR}")
    print("=" * 80)


if __name__ == '__main__':
    main()
