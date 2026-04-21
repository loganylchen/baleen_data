#!/usr/bin/env python3
"""
Per-(cell_line, modification) overlap analysis across ground-truth sources.

For every `{cell_line}/{mod}/` directory:
  - computes site counts per source
  - (if ≥2 sources) computes pairwise overlap counts and Jaccard
  - writes PNG figures: site-count bar chart, Jaccard heatmap, UpSet plot
  - writes/updates `{cell_line}/{mod}/readme.md` with counts, descriptions,
    overlap tables, and embedded figure links

Figures go to `output/overlap/{cell}_{mod}_{kind}.png`.
Combined markdown summary → `output/overlap_summary.md`.
"""

import os
from pathlib import Path
from itertools import combinations
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from upsetplot import UpSet

ROOT = Path(__file__).resolve().parent
OUT_MD = ROOT / 'overlap_summary.md'  # repo-root summary (tracked)
# Per-mod figures live inside each `{cell}/{mod}/figures/` directory.

CELL_LINES = ['hek293t', 'hela', 'a549']
MODS = ['m6a', 'm5c', 'm7g', 'm1a', 'nm', 'pseudo-u']

KEY_COLS = ['chr', 'start', 'end', 'strand']

# --- metadata per source ---------------------------------------------------
SOURCE_META = {
    # hek293t/m6a
    ('hek293t', 'm6a', 'm6aceseq'):   ('m6ACE-seq (raw, author supplementary)', 'WT/KO RML ratio ≥ 4 → 1 else 0'),
    ('hek293t', 'm6a', 'xpore'):      ('Paper supplementary 4 / xPore zenodo', 'WT/KO RML ratio ≥ 4 → 1 else 0'),
    ('hek293t', 'm6a', 'glory'):      ('GLORY cluster annotations', 'NA (positive-only)'),
    ('hek293t', 'm6a', 'gt'):         ('Consensus ground-truth derived from m6ACE-seq + xPore', 'mettl3-m6a (0/1)'),
    ('hek293t', 'm6a', 'gt_positive'):('gt positives only', 'mettl3-m6a = 1'),
    ('hek293t', 'm6a', 'gt_negative'):('gt negatives only', 'mettl3-m6a = 0'),
    ('hek293t', 'm6a', 'gt_extend'):  ('Extended ground-truth (gt ∪ m6ACE-seq)', 'mettl3-m6a (0/1)'),
    ('hek293t', 'm6a', 'rmbase'):     ('RMBase v3 HEK293T-filtered m6A sites', 'NA (positive-only)'),
    # hek293t/m5c
    ('hek293t', 'm5c', 'gse122254'):  ('GSE122254 HEK293T m5C', '`if m5C site` → 0/1'),
    ('hek293t', 'm5c', 'gse225614'):  ('GSE225614 HEK293T-WT sites (rRNA dropped)', 'NA (positive-only)'),
    ('hek293t', 'm5c', 'rmbase'):     ('RMBase v3 HEK293T-filtered m5C sites', 'NA (positive-only)'),
    # hek293t/m7g
    ('hek293t', 'm7g', 'rmbase'):     ('RMBase v3 HEK293T-filtered m7G sites', 'NA (positive-only)'),
    # hek293t/m1a
    ('hek293t', 'm1a', 'rmbase'):     ('RMBase v3 HEK293T-filtered m1A sites', 'NA (positive-only)'),
    # hek293t/pseudo-u
    ('hek293t', 'pseudo-u', 'bidseq'):('GSE179798 HEK293T mRNA WT BID-seq', 'NA (positive-only)'),
    # hela/m6a
    ('hela', 'm6a', 'ftom_ftop'):     ('GSE211303 FTOm vs FTO+ control, rep1', 'NA (called positives)'),
    ('hela', 'm6a', 'ftom_ivt'):      ('GSE211303 FTOm vs IVT control, rep1', 'NA (called positives)'),
    ('hela', 'm6a', 'rmbase'):        ('RMBase v3 HeLa-filtered m6A sites', 'NA (positive-only)'),
    # hela/m5c
    ('hela', 'm5c', 'gse140995'):     ('GSE140995 transcriptome-wide m5C', 'NA (positive-only)'),
    ('hela', 'm5c', 'gse225614'):     ('GSE225614 HeLa-WT sites (rRNA dropped)', 'NA (positive-only)'),
    ('hela', 'm5c', 'gse93749'):      ('GSE93749 (hg19 → hg38 via pyliftover)', 'siCTRL Status_rep1: "1"→1, "/"→0'),
    ('hela', 'm5c', 'rmbase'):        ('RMBase v3 HeLa-filtered m5C sites', 'NA (positive-only)'),
    # hela/m7g, nm, pseudo-u
    ('hela', 'm7g', 'rmbase'):        ('RMBase v3 HeLa-filtered m7G sites', 'NA (positive-only)'),
    ('hela', 'nm', 'rmbase'):         ('RMBase v3 HeLa-filtered 2′-O-methyl sites', 'NA (positive-only)'),
    ('hela', 'pseudo-u', 'bidseq'):   ('GSE179798 HeLa mRNA WT BID-seq', 'NA (positive-only)'),
    # a549
    ('a549', 'm6a', 'rmbase'):        ('RMBase v3 A549-filtered m6A sites', 'NA (positive-only)'),
    ('a549', 'm7g', 'rmbase'):        ('RMBase v3 A549-filtered m7G sites', 'NA (positive-only)'),
    ('a549', 'pseudo-u', 'bidseq'):   ('GSE179798 A549 mRNA WT BID-seq', 'NA (positive-only)'),
}

MOD_TITLE = {
    'm6a': 'N6-methyladenosine (m6A)',
    'm5c': '5-methylcytosine (m5C)',
    'm7g': 'N7-methylguanosine (m7G)',
    'm1a': 'N1-methyladenosine (m1A)',
    'nm': "2′-O-methylation (Nm)",
    'pseudo-u': 'Pseudouridine (Ψ)',
}


# --- helpers ---------------------------------------------------------------

def load_sites(path: Path) -> set:
    df = pd.read_csv(path, sep='\t', usecols=KEY_COLS, dtype={'chr': str})
    df = df.dropna(subset=KEY_COLS)
    return set(zip(df['chr'], df['start'].astype(int), df['end'].astype(int), df['strand']))


def md_table(rows, headers):
    out = ['| ' + ' | '.join(headers) + ' |',
           '|' + '|'.join(['---'] * len(headers)) + '|']
    for r in rows:
        out.append('| ' + ' | '.join(str(x) for x in r) + ' |')
    return '\n'.join(out)


def plot_counts(sources, cell, mod, outpath):
    names = list(sources.keys())
    counts = [len(sources[n]) for n in names]
    fig, ax = plt.subplots(figsize=(max(5, 0.8 * len(names) + 2), 4))
    bars = ax.bar(names, counts, color='steelblue')
    ax.set_ylabel('Unique sites')
    ax.set_title(f'{cell}/{mod}: site counts per source')
    ax.tick_params(axis='x', rotation=30)
    for b, c in zip(bars, counts):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height(), f'{c:,}',
                ha='center', va='bottom', fontsize=9)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()


def plot_jaccard(sources, cell, mod, outpath):
    names = list(sources.keys())
    n = len(names)
    mat = np.zeros((n, n))
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            A, B = sources[a], sources[b]
            union = len(A | B)
            mat[i, j] = len(A & B) / union if union else 0.0

    fig, ax = plt.subplots(figsize=(0.7 * n + 3, 0.7 * n + 2))
    im = ax.imshow(mat, cmap='viridis', vmin=0, vmax=1)
    ax.set_xticks(range(n)); ax.set_xticklabels(names, rotation=45, ha='right')
    ax.set_yticks(range(n)); ax.set_yticklabels(names)
    for i in range(n):
        for j in range(n):
            color = 'white' if mat[i, j] < 0.5 else 'black'
            ax.text(j, i, f'{mat[i, j]:.2f}', ha='center', va='center',
                    fontsize=8, color=color)
    ax.set_title(f'{cell}/{mod}: pairwise Jaccard')
    fig.colorbar(im, ax=ax, shrink=0.7)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()


def plot_upset(sources, cell, mod, outpath):
    names = list(sources.keys())
    all_sites = set().union(*sources.values())
    rows = [{n: s in sources[n] for n in names} for s in all_sites]
    df = pd.DataFrame(rows)
    counts = df.groupby(names).size()

    fig = plt.figure(figsize=(max(8, 0.8 * len(names) + 6), 5))
    up = UpSet(counts, show_counts=True, sort_by='cardinality',
               min_subset_size=1)
    up.plot(fig=fig)
    plt.suptitle(f'{cell}/{mod}: UpSet of source intersections', y=1.02)
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()


# --- readme generation ------------------------------------------------------

def build_readme(cell, mod, sources, pairwise_rows, fig_rel_dir):
    """Build full readme.md content for {cell}/{mod}/."""
    title = MOD_TITLE.get(mod, mod)
    lines = [
        f'# {cell}/{mod}',
        '',
        f'{cell.upper()} {title} benchmark datasets.',
        'All TSV files share standardized first 5 columns: `chr`, `start`, `end`, `strand`, `label`.',
        '',
        '## Sources',
        '',
    ]
    src_rows = []
    for name, sites in sources.items():
        desc, label = SOURCE_META.get((cell, mod, name), ('—', '—'))
        src_rows.append((f'`{name}_genome.tsv`', f'{len(sites):,}', label, desc))
    lines.append(md_table(src_rows, ['File', 'Sites', 'Label', 'Description']))
    lines.append('')

    # figures
    if len(sources) >= 1:
        lines += ['## Figures', '',
                  f'![counts]({fig_rel_dir}/counts.png)',
                  '']
    if len(sources) >= 2:
        lines += [f'![jaccard]({fig_rel_dir}/jaccard.png)',
                  '',
                  f'![upset]({fig_rel_dir}/upset.png)',
                  '']

    # overlap table
    if pairwise_rows:
        lines += ['## Pairwise overlap',
                  '',
                  'Site key: `(chr, start, end, strand)`. Jaccard = |A ∩ B| / |A ∪ B|.',
                  '',
                  md_table(pairwise_rows,
                           ['A', 'B', '|A∩B|', 'Jaccard', '|A∩B|/|A|', '|A∩B|/|B|']),
                  '']
    else:
        lines += ['_Only one source — no pairwise overlap._', '']

    # regeneration hint
    lines += ['## Regenerating',
              '',
              '```bash',
              'python analyze_overlap.py   # from repo root',
              '```',
              '']
    return '\n'.join(lines)


# --- main ------------------------------------------------------------------

def main():
    md_parts = ['# Ground-truth overlap summary', '',
                'Site key: `(chr, start, end, strand)`. Jaccard = |A ∩ B| / |A ∪ B|.',
                '']

    for cell in CELL_LINES:
        for mod in MODS:
            d = ROOT / cell / mod
            if not d.is_dir():
                continue
            files = sorted(d.glob('*_genome.tsv'))
            if not files:
                continue

            header = f'{cell}/{mod}'
            print('\n' + '=' * 70)
            print(header)
            print('=' * 70)
            md_parts += [f'## {header}', '']

            sources = {}
            for f in files:
                name = f.stem.replace('_genome', '')
                sources[name] = load_sites(f)
                print(f'  {name:20s} {len(sources[name]):>10,} sites')

            # figures → {cell}/{mod}/figures/
            fig_dir = d / 'figures'
            fig_dir.mkdir(exist_ok=True)
            fig_rel = 'figures'
            plot_counts(sources, cell, mod, fig_dir / 'counts.png')
            if len(sources) >= 2:
                plot_jaccard(sources, cell, mod, fig_dir / 'jaccard.png')
                plot_upset(sources, cell, mod, fig_dir / 'upset.png')

            # pairwise stats
            pairwise_rows = []
            for a, b in combinations(sources.keys(), 2):
                A, B = sources[a], sources[b]
                inter = len(A & B)
                union = len(A | B)
                jac = inter / union if union else 0.0
                pa = inter / len(A) if A else 0.0
                pb = inter / len(B) if B else 0.0
                pairwise_rows.append((a, b, f'{inter:,}', f'{jac:.4f}',
                                       f'{pa:.3f}', f'{pb:.3f}'))

            # per-dir readme
            readme = build_readme(cell, mod, sources, pairwise_rows, fig_rel)
            (d / 'readme.md').write_text(readme)
            print(f'  wrote {d / "readme.md"}')

            # summary md
            count_rows = [(n, f'{len(s):,}') for n, s in sources.items()]
            md_parts.append(md_table(count_rows, ['source', 'sites']))
            md_parts.append('')
            if pairwise_rows:
                md_parts.append(md_table(pairwise_rows,
                    ['A', 'B', '|A∩B|', 'Jaccard', '|A∩B|/|A|', '|A∩B|/|B|']))
                md_parts.append('')

    OUT_MD.write_text('\n'.join(md_parts) + '\n')
    print(f'\nSummary  → {OUT_MD.relative_to(ROOT)}')
    print('Figures  → {cell}/{mod}/figures/')


if __name__ == '__main__':
    main()
