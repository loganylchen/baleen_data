# baleen_data

Benchmark ground truth datasets for m6A modification detection.

## Data Layout

```
hek293t/m6a/
├── m6aceseq_genome.tsv        # m6ACE-seq raw from author (48,027 sites)
├── xpore_genome.tsv           # Paper supp4 / xPore zenodo (15,073 sites)
├── m6anet_transcriptome.tsv   # m6anet ML predictions, transcript coords (15,790 sites)
├── glory_genome.tsv           # GLORY (170,240 sites)
├── gt_genome.tsv              # Consensus ground truth (44,186 sites)
├── gt_positive_genome.tsv     # GT positive only (13,768 sites)
├── gt_negative_genome.tsv     # GT negative only (30,418 sites)
└── gt_extend_genome.tsv       # GT extended (48,802 sites)
```

All TSVs share a standardized first 5 columns: `chr`, `start`, `end`, `strand`, `label`, followed by source-specific columns.

### Label definitions

| File | label |
|---|---|
| m6aceseq | `1` if WT/KO RML ratio >= 4, else `0` |
| xpore | `1` if WT/KO RML ratio >= 4, else `0` |
| m6anet | `1` if y_true=True, else `0` |
| glory | NA (no m6A label) |
| gt_* | `mettl3-m6a` (1=positive, 0=negative) |

## Regenerating

```bash
docker build -f Dockerfile.analysis -t baleen-analysis .
docker run --rm -u 1000:1000 -v "$(pwd):/data" baleen-analysis python /data/reorganize.py
```

## Analysis

```bash
docker run --rm -u 1000:1000 -v "$(pwd):/data" baleen-analysis python /data/analyze_groundtruth.py
```

Output plots are saved to `output/`.
