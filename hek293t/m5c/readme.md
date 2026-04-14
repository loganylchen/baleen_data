# hek293t/m5c

HEK293T m5C modification benchmark datasets. All TSV files share standardized first 5 columns: `chr`, `start`, `end`, `strand`, `label`.

## Source Data

| File | Sites | Coordinates | Label | Description |
|---|---|---|---|---|
| `gse122254_genome.tsv` | 260 | Genomic | 1 (all positive) | GSE122254 bisulfite sequencing m5C sites |
| `gse225614_genome.tsv` | ~2,400 | Genomic | NA | GSE225614 HEK293T-WT sites (rRNA rows dropped) |

## Raw Files

- `gse122254.tsv` — 1-based positions, no chr prefix. All rows TRUE for m5C.
- `GSE225614_HEK293T-WT_sites.tsv` — Rows with chromosome="." are rRNA (dropped).
