# hek293t/pseudo-u

HEK293T Pseudouridine (Ψ) benchmark datasets.
All TSV files share standardized first 5 columns: `chr`, `start`, `end`, `strand`, `label`.

## Sources

| File | Sites | Label | Description |
|---|---|---|---|
| `bidseq_genome.tsv` | 543 | NA (positive-only) | GSE179798 HEK293T mRNA WT BID-seq |

## Figures

![counts](figures/counts.png)

_Only one source — no pairwise overlap._

## Regenerating

```bash
python analyze_overlap.py   # from repo root
```
