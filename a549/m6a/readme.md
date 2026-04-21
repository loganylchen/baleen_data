# a549/m6a

A549 N6-methyladenosine (m6A) benchmark datasets.
All TSV files share standardized first 5 columns: `chr`, `start`, `end`, `strand`, `label`.

## Sources

| File | Sites | Label | Description |
|---|---|---|---|
| `rmbase_genome.tsv` | 182,560 | NA (positive-only) | RMBase v3 A549-filtered m6A sites |

## Figures

![counts](figures/counts.png)

_Only one source — no pairwise overlap._

## Regenerating

```bash
python analyze_overlap.py   # from repo root
```
