# hek293t/m7g

HEK293T N7-methylguanosine (m7G) benchmark datasets.
All TSV files share standardized first 5 columns: `chr`, `start`, `end`, `strand`, `label`.

## Sources

| File | Sites | Label | Description |
|---|---|---|---|
| `rmbase_genome.tsv` | 21 | NA (positive-only) | RMBase v3 HEK293T-filtered m7G sites |

## Figures

![counts](figures/counts.png)

_Only one source — no pairwise overlap._

## Regenerating

```bash
python analyze_overlap.py   # from repo root
```
