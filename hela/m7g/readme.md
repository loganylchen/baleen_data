# hela/m7g

HELA N7-methylguanosine (m7G) benchmark datasets.
All TSV files share standardized first 5 columns: `chr`, `start`, `end`, `strand`, `label`.

## Sources

| File | Sites | Label | Description |
|---|---|---|---|
| `rmbase_genome.tsv` | 824 | NA (positive-only) | RMBase v3 HeLa-filtered m7G sites |

## Figures

![counts](figures/counts.png)

_Only one source — no pairwise overlap._

## Regenerating

```bash
python analyze_overlap.py   # from repo root
```
