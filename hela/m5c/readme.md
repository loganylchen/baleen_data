# hela/m5c

HeLa m5C modification benchmark datasets. All TSV files share standardized first 5 columns: `chr`, `start`, `end`, `strand`, `label`.

## Source Data

| File | Sites | Coordinates | Label | Description |
|---|---|---|---|---|
| `gse225614_genome.tsv` | ~2,700 | Genomic | NA | GSE225614 HeLa-WT sites (rRNA rows dropped) |
| `gse140995_genome.tsv` | 1,034 | Genomic | NA | GSE140995 transcriptome-wide m5C sites |
| `gse93749_genome.tsv` | ~16,800 | Genomic | 0/1 | GSE93749 HeLa sites (hg19→hg38 via pyliftover) |

## Raw Files

- `GSE225614_HeLa-WT_sites.tsv` — Rows with chromosome="." are rRNA (dropped).
- `GSE140995_transcriptome-wide_sites.xlsx` — Already 0-based BED coordinates with chr prefix.
- `GSE93749_hg19_Human_m5C_sites_information.xls` — hg19 coordinates, lifted over to hg38. Label from siCTRL Status_rep1: "1"→1, "/"→0.
