# hek293t/m6a

HEK293T m6A modification benchmark datasets. All TSV files share standardized first 5 columns: `chr`, `start`, `end`, `strand`, `label`.

## Source Data

| File | Sites | Coordinates | Label | Description |
|---|---|---|---|---|
| `m6aceseq_genome.tsv` | 48,027 | Genomic | WT/KO ratio >= 4 | m6ACE-seq raw data from author |
| `xpore_genome.tsv` | 15,073 | Genomic | WT/KO ratio >= 4 | Paper supplementary 4 (= xPore zenodo) |
| `m6anet_transcriptome.tsv` | 15,790 | Transcriptomic | y_true | m6anet ML predictions |
| `glory_genome.tsv` | 170,240 | Genomic | NA | GLORY cluster annotations |

## Ground Truth

| File | Sites | Positive | Negative | Description |
|---|---|---|---|---|
| `gt_genome.tsv` | 44,186 | 13,768 | 30,418 | Consensus ground truth |
| `gt_positive_genome.tsv` | 13,768 | 13,768 | 0 | Positive sites only |
| `gt_negative_genome.tsv` | 30,418 | 0 | 30,418 | Negative sites only |
| `gt_extend_genome.tsv` | 48,802 | — | — | Extended ground truth |

## Column Details

**Genome files** (`chr` = chromosome without "chr" prefix, 0-based coordinates):

- **m6aceseq**: + WT/KO rep RML values, Gene, Annotation, Motif, DESeq2 padj, WT-KO ratio
- **xpore**: + WT/KO rep RML values, Gene, Annotation, Motif, DESeq2 padj, WT-KO ratio
- **glory**: + Gene, AGCov_rep1/2, m6A_level_rep1/2, Cluster_info
- **gt_***: + Motif, mettl3-m6a

**Transcriptome file** (`chr` = transcript_id, transcript-relative coordinates):

- **m6anet**: + gene_id, genomic_position, HCT116/HEK293T/Arabidopsis model scores, n_reads
