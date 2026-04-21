# Ground-truth overlap summary

Site key: `(chr, start, end, strand)`. Jaccard = |A ∩ B| / |A ∪ B|.

## hek293t/m6a

| source | sites |
|---|---|
| glory | 170,240 |
| gt_extend | 48,802 |
| gt | 44,186 |
| gt_negative | 30,418 |
| gt_positive | 13,768 |
| m6aceseq | 48,027 |
| rmbase | 359,318 |
| xpore | 15,073 |

| A | B | |A∩B| | Jaccard | |A∩B|/|A| | |A∩B|/|B| |
|---|---|---|---|---|---|
| glory | gt_extend | 21,408 | 0.1083 | 0.126 | 0.439 |
| glory | gt | 18,915 | 0.0967 | 0.111 | 0.428 |
| glory | gt_negative | 7,995 | 0.0415 | 0.047 | 0.263 |
| glory | gt_positive | 10,920 | 0.0631 | 0.064 | 0.793 |
| glory | m6aceseq | 21,040 | 0.1067 | 0.124 | 0.438 |
| glory | rmbase | 79,856 | 0.1776 | 0.469 | 0.222 |
| glory | xpore | 11,640 | 0.0670 | 0.068 | 0.772 |
| gt_extend | gt | 44,186 | 0.9054 | 0.905 | 1.000 |
| gt_extend | gt_negative | 30,418 | 0.6233 | 0.623 | 1.000 |
| gt_extend | gt_positive | 13,768 | 0.2821 | 0.282 | 1.000 |
| gt_extend | m6aceseq | 48,027 | 0.9841 | 0.984 | 1.000 |
| gt_extend | rmbase | 16,496 | 0.0421 | 0.338 | 0.046 |
| gt_extend | xpore | 15,073 | 0.3089 | 0.309 | 1.000 |
| gt | gt_negative | 30,418 | 0.6884 | 0.688 | 1.000 |
| gt | gt_positive | 13,768 | 0.3116 | 0.312 | 1.000 |
| gt | m6aceseq | 44,186 | 0.9200 | 1.000 | 0.920 |
| gt | rmbase | 14,625 | 0.0376 | 0.331 | 0.041 |
| gt | xpore | 13,768 | 0.3027 | 0.312 | 0.913 |
| gt_negative | gt_positive | 0 | 0.0000 | 0.000 | 0.000 |
| gt_negative | m6aceseq | 30,418 | 0.6334 | 1.000 | 0.633 |
| gt_negative | rmbase | 6,097 | 0.0159 | 0.200 | 0.017 |
| gt_negative | xpore | 0 | 0.0000 | 0.000 | 0.000 |
| gt_positive | m6aceseq | 13,768 | 0.2867 | 1.000 | 0.287 |
| gt_positive | rmbase | 8,528 | 0.0234 | 0.619 | 0.024 |
| gt_positive | xpore | 13,768 | 0.9134 | 1.000 | 0.913 |
| m6aceseq | rmbase | 16,221 | 0.0415 | 0.338 | 0.045 |
| m6aceseq | xpore | 14,298 | 0.2930 | 0.298 | 0.949 |
| rmbase | xpore | 9,086 | 0.0249 | 0.025 | 0.603 |

## hek293t/m5c

| source | sites |
|---|---|
| gse122254 | 260 |
| gse225614 | 2,191 |
| rmbase | 258 |

| A | B | |A∩B| | Jaccard | |A∩B|/|A| | |A∩B|/|B| |
|---|---|---|---|---|---|
| gse122254 | gse225614 | 1 | 0.0004 | 0.004 | 0.000 |
| gse122254 | rmbase | 3 | 0.0058 | 0.012 | 0.012 |
| gse225614 | rmbase | 193 | 0.0855 | 0.088 | 0.748 |

## hek293t/m7g

| source | sites |
|---|---|
| rmbase | 21 |

## hek293t/m1a

| source | sites |
|---|---|
| rmbase | 1,014 |

## hek293t/pseudo-u

| source | sites |
|---|---|
| bidseq | 543 |

## hela/m6a

| source | sites |
|---|---|
| ftom_ftop | 40,096 |
| ftom_ivt | 69,834 |
| rmbase | 458,320 |

| A | B | |A∩B| | Jaccard | |A∩B|/|A| | |A∩B|/|B| |
|---|---|---|---|---|---|
| ftom_ftop | ftom_ivt | 35,352 | 0.4740 | 0.882 | 0.506 |
| ftom_ftop | rmbase | 24,731 | 0.0522 | 0.617 | 0.054 |
| ftom_ivt | rmbase | 36,397 | 0.0740 | 0.521 | 0.079 |

## hela/m5c

| source | sites |
|---|---|
| gse140995 | 899 |
| gse225614 | 2,441 |
| gse93749 | 16,662 |
| rmbase | 686 |

| A | B | |A∩B| | Jaccard | |A∩B|/|A| | |A∩B|/|B| |
|---|---|---|---|---|---|
| gse140995 | gse225614 | 458 | 0.1589 | 0.509 | 0.188 |
| gse140995 | gse93749 | 423 | 0.0247 | 0.471 | 0.025 |
| gse140995 | rmbase | 311 | 0.2441 | 0.346 | 0.453 |
| gse225614 | gse93749 | 703 | 0.0382 | 0.288 | 0.042 |
| gse225614 | rmbase | 452 | 0.1690 | 0.185 | 0.659 |
| gse93749 | rmbase | 448 | 0.0265 | 0.027 | 0.653 |

## hela/m7g

| source | sites |
|---|---|
| rmbase | 824 |

## hela/nm

| source | sites |
|---|---|
| rmbase | 2,887 |

## hela/pseudo-u

| source | sites |
|---|---|
| bidseq | 575 |

## a549/m6a

| source | sites |
|---|---|
| rmbase | 182,560 |

## a549/m7g

| source | sites |
|---|---|
| rmbase | 803 |

## a549/pseudo-u

| source | sites |
|---|---|
| bidseq | 922 |

