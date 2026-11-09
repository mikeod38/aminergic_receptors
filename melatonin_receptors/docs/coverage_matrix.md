# Melatonin receptors in *C. elegans* — lineage x method coverage

One row per lineage, one column per independent axis. A cell is only a negative if
the method was shown to be able to find a true positive in a comparable lineage.

| lineage | proteome complete | BLASTP RBH | MTNR profile (35 seqs) | verified-clade profile (8 seqs) | tree placement, 178 taxa | structural |
|---|---|---|---|---|---|---|
| *C. elegans* (Nematoda) | yes, 26,629/26,629 | **negative** — 0 of 8 candidates reciprocate to a melatonin receptor | **inconclusive** — ceiling 232 is inside the method's own blind spot | — | **negative** — neither PCDR-1 nor any of the 71 competing worm proteins falls inside the melatonin clade (100/100) | not run |
| *D. melanogaster* (Arthropoda) | yes | 3 RBH | held-out ortholog rank 12, below background | 134–189 band (background) | **RBH calls not supported** — none joins the melatonin clade | not run |
| *P. caudatus* (Priapulida) | yes | 2 RBH | held-out ortholog rank 3, below background | background | **RBH calls not supported** | not run |
| *N. vectensis* (Cnidaria) | yes | 8 RBH | held-out ortholog rank 2, below background | 168 (background) | **RBH calls not supported** | not run |
| *M. yessoensis*, *L. anatina* (Lophotrochozoa) | yes | 3, 2 RBH | held-out ortholog rank 1 | 148–181 (background) | **RBH calls not supported** | not run |
| *S. purpuratus*, *C. intestinalis*, *B. floridae* | yes | 5, 3, 1 RBH | held-out ortholog rank 1 | 171–185, 169, 179 (background) | **RBH calls not supported** | not run |
| *P. marinus* (lamprey, Vertebrata) | yes | 2 RBH to MTNR1A (429, 377 bits) | — | **positive, 544 and 451** — both records self-named melatonin receptors | not in panel | not run |
| *D. rerio* (Vertebrata) | yes | 5 RBH | — (reference lineage) | — (reference lineage) | **positive** — 5 of 5 inside the melatonin clade | not run |

Only 5 of the 32 invertebrate RBH calls — all of them zebrafish — survived tree
placement. That is what resolved the circularity recorded in the README: the RBH screen
was wrong wholesale, and the profile's invertebrate negative was the profile behaving
correctly rather than a blind spot.

**Not yet tested at all**: other nematodes (*C. briggsae*, *P. pacificus*, *Trichinella* —
the panel used for the dop-6 family history in `../dop6_c24a8_6` §7 already has these);
structural search (foldseek is available and the AlphaFold models are not) — the whole
right-hand column is empty, and it is the only axis independent of sequence similarity.

## Correction, 2026-11-03

The first tree panel (`output/tree/`) held only **20 *C. elegans* tips** — the 19
canonical aminergic receptors from the b120 backbone, plus PCDR-1. It therefore could
not address whether some *other* worm gene is closer to the melatonin receptors than
PCDR-1 is; no competing worm gene was in it. 71 worm proteins outscore PCDR-1 against
at least one melatonin query, and `output/tree2/` adds all of them.

"PCDR-1 family" in the first figure is also a label that claims more than it earned:
it is the single best BLASTP hit to PCDR-1 in each of ten proteomes — nearest
neighbour per species, not a curated family. Only four of the ten join PCDR-1's clade;
the zebrafish member falls inside the aminergic backbone. Renamed to "PCDR-1 nearest
neighbours".
