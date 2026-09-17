# Tree bundle: v2_b120_boundary120

Built 2026-09-16.  Boundary at monoamine HMM >= 120; ador-1/pcdr-1/F35H10.10 excluded with reasons; 12-gene boundary test included. Current.

| | |
|---|---|
| tree | `tree.newick` |
| source | `output/tree/b120/trees/gt0.1.treefile` |
| sha256 (first 16) | `90cf118cdefa5e06` |
| tips | 173 |
| model | LG+F+I+R7 |
| alignment | `output/tree/b120/panel.gt0.1.faa` |

## Annotation datasets, validated against THIS tree

| file | rows | IDs matching a tip | status |
|---|---|---|---|
| `ann_itol_coupling.txt` | 52 | 100 % | ok |
| `ann_itol_ligand.txt` | 71 | 100 % | ok |
| `ann_itol_res332.txt` | 161 | 100 % | ok |
| `ann_itol_testset.txt` | 14 | 100 % | ok |
| `ann.figtree.nexus` | 173 tips annotated | 100 % | ok |

## Use

**iTOL** — upload `tree.newick`, then drag the `ann_itol_*.txt` files onto the rendered tree.
**FigTree** — open `ann.figtree.nexus`; annotations appear on the Tip Labels panel.

Do not mix these files with another tree. Every dataset here was generated from the tree whose
sha256 is recorded above; loading them onto a different taxon set fails silently, showing an
empty or partial overlay with no error.

