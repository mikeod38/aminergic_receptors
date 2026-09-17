# Tree bundle: v1_163taxa

Built 2026-09-16.  First full panel; still contains ador-1 and pcdr-1, no boundary test. Superseded by b120.

| | |
|---|---|
| tree | `tree.newick` |
| source | `output/tree/trees/gt0.1.treefile` |
| sha256 (first 16) | `96cf800058479c57` |
| tips | 163 |
| model | LG+F+R7 |
| alignment | `output/tree/aln/panel.gt0.1.faa` |

## Annotation datasets, validated against THIS tree

| file | rows | IDs matching a tip | status |
|---|---|---|---|
| `ann_itol_coupling.txt` | 52 | 100 % | ok |
| `ann_itol_ligand.txt` | 71 | 100 % | ok |
| `ann_itol_res332.txt` | 163 | 100 % | ok |
| `ann_itol_testset.txt` | 4 | 100 % | ok |
| `ann.figtree.nexus` | 163 tips annotated | 100 % | ok |

## Use

**iTOL** — upload `tree.newick`, then drag the `ann_itol_*.txt` files onto the rendered tree.
**FigTree** — open `ann.figtree.nexus`; annotations appear on the Tip Labels panel.

Do not mix these files with another tree. Every dataset here was generated from the tree whose
sha256 is recorded above; loading them onto a different taxon set fails silently, showing an
empty or partial overlay with no error.

