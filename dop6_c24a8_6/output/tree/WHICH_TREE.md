# Which tree to look at

Three panels exist. Only the last two are current.

| panel | taxa | trees in | status |
|---|---|---|---|
| v1 | 163 | `output/tree/trees/` | **SUPERSEDED.** Boundary drawn on WormBase annotation labels, which the referee panel showed was wrong — 36 excluded genes scored higher on the monoamine HMM than an included one. Figures in `output/tree/*.pdf` are from this panel. |
| **b120** | 173 | `output/tree/b120/trees/` | **CURRENT — use this for where DOP-6 sits.** Boundary at monoamine HMM ≥ 120, plus 12 neuropeptide-receptor tips that test the boundary from outside. |
| **placement** | 268 | `output/tree/placement/trees/` | **CURRENT — use this for the family's history.** b120 plus every nematode monoamine-HMM candidate ≥ 120 from 8 proteomes spanning Blaxter clades I, III, IV and V. No orthology call made; the topology assigns paralogy. |

## To view

**Interactive (collapsible, searchable, no software needed):**
- `output/tree/b120/viewer_gt0.1.html`
- `output/tree/placement/viewer_gt0.1.html`

Open in a browser. Click internal nodes to collapse/expand; support values on nodes.

**Publication figures:**
- `output/tree/b120/tree_final_gt0.1.pdf` (and `.png`)
- `output/tree/placement/tree_final_gt0.1.pdf`

Tips are coloured by experimentally determined ligand where one exists. Proteome-scan tips carry
only an accession and are deliberately left uncoloured — they are breadth, not evidence about ligand.

**iTOL:** `output/tree/bundles/v2_b120_boundary120/` holds the b120 treefile with matched
annotation files (ligand, coupling, residue 3.32, boundary test set).

## What each tree answers

**b120** — DOP-6 sits in the D2-like radiation. Nesting outward at gt0.1 (UFboot):
dop-6 orthologues 100 → +dop-5 95 → +dop-3 90 → +dop-2 90 → +hDRD2/3/4 95 → +5-HT2 group 100.
Identical topology at all four trim levels; the 13-tip clade falls to 60 at gappyout, which must be
quoted alongside the 90. DOP-6 groups with no octopamine receptor — but see README §6, an AU test
puts an OA placement at 23–28 logL worse and ~12–15 % of the likelihood weight, i.e. disfavoured,
NOT excluded.

**placement** — the family is old. dop-2 is pan-nematode (clade I); dop-6 reaches clade IV
(*Strongyloides*) in an 8-tip orthologue clade at **100/100 on all four trim levels**; dop-3 is
clade V only. The dop-5/dop-6 split predates the clade III divergence. The "Rhabditina-restricted"
reading from the Compara-built panel was a sampling artefact.
