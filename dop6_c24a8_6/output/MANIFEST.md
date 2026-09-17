# Where everything is

Two folders, one of which you can ignore unless you are re-deriving the panel.

## `tree/` — the tree and everything read off it

| file | what |
|---|---|
| **`tree_final_gt0.1.png` / `.pdf`** | **the full 163-taxon figure.** PDF is vector, zoom freely |
| `tree_final_gt0.1_clade25.png` / `.pdf` | DOP-6's immediate clade, readable at normal size |
| `tree_final_gt0.1_clade150.png` / `.pdf` | the whole aminergic ingroup, muscarinic outgroup dropped |
| `panel.faa` | the 164 full-length sequences |
| `panel_domains.faa` | the 163 7TM excerpts the tree is actually built from |
| `basis.tsv` | per-taxon inclusion basis: anchor / ortholog / scan / subject |
| `aln/` | alignment + the four trimAl levels |
| `trees/` | IQ-TREE output per trim level (`.treefile`, `.iqtree`, `.log`) |
| `pocket.tsv` | residue at 3.32 and the 20-position pocket string, per taxon |
| `resolution.txt` | node support by clade depth — does the matrix resolve anything |
| `domain_report.txt` | what the 7TM excerpting kept and cut |
| `build_report.txt` | alignment/trim/model record |

## `panel_build/` — how the panel was assembled (intermediates)

| folder | what |
|---|---|
| `v2_anchor_panel/` | verified human/fly/worm reference receptors + accession provenance |
| `v3_breadth_panel/` | nematode Compara orthologs added |
| `proteome_scan/` | monoamine-HMM scan of 31 metazoan proteomes → `hits.tsv`, `candidates.faa` |

## Top level of `output/` — the non-tree analyses

Locus and duplication: `repeat_structure.txt`, `assembly_survey.tsv`, `wild_locus_structure.txt`,
`haplotype_context.tsv`, `haplotype_extent.tsv`, `wild_isolate_variants.tsv`,
`dop6_polymorphism.tsv`, `crossspecies.txt`, `ortholog_divergence.txt`,
`dnds_vs_briggsae.tsv`, `mirror_bias_check.txt`

Receptor biology: `pocket_332_v2.tsv` (3.32 calls), `gprotein_interface_multi.tsv` (Gα interface
across 11 complexes), `gprotein_survival.txt` (does excerpting keep it), `blast_reciprocity.txt`,
`icl3_information.txt`, `accession_audit.tsv`

Reagents: `dop6_guides.tsv`, `specific_tools.txt`

## Removed

`output/tree/` previously held a superseded 49-taxon tree built on a panel that contained a squid
rhodopsin, a bacterial decarboxylase and a ubiquitin ligase. It is deleted rather than kept, so
nobody reads a figure off it; it survives in git commit 65c7db2 if ever needed.
