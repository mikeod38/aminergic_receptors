# Aminergic receptor evolution in nematodes

Two sequence-evolution sub-projects on aminergic and indoleamine GPCRs in *C. elegans* and
related nematodes. Both were split out of
[Malaiwong_Oglu_2026](https://github.com/ODonnellLab/Malaiwong_Oglu_2026), with their history,
because they are side threads rather than part of that manuscript. Commit messages and some
internal path references from before 2026-11 still refer to that repository.

```
dop6_c24a8_6/          dop-6 / C24A8.6: an inverted duplication on X, and whether DOP-6
                       is an octopamine receptor
melatonin_receptors/   is there a C. elegans melatonin receptor?
```

Each sub-project has its own README with the full argument, the data behind each claim, and the
scripts that produced it. What follows is only an index.

## `dop6_c24a8_6/`

**`C24A8.6` is not an independent gene.** It lies entirely within arm B of a perfect inverted
repeat on chromosome X — arm A 4,324,178–4,328,106 and arm B 4,329,303–4,333,231, 3,929 bp each,
**3,929/3,929 = 100.000 % identity, zero mismatches**. None of C24A8.6's 2,330 bp is unique, and
C24A8.6a (287 aa) is DOP-6 residues 406–692 exactly, 287/287 identical, with C24A8.6a's start
codon being DOP-6's own internal Met406. That segment carries no ligand pocket: the amine-binding
D3.32 (D103) and the DRL motif are both in dop-6's single-copy half.

**A second perfect palindrome abuts the first.** Arms C and D, 4,489 bp each, also at 100.000 %,
beginning one base after arm B ends — so X:4,324,178–4,342,446, **18,269 bp, is inverted repeat
almost end to end**. `cst-2` is duplicated across it, and its duplicate is `cst-1`. Both
palindromes are polymorphic in wild isolates and rise and fall together.

**The strains are not what their names imply, and this is the most consequential result here.**
`dop-6` and `C24A8.6` were never targeted independently: Pu et al.'s own supplementary data give
*the same* guide, PAM and ssODN for both, and that guide has exactly two perfect full-length
matches in WBcel235 — one in each arm, as it must, since the arms are identical. One reagent made
two double-strand breaks 6.2 kb apart in inverted orientation, so `yum2349`/`yum2350` are the two
ends of a single event rather than two alleles, and CHS1223 is a complex rearrangement carrying a
75 kb duplication, not a double null.

**Consequences.** C24A8.6 cannot be measured apart from `dop-6` by any short-read method, cannot
be mutated without also cutting `dop-6`, and has no diagnostic base anywhere. `dop-6` itself
remains fully tractable through exons 1–9.

**What this work does not establish**, recorded because earlier drafts of it overreached. DOP-6
is an orphan with no pharmacology, and **its ligand cannot be resolved from sequence** — the
binding pocket separates dopamine from serotonin in zero alignment columns, so two pocket-based
ligand predictions made earlier were withdrawn. The dN of 0.0252 is the 4th lowest of 15, not the
lowest; `octr-1` is lower. A claimed affinity with Gq-coupled 5-HT neighbours was also withdrawn.

`reviews/2026-09-16_dop6_conclusions/` is an **internal red-team audit** of this sub-project's own
conclusions, in which two adversarial reviewers re-derived every number from the treefiles and
alignments. Despite the filenames, it has nothing to do with journal peer review — it is
self-criticism, and it is where several of the withdrawals above came from.

## `melatonin_receptors/`

There is **no *C. elegans* melatonin receptor ortholog**. PCDR-1, suggested as a candidate, is a
protostome receptor outside the melatonin receptor clade.

## Reproducing

Scripts fetch from public sources (Ensembl and Ensembl Metazoa REST, UniProt, WormBase, CaeNDR)
and write into each sub-project's `output/`. Large downloaded reference data is not tracked —
see `.gitignore` for what is excluded and how to re-fetch it. Alignments use MAFFT, tree
inference IQ-TREE. Coordinates are WBcel235 unless stated otherwise.
