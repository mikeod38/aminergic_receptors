# C35 audit — Harris et al. 2009, now that the full text is in hand

Closes the one BLOCKED item in the panel. PDF supplied by MOD 2026-09-17 at
`dop6_c24a8_6/data/1446.full.pdf`; text at `data/harris2009.txt`.
Harris GP, Hapiak VM, Wragg RT, Miller SB, Hughes LJ, Hobson RJ, Steven R,
Bamber B, Komuniecki RW. *J Neurosci* 29(5):1446–1456, 2009. PMID 19193891.

**Verdict: my critique was directionally right about the method and wrong on two
specifics, and it omitted that SER-5's assignment has independent support. Net
effect — the paper is better than I said, and my own method is closer to theirs
than I claimed.**

## Confirmed, verbatim from Methods (p. 1448) and the Fig. 5 legend

| my claim (C35) | verdict | evidence |
|---|---|---|
| 5-HT receptors only | **confirmed** | Fig. 5B legend lists 27 sequences: 12 human 5-HT, 4 *C. elegans* SER, 4 *C. briggsae* SER, 4 *Drosophila*, 1 *Aedes*, 1 *Haemonchus*, 1 *Ascaris*. No non-5-HT receptor in the panel. |
| rooted on GAR-1 | **confirmed** | "The *Caenorhabditis elegans* GAR-1 sequence was used as an outgroup." |
| ClustalW in MegAlign/DNASTAR, hand-tuned | **confirmed** | "aligned using MegAlign in DNAStar with ClustalW … and fine-tuned by hand." |
| NJ distance in PAUP 4.0b10, DNASTAR bootstrapping | **confirmed** | "Bootstrapping was performed with DNASTAR (1000 replicates, random seed), and the resulting tree was then imported into … PAUP, version 4.0b10, where it was analyzed under the distance, neighbor-joining method." |
| no ML, no model selection | **confirmed** by absence | no substitution model is named anywhere. |
| the validation is one eyeball sentence, no statistical test | **confirmed, verbatim** | "all of the previously characterized 5-HT receptors clustered according to their documented G-protein coupling, regardless of phylogenetic origin, validating the utility of the tree and supporting the identification of SER-5 as a Gαs-coupled 5-HT receptor." No test statistic accompanies it. |

## Correction 1 — they did not delete all of ICL3, and what they kept is the part that matters

My §6 said *"hypervariable regions were deleted, including the third intracellular
loop — which is where part of the Gα interface actually lives … so coupling
clustering was assessed on an alignment stripped of some of its own determinants."*

**That overstates it.** The actual operation (Methods, and Fig. 5B legend):

> "…the third intracellular loops, **8 aa after and 6 aa before predicted
> transmembrane domains 5 and 6, respectively**…"

They excised the **middle** of ICL3 and **kept 8 residues after TM5 and 6 before
TM6** — i.e. they kept the loop's ends. That is precisely where this project's own
§5 analysis found the Gα contacts sit: the consensus TM5/ICL3 contacts (I206,
L210, R213 in DOP-6 numbering) and the TM6 entry contacts. An 8-residue tail after
TM5 would capture I206 and L210 but miss R213, so the excision loses *some* Gα
signal — but "stripped of its determinants" is not a fair description of a
procedure that deliberately retained the loop ends.

**The irony is worth recording:** this project's §6b describes its own
concatenated-HMM-block excerpting as "the same operation Harris et al. 2009
performed by hand". That sentence is *more* accurate than the criticism in §6, and
the two should not have coexisted in the same document.

## Correction 2 — they tested the full-length alignment

Immediately after the validation sentence:

> "Although this tree was generated using truncated amino acid sequences with the
> more variable N and C termini and third intracellular loops deleted, **similar
> results were obtained using the full-length receptors (data not shown)**."

So the obvious control was run. It is "data not shown", which is weak by modern
standards and cannot be checked — but my critique was written as though the
question had never been asked, and that is not the case.

## Correction 3 — SER-5's Gs call does not rest on the tree alone

I presented the Gs assignment as a tree-clustering inference. The paper gives the
tree as one of **five** lines, four of which are sequence-level and independent of
topology (p. 1453):

1. 35 % identity to human 5-HT6 within the predicted ligand-binding TMs;
2. short ICL3 — 67 aa in SER-5 vs 60 aa in 5-HT6;
3. a **BBXXB** motif at the C-terminus of ICL3 (B = basic), associated with
   constitutive activity in other GPCRs;
4. a **PLRYK** motif in ICL2 "that is not found in any other 5-HT receptors";
5. **two intron positions conserved with the mammalian 5-HT6 genes**, in the
   sequences encoding ICL3 and ECL3 (Hapiak et al. 2009).

Point 5 is a discrete shared-derived character, which is stronger evidence of
common descent than any distance clustering, and I gave it no weight because I
had not read the paper. **My §6 sentence that the coupling call rests on an
eyeball assessment describes the *tree's* validation correctly but
mischaracterises the SER-5 assignment as a whole.**

## Minor: an internal inconsistency in the paper

Methods say GAR-1 was the outgroup; the Fig. 5B legend calls the tree
**"Unrooted"** and does not list GAR-1 among the accessions. One of the two is
wrong. It does not affect anything this project uses, but a reader reconstructing
the panel from the legend would not recover the outgroup.

## What survives, and what this changes for the project

**Survives.** The method is NJ on a hand-tuned ClustalW alignment with no model
selection and no ML, and the coupling-clustering claim is asserted from
inspection rather than tested. Judged as a basis for *inferring coupling from tree
position*, it remains weaker than its citation implies, and §6c's formal
replacement — Fitch parsimony against 10,000 permutations across the whole
aminergic family, which supports the intuition — is still the right way to have
settled it.

**Changes.** Two of the four criticisms in §6 are withdrawn (ICL3 deletion,
untested truncation), and the SER-5 assignment must be described as
multi-line rather than tree-based. README §6 corrected accordingly.

**New, and useful to this project.** Fig. 5A annotates SER-5's **ligand-binding
residues and G-protein-coupling residues** from a named primary literature set —
Strader 1989, Choudhary 1993, Moro 1993, Barak 1994, Almaula 1996, Roth 1997.
That is a **sequence-level, non-phylogenetic** route to exactly the two questions
still open for DOP-6 (ligand identity and coupling class), and it is the route
this project has not yet taken. Followed up in `output/motif_prediction.txt`.
