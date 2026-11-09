# Melatonin receptors in *C. elegans* — is there any evidence of orthology?

Sub-project started 2026-11-02, extending the receptor identification work in
`../dop6_c24a8_6` (where DOP-6's placement among the D2-like dopamine receptors
was settled, and whether it could be an octopamine receptor was tested formally).

**Status (2026-11-03): no evidence that *C. elegans* has a melatonin receptor
ortholog, and specifically no support for PCDR-1 / F59D12.1 being one.** Tree
placement across 178 taxa puts no worm sequence inside the melatonin clade; PCDR-1
is not even the nearest worm gene to it; forcing it in is rejected formally, but at
exactly the cost a dopamine receptor pays, so the right statement is that PCDR-1 is
no more melatonin-related than DOP-1 is. The binding pocket carries neither the
aminergic nor the melatonin signature. None of this contradicts the 2020 functional
result that melatonin acts through PCDR-1 — ligand and ancestry are separate claims.

**PCDR-1's actual family is a PROTOSTOME one, and the human-orphan version of that
claim is formally rejected.** In a 258-taxon panel that samples the family properly
(rather than one member per species), PCDR-1 sits in an 8-member protostome clade at
99.7/100 and 99.8/100 across two trim levels: PCDR-1, *Drosophila* DopEcR, *Priapulus*,
oyster, brachiopod and three *Limulus* receptors. Forcing that clade costs **zero** log
units. Forcing the earlier published-style claim - PCDR-1 together with human GPR52 and
GPR21 - costs **856** (p-AU 1.4e-90), against **0.58** for the melatonin clade that
independent evidence says is real. See the sections below; the 2026-11-03 claim that
this family runs from cnidarian to human was wrong and is retracted there.

**The melatonin family reaches the base of vertebrates and no further detectably**:
lamprey carries two clear orthologs; no pre-vertebrate lineage scores above class-A
background.

**GPR88 is unrelated to either family** and reaches back only to jawed vertebrates;
there is no *C. elegans* ortholog (two methods agree). See `output/gpr88/result.txt`.

**Still genuinely open**: whether melatonin receptors are absent from *C. elegans*
as a lineage (only one species and one method-pair has been tested properly), whether
a structural search finds what the vertebrate-anchored profiles cannot, and what
Keating 2003's supplementary table actually says.

---

## Why the existing tree does not already answer this

The b120 panel (173 taxa) in `../dop6_c24a8_6` contains `hMTNR1A`, `hMTNR1B` and
`hGPR50`, and they form a 3-tip clade whose sister is three *Mnemiopsis* sequences,
with no nematode near them. That is an artefact of panel construction, not a result.

Panel membership required monoamine HMM ≥ 120. **The melatonin receptors score
95.3–102.5** — they are in the panel as declared exceptions (see `build_strict.py`,
which removes them for exactly this reason). Any worm melatonin-receptor ortholog
would score in that same band or below, and 21 *C. elegans* class-A GPCRs do. They
could not have entered the panel. The cell is untested, not negative.

## What was run

Full log in `docs/search_registry.tsv`; coverage by lineage in
`docs/coverage_matrix.md`.

**1. Proteome completeness.** 12 UniProt reference proteomes checked against each
one's own `proteinCount`. All complete; *C. elegans* 26,629/26,629, which is what
licenses an absence call from it at all. (`data/proteome_manifest.tsv`. The human
set is 20,416/147,520 because it is the canonical one-per-gene set, not truncation —
it is used only as a reciprocal target.)

**2. BLASTP, three human queries against 11 proteomes.** Every species returns hits
at 24–30 % identity. At class-A GPCR distance that is what any 7TM protein returns
for any other, so the top hit carries no information on its own.

**3. Reciprocal best hit back against the human canonical proteome.** This separates
cleanly. 32 sequences across 9 species reciprocate to a melatonin receptor —
cnidarian 8, echinoderm 5, fish 5, tunicate 3, mollusc 3, insect 3, brachiopod 2,
priapulid 2, amphioxus 1. **Of the 8 *C. elegans* candidates, 0 do**; they reciprocate
to neuropeptide Y, neuropeptide FF, neuromedin-K and dopamine D3 receptors instead.
Both other ecdysozoans sampled (*Drosophila*, *Priapulus*) return hits, so the
negative is not a clade-wide detection failure.

**4. A melatonin profile HMM, built from those 32 plus the 3 human references**
(368 positions, `data/MTNR.hmm`) and scanned against the same proteomes. The best
*C. elegans* hit is **SER-2 at 232.1** — a tyramine receptor, i.e. the ordinary
aminergic background, with nothing standing out above it.

**5. The control that undercuts (4).** Leave-one-species-out: rebuild the profile
without a species, then ask whether it recovers that species' own ortholog.

| held out | its ortholog scores | rank | best non-ortholog in that species |
|---|---|---|---|
| *D. melanogaster* | 215 | **12** | 232 |
| *N. vectensis* | 231 | 2 | 275 |
| *P. caudatus* | 266 | 3 | 293 |
| *M. yessoensis* | 364 | 1 | 268 |
| *S. purpuratus* | 313 | 1 | 303 |
| *C. intestinalis* | 286 | 1 | 193 |
| *L. anatina* | 293 | 1 | 260 |
| *B. floridae* | 364 | 1 | 334 |

In three of eight species the profile cannot pull a known ortholog above that
species' own aminergic background, and *Drosophila*'s sits at 215 — **below the
*C. elegans* ceiling of 232**. The profile's detection limit is inside the range
where a worm ortholog would have to be found. It therefore cannot report an absence.

## PCDR-1 / F59D12.1 — the proposed candidate, and why it is not a clear ortholog

The melatonin-receptor assignment for this gene originates in **Keating et al.,
*Curr Biol* 13:1715–20 (2003), PMID 14521838** (verified via esummary), an RNAi
knockout survey of 60 GPCRs that assigned families by whole-genome homology. UniProt
O17899 now carries the name "Melatonin receptor pdcr-1"; its other name is
**Pathogen clearance-defective receptor 1**, from the orphan-GPCR screen that
rediscovered it (PMID 30692178), and the melatonin association was revived
functionally by PMID 32958651 (melatonin promotes sleep via the BK channel).

**The same gene is `cel_F59D12.1` in `../dop6_c24a8_6`**, where it was already
characterised for unrelated reasons (§6c, "Long-branch orphans") — so part of the
evidence below was in hand before this sub-project existed.

**For the assignment:**

* Its single best hit in the human canonical proteome **is MTNR1A** (71.6 bits),
  with GPR50 third and MTNR1B seventh — three of its top eight hits are the
  melatonin family.
* Direct functional evidence exists that melatonin acts through it (PMID 32958651).

**Against it, on the orthology question specifically:**

1. **It is not a reciprocal best hit.** In the forward direction it ranks **35th of
   218** for MTNR1A, 65th for MTNR1B, 43rd for GPR50. Thirty-four *C. elegans*
   proteins are better MTNR1A hits than it is.
2. **Its margin over a competing family is 3.8 bits** — 71.6 to MTNR1A against 67.8
   to neuropeptide Y receptor 2, a ratio of 1.06×. The membership rule the
   methodology recommends asks for ≥ 2×.
3. **Its own family does not point at melatonin.** PCDR-1 has clear, much closer
   relatives across ten species (*Drosophila* Q9VZD1 at 177 bits, *Priapulus*
   A0ABM1EQA9 at 185 — both far closer to it than any melatonin receptor at 71.6).
   **Not one of those ten reciprocates to a melatonin receptor**; they return 5-HT4,
   β2-adrenergic, NPY2, NPFF, GPR52, GPR21 and α1B. A melatonin top-hit that no
   member of the gene's own family shares is a property of one sequence, not of a
   lineage — which is what orthology would require.
4. **It is extremely divergent relative to an aminergic panel — though not relative
   to its own relatives.** In the dop-6 aminergic tree it carries the longest branch
   of all: 5.98 subs/site, robust z = 17.4 against a median of 0.348
   (`../dop6_c24a8_6` §6c). A sequence that divergent is attracted to whichever panel
   members are themselves most divergent, and in an aminergic panel that is the
   melatonin receptors — a sufficient explanation for a 2003 whole-genome homology
   call without any shared ancestry. **But that outlier status is a property of that
   panel, not of the sequence.** In the 178-taxon panel here, which contains PCDR-1's
   actual relatives, its branch is 2.652 subs/site, robust z = 1.5, **rank 10 of 178**
   against a median of 0.852; the longest in this panel is `PCDRfam_Cint_H2XQ43` at
   4.334. The divergence is real but unremarkable once the panel is not restricted to
   aminergics, so the long-branch explanation is not needed — and, per the controls in
   §"PCDR-1's real family", is not what is driving the placement.
5. **Position 3.32 is Q** — neither the aspartate that marks an amine receptor nor
   the M/T that the melatonin receptors carry (`../dop6_c24a8_6/output/tree/pocket.tsv`).
6. **Monoamine HMM 52.8**, the lowest of the *C. elegans* set, and only 318 aa.

### What the three papers actually claim (abstracts; full texts are paywalled)

Checked 2026-11-03 via PubMed efetch, since the melatonin assignment had to come from
somewhere.

* **Keating 2003 (PMID 14521838)** — *does not call F59D12.1 a melatonin receptor
  anywhere in the abstract*. It screened "all GPCRs predicted to bind either
  small-molecule neurotransmitters or neuropeptides" by RNAi, and reports F59D12.1
  only among genes with "a role in locomotion". Its own summary line is that
  "most phenotypes were observed for NPY-like receptors that may bind neuropeptides"
  — which is where this project's tree puts PCDR-1.
* **Anderson 2019 (PMID 30692178)** — names the gene PCDR-1 and calls it an
  **orphan** GPCR, required for pathogen clearance. No ligand assignment.
* **Niu 2020 (PMID 32958651)** — "the Mel receptor PCDR-1", on functional grounds:
  knockouts of *pcdr-1*, *slo-1* and *homt-1* are non-additive, and exogenous
  melatonin fails to act in *pcdr-1* mutants. A ligand result, not a homology one.

So the UniProt protein name traces to a 2020 functional paper, and the family
assignment it echoes is not visibly argued in any of the three abstracts.

**OPEN — needs Yale access**: Keating 2003's supplementary classification table is
where the per-gene family assignment for all 60 receptors will live, and it is the
one place a melatonin-homology claim could actually be stated. Worth fetching before
this question is called settled.

### The binding pocket says neither aminergic nor melatonergic

`pocket_compare.py`, output in `output/pocket_compare.txt`. Position 3.32 is anchored
on the alignment column that all four ligand-bound reference structures and 35/35
human aminergic receptors occupy with the canonical aspartate; the call is then
cross-checked against the structurally anchored run in `../dop6_c24a8_6`, where
**71/71 shared taxa return the same residue**. Group consensus is a majority vote,
and the melatonin group is the 8-tip clade this project recovered at 100/100 rather
than a name-based set.

| group | 3.32 |
|---|---|
| human aminergic receptors | **D** (35/35) |
| melatonin clade (verified) | **M** (7/8; the exception is GPR50 at T) |
| human neuropeptide receptors | Q (3/6) |
| PCDR-1 / F59D12.1 | **Q** |

Across the 15-column TM3 window, PCDR-1 matches the melatonin consensus at **1 of 14**
columns and the aminergic consensus at **0 of 7**. It carries neither signature.

*Caveat on how much this proves.* PCDR-1 also matches its own nearest neighbours at
0 of 1 informative columns — it is divergent enough to match nothing. So "the pocket
is not melatonergic" is sound, but "the pocket is diagnostic of something else" is
not supported; the result is an absence of the melatonin signature, not positive
evidence for another ligand class.

**The distinction that matters here**: melatonin activating PCDR-1 is a claim about
ligand, and it may well be true. Orthology to MTNR1A/1B is a claim about descent.
A receptor can acquire a ligand independently — that is ordinary, not exotic — and
the evidence above is consistent with melatonin signalling in *C. elegans* running
through a receptor that is **not** a melatonin-receptor orthologue. Nothing here
argues against the 2020 functional result; it argues against the 2003 family
assignment that the gene's name still carries.

## The trees

Two rounds. Both at MAFFT L-INS-i, LG+F+I+R7 (the model the b120 run converged on;
fixed rather than re-selected, per the two-tier workflow), trimmed at gt0.1 and
gappyout, SH-aLRT and UFBoot from 1000 replicates each.

**Round 1, `output/tree/`, 120 taxa.** Melatonin clade recovered at 100/100 and
containing only the 3 human receptors plus 5 zebrafish sequences. PCDR-1 outside it,
in a 5-tip clade at 99.7/100 that attaches to the neuropeptide receptors. Identical
at both trim levels.

**Round 1 had a gap** that round 2 exists to close: it held only 20 *C. elegans*
tips — the 19 aminergic reference receptors plus PCDR-1 — so it could place PCDR-1
but could not ask whether some *other* worm gene is closer to the melatonin
receptors. 71 worm proteins outscore PCDR-1 against at least one melatonin query.

**Round 2, `output/tree2/`, 178 taxa, 78 of them *C. elegans*** (round 1 plus all 71
competitors). Results in `output/tree2/placement.txt`:

| | gt0.1 | gappyout |
|---|---|---|
| melatonin clade | 8 tips, **100/100** | 7 tips, **43.5/30** |
| worm sequences inside it | none | none |
| PCDR-1's own 5-tip clade | 99.8/100 | 100/100 |
| nearest clade out holding worm genes | 4 nodes out, 48 worm members | 6 nodes out, 48 worm members |
| PCDR-1 among those | no | no |

Three things to read off this:

1. **No *C. elegans* sequence falls inside the melatonin clade** at either trim level,
   now that every worm gene that outscores PCDR-1 is in the panel.
2. **PCDR-1 is not even the nearest worm gene to the melatonin receptors.** 48 other
   worm sequences sit in a clade closer to them than PCDR-1's does. PCDR-1's own
   placement, by contrast, is one of the most stable things in the tree (99.8–100/100
   for its 5-tip clade at both trim levels).
3. **The melatonin clade itself is trim-sensitive in the expanded panel** — 100/100 at
   gt0.1 but 43.5/30 at gappyout, where it also loses a zebrafish member. Adding 58
   divergent worm sequences destabilised it. The conclusions above do not depend on
   that node's support (PCDR-1 is far outside it either way), but any claim that
   *rests* on melatonin-clade membership should use gt0.1 and say so.

### Formal test

Forcing PCDR-1 into exclusive monophyly with the melatonin clade, on the 120-taxon
panel. The constraint is well formed: the ML tree already holds those 8 taxa as an
exclusive clade, so nothing is ejected and the cost is attributable to PCDR-1 alone.

| trim | ΔlogL | p-AU | p-KH | p-SH | c-ELW |
|---|---|---|---|---|---|
| gt0.1 | 102.8 | 3.2e-08 | 0 | 0 | 1.6e-15 |
| gappyout | 98.7 | 2.8e-04 | 0 | 0 | 1.1e-16 |

**Rejected at both trim levels by all four tests.** Read the ΔlogL, not the p-value
exponents — those move four orders of magnitude between trim levels while ΔlogL moves
by 4, which is what quantities computed against an empty bootstrap tail look like.

A first attempt used a *partial* constraint and is recorded in the registry as
uninformative: the constrained search placed 45 tips inside the forced group, so it
tested only that the taxa share some clade (ΔlogL 3.2, not rejected), and at gappyout
it returned the ML tree unchanged. That is a constraint-construction artefact, not
evidence for the hypothesis.

### The calibration control, and what it does to the headline

~100 log units means nothing without a denominator, so the same exclusive constraint
was applied to two ordinary worm aminergic receptors that nobody has ever proposed as
melatonin receptors (gt0.1, `output/tree/calib2.iqtree`):

| taxon forced into the melatonin clade | ΔlogL | p-AU |
|---|---|---|
| **PCDR-1 / F59D12.1** | **104.8** | 2.0e-07 |
| *dop-1*, a dopamine receptor | 105.6 | 4.8e-08 |
| *ser-2*, a tyramine receptor | 111.0 | 8.6e-52 |

All three are rejected, and **PCDR-1's cost is indistinguishable from DOP-1's** — 0.8
log units apart, inside the 0–10 unit seed-to-seed spread this alignment family was
already documented to have (`../dop6_c24a8_6` §6).

**This changes how the rejection should be stated.** It is not that PCDR-1 specifically
fails a melatonin test; it is that PCDR-1 pays exactly the ordinary price any aminergic
receptor pays to be dragged into that clade. The ~100 units measures the distance from
the melatonin receptors to the rest of the panel, not anything about PCDR-1. The
defensible sentence is: *PCDR-1 is no more closely related to the melatonin receptors
than a dopamine receptor is.*

The `ser-2` row also illustrates why the p-value column must not be ranked: 8.6e-52
against *dop-1*'s 4.8e-08 for trees four log units apart. Those are extrapolations from
an empty bootstrap tail and carry no ordering.

*Not done*: a KH standard error on the ΔlogL (needs `-wsl` per-site output), and the
positive control the methodology also asks for — a group an independent line of evidence
says must be monophyletic — which would show whether this alignment can adjudicate
anything at this depth.

## PCDR-1's real family

> **Superseded, 2026-11-05.** Everything below was derived from the 178-taxon panel
> (`output/tree2`), which carried only ONE member per species for this family. With
> the family sampled properly - 258 taxa, 88 above-floor hits, the full aminergic
> backbone - the GPR52/GPR21 grouping does not survive, and an exclusive-constraint
> AU test rejects it at dLogL 856 / p-AU 1.4e-90 (gt0.1) and 709 / 5.1e-31
> (gappyout), calibrated against 0.58 for a group known to be real. The claim that
> the family runs from cnidarian to human is retracted: both cnidarian hits fall
> OUTSIDE the clade. What replaces it is in "The protostome clade" below. The
> long-branch controls remain valid as far as they go - they tested whether the
> tree2 topology was an attraction artifact, not whether tree2 was adequately
> sampled, and inadequate sampling is what was actually wrong.

### The original tree2 analysis (kept for the record)

In the 178-taxon tree PCDR-1 does not sit near the melatonin receptors. It sits inside a
small, well-supported clade of orphan receptors: **human GPR52 and GPR21** plus single
members from sea urchin (`Spur_A0A7M7GHW5`), *Drosophila* (`Dmel_Q9VZD1`), *Priapulus*
(`Pcau_A0ABM1EQA9`), oyster (`Cvir_A0A8B8D717`) and lamp shell (`Lana_A0A1S3KFU7`).
That is a positive placement, not a leftover, and it needed testing on its own terms —
a divergent sequence can be dragged into any clade by long-branch attraction.

**Membership evidence that does not come from the tree.** Running BLAST *forward* from
hGPR52 and hGPR21 ranks the other clade members near the top of their own proteomes:
sea urchin 1st against both (170 / 179 bits, so reciprocal best hits), *Drosophila* 1st
for GPR52, lamp shell 2nd, *Priapulus* 4th–5th, oyster 12th–34th. **PCDR-1 itself ranks
only 15th–16th** — it is not an RBH for this family either, consistent with a sequence
that is weakly similar to everything. So the family is corroborated independently of the
tree, but PCDR-1's membership in it rests on the tree alone.

**Two controls, both in `output/lba/`.**

| run | model | orphan-family clade | PCDR-1 inside it? | melatonin clade |
|---|---|---|---|---|
| ML, PCDR-1 present (`output/tree2/gt0.1`) | LG+F+I+R7 | 8 tips, **98.3/98** | yes | 8 tips, 100/100; PCDR-1 outside |
| PCDR-1 deleted, tree rebuilt (`no_pcdr1`) | LG+F+I+R7 | 7 tips, **97.7/93**, exactly monophyletic | n/a | 8 tips, 100/100 |
| site-heterogeneous (`pmsf`) | LG+C20+F+G, PMSF guide = `tree2/gt0.1` | 8 tips, **96/99** | yes | 8 tips, 100/100; PCDR-1 outside |

Both controls support a genuine placement:

1. **The family does not depend on PCDR-1.** Delete it and rebuild, and the remaining
   seven members come out *exactly* monophyletic at 97.7/93 — nothing else falls inside.
   PCDR-1 was therefore not holding together or destabilising a group that would
   otherwise dissolve.
2. **The placement survives the LBA-resistant model.** `LG+C20+F+G` via PMSF is the
   standard defence against long-branch attraction, and it is where an attraction-driven
   placement typically collapses. PCDR-1 stays inside the family at 96/99.

Taken with the branch-length correction in §4 above — z = 1.5, rank 10 of 178, not an
outlier in this panel — the consistent reading is that **PCDR-1 is a divergent but real
member of a GPR21/GPR52-related orphan family with invertebrate representatives**, and
that its apparent rogue-long-branch behaviour in the dop-6 tree was an artifact of a
panel that did not contain its relatives.

*Caveats.* The two controls are single runs, not replicates over seeds. The sequence-level
evidence for PCDR-1's own membership is weak (rank 15–16), so the tree is carrying that
specific claim; the family's existence is better supported than PCDR-1's place in it.
Neither point changes the melatonin conclusion in either direction — the melatonin clade
is 8 tips at 100/100 with PCDR-1 outside it in all three runs.

## The protostome clade — what the 258-taxon panel actually supports

Panel built by `build_orphan_tree.py`: the b120 aminergic backbone minus its nematode
neuropeptide expansion, all 88 above-floor hits from `orphan_family.py` unfiltered
(membership decided by the tree, not by what a record is named), PCDR-1, the melatonin
family including the two lamprey orthologs, and 3 human neuropeptide receptors to root
on. Read with `orphan_clade.py`.

**Positive control first** (methodology 7.2). The melatonin family is monophyletic on
evidence independent of this panel, so if it fragments here nothing else is
interpretable. It comes out **exactly monophyletic, 10 tips, 100/100, at both trim
levels**. The alignment resolves this depth.

**The supported group.** PCDR-1 sits in an 8-member clade at **99.7/100** (gt0.1) and
**99.8/100** (gappyout), identical membership at both:

| | |
|---|---|
| *C. elegans* | PCDR-1 / F59D12.1 |
| *Drosophila* | **DopEcR** (Q9VZD1) |
| *Priapulus* | A0ABM1EQA9 — PCDR-1's single closest relative, 82/100 |
| *Crassostrea* | A0A8B8D717 |
| *Lingula* | A0A1S3KFU7 |
| *Limulus* | A0ABM1BHV8, A0ABM1BZC6, A0ABM1S0K4 |

All protostomes. Adding the two *Saccoglossus* (hemichordate) sequences gives 10 tips
at 94.7-95.9 / 98-100.

**What is NOT supported.** The node that adds hGPR52 and hGPR21 is **36.2/91** at
gt0.1 and **3.2/52** at gappyout - it collapses. The exclusive-constraint AU test
rejects that grouping outright (table above). So PCDR-1 belongs to a deeply conserved
protostome orthology group; whether that group is specifically sister to the human
orphans is unresolved, and the version of the claim that named them is rejected.

**The one member with measured ligand affinity is DopEcR**, which binds
20-hydroxyecdysone and dopamine, and it is inside PCDR-1's 6-tip subclade at 97.3/100.
Every invertebrate member named "Probable G-protein coupled receptor 52" in UniProt is
`proteinExistence=3, Inferred from homology` with no FUNCTION annotated, so those names
carry no experimental weight in either direction. GPR52 itself has synthetic-agonist
pharmacology but no endogenous ligand; GPR21 is an orphan outright. The family is
conserved in sequence and NOT in ligand - which is the same pattern as PCDR-1 acquiring
melatonin without inheriting it.

**Bleed-through, quantified.** 63 of the 88 above-floor hits fall outside the family
clade. The ADRB2-based floor in `orphan_family.py` was porous, so the "93 members"
figure was mostly aminergic receptors from other species clearing a human threshold.
The tree was the right arbiter and this is what it was for.

*Caveats.* Single runs, not replicated over seeds. PCDR-1's own membership ranks only
15th-16th by forward BLAST from hGPR52/hGPR21, so the tree carries that specific claim
while the family's existence is better supported. The `calib_proto` constrained search
returned a HIGHER likelihood than the unconstrained ML search, which means the ML
search settled in a local optimum - a reason to treat the backbone outside these
clades as unresolved rather than merely weakly supported.

## How far back does the melatonin clade go?

A profile built from **only the 8 verified melatonin-clade members** (`data/MEL_verified.hmm`,
354 positions) was scanned across the pre-vertebrate proteome set.

* **Lamprey is a clean positive**: 544 and 451 bits, where every other lineage tops out
  at 134–189. The two hits are `A0AAJ7T8G6` ("Melatonin receptor type 1A-A-like") and
  `A0ACM8ASR7` ("Melatonin receptor type 1B-A-like"), both reciprocal best hits to human
  MTNR1A at 429 and 377 bits. The family therefore predates the jawed/jawless split.
* **Nothing below that.** Tunicate 169, amphioxus 179, acorn worm 183, echinoderms
  171–185, cnidarian 168, molluscs 148–181 — indistinguishable from the general class-A
  background.
* **Profile-walking through lamprey did not change this.** Adding the two lamprey
  sequences and rebuilding moved invertebrate scores by a few bits only (169 → 172,
  179 → 185), so the invertebrate negative is not simply a reach problem.

*Caveat, and it is the limiting one*: the profile is anchored on ten vertebrate sequences,
and lamprey is also a vertebrate, so the walk never gained a foothold outside the subphylum.
A genuinely divergent invertebrate member could still be invisible to it. Discriminating
that needs the structural axis (foldseek against the pre-vertebrate chordates), the one
method here that does not inherit sequence-level bias. Not run.

## The circularity to resolve before anything is concluded

The leave-one-out table has two readings and they lead opposite ways:

* **The profile is insensitive at ecdysozoan distance**, in which case the *C.
  elegans* negative in (3) is not trustworthy either.
* **Or the weak RBH calls are not real orthologs.** *Drosophila*'s three reciprocate
  at only 108–113 bits, and *Drosophila* is not generally held to have melatonin
  receptors. If those are false positives, the "truth set" in (5) is wrong, and the
  low ranks are the profile behaving correctly.

Nothing run so far distinguishes these. Both are resolved the same way: by tree
placement, which is the arbiter the methodology calls for and the one axis not yet
run.

## Next

*Done since this list was first written*: the tree round (178 taxa, §"The trees" — and
it settled the circularity: only 5 of 32 invertebrate RBH calls joined the melatonin
clade, so the RBH screen was wrong wholesale and the profile was behaving correctly);
the literature check (§"What the three papers actually claim").

1. **Structural axis.** `foldseek` is installed and the AlphaFold models are not yet
   used. This is now the highest-value remaining item, and for two separate questions:
   it is the genuinely orthogonal second method the *C. elegans* negative needs
   (sequence search and profile search are not independent of each other), and it is
   the only axis that does not inherit the vertebrate anchoring that limits the
   pre-vertebrate scan.
2. **Other nematodes.** *C. briggsae*, *P. pacificus* and an outgroup nematode
   (*Trichinella*) — the dop-6 family-history panel (`../dop6_c24a8_6` §7) already
   has these assembled. A gene lost in *C. elegans* may survive elsewhere in the
   phylum, and a phylum-wide absence is a much stronger claim than a species one.
3. **Keating 2003's supplementary classification table** (needs Yale library access).
   The one place the original homology argument could be written down.
4. **Replicate the two LBA controls over seeds**, and add the §7.2 positive control
   the methodology asks for — a group an independent line of evidence says must be
   monophyletic — to show whether this alignment can adjudicate anything at this depth.
   Also the KH standard error on the ΔlogL (needs `-wsl`).

## Layout

```
data/     queries, profile seed set, MTNR.hmm, blast DBs, proteome manifest
output/   blast/ — raw search output, kept (not /tmp)
docs/     search_registry.tsv, coverage_matrix.md
```
