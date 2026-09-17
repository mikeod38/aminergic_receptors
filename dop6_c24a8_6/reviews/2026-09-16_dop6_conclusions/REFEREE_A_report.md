# Referee A — inferential scrutiny, dop-6 / C24A8.6 sub-project

Report date 2026-09-17. Charge items A1–A6 in `CLAIMS.md` (including the ADDENDUM).
All tree work below was re-derived from the treefiles and alignments, not read off the
prose. Scratch work in `/private/tmp/.../scratchpad/au/`; nothing in `README.md` or
`output/` was modified.

## Nomenclature

| term | meaning |
|---|---|
| OA / TA / DA / 5-HT | octopamine / tyramine / dopamine / serotonin |
| p-AU | p-value of the approximately-unbiased topology test |
| c-ELW | expected-likelihood weight; the share of total likelihood weight a topology holds |
| ΔlogL | log-likelihood cost of a constrained topology relative to the best tree in the set |
| PI | parsimony-informative (alignment columns) |
| Rhabditina | the nematode suborder containing *Caenorhabditis* and *Pristionchus* |
| v1 / b120 | the 163-taxon and 173-taxon panels |
| h1–h4 | my four constrained OA topologies (see MAJOR-1) |

---

## Headline finding

**C31 reaches a defensible conclusion by an indefensible route, and the route can be
replaced with a real test that the project already owns.**

C31 argues the tree is "close to uninformative *against*" an OA assignment because "a
nematode-specific lineage is exactly where a phylum-specific ligand could sit, and such a
lineage would not be expected to group with insect OA receptors." Two things are wrong
with that as an argument, and one thing is right about the conclusion.

**Wrong 1 — the premise is refuted by the panel's own data.** Octopamine is not a
phylum-specific ligand, and in this very tree every characterised OA receptor sits in a
clade that crosses phyla at high support. Re-derived from `b120/trees/gt0.1.treefile`:

| OA receptor | smallest clade containing it and a non-nematode | SH-aLRT/UFboot |
|---|---|---|
| CEL SER-3 / SER-6 | 11 tips: + dmOAMB, *Limulus*, *Mizuhopecten*, *Lingula*, *Saccoglossus* | **96.4/98** |
| CEL OCTR-1 | 14 tips: + hADRA2A/2B, *Danio* ADA2A, dmOctα2R, *Branchiostoma*, *Crassostrea*, *Apostichopus* | **98.1/100** |
| dmOctβ1/2/3R | 6 tips: + *Mizuhopecten*, hHTR4, hTAAR1 | 80.1/97 |
| bmOAR1 / dmOctTyrR | 9 tips: + CEL SER-2, *Sra* SER-2, *Crassostrea*, *Priapulus* | 92.9/99 |

So the panel contains four internal positive controls demonstrating that OA-receptor
identity *is* recoverable from phylogenetic position across Bilateria. C31 asserts the
test is uninformative while the tree contains four demonstrations that the test works.
The premise would be legitimate for a genuinely nematode-restricted ligand (an
ascaroside, betaine); for octopamine it is not.

**Wrong 2 — it is deployed asymmetrically.** "Nematode-specific lineage ⇒ tree
uninformative about ligand" is ligand-agnostic. If accepted, it equally voids the project's
own "D2-like by descent is defensible" and the local Gi/o inference, which are ligand-and-
coupling inferences drawn from the same clade by the same logic. As written, C31 is an
escape hatch applied only where the answer is unwelcome.

**Right — but for a reason never established.** The conclusion is close to correct, and it
is testable. `topology_test.py` already implements exactly the machinery (constrained
search + AU/SH), including a docstring documenting a previous failure mode; it was never
pointed at the OA question and has produced no output file. I ran it. Four hard constraints
forcing the five dop-6 orthologues to be monophyletic with each OA group **to the exclusion
of** dop-5, dop-3, dop-2, dmDop2R and hDRD2/3/4, under the published model:

AU test, b120 gt0.1 (LG+F+I+R7, 10,000 RELL replicates); gappyout in brackets (LG+F+R8):

| topology | ΔlogL | c-ELW | **p-AU** | rejected at 0.05? |
|---|---|---|---|---|
| published ML tree | 0 [0] | 0.876 [0.847] | 1.0 [1.0] | — |
| h1  dop-6 + SER-3/SER-6/dmOAMB | 26.9 [25.5] | 0.069 [0.071] | **0.132** [0.157] | no |
| h2  dop-6 + OCTR-1/α2-adrenergic | 25.7 [23.3] | 0.021 [0.023] | **0.063** [0.103] | no |
| h3  dop-6 + dmOctβ1/2/3R | 27.9 [24.5] | 0.028 [0.039] | **0.096** [0.133] | no |
| h4  dop-6 + SER-2/bmOAR1/dmOctTyrR | 26.0 [24.3] | 0.006 [0.020] | **0.060** [0.091] | no |

**No OA placement is rejected at α = 0.05, at either trim level.** But none is close to
neutral either: each costs 23–28 log-likelihood units, and the four together hold ~12 %
(gt0.1) to ~15 % (gappyout) of the expected-likelihood weight against ~85–88 % for the
published topology. So the correct statement is neither C30's bare negative nor C31's
"close to uninformative": **the data disfavour every OA placement by a consistent and
substantial margin, but do not exclude any of them.** That is a quantitative claim a
reviewer can check, and it replaces a rhetorical move with a number.

There is also a *legitimate* version of C31's argument that the project does not make, and
should: **ligand is not a topological property.** The project's own Fitch reconstruction
(C33) counts 15 ligand changes across 56 labelled tips — roughly one transmitter switch per
3.7 labelled tips. A ligand switch on the dop-6 stem branch is therefore common-rate for
this family and is **invisible to any topology test**, including mine. That, not
"phylum-specific ligand in a phylum-specific lineage", is why the tree cannot exclude OA.

---

## MAJOR findings

### MAJOR-1 — C31 (and C30): the stated reason for the negative result is refuted by the panel; the quantitative test was available and unrun

Detail as in the Headline. Supplementary point against C31: *C. elegans* already has three
characterised OA receptors (SER-3, SER-6, OCTR-1) and C30 places all three elsewhere. An OA
assignment for DOP-6 therefore requires a *fourth*, lineage-restricted OA receptor — a
possible but unparsimonious auxiliary hypothesis that C31 never prices in.

C30 itself I **confirm as a description of the tree**: DOP-6's orthologue clade shares no
node with any OA receptor short of the deep backbone, on both panels and all trim levels.

**Required action (compute-executable now).** Run `topology_test.py` retargeted at the four
OA constraints (constraint files and results are in the scratch directory and can be copied
in). Report p-AU and c-ELW, not a verbal characterisation. Then reword C31 as given in the
disposition list. Owner: whoever next touches the tree section; no external data needed.

### MAJOR-2 — C16: "the lowest of the aminergic receptors tested" is contradicted by the file it cites

`output/dnds_vs_briggsae.tsv` contains **15 genes with a whole-gene dN**, not the 9 that
C16 and README §4 list. Recomputed from the table's own `Nd`/`N` columns:

| rank | gene | dN | 95 % CI (Poisson on Nd) |
|---|---|---|---|
| 1 | **tyra-2** | **0.0116** | 0.0061–0.0201 |
| 2 | **octr-1** | **0.0236** | 0.0147–0.0360 |
| 3 | **ser-5** | **0.0246** | 0.0158–0.0365 |
| 4 | dop-6 | 0.0253 | 0.0179–0.0347 |
| 5 | dop-5 | 0.0261 | 0.0190–0.0351 |

dop-6 is **4th of 15**, not 1st. Six comparators are omitted from the claim (ser-5, ser-6,
octr-1, tyra-2, tyra-3, gar-3) and **the three that beat dop-6 are all in the omitted set**.
I make no inference about intent, but the omission is not random with respect to the claim.
`README.md` line 266 states the script covers 17 receptors, so the fuller set was computed.

Two further problems. (i) No uncertainty is reported anywhere; dop-6's CI overlaps ranks
1–8, so the ordering among the top half is not resolved by these data. (ii) The claim is
load-bearing: §4's heading is "an orphan, and a strongly constrained one", and the memory
index records "most constrained aminergic receptor by dN". With octr-1 and tyra-2 below it,
dop-6 is unremarkable within the family.

**Required action (compute-executable now).** Reprint the full 15-gene ranking with CIs and
restate the claim (wording in the disposition list). Also correct the memory file
`dop6_is_an_orphan_receptor.md`, which propagates the superlative.

### MAJOR-3 — C17 / charge A6: the Fisher test has 20 % power at the observed effect, and the same difference it calls non-significant is what produces C16's ranking

Recomputed from `output/mirror_bias_check.txt`: N-half 27/929 vs C-half 13/668 nonsynonymous
differences; rates 0.0291 vs 0.0195, **ratio 1.49 with the C-half lower** — i.e. the point
estimate lies in exactly the direction the bias hypothesis predicts. Fisher p = 0.258.

Power, by simulation (8,000 replicates per point, α = 0.05 two-sided, N = 929 vs 668 sites,
baseline 0.0291):

| true C-half deflation | 1.25× | **1.49× (observed)** | 2.0× | 2.5× | 3.0× |
|---|---|---|---|---|---|
| power | 0.09 | **0.20** | 0.47 | 0.65 | 0.78 |

The test attains 80 % power only at roughly a **three-fold** suppression. Woolf 95 % CI on
the rate ratio is **0.77–2.95**: the data are compatible with the C-half being 1.3× *faster*
or 2.9× *slower*. "Not significant" here carries almost no information, and C17's wording
("the difference is not significant **anyway**", README line 172) presents it as a third
supporting check when it is an uninformative one whose point estimate runs the other way.

**The entanglement is the real problem.** If the C-half rate is deflated relative to the
N-half by a factor the test cannot detect, dop-6's whole-gene dN rises and its rank changes:

| assumed deflation | 1.00× | 1.25× | 1.49× | 2.00× |
|---|---|---|---|---|
| dop-6 whole-gene dN | 0.0253 | 0.0274 | 0.0294 | 0.0337 |
| rank among the 9 genes C16 lists | 1st | 2nd | 3rd | 4th |

A deflation of 1.25× — inside the CI, and one the test detects 9 % of the time — is enough to
lose the headline. The project cannot simultaneously treat the N/C difference as real enough
to deliver the lowest dN in the family and as too small to be worth worrying about.

Second problem with C17: the check "recomputing dN from four palindrome-free isolates' own
CDSs gives identical values" is near-tautological. Test 1 of the same file shows those
isolates differ from N2 by **0–2 nt in 2,079**, so the NG86 counts are arithmetically
forced to be identical (they are, to the decimal, for all four). It bounds N2-specific gene
conversion at ≈2/864 sites — which is worth saying — but it is not an independent
confirmation and should not be counted as one of "three ways".

**What actually carries C17 is the timing argument**, and it is sound: the duplication is
absent from CB4856, XZ1516, DL238 and MY2147, so it postdates the common ancestor of those
*C. elegans* isolates and occupies a negligible fraction of the *Cel*/*Cbr* branch on which
dN was measured. That argument alone is close to decisive and does not need the other two.

**Required action (compute-executable now).** Delete the Fisher test from the supporting
list or report it with its power and CI; demote the four-isolate recomputation to a bound on
gene conversion; lead C17 with the timing argument. Add CIs to C16 (MAJOR-2).

### MAJOR-4 — C29 / charge A2: "exclusively nematode" is exclusively *Rhabditina*, and the project's own output prints the evidence

The 13-tip clade, re-derived from both b120 trim levels, contains five species: *C. elegans*,
*C. brenneri*, *C. remanei*, *C. briggsae*, *P. pacificus*. **All five are Rhabditina** —
`output/tree/b120/nematode_expansion.txt` prints exactly that in its own `clade` column
("Rhabditina (Caenorhabditis)", "Rhabditina (Diplogastero.)") and then concludes "a radiation
inside the nematodes ... no vertebrate or insect orthologue anywhere".

The panel is not silent on the rest of the phylum. It contains four non-Rhabditina or
non-*Caenorhabditis* nematodes, and **every one of them falls outside the 13-tip clade**:

| species | nematode clade | genes in panel | where they fall |
|---|---|---|---|
| *Trichinella spiralis* | I | dop-2 | enters only at 16 tips |
| *Brugia malayi* | III | tyra-2, tyra-3 | tyra-3 clade; 57-tip clade |
| *Strongyloides ratti* | IV | ser-2, ser-6 | ser-2 clade; ser-7 clade |
| *Necator americanus* | V | ser-3, ser-7 | ser-3 clade; ser-7 clade |

So the supported statement is that the clade is restricted to Rhabditina — one suborder —
not to Nematoda. The difference matters for the inference drawn from it. "Nothing outside
the phylum is orthologous, so nothing outside can lend it a ligand or a coupling" becomes
"nothing we sampled outside Rhabditina is orthologous", and Compara returned no
dop-3/dop-5/dop-6 orthologue for *Trichinella*, *Brugia*, *Strongyloides* or *Necator* —
which is precisely the RBH-into-a-recent-radiation artefact the project itself withdraws as
W5. The absence is unsearched, not established.

**Robustness to taxon sampling (A2, direct answer): it is not robust, and the two known
gaps are already documented in the project's own files.**

1. `Ppa_dop-5` is **absent from `panel_b120.faa`**, although `output/ppa_dop5_search.txt`
   (written 21:39) establishes `A0A2A6CA79` as a clean reciprocal best hit at 764 bits.
   The b120 panel was assembled at 22:30 — 51 minutes later — and does not contain it.
   **W4's correction did not propagate into the panel it corrects.**
2. `Cbr_dop-3` is absent likewise.

Both gaps sit inside the very clade whose composition is the claim. A 13-tip clade missing
two of its expected members is not a safe basis for a phylum-level exclusivity statement.

**Support spread (per the coordinator's point 1).** The 13-tip clade takes UFboot 90
(gt0.1), 79 (gt0.2), 82 (gt0.3) and **60 (gappyout)** on b120, and 69 on v1 gt0.1. Any
statement about this clade must carry 60, not 90. At gappyout it is below the conventional
70 threshold — i.e. at the trim level with the least gap-driven signal, the nematode-
exclusive clade is not supported at all.

**Required action.** (a) *compute-executable now*: add `Ppa_dop-5` (A0A2A6CA79) and search
for `Cbr_dop-3` directly in the *C. briggsae* proteome; rebuild b120; re-report the clade.
(b) *PI-only*: decide whether to extend the panel to Clade I/III/IV nematode proteomes
(*Trichinella*, *Brugia*, *Strongyloides*) by direct HMM+phylogenetic placement rather than
Compara, which is the only way to convert C29 from a panel statement into a searched
absence. (c) *compute-executable now*: reword C29 to "Rhabditina" and carry the gappyout
value. Wording in the disposition list.

### MAJOR-5 — charge A5: the Gi/o inference does depend on a blocking claim, and the leakage is at the choice of clade level

Charge A5 asks whether any inference depends on the monophyly C33 denies. **Yes, in two
places.**

**(i) README line 419–420** states the Gi/o prediction "really says 'DOP-6 falls inside the
D2-like radiation, **which is Gi/o throughout**'". That is a blocking claim at clade scale.
C33 licenses *local* claims, so it is not forbidden — but it is asserted, never tested, and
the level at which it holds is chosen because it holds there. Re-derived from b120 gt0.1:

| clade | UFboot | characterised members | coupling |
|---|---|---|---|
| 13 tips | 90 (**60** gappyout) | 1 (dop-3) | Go — and C36 weakens it to IMP-only |
| 16 tips | 90 | 3 | 3/3 Gi/o-family |
| **19 tips** | **95** | 6 | **6/6 Gi/o-family** ← the level the claim uses |
| **24 tips** | **100** | 10 | **6 Gi/o vs 4 Gq** (ser-1, hHTR2A/B/C) |

The clade with the *highest* support containing DOP-6 is coupling-**mixed**, 6:4. The Gi/o
inference requires stopping at 19 tips, and the reason to stop there is that 24 tips gives a
mixed answer. Additionally, 13 of the 19 tips are uncharacterised, and C33's Fitch count (8
coupling changes over 52 labelled tips) gives no basis for assuming zero changes on those 13
branches. Nothing in the permutation test provides a per-clade probability.

**(ii) C30's own wording carries the withdrawn W3 framing.** README line 247 describes the
three OA positions as "the α2-adrenergic/Gi-o **block**", "the Gq/Gs **block**", "the
β-adrenergic/Gs **block**". W3 — "the canonical G-protein blocks hold" — is listed as
withdrawn 214 lines later. The sentence stating C30 is written in the vocabulary of the
claim the document withdraws.

**Required action.** (a) *compute-executable now*: replace the assertion with a marginal
ancestral-state reconstruction of coupling at the dop-6 stem node (Mk or parsimony with
bootstrap), reporting the posterior/support for Gi/o rather than asserting clade uniformity;
report it for the 19- **and** 24-tip levels so the level-dependence is visible.
(b) *compute-executable now*: rewrite README line 247 without "block".

### MAJOR-6 — README §6 is a superseded analysis still asserted in the present tense, including the withdrawn W2

`README.md` §6 (lines 188–256) describes a **50-sequence panel**, "Best-fit model
Q.PFAM+F+R6 (gt0.1, 1,034 columns, 50 taxa)", a trait table with 36 labelled tips and 12–13
ligand changes, and support values "SH-aLRT 73–79, UFboot 37–62". Every one of those is
superseded by §6b/§6c and by b120 (173 taxa, LG+F+I+R7, 367 columns, 56/52 labelled tips,
15/8 changes, UFboot 90/95). The two sections are not marked as sequential; a reader
encountering §6 first will take its numbers as current.

Most seriously, **line 249 states W2 verbatim**: "Its nearest characterised relatives are two
**5-HT receptors, both Gq**." Line 459 withdraws that exact sentence. The document asserts
and withdraws the same claim 210 lines apart.

**Required action (compute-executable now).** Mark §6 as superseded at its head, or fold its
surviving content into §6c and delete the stale numbers; delete line 249. This is
Referee B's charge B5 as well, but it bears on A1 because C30's wording is taken from §6.

---

## MINOR findings

**MINOR-1 — C3: D3.32 and the DRY aspartate are conflated.** C3 and README line 27 state
"D3.32 and the `DRL` motif (D120-R121-L122) sit at ~residue 120". They are different
residues. In DOP-6 (`output/dop6_c24a8_6_proteins.faa`) the only aspartates in 95–135 are
**D103** and **D120**; the pocket window in `output/tree/pocket.tsv`
(`YCDVST.YYMSSAWFFNFTY`) resolves 3.32 to **D103**, and D120-R121-L122 is the D3.49-R3.50
motif at the *cytoplasmic* end of TM3. The 103→120 spacing of 17 is canonical. C19 correctly
treats R121 as the DRY arginine, so the two claims are mutually inconsistent. **C3's
conclusion is unaffected** — both residues are below 406 and so in the single-copy half.
*Required action (compute-executable now):* change C3/README to "D3.32 (D103) and the DRL
motif (D120-R121-L122) both lie in the single-copy N-half".

**MINOR-2 — C18: the class breadth is thinner than "spanning five classes" implies, and one
structure is misclassified.** `output/gprotein_interface_multi.tsv` is n = 4 Gi, n = 4 Gs,
and **n = 1 each for Go, Gq and G11**. "Per class 29–41 %, no class rescues it" rests on a
single structure for three of the five classes. Separately, **6WHA is assigned `Gs` at 88.6 %
identity while its own title in the same file reads "IN COMPLEX WITH A MINI-GALPHA-Q
PROTEIN"** — a mini-Gs/q chimera, Gs scaffold with the Gq α5 helix, and α5 is the element
that makes the receptor contacts. Reclassifying it moves Gs to 25 % (n=3) and Gq to 40 %
(n=2). Also, two of the four "Gs" structures (6GDG, 6WHA) use engineered mini-G constructs.
Finally the per-**structure** range is **15–46 %**, not 29–41 %; 3SN6 — the β2AR–Gs complex,
the best-characterised aminergic–G protein structure in the set — gives the lowest value.
*Required action (compute-executable now):* state n per class, reclassify 6WHA as Gq or drop
it, and report the per-structure range alongside the class means.

**MINOR-3 — C19: "~35 %" is a position count and should not be read as a functional
fraction.** Charge A4 asks whether the "no TM3, no cavity" argument is real. **It is real**:
receptor activation opens the intracellular cavity by moving the TM6 cytoplasmic end away
from TM3, and a fragment beginning at Met406 (mid-ICL3) has neither TM3 nor the TM1–TM5
bundle that positions TM6 and TM7 relative to each other. The seven retained positions are
not held in any defined geometry. So 35 % is a **generous upper bound on an upper bound**,
and C19's conclusion is safer than its number. The wording should not invite a reader to
infer 35 % residual coupling. What would settle it structurally is an AlphaFold model of the
287-aa C24A8.6a product, which has not been done.
*Required action:* (a) *compute-executable now*: reword to "retains 7 of 20 contact positions
and none of the geometry that presents them". (b) *deferred-with-owner*: AlphaFold2/3 model
of C24A8.6a — owner MOD, low priority, since C5 already makes the gene experimentally
inseparable.

**MINOR-4 — C23: the README's headline comparison is not controlled.** The table at README
line 338 compares full-length **49 × 1034** against excerpts **163 × 362**. Gap fraction is
strongly taxon-number dependent, so 48.5 % → 11.4 % conflates the excerpting with a 3.3×
change in panel size. C23 *as worded in CLAIMS.md* is fine — it cites the same-panel figure
(84.4 % gaps, 164 sequences). It is the README table that is uncontrolled.
*Required action (compute-executable now):* align the same 163/173 full-length sequences,
trim at gt0.1, and report the same-panel comparison; or label the table's columns with their
differing taxon counts.

**MINOR-5 — C25 should be restated for b120.** `b120/trees/gt0.1.iqtree` reports 361 PI sites
and 375 free parameters: **0.96 sites per free parameter**, slightly *worse* than the 1.1
quoted for v1. Since the ADDENDUM makes b120 the superseding panel, C25's number is stale in
the direction that understates the limit.
*Required action (compute-executable now):* restate C25 as 361/375 = 0.96 for b120, retaining
1.1 for v1.

**MINOR-6 — C28: "identical tip counts at every level (5→9→13→16)" is true but not
reassuring.** I confirm the counts and the topology on both panels. But the counts embed the
two missing orthologues of MAJOR-4; they are identical across trim levels because the same
incomplete orthologue set was used at every level, which is stability of input, not of
inference. *Required action (compute-executable now):* qualify with "at the current
orthologue sampling", and re-check after the W4 fix.

**MINOR-7 — C34: "roughly 20 points" is the largest of three values.** `rogue_test.txt` gives
the 9-tip clade at 71 → **91** (drop_F59), **80** (drop_T21), **81** (drop_both). Dropping
both rogues gives *less* lift than dropping F59D12.1 alone, which is not the behaviour of a
clean rogue effect and suggests a component of search/support noise. C34's "roughly 20
points" quotes only drop_F59. *Required action (compute-executable now):* quote the range
+9 to +20, and note the non-monotonicity.

**MINOR-8 — C7: the cross-species survey stops where a positive result begins.**
`output/crossspecies.txt` prints "NO internal repeat ... single-copy in this species" for
*C. briggsae*, *C. brenneri*, *C. remanei* and *C. japonica*, and then for *P. pacificus*
reports **two inverted repeats > 200 bp** (611 bp at 89.4 %, 524 bp at 86.6 %) at the
orthologous locus — with no verdict line. These are nothing like the 3,929 bp / 100 %
structure and do not threaten C7, but an unremarked inverted-repeat signal at the
orthologous locus in a distant nematode is relevant to whether this locus is repeat-prone,
and it is the kind of silence a reviewer will notice. *Required action (compute-executable
now):* print a verdict line for *P. pacificus* and state explicitly that its repeats are not
homologous to the *C. elegans* palindrome.

---

## Confirmations

Stated explicitly, per the audit requirement.

- **C26 — confirmed, and stronger than claimed.** The 12 `TEST_*` tips form one exclusive
  12-tip clade on b120 gt0.1 at **SH-aLRT 98.3 / UFboot 100**, with zero aminergic receptors
  inside. I additionally recomputed the full monoamine-HMM ranking (142 genes,
  `cel_gpcr_monoamine_scores.tsv`) and confirm that the **largest single gap in the entire
  ranking is 27.0 points, between gar-1 (124.4, included) and npr-35 (97.4, excluded)** —
  exactly where the 120 cut falls. Stronger still: **the b120 panel is invariant to any cut
  placed anywhere in (97.4, 124.4]**, a 27-point interval, so the choice of "120" carries no
  information. That invariance is a better defence of the boundary than anything currently
  written down and should be stated in C26.
- **Boundary-decision ordering (coordinator's item 3) — confirmed mechanically, with one
  correction to the wording.** `rebuild_boundary120.py` reads only the HMM score file and
  applies `CUT=120.0`; it never opens a `.treefile`, `.contree`, `.iqtree`, `.splits.nex` or
  any support value, and neither does `add_boundary_testset.py`. The only tree-touching
  script, `build_b120.py`, runs afterwards. Timestamps: panel written **22:30:12**, first
  b120 tree **22:49:31**. So no b120 support value could have influenced the drop.
  **However**, `output/tree/rogue_test.txt` — which reports that dropping F59D12.1 and
  T21B4.4 lifts the 9-tip clade from 71 to 91 — was written at **22:21:34**, *8 min 38 s
  before* the panel rebuild, and two of the three dropped genes are exactly those two. So
  "the drop was decided before any tree was built" is **false as written** and true only as
  "before any *b120* tree was built"; the support consequence was demonstrably known. This
  does not make the procedure illegitimate — the criterion is a single support-blind
  threshold, it also drops F35H10.10 which was never a rogue-test subject, and it adds 12
  test tips that could have falsified the boundary. But the defence should rest on the
  27-point gap and the (97.4, 124.4] invariance, which are checkable facts, not on a timeline
  that does not quite say what it claims. *Required action (compute-executable now):* reword
  the ADDENDUM's point 1 accordingly.
- **C32 — confirmed, and MOD has not blurred it.** The distinction MOD asked me to check
  holds. v1 and b120 both give the 13-tip clade **4 independent genes, 1 characterised
  (dop-3)**, while UFboot rose 69 → 90. Support rose; independent evidence did not. One
  addition: at 24 tips, where b120 support reaches 100, the characterised set becomes 10 but
  splits **6 Gi/o vs 4 Gq**, so the extra evidence bought by higher support is internally
  contradictory about coupling (see MAJOR-5).
- **C28 — confirmed exactly.** Re-derived from `b120/trees/gt0.1.treefile`: 5 tips 100/100,
  +dop-5 **95/95**, +dop-3 **85.7/90**, +dop-2 **88/90**. Matches C28 to the decimal. The
  5→9→13→16 nesting and composition are identical on gt0.1, gt0.2, gt0.3 and gappyout.
- **Topological robustness — confirmed beyond what the project did.** I ran six independent
  unconstrained ML searches on `b120/panel.gt0.1.faa` (different seeds, one with
  `--ninit 200 --nstop 500 -pers 0.2`). **All six recover the identical dop-6 nesting**
  5→9→13→16→19→24 with identical membership. Caveat worth recording: per-run log-likelihoods
  spanned −77505.7 to −77519.0 because each run independently optimises the LG+F+I+R7
  free-rate parameters, whose surface is multimodal at this taxon:site ratio; when all trees
  are scored under one shared parameterisation the published topology is best. **Raw logL
  values from separate IQ-TREE runs in this project are therefore not comparable**, and a
  reader could easily mistake the spread for topological disagreement. It is not.
- **C33 — confirmed and slightly understated on b120.** Recomputed the smallest clade
  containing all tips of each coupling class on `b120/trees/gt0.1.treefile`: Gs **136/173**
  (79 %), Gi/o **170/173** (98 %), Gq **173/173** (100 %). The v1 values in C33 (125/163,
  160, 163) are correct and the b120 values are, if anything, worse for monophyly. One
  bookkeeping note: the annotation file uses four coupling states including both `Gi/o` and
  `Go` as separate labels, so `confirm.txt` item 3's "3/3 are Gi/o" silently merges them.
- **C20 — confirmed as reported.** `pocket_332_v2.txt` records all three ligand-bound
  structures mapping to one alignment column and the 35/35 positive control passing.
- **C24(a) — confirmed.** All 20 consensus Gα-contact positions fall inside DOP-6's retained
  blocks (12–239, 451–561). The structural reason given (ICL3 contacts Gα at its ends, which
  are the conserved helix boundaries the HMM recognises) is correct.
- **C12 — the genome-wide basis is real.** `design_specific_tools.py` BLASTs against a
  WBcel235 database counting 100 %-identity full-length hits on both strands, so "exactly two
  perfect full-length matches" is a searched result, not an assumption.
- **C15, C5, C8 — I have no criticism.** These are the strongest claims in the document. C15
  in particular follows from C12–C14 by a clean argument and is correctly stated as an
  inability to attribute, not as a negative result.

---

## Disposition

### Let stand as worded
C1, C2, C4, C5, C6, C8, C10, C11, C12, C13, C14, **C15**, C18 *(subject to MINOR-2 on the
class-n footnote)*, C20, C21, C22, C24, C27, C30, C32, C33, C35, C36, C37, C38, C39.
W1–W5 are correctly withdrawn; W2 and W4 have **not** propagated (MAJOR-6, MAJOR-4).

### Need rewording

**C3** → "The duplicated segment is end-of-ICL3 + TM6 + ECL3 + TM7 + C-tail. The ligand
pocket is not in it: D3.32 (D103) and the DRL motif (D120-R121-L122, the D3.49-R3.50 of the
TM3 cytoplasmic end) both lie in the single-copy half."

**C7** → add: "*P. pacificus* carries two short inverted repeats (611 bp at 89 %, 524 bp at
87 %) at the orthologous locus; these are unrelated in length and identity to the
*C. elegans* palindrome and are not a second instance of it."

**C16** → "Of the 15 aminergic and cholinergic receptors with a whole-gene dN vs
*C. briggsae*, dop-6 (0.0253, 95 % CI 0.018–0.035) is **fourth lowest**, behind tyra-2
(0.0116), octr-1 (0.0236) and ser-5 (0.0246) and marginally below dop-5 (0.0261). The
confidence intervals of the top eight overlap, so the ordering within that group is not
resolved by these data. dop-6 is under strong purifying selection; it is not exceptional
within the family. dS is saturated throughout, so dN is the usable column."

**C17** → "The palindrome does not detectably bias C16. The decisive argument is timing: the
duplication is absent from CB4856, XZ1516, DL238 and MY2147, so it postdates their common
ancestor with N2 and occupies a negligible fraction of the *Cel*/*Cbr* branch on which dN was
measured. Two supporting bounds: N2-vs-wild divergence is equal or higher in the C-half, and
N2's C-half differs from the four palindrome-free isolates at only 1–2 of 864 nt, bounding any
N2-specific gene conversion. **The Fisher test of the N-half/C-half dN difference (p = 0.258)
is not evidence of no difference**: the point estimate is a 1.49× deflation in the C-half,
the 95 % CI on the rate ratio is 0.77–2.95, and the test has 20 % power at the observed
effect size and 80 % power only at a three-fold deflation."

**C25** → "361 parsimony-informative sites for 173 taxa on b120 is **0.96** sites per free
parameter (375 free parameters); 356 for 163 taxa on v1 is 1.1. Deep nodes should be
expected to carry modest support, and free-rate model parameters are multimodal at this
ratio, so log-likelihoods from independent runs are not directly comparable."

**C26** → append: "The 120 cut sits in the largest gap in the entire 142-gene monoamine-HMM
ranking (27.0 points, gar-1 124.4 → npr-35 97.4). The resulting panel is **invariant to any
cut in (97.4, 124.4]**, so the specific value 120 carries no information. The drop of
F59D12.1, T21B4.4 and F35H10.10 was made by `rebuild_boundary120.py`, which reads only the
score file and no tree or support value; it preceded every b120 tree, though the v1 rogue
test showing that two of those three depress support was already in hand."

**C28** → append: "at the current orthologue sampling, which is missing `Ppa_dop-5`
(A0A2A6CA79, W4) and `Cbr_dop-3`."

**C29** → "The clade containing DOP-6 with no non-nematode member is a 13-tip group spanning
five species, **all of them Rhabditina** (*C. elegans*, *C. brenneri*, *C. remanei*,
*C. briggsae*, *P. pacificus*); the first outside member, *Drosophila* Dop2R, appears at 16
tips. Support ranges **60 (gappyout) to 90 (gt0.1)** on b120. The panel's four other
nematodes — *Trichinella* (I), *Brugia* (III), *Strongyloides* (IV), *Necator* (V) — all fall
outside it, but Compara returned no dop-3/dop-5/dop-6 orthologue for any of them, which under
W5 is an unsearched absence rather than an established one. The clade is therefore
**Rhabditina-restricted in this panel**; whether the radiation predates Nematoda is not
determined, and cannot be until non-Rhabditina proteomes are searched by phylogenetic
placement rather than best hit."

**C31** → "The tree gives no positive support for an OA assignment, and disfavours it without
excluding it. Constrained AU tests (b120 gt0.1 and gappyout, 10,000 RELL replicates) forcing
the dop-6 orthologue clade to be monophyletic with each OA receptor group cost 23–28
log-likelihood units and return **p-AU 0.06–0.16 — none rejected at α = 0.05** — while
holding only ~12–15 % of the expected-likelihood weight. The reason the tree cannot exclude
OA is **not** that a nematode-specific lineage is where a phylum-specific ligand would sit:
octopamine is not phylum-specific, and this panel's four characterised OA-receptor lineages
all sit in cross-phylum clades at UFboot 97–100 (SH-aLRT 80–98). The reason is that **ligand is not a
topological property** — the project's own Fitch reconstruction counts 15 transmitter
changes across 56 labelled tips, so a ligand switch on the dop-6 stem is common-rate for this
family and invisible to any topology test."

**C34** → "...lift support by +9 to +20 points depending on which rogue is dropped (9-tip
clade 71 → 91 dropping F59D12.1, → 80 dropping T21B4.4, → 81 dropping both). The
non-monotonicity indicates part of the lift is search noise. C34 is a **v1-only** result."

### Should be withdrawn
- **The sentence at README line 249** ("Its nearest characterised relatives are two 5-HT
  receptors, both Gq") — this *is* W2, already withdrawn, still asserted.
- **The "block" vocabulary at README line 247** — this is W3's framing, already withdrawn.
- **README §6 as a whole** should be withdrawn or explicitly marked superseded: it presents a
  50-taxon panel, a different model, different column counts and different trait-signal
  numbers as current.
- **The reading paragraph at the foot of both `nematode_expansion.txt` files** ("they are a
  radiation inside the nematodes from a D2-like ancestor ... nothing outside the phylum can
  lend them a ligand or a coupling") — two independent overreaches: "nematodes" should be
  "Rhabditina", and "from a D2-like ancestor" is an ancestral-state claim that no ancestral-
  state reconstruction supports.

### Open, not reviewable by me
C39 (expression) is correctly flagged as unaddressed scope. The *P. pacificus* / Clade I–IV
orthologue search (MAJOR-4b) and the AlphaFold model of C24A8.6a (MINOR-3b) are **PI-only**
and **deferred-with-owner** respectively.
