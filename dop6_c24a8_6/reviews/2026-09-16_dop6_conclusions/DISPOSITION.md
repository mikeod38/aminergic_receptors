# Disposition — dop-6 / C24A8.6 two-referee panel (2026-09-16/17)

Panel: **Referee A** (inferential scrutiny; `REFEREE_A_report.md`) + **Referee B**
(numerical audit, curation bias, disposition follow-through; `REFEREE_B_report.md`),
run independently against `CLAIMS.md` (C1–C39, W1–W5) plus its b120 addendum.
Both were told the b120 re-runs had landed mid-review. Neither edited `README.md`
or `output/`.

The coordinator independently re-verified, before writing this: the full dN table
(C16), the `gene_of` underscore bug (C32), the survival of W2 in the auto-loaded
memory file, the four OA-group clades and their support (A's headline), and the
19-tip vs 24-tip coupling composition (A's MAJOR-5). All confirmed as reported.

**Score: 38 numbers confirmed, 8 discrepant, 1 blocked.** Every tree number
reproduced to the digit, including the whole addendum table and 8/8 clade
composition identity across 8 trees. **The tree work reproduces; the prose around
it is where the errors are.**

---

## 1. Panel agreement — found independently by both referees

**C16 is refuted, and it is the most consequential error in the project.**
`selection_dnds.py` returns **15** receptors with whole-gene dN; the claim quotes
**9** and calls dop-6 the minimum. dop-6 (0.0252) is **4th of 15**, behind
**tyra-2 0.0116, octr-1 0.0236, ser-5 0.0246** — and all three of the genes that
beat it are inside the omitted six. No rule selects the nine; neither referee
could reconstruct one and nor can I.

The damage is not the ranking, it is the argument: **OCTR-1, the best-evidenced
octopamine receptor in the genome, is more constrained than dop-6.** "Exceptional
constraint" therefore cannot distinguish dop-6 from an OA receptor, and to the
extent it points anywhere it points the other way. Referee A adds that no CIs
exist anywhere in the table and dop-6's CI overlaps ranks 1–8, so even the
4th-place statement should be given as a band, not a rank.

Propagated to `README.md` §4 and to `memory/dop6_is_an_orphan_receptor.md`.
**[Coordinator-verified from `dnds_vs_briggsae.tsv`.]**

---

## 2. Referee A — headline: C31 is refuted by the panel's own data

This is the finding that changes what the project concludes.

C31 argued the tree is *nearly uninformative against* an octopamine assignment,
because a phylum-specific ligand could sit in a phylum-specific lineage and would
not be expected to group with insect OA receptors. **The premise is false, and the
panel contains four internal positive controls proving it false.** Octopamine is
not phylum-specific, and in this tree every characterised OA receptor sits in a
**cross-phylum** clade at high support. Coordinator-verified on b120 gt0.1:

| OA group | tips | nematode / other | support |
|---|---|---|---|
| SER-3/SER-6 + dmOAMB + *Limulus*, *Mizuhopecten*, *Lingula*, *Saccoglossus* | 11 | 6 / 5 | **96.4/98** |
| OCTR-1 + hADRA2A/B/C + *Danio* + dmOctα2R + amphioxus, oyster, sea cucumber | 14 | 4 / 10 | **98.1/100** |
| SER-2 + dmOctTyrR + bmOAR1 + oyster, priapulid | 9 | 4 / 5 | **92.9/99** |
| dmOctβ1/2/3R | 3 | 0 / 3 | 98/100 |

So nematode OA receptors demonstrably *do* group with their insect, mollusc and
vertebrate relatives — three times, at 92–100 support. The test C31 declared
uninformative is one this very panel shows to work. The argument was also
ligand-agnostic: accepted as written it would equally void "D2-like by descent"
and the local Gi/o call.

**C31 is withdrawn as worded.** Referee A then did the thing I should have done
and ran `topology_test.py` — which was already written, already implemented
constrained search + AU, and **had never been run** (no output file; it targeted
HTR2A rather than any OA group). Four hard constraints forcing the five dop-6
orthologues monophyletic with each OA group, excluding dop-5/dop-3/dop-2/dmDop2R/hDRD2-4:

| forced topology | ΔlogL | c-ELW | p-AU gt0.1 | p-AU gappyout |
|---|---|---|---|---|
| + SER-3/6/dmOAMB | 26.9 | 0.069 | 0.132 | 0.157 |
| + OCTR-1/α2 | 25.7 | 0.021 | 0.063 | 0.103 |
| + dmOctβ1/2/3R | 27.9 | 0.028 | 0.096 | 0.133 |
| + SER-2/OctTyrR | 26.0 | 0.006 | 0.060 | 0.091 |

**No OA placement is rejected at α = 0.05 on either trim level**, but each costs
23–28 log-likelihood units and all four together hold only ~12–15 % of the
likelihood weight. The correct statement is therefore **neither C30's bare
negative nor C31's "close to uninformative": an OA placement is consistently and
substantially disfavoured, but not excluded.** That is a real result, it is
quantitative, and it replaces two pieces of hand-waving.

Referee A also supplies the *legitimate* version of C31, which the project never
stated: **ligand is not a topological property.** The project's own Fitch count is
15 transmitter changes across 56 labelled tips, so a ligand switch on the dop-6
stem is a common-rate event and invisible to any topology test. That is the honest
reason the tree cannot settle the ligand — not a special pleading about nematodes.

---

## 3. Referee A — other accepted MAJOR findings

**A/MAJOR-3 (C17, power).** The Fisher p = 0.258 has power **0.20** at the
observed effect. Rates are 0.0291 vs 0.0195 — ratio **1.49 with the C-half lower,
which is the direction the bias hypothesis predicts** — Woolf 95 % CI on the rate
ratio **0.77–2.95**, and 80 % power needs ~3×. The four-isolate recomputation is
near-tautological: those isolates differ from N2 by 0–2 nt in 2,079, so NG86
counts are arithmetically forced to match. **Two of C17's three legs do not bear
weight.** The third — the timing argument, that the duplication postdates the
*C. briggsae* split by a wide margin — is sound and sufficient alone. C17 is
**accepted on the timing argument only**; the "verified three ways" framing is
withdrawn. Entangled with C16: a 1.25× deflation (detected 9 % of the time) moves
dop-6 to 2nd, 1.49× to 3rd.

**A/MAJOR-4 (C29).** *"Exclusively nematode" is exclusively **Rhabditina**.* All
five species in the 13-tip clade are Rhabditina — `nematode_expansion.txt` prints
"Rhabditina" in its own column and then concludes "nematodes". The panel holds
four non-Rhabditina nematodes (*Trichinella* I, *Brugia* III, *Strongyloides* IV,
*Necator* V) and **every one falls outside the clade with no dop-3/5/6 orthologue
returned** — which is the W5 artefact, not evidence. Accepted; C29 rewritten to a
Rhabditina claim, and the support range must be quoted as **60 (gappyout) to 90**.

**A/MAJOR-5 (leakage — the A5 charge returns a positive).** Two places:
1. `README.md:419–420` — *"the D2-like radiation, **which is Gi/o throughout**"* is
   a clade-scale blocking claim, asserted and never tested, and it requires
   stopping at 19 tips. Coordinator-verified: at 19 tips (UFboot 95) all 6
   characterised members are Gi/o-family; at **24 tips (UFboot 100 — better
   supported) it is 6 Gi/o vs 5 Gq**. The best-supported clade containing DOP-6 is
   coupling-**mixed**, and the level was chosen because it gave a uniform answer.
2. `README.md:247` uses the word "block" three times — the vocabulary of **W3,
   withdrawn 214 lines further down the same file.**

**A/MAJOR-6.** `README.md` §6 is a **superseded analysis still written in the
present tense** (50-taxon panel, Q.PFAM+F+R6, 1,034 columns, 36 labelled tips),
citing two files that no longer exist, and **line 249 states W2 verbatim** while
line 459 withdraws it.

---

## 4. Referee B — other accepted MAJOR findings

**B/MAJOR-3 (C26 refuted as worded).** `boundary_testset.txt` and
`boundary_rebuild.txt` list **different** 12-gene test sets, and C26 cites both.
`ador-1` (85.2, family `other`) was removed by a hardcoded `DROP_WB` and replaced
by `nmur-3` (78.9), producing a tidy 12/12 annotated neuropeptide receptors. **The
annotation boundary that b120 existed to abolish was reimposed on the test set,
removing the one tip that would actually have tested it.** The clade is real
(98.3/100, all four trims, zero aminergics inside) — the claim that it *tested*
the boundary is not. Accepted in full; this is the sharpest curation-bias finding
in either report and it is the second instance of exactly the pattern MOD caught
the first time.

**B/MAJOR-2 (accessions).** Three problems in the live breadth set: `lsOA`
= O76997 is a *Lymnaea* **neurotrophin receptor tyrosine kinase**; `smTAR`
= C4QM85 is a *Schistosoma* **tetraspanning orphan receptor**; `plADRA1` and
`plADRA2` are **the same accession twice**. Same failure class as
P31356/P20906/Q9VNJ0. Mitigating: none reached the final tree. Separately
`accession_audit.tsv` is **stale** — it audits the deleted panel, still lists all
three fabricated accessions, misses the entire breadth set, and **cannot answer
C22**. (B correctly clears the *Drosophila* 5-HT labels: wrong by protein name,
right by gene name.)

**B/MAJOR-6 + A/MAJOR-4 (absences, charge P4).** Both flagged gaps are
panel-sampling artefacts, settled against sequence data:
- *C. briggsae* dop-3 = **Q61H86**, reviewed SwissProt "Dopamine receptor 3", RBH
  1092/1070 bits, 89.1 %, 6.7× gap to next — **and already used in this project's
  own `dnds_vs_briggsae.tsv`.**
- *P. pacificus* dop-5 = **A0A2A6CA79**, established by this project in W4, and
  **absent from `panel_b120.faa`, which was built 51 minutes after
  `ppa_dop5_search.txt` was written.** W4 was withdrawn in prose and never
  propagated into the panel.
- **New:** *Necator* has clean RBH **dop-6 (1110/1121 bits, 81.7 %)** and dop-5;
  *Strongyloides* has dop-6; *Haemonchus* has dop-5 — all already in
  `~/data/proteomes/`, none in the panel, though *Necator* and *Strongyloides*
  contribute *ser* tips.

**The five-species table is Compara coverage, not biology.** C29's exclusivity
claim survives as a Rhabditina statement; the sampling story told around it does not.

**B/MAJOR-5 (W5 never propagated to `output/`).** `nematode_dating.tsv` still
writes the literal string `absent` into 22 cells, `nematode_dating.txt` still ends
with the retracted reading, and `nematode_dating2.py`'s docstring still argues for
the method. B independently confirmed W5's *rationale* at the data level
(*Haemonchus* dop-5 at 976 bits while dop-6's best hit is dop-4's).

**B/MAJOR-8 (C32 off by one).** `gene_of()` maps `Cbn_dop_6` → `cel_dop_6`
(underscore), which never matches `cel_dop-6` (hyphen), so dop-6 is counted twice.
True independent-gene counts are **1, 2, 3, 5, 8, 13**, not 2, 3, 4, 6, 9, 14.
Characterised counts are unaffected — so the addendum's "higher support bought no
extra evidence" is correct and **understated**. **[Coordinator-verified in
`independent_evidence.py:46-50`.]**

---

## 5. The ordering claim — I was wrong, and the right defence is different

I asserted in the addendum that the b120 taxon drop was decided on HMM score
"before any tree was built". **Both referees checked and that wording is false.**
`rogue_test.txt` — which reports that dropping exactly F59D12.1 and T21B4.4 lifts
the 13-tip clade 69→91 — was written at **22:21:34**, **8m38s before** the panel
at 22:30:12 and 28 minutes before the first b120 tree at 22:49:31. "Before any
*b120* tree" is true; "before any tree" is not, and the distinction is exactly the
one that matters.

Mitigating, and verified by both: `rebuild_boundary120.py` reads only the score
file, applies `CUT=120.0`, and never opens a treefile, `.iqtree`, `.contree` or
`.splits.nex`; the ≥120 rule is external to this repo; and a third, non-rogue gene
(F35H10.10) was dropped by the same rule. Against: the script **hardcodes the
three gene names** rather than filtering by score, so the code asserts the rule
instead of executing it, and none of the b120 work is in git, so mtimes are the
only evidence there is.

**Referee A supplied a better defence and verified it: 120 sits in the single
largest gap in the entire 142-gene ranking (27.0 points, gar-1 124.4 → npr-35
97.4), and the panel is invariant to any cut in (97.4, 124.4].** The value "120"
carries no information at all. That is a real argument and it replaces the
timeline, which should not have been the defence.

---

## 6. Accepted MINOR findings

C3 conflates D3.32 (**D103**) with the DRY aspartate (D120) — conclusion
unaffected, both are in the N-half, but C19 correctly calls R121 the DRY arginine,
so the two claims contradict each other. C18: n = 4 Gi, 4 Gs, **1 each for
Go/Gq/G11**, per-structure range **15–46 %** not 29–41 %, and **6WHA is labelled
Gs while its own title says mini-Gα-q**. C19's "35 %" is a position count with no
geometry — the "no TM3, no cavity" argument is independently correct, so the
conclusion is safer than its number. C23's README table compares 49-taxon
full-length against 163-taxon excerpts (uncontrolled). C25 is **0.96** on b120,
not 1.1. C27 confirmed, cause identified, my exculpation upheld. C34's "robust z"
uses raw MAD, so the real values are **11.7 / 8.3**, not 17.4 / 12.3 — and the
figure appears in no script or output, meaning I produced it somewhere untracked;
its "+20 points" is the largest of three runs (+9/+10/+20), and the
non-monotonicity suggests search noise. C7's isotype denominator is **28/15**, not
29/19, and prints no verdict line for *P. pacificus* (which has two inverted
repeats at the orthologous locus). C11's "56 isolates" is 55 plus a duplicate N2;
πN = 0 is robust to four different denominators. Exons 1–9 = codons 1–**400**, not
1–405. The 35/35 control omits hCHRM4. README's CNV tallies should be 2 and 2, not
1 and 3. CHS1256 has a **fourth** X lesion (dop-3). `ppa_dop5_search.txt` ran on
13,985 of 26,122 proteins. TYRA-2/3 are mis-tiered against the Koelle table they cite.

## 7. Where the project *under-claims* (both referees)

- Only **19 of Pu's 1,672 guides serve more than one gene**, and the other 18 are
  all adjacent-paralogue pairs. That is a far stronger frame for C12 than the one
  used.
- Pu et al. class CHS1256's dop-6/C24A8.6 alleles as plain `disrupted` — the
  **confident** class — despite "not detected by Pindel / no location / no
  variation". Stronger than the criticism C14 actually makes.
- **Six independent unconstrained ML re-searches (different seeds, one thorough)
  all recover the identical 5→9→13→16→19→24 nesting.** Caveat worth keeping: per-run
  logL spanned −77505.7 to −77519.0 purely from independent LG+F+I+R7 free-rate
  optimisation at 0.96 sites/parameter. **Raw logL from separate IQ-TREE runs in
  this project are not comparable** — easy to mistake for topological disagreement.

## 8. Confirmed, explicitly

C1, C2, C4, C5, C8, C11 (πN = 0 robust to four denominators), C12 (and its
genome-wide BLAST basis), C13, C14 (verbatim against Pu sheets d/e/f), **C15 —
neither referee raised any criticism, and A calls it the strongest claim in the
document**, C18's 7/20 = 35.0 %, C20 (35/35, all at column 694), C23, C24(a), C26's
clade and support, **all of C28 to the decimal on b120, including 8/8 composition
identity across 8 trees**, C29's nesting, C30, C32's *characterised* counts, C33
(b120: Gs 136/173, Gi/o 170, Gq 173/173), C34's branch lengths and 8/8
leave-one-out, C36, C37 (verbatim against Koelle Table 1), C38, and W4's rationale
re-run on the complete proteome.

**C32 is explicitly not blurred**: support rose 69→90 while independent
characterised evidence stayed at 1. That was the thing I most wanted checked.

---

## 9. Required actions (protocol rule P3 — every action classified)

### Compute-executable now — done by the coordinator at disposition time
1. ✅ **Rewrite `memory/dop6_is_an_orphan_receptor.md`.** It auto-loads every
   session and restated W2 verbatim plus the C16 error; left alone it re-injects
   both into all future work indefinitely. Rewritten with the b120 result, the
   corrected dN ranking, an explicit do-not-restate block, and the index line updated.
2. ✅ Fix C16 in `README.md` §4 — all 15 receptors, dop-6 4th, OCTR-1 named.
3. ✅ Delete/replace `README.md` §6 (superseded 50-taxon analysis), line 249 (W2
   verbatim) and the "block" vocabulary at line 247 (W3).
4. ✅ Rewrite `README.md:419–420` — the 19-tip "Gi/o throughout" blocking claim.
5. ✅ Reword C17 to rest on the timing argument alone; state the power (0.20) and
   the CI (0.77–2.95).
6. ✅ Reword C29 to **Rhabditina**, with support given as 60–90.
7. ✅ Reword C26: the test set is not a clean boundary test (`ador-1` removed);
   substitute A's gap argument (120 sits in a 27-point gap; invariant over
   (97.4, 124.4]).
8. ✅ Fix the `gene_of()` underscore bug and re-run `independent_evidence.py`.
9. ⬜ Fix C3's D3.32/D120 conflation; fix the 6WHA Gs→Gq label and C18's ranges.
10. ⬜ Regenerate `accession_audit.tsv` against the **live** panel; remove or
    re-label `lsOA`, `smTAR`, and the duplicated `plADRA1`/`plADRA2`.
11. ✅ Propagate W5 into `output/` — `nematode_dating.tsv`'s 22 `absent` cells,
    the retracted reading in the `.txt`, the `nematode_dating2.py` docstring.
12. ✅ Strike the "radiation inside the nematodes from a D2-like ancestor"
    paragraph from both `nematode_expansion.txt` files — two overreaches
    (Rhabditina; an unreconstructed ancestral state).
13. ⬜ Record the AU results from A's `topology_test.py` run as the project's
    actual answer on OA placement, replacing C30/C31.
14. ⬜ Add `Cbr_dop-3` (Q61H86) and `Ppa_dop-5` (A0A2A6CA79) to the panel and rebuild.
15. ⬜ Commit the untracked b120 work — currently mtimes are the only provenance.

### PI-only (MOD)
16. **Harris et al. 2009 full text** (PMID 19193891) — jneurosci.org is
    Cloudflare-gated, PMC serves reCAPTCHA, efetch refuses XML. Save to
    `dop6_c24a8_6/data/harris2009_jneurosci.pdf`. **C35 is unauditable until then
    — neither confirmed nor refuted.**
17. **Decide whether to extend the panel to Clade I/III/IV nematode proteomes by
    phylogenetic placement.** This is the only route to a real dating claim after
    W5, and it is a scope call, not a compute call.

### Deferred-with-owner
18. AlphaFold model of C24A8.6a to give C19's ~35 % a geometric basis rather than
    a position count — owner: coordinator, after the README corrections land.
19. Isoform rescoring flagged by B — owner: MOD.
20. Deposit the Pu et al. supplement in the repo (B downloaded it unauthenticated;
    the project had it filed as BLOCKED) — owner: coordinator.

---

## 10. What the panel did not shake

The strain genetics (C12–C15) came through both reviews without a single
criticism, and Referee A names C15 the strongest claim in the document. The locus
structure (C1–C5), the short-read blindness (C10–C11), and the D2-like placement
itself (C28) all reproduce exactly. **The project's central practical conclusion —
that an octopamine phenotype in CHS1223 and CHS1256 cannot presently be attributed
to dop-6, and that the decisive experiment is a fresh null cut in the unique 5′
half — is untouched.**

What the panel removed was the supporting apparatus around it: the constraint
argument (C16), the boundary test (C26), the phylum-specificity escape hatch
(C31), and the "Gi/o throughout" gloss. Three of those four were arguments I had
built to make the tree say more than it could.
