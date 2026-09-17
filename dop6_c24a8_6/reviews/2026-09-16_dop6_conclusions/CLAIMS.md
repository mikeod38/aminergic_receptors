# Claims under review — dop-6 / C24A8.6 sub-project (2026-09-16)

This is the claim set put to the two-referee panel. Every claim is numbered so
reports can address them by ID. Each carries a pointer to the file that produced
it. Where I already know a claim is weak, that is stated — referees should not
have to rediscover a weakness I am aware of, and should spend their budget on the
ones I believe are solid.

**Repo root for all paths:** `/Users/mikeodonnell/git/Malaiwong_Oglu_2026/dop6_c24a8_6/`

## Background in one paragraph

Two CGC strains from Pu et al. 2023 (*Nat Commun* 14:8410, PMID 38110404),
CHS1223 and CHS1256, show an octopamine-related behavioural phenotype. The
project began from the observation that the only genes those two strains share
are `dop-6` (C24A8.1) and `C24A8.6`. Three goals: (a) is DOP-6 an octopamine
receptor; (b) where is it expressed; (c) how can dop-6 be manipulated without
also altering C24A8.6.

## Panel provenance — read this before checking any tree number

Two tree panels exist and their numbers are **not** interchangeable:

| panel | taxa | built by | trees in | boundary |
|---|---|---|---|---|
| **v1 rebuild** | 163 | `assemble_final_panel.py` + `build_tree_final.py` | `output/tree/trees/` | annotation-based |
| **b120** | 173 | `rebuild_boundary120.py` + `build_b120.py` | `output/tree/b120/trees/` | monoamine HMM ≥ 120, plus 12 neuropeptide-receptor test tips |

b120 exists because the v1 boundary was drawn on WormBase annotation labels
rather than on a score, and 36 excluded genes scored higher on the monoamine HMM
than `F59D12.1`, which was inside. b120 drops `F59D12.1`, `T21B4.4` and
`F35H10.10` and adds 12 `TEST_*` neuropeptide receptors to test the new
boundary from the outside.

**Known staleness at the time of writing:** `output/tree/confirm.txt`,
`independent_evidence.txt`, `nematode_expansion.txt`, `rogue_test.txt` and
`resolution.txt` were all generated on the **v1 163-taxon** panel. The b120
re-runs of the first three are regenerating into `output/tree/b120/` as the
panel convenes. Referees should treat any conflict between a v1 file and a b120
file as expected, and should flag any claim below whose support depends on which
panel it was computed from.

---

# A. Locus structure

**C1.** The dop-6 locus contains a perfect inverted repeat on chromosome X:
arm A 4,324,178–4,328,106, spacer 4,328,107–4,329,302 (1,196 bp, contains
C24A8.5), arm B 4,329,303–4,333,231. Arm B is the reverse complement of arm A at
**3,929/3,929 = 100.000 %, zero mismatches**, with sharp junctions.
→ `output/repeat_structure.txt`, `analyze_duplication.py`

**C2.** `C24A8.6a` (287 aa) is `DOP-6` residues 406–692 **exactly, 287/287
identical**; C24A8.6 has **0 bp of unique sequence** in its 2,330 bp. Its start
codon is DOP-6's internal Met406.
→ `output/dop6_c24a8_6_proteins.faa`, `analyze_duplication.py`

**C3.** The duplicated segment is end-of-ICL3 + TM6 + ECL3 + TM7 + C-tail. The
ligand pocket is **not** in it: D3.32 and the DRL motif (D120-R121-L122) sit at
~residue 120, in the single-copy half. So C24A8.6 can encode at most a two-TM
fragment with no binding pocket.

**C4.** The two mRNAs carry the **identical** 864 nt sequence, not reverse
complements — so no dsRNA/endo-siRNA mechanism is available and reads from either
gene are the same string.

**C5 (operational consequence).** C24A8.6 cannot be measured apart from dop-6 by
any short-read method, cannot be mutated without also cutting dop-6, and has no
diagnostic base anywhere. dop-6 *is* fully tractable through exons 1–9.
→ `output/specific_tools.txt`, `output/dop6_guides.tsv`

# B. Is it real, and how common

**C6.** Not an assembly artefact. Present in three independent N2 assemblies
(Sanger cosmid, U. Oregon long-read, Wellcome Sanger). Not an assembler collapse
in non-carriers either: the U. Oregon N2 and CB4856 assemblies share BioProject,
pipeline and release date and differ in copy number.
→ `output/assembly_copy_number.txt`

**C7.** Across **56 assemblies**: 29 carriers / 19 non-carriers among natural
isolates; all 8 N2-derived lab strains carry it. Absent from *C. briggsae*,
*C. brenneri*, *C. remanei*, *C. japonica*. Replicate assemblies of the same
isotype agree (CB4856 1 copy ×3, DL238 1 copy ×2, ECA36 1 copy ×2, EG4725 2
copies ×2).
→ `output/assembly_survey.tsv`, `survey_assemblies.py`

**C8.** Carriage is **identity-by-descent, not independent**. Carriers are
99.977 % identical to N2 in the flanking single-copy windows vs 99.477 % for
non-carriers, while both look like ordinary wild isolates on chromosomes II and
III. The carrier-specific N2-like block runs about −1 Mb to +1 Mb, has decayed by
±3 Mb, gone by ±6 Mb. **The 29 carriers are essentially one haplotype observed 29
times.** Whether the duplication is under selection or hitchhiking cannot be
separated with these data. Block length must not be converted into a date
(*C. elegans* selfs).
→ `output/haplotype_context.tsv`, `output/haplotype_extent.tsv`

**C9 — previously withdrawn, kept here for audit.** An earlier claim that the
palindrome is N2-lineage-specific, made from 4 wild isolates, was **withdrawn**
when the panel was extended to 56.

# C. Short reads are blind here

**C10.** CaeNDR 20250625, 684 isotypes, variants per kb: unique flank 24.75,
spacer 34.28, arm A 3.82, arm B **0.25**. This is a MAPQ-0 mapping-ambiguity
desert, not a conserved region.
→ `output/wild_isolate_variants.tsv`

**C11.** Rebuilding the dop-6 CDS exon-by-exon from the assemblies recovers what
mapping cannot: segregating sites in the C-half exons are **0** from 684 short-read
isotypes but **7** from 56 assemblies. **πN = 0 in both halves** — all 10
segregating sites are synonymous across 56 isolates. Carriers 1 site vs
non-carriers 8, which is the C8 haplotype signal, not extra constraint.
→ `output/dop6_polymorphism.tsv`

# D. The strains (the most consequential result)

**C12.** dop-6 and C24A8.6 were **never targeted independently**. Pu et al.
Supplementary Data 1d gives *the same* guide, PAM and ssODN for both
(`CAAGTTGGTTGCAACTGTAA`, PAM TGG). That guide has **exactly two perfect
full-length matches in WBcel235**, X:4,325,604–4,325,623 (arm A, in dop-6) and
X:4,331,786–4,331,805 (arm B, in C24A8.6, minus strand), as it must, since the
arms are identical. Only the genotyping primers differ, because those sit in the
unique flanks. So `yum2349`/`yum2350` are the two ends of one event, not two
alleles.

**C13.** CHS1223 is a complex rearrangement, not a double null. Sheet f, both
alleles class *"disrupted larger CNV"*; `yum2349` = inversion with 50 bp
insertion at X:4,325,629–4,330,411; shared comment ***"4kb deletion into dop-6
and C24A8.5, plus 75kb TD spanning multiple genes"***.

**C14.** CHS1256's dop-6 and C24A8.6 alleles were **never actually detected** —
sheet f records both as *"not detected by Pindel"*, *"manually curated via IGV"*,
no mutation location, no detected variation, in exactly the MAPQ-0 interval of
C10. CHS1256 also carries `dop-1(yum2598)` at X:4,415,971, ~90 kb from the dop-6
cut sites — three breaks inside 90 kb of X in one injection — and two further
CNV-class alleles (`ser-6`, `dop-4`).

**C15 (the conclusion that matters).** The shared-gene argument that motivated
the project does not survive C12–C14. **An octopamine phenotype in these strains
cannot presently be attributed to dop-6.** The decisive experiment is a fresh
dop-6 null cut in the unique 5′ half with a single-target guide, outcrossed, plus
cDNA rescue.

# E. Selection

**C16.** dN vs *C. briggsae* (NG86): dop-6 **0.0252**, the lowest of the
aminergic receptors tested (dop-5 0.0261, ser-2 0.0290, dop-1 0.0322, ser-3
0.0500, dop-2 0.0508, ser-4 0.0609, dop-3 0.0618, dop-4 0.0733). dS is saturated
throughout, so dN is the usable column.
→ `output/dnds_vs_briggsae.tsv`

**C17.** The palindrome does not bias C16, checked three ways: N2-vs-wild
divergence is equal or higher in the C-half; recomputing dN from four
palindrome-free isolates' own CDSs gives identical values; the N-half/C-half dN
difference is not significant (Fisher p = 0.258). The duplication postdates the
*C. briggsae* split by a wide margin.
→ `output/mirror_bias_check.txt`

# F. Could the duplicated C-terminus bind a G protein alone?

**C18.** Interface defined from ≤4.5 Å contacts in **11 receptor–G protein
complexes** spanning Gs, Gi, Go, Gq, G11, with receptor and Gα chains identified
by sequence rather than chain label. Consensus contacts (≥6 of 11) = 20
positions; the fragment retains 7, loses 13 → **~35 %**; per class 29–41 %, no
class rescues it.
→ `output/gprotein_interface_multi.tsv`

**C19.** What is lost is decisive: **R121, the DRY arginine, contacts Gα in 11 of
11 structures**, plus the whole ICL2 patch and the TM5/ICL3 contacts. The α5
helix docks into a cavity opened by TM6 swinging away from TM3; with no TM3
there is no cavity. The live hypothesis is therefore not autonomous signalling
but complementation/competition with full-length DOP-6, or GRK/arrestin
sequestration via the serine-rich tail.

# G. The ligand pocket

**C20.** Position 3.32 was located independently from three ligand-bound
structures (2RH1/CAU/hADRB2 → Asp113; 4IB4/ERM/hHTR2B → Asp135; 3PBL/ETQ/hDRD3 →
Asp110) which **all three map to one alignment column (591)**. Positive control:
**35/35 human aminergic receptors return D**.
→ `output/pocket_332_v2.tsv`, `pocket_positions_final.py`

**C21.** Calls: dop-6 = **D** in all five orthologues; **dop-5 = G** in all five
including *P. pacificus*; F59D12.1 = Q; T21B4.4 = L; melatonin receptors M/M/T.

**C22 — previously withdrawn, kept for audit.** Three accessions in an earlier
version of the pocket analysis were hardcoded from memory and were a squid
rhodopsin (P31356), a *Pseudomonas* decarboxylase (P20906) and a *Drosophila*
ubiquitin ligase (Q9VNJ0). Caught by the 35/35 positive control after MOD said
"the binding residue of the 5-HT receptors should be known". `verify_accessions.py`
was written in response. **Referees should check whether any accession anywhere in
the current panel is still unverified.**

# H. Panel construction and the alignment

**C23.** The alignment is built from **7TM excerpts, not full-length sequences**.
Each sequence is scanned with the monoamine HMM and its HMM blocks are
concatenated with a 15-residue margin; ICL3's middle is excised. Justification:
ICL3 ranges 14–821 residues in this panel (DOP-6's is ~240), and full-length
alignment gave 84.4 % gaps with all 164 sequences above 50 %. Sequences of 413,
692 and 1256 aa converge on 254–263 aa excerpts. Effect: gaps 48.5 % → 11.4 %,
parsimony-informative 76 % → 98 %.
→ `extract_domains.py`, `output/tree/domain_report.txt`

**C24.** Three checks that the excision is safe: (a) **all 20 consensus Gα-contact
positions fall inside retained blocks** in DOP-6, because ICL3 contacts Gα at its
ends where it leaves TM5 and enters TM6, and those ends are the conserved helix
boundaries the HMM recognises; (b) D3.32 survives in every receptor known to carry
it; (c) ICL3's own signal was measured — its **length** is phylogenetically
conserved (mean sister-tip |Δlength| 52.7 aa vs 110.8 shuffled, p = 0.0001) but its
**sequence** stops being comparable with depth (76.8 % between orthologues, 28.4 %
between *C. elegans* paralogues, 26.4 % worm vs human).
→ `output/gprotein_survival.txt`, `output/tree/icl3_information.txt`

**C25 (stated limit).** 356 parsimony-informative sites for 163 taxa is **~1.1
sites per free parameter** (323 branches on an unrooted 163-taxon tree). Deep
nodes should be expected to carry modest support.

**C26 — the boundary test.** The b120 boundary (monoamine HMM ≥ 120) was tested
from the outside by including the 12 highest-scoring excluded genes as `TEST_*`
tips (npr-35 97.4 down to nmur-3 78.9, all annotated neuropeptide receptors).
**All 12 form one 12-tip clade at 98/100 support with 0 aminergic receptors
inside it.** Claim: the boundary is tested, not assumed, and there is no
unannotated OA receptor sitting in the top of the excluded band.
→ `output/tree/boundary_testset.txt`, `output/tree/boundary_rebuild.txt`

**C27 — known defect, flagged by me.** `output/tree/boundary_rebuild.txt` prints
the three dropped genes with `HMM nan` although `boundary_testset.txt` gives
F59D12.1 = 52.8. The score lookup in `rebuild_boundary120.py` is failing for the
dropped set. This does not change which genes were dropped (that was decided on
the scores in `boundary_testset.txt`) but **the referee should confirm that and
not take my word for it.**

# I. Where DOP-6 falls

**C28.** dop-6 orthologue clade (Cel, Cbr, Cbn, Cre, Ppa) is recovered at
**UFboot 100 at every trim level on both panels**. Nesting outward on b120
(SH-aLRT/UFboot, gt0.1): +dop-5 95/95, +dop-3 85.7/90, +dop-2 88/90. The same
nesting on gappyout is 96.1/84, 83.3/60, 89/62 — **identical topology and tip
counts at every level (5→9→13→16), trim-dependent support**.

**C29.** DOP-6 sits in the **D2-like radiation**. The dop-3/dop-5/dop-6 clade is
**exclusively nematode**: the deepest clade containing DOP-6 with no non-nematode
member is the 13-tip group, and the first outside member (*Drosophila* Dop2R)
appears only at 16 tips.
→ `output/tree/nematode_expansion.txt`

**C30.** **DOP-6 does not group with any octopamine receptor.** The OA receptors
sit in three separate places (OCTR-1 with α2-adrenergic/Gi-o; SER-3/SER-6 with
dmOAMB; dmOctβ1/2/3R with β-adrenergic/Gs) and DOP-6 is in none.

**C31 — the honest limit on C30, and I believe this is the single most
review-worthy claim in the document.** The tree's failure to place DOP-6 with the
OA receptors is **weaker evidence against an OA assignment than it appears**: a
nematode-specific lineage is exactly where a phylum-specific ligand could sit, and
such a lineage would not be expected to group with insect OA receptors. So the
tree gives *no positive support* for OA, but is close to uninformative *against* it.

**C32 — evidence is not clade size.** About 40 of 163 tips are Compara orthologues
that inherit ligand and coupling by descent. Recounted in independent genes: the
5-, 9- and 13-tip clades contain **0, 0 and 1** characterised receptors
respectively. DOP-6's first characterised neighbour is DOP-3, four genes out.
→ `output/tree/independent_evidence.txt`

**C33 — trait signal is real but is clustering, not blocking.** Fitch parsimony vs
10,000 permutations, all trim levels: ligand 15 changes vs 35.2 random, coupling 8
vs 23.7, both p = 0.0001. But the smallest clade containing all Gs receptors holds
**125 of 163 tips**, Gi/o 160, Gq all 163. Receptor *subtype* groups are
monophyletic; *transmitter* classes are not. Any Gi/o inference for DOP-6 is a
**local** claim about its own clade.

**C34.** `cel_F59D12.1` has the longest terminal branch in the v1 tree, 5.98
subs/site, robust z = **17.4** against a median of 0.348; `cel_T21B4.4` second at
4.35 (z = 12.3); DOP-6's own is 0.020. Leave-one-out rebuilds (drop F59 / drop T21
/ drop both, same alignment, model fixed) change **no topology — 8/8 clades
identical** — but lift support by roughly 20 points (9-tip clade 71 → 91).
F59D12.1 has three further strikes: lowest monoamine HMM score of the *C. elegans*
set (52.8), **Q** not D at 3.32, and only 318 aa.
→ `output/tree/rogue_test.txt`

**C35.** The Harris et al. 2009 precedent for reading coupling off a tree is
weaker than its citation implies: 5-HT receptors only, rooted on GAR-1, ClustalW
in MegAlign hand-tuned, neighbour-joining in PAUP, no ML, no model selection, and
**hypervariable regions including ICL3 deleted** — which is where part of the Gα
interface lives. Its validation is one sentence of eyeball assessment with no
statistical test. C33 is the formal version of the same intuition and supports it.

# J. Literature re-reads

**C36.** Sugiura et al. 2005 (PMID 16001968) is **not a cAMP measurement** — it is
a CRE-luciferase transcriptional reporter in transiently transfected CHO cells, 3 h
with 1 µM forskolin. In DOP-3 cells **5-HT beat dopamine on both axes**: pEC50
8.90 vs 7.57 and ~48 % vs ~40 % maximal inhibition. The authors discount 5-HT
because untransfected CHO cells also respond (pEC50 8.96, n.s.) via endogenous
5-HT1B — but that compares **potency only**; in their own inset, control cells
bottom out near 90 % (~10 % reduction) against ~52 % (~48 %) in DOP-3 cells, a
~5-fold efficacy difference in the direction of the transfected receptor. The D87N
(D3.32) and S167A mutants that carry dopamine's specificity claim were **never
tested against 5-HT**. Conclusion: "DA > TA > OA" should not be cited as a
selectivity order. DOP-3's Gαo coupling is separately weak — IMP from a mutant
phenotype (Chase, Pepper & Koelle 2004), not a coupling assay.

**C37.** Koelle 2018 WormBook (NBK535514) Table 1 lists `C24A8.1` under
"Poorly-studied putative small-molecule GPCRs" with every ligand, coupling and
assay column empty, **and lists `C24A8.6` separately** — propagating the reference
artefact into the canonical orphan list. Of Koelle's six orphans, four are now
resolved: `ador-1` is an adenosine receptor *homolog* (a name, not a homology
result), `pcdr-1` is a phenotype name, `F35H10.10` is **class C** (7tm_3 141.6 +
ANF/VFT 57.9, no class-A signal), and `C24A8.6` is not a gene.

**C38.** DOP-6 has **no pharmacology of any kind**. It was named by homology in
Chase & Koelle 2004.

# K. Expression (goal b — largely unaddressed)

**C39.** Reported expression from a transgenic reporter (PLoS Comput Biol
2016, `10.1371/journal.pcbi.1005283`, Expr11993): OLL, plus AUA, RID, RMD, RIB,
ASI, PHA, IL2, PVQ, URA, AVF, ADF, RIH, URX. MOD notes the published images show
**gut expression** as well. **No analysis of my own has been done on goal (b).**
Referees should note this as an open scope item rather than review it.

# L. Already withdrawn — do not re-derive, but do check the withdrawal stuck

**W1.** "The palindrome is N2-lineage-specific" (from 4 isolates) — withdrawn, see C9.

**W2.** "DOP-6's nearest characterised relatives are Gq-coupled 5-HT receptors" —
withdrawn; rested on UFboot 37–62 in a 48 %-gap alignment with 2 of 13 human
serotonin receptors present.

**W3.** "The canonical G-protein blocks hold" — withdrawn; asserted from a
thumbnail figure, refuted by the monophyly test in C33.

**W4.** "dop-5 is absent from *P. pacificus*" — withdrawn; direct proteome search
finds a clean reciprocal best hit (`A0A2A6CA79`, bits 764). Compara missed it. The
first run of that search used a **truncated download** (15,999 of 26,122 proteins),
so any negative from it is unusable.

**W5 — the nematode dating matrix — withdrawn as a method artefact.** RBH cannot
resolve members of a recent paralogue radiation: in *Haemonchus* one protein
reciprocates with dop-5 at bits 969 while reverse hits pile onto dop-4 and dop-1
indiscriminately. Every "absent" in that matrix means "reciprocated to a different
paralogue". Only the dop-2 rows survive (present in all six species, consistent
with dop-2 being the ancestral lineage shared with *Drosophila*). **Dating the
dop-3/dop-5/dop-6 radiation requires phylogenetic placement, not best-hit
assignment, and has not been done.**
→ `output/nematode_dating.txt`, `nematode_dating2.py`

---

# Charges

## Referee A — inferential scrutiny

A1. For each of C1–C38, does the cited evidence support the claim **as worded**?
Name every claim where the wording outruns the evidence, and give the wording you
would accept instead.

A2. **C28–C34 are the core.** Is "DOP-6 sits in the D2-like radiation" supportable
at the stated support values, given C25 (1.1 informative sites per free parameter)?
Is the 13-tip nematode-exclusivity claim (C29) robust to the taxon sampling that
produced it?

A3. **C31 is the claim I most want attacked.** I have argued that the tree is
nearly uninformative *against* an OA assignment because a phylum-specific ligand
could sit in a phylum-specific lineage. Is that reasoning sound, or is it an
unfalsifiable escape hatch that would excuse any negative result?

A4. C18–C19: is a 20-position consensus interface from 11 structures an adequate
basis for a ~35 % retention figure? Is the "no TM3, no cavity" argument a real
structural argument or a plausible-sounding one?

A5. C33 vs C28: the trait-signal permutation test establishes clustering. Does any
inference actually drawn in this project depend on the stronger blocking claim
that C33 explicitly denies? Check for leakage.

A6. Statistical specifics: the Fisher p = 0.258 in C17 — what is its power? A
non-significant difference in a low-count comparison is not evidence of no
difference, and C17 is worded as if it were.

## Referee B — audit

B1. **Reproduce the numbers.** Independently recompute, from the output files, at
least: C1 (3,929/3,929), C2 (287/287), C7 (29/19 of 56), C11 (0 vs 7 segregating
sites; πN = 0), C16 (the dN ordering), C18 (7 of 20 retained), C20 (35/35), C26
(the 12-tip clade and its support), C28 (the nesting and support on both panels),
C34 (z = 17.4, 8/8 clades). Report every discrepancy, however small.

B2. **C27 is a known defect — confirm or refute it**, and determine whether the
`nan` affected which genes were dropped.

B3. **Curation-bias hunt.** Why is each taxon in the panel, and each excluded one
out? The boundary moved once already (annotation → HMM ≥ 120) because MOD caught
that 36 higher-scoring genes were excluded while a lower-scoring one was in. Is
there a second instance of the same pattern anywhere — in the nematode orthologue
set, the invertebrate breadth set, or the reference receptors?

B4. **Accession verification (C22).** Is every accession in the current panel
traceable to a real record for the protein it is labelled as? The earlier failure
was three fabricated-from-memory accessions that survived until a positive control
caught them.

B5. **Disposition follow-through.** For each of W1–W5, verify the withdrawal
actually propagated: check `README.md`, the memory files, the figure annotation
tables and the iTOL/ggtree bundles for surviving statements of the withdrawn
claim. W2 in particular was a headline claim.

B6. **P4 charge (protocol rule from `receptor_phylogeny/analysis/reviews/
REFEREE_PROTOCOL_POSTMORTEM_2026-09-10.md`): any claim of the form "X is
absent / root-state / first-in-lineage" gets a charge item.** In this document
that covers C7 (absent from four *Caenorhabditis* species), C29 (no non-nematode
member; no vertebrate or insect orthologue anywhere), C30 (DOP-6 in none of the
three OA groups), C37 (C24A8.6 is not a gene), C38 (no pharmacology of any kind),
and all of W5. For each: enumerate the published datasets that could contain a
counterexample and state which were checked **at the data level**, not cited.

B7. **P1 dataset diff.** Produce a table, not prose: for Pu et al. 2023, Koelle
2018 Table 1, and the Harris 2009 panel — what is in their dataset vs ours, both
directions. If a supplement is inaccessible to you, write `BLOCKED: <what to
fetch>` at the **top** of your report, not in a reading list.

## Both referees

- Write to `reviews/2026-09-16_dop6_conclusions/REFEREE_[A|B]_report.md`.
- Mark each finding MAJOR / MINOR and give each a required action.
- A required action must be classified **compute-executable now**, **PI-only**, or
  **deferred-with-owner** (protocol rule P3). No action may close as a bare
  recommendation.
- Where you disagree with a claim, say what evidence would settle it.
- Do not edit `README.md` or any analysis output. Reports only.

---

# ADDENDUM (added after the panel convened) — b120 re-runs are now complete

The staleness noted under "Panel provenance" is resolved. `confirm_tree.py`,
`independent_evidence.py` and `nematode_expansion.py` have been re-run on the
**b120 173-taxon panel** and written to `output/tree/b120/confirm.txt`,
`independent_evidence.txt`, `nematode_expansion.txt`. The v1 copies in
`output/tree/` are untouched, so both are available for comparison.

`rogue_test.py` was **not** re-run on b120 and will not be: its two rogue taxa
(`cel_F59D12.1`, `cel_T21B4.4`) are not in the b120 panel, so the leave-one-out
has no meaning there. **C34 therefore remains a v1-only result** — referees
should treat it as evidence about the v1 panel and about why b120 exists, not as
a b120 claim.

## Revised numbers — these supersede the v1 values in C28, C29 and C32

**Nesting and support (UFboot), b120:**

| clade | gt0.1 | gt0.2 | gt0.3 | gappyout | v1 gt0.1 (superseded) |
|---|---|---|---|---|---|
| 5 tips — dop-6 orthologues | 100 | 100 | 100 | 100 | 100 |
| 9 tips — + dop-5 | **95** | 87 | 94 | 84 | 71 |
| 13 tips — + dop-3 | **90** | 79 | 82 | **60** | 69 |
| 16 tips — + dop-2 | **90** | 79 | 82 | 62 | 71 |
| 19 tips — + hDRD2/3/4 | **95** | — | — | 65 | 88 |
| 24 tips — + 5-HT2 group | **100** | — | — | 99 | — (25 tips, 90) |

Clade composition agreement with gt0.1 is now **8/8 at all three other trim
levels** (v1 gappyout was 7/8).

**Two things referees should weigh, and I am flagging both against myself:**

1. Support rose substantially at the 9-, 13- and 16-tip clades when the three
   sub-boundary genes were dropped (69 → 90 at 13 tips). That is the expected
   consequence of removing long-branch rogues — but it is also exactly the
   pattern that makes "drop the taxon, support improves" a tempting and
   illegitimate procedure. The b120 drop was decided on **HMM score**, before
   any tree was built, not on the support it produced. The referees should
   check that ordering independently rather than take my word for it.
2. **gappyout moves the opposite way**: the 13-tip clade falls from 82 (v1) to
   60 (b120), and 60 is the lowest support the dop-3/dop-5/dop-6 clade takes
   anywhere. Any statement of the form "the nematode-exclusive clade is
   well-supported" must carry that number, not just the gt0.1 value of 90.

**C32 is unchanged in substance.** Recounted in independent genes on b120:

| tips | UFboot | independent genes | characterised | which |
|---|---|---|---|---|
| 5 | 100 | 2 | **0** | — |
| 9 | 95 | 3 | **0** | — |
| 13 | 90 | 4 | **1** | dop-3 (Go, and see C36) |
| 16 | 90 | 6 | 3 | + dop-2, dmDop2R |
| 19 | 95 | 9 | 6 | + hDRD2/3/4 |
| 24 | 100 | 14 | 10 | + ser-1, hHTR2A/B/C |

DOP-6's first characterised neighbour is still DOP-3, four genes out, and DOP-3
is the claim weakened by C36. Higher support has **not** bought more independent
evidence — it has only made the same thin neighbourhood better resolved. That
distinction is the one I most want Referee A to check I have not blurred.

**C29 (nematode exclusivity) on b120**, with the species composition now printed:

| species | genes in the 13-tip clade |
|---|---|
| *C. elegans* | dop-3, dop-5, dop-6 |
| *C. brenneri* | dop-3, dop-5, dop-6 |
| *C. remanei* | dop-3, dop-5, dop-6 |
| *C. briggsae* | dop-5, dop-6 — **no dop-3** |
| *P. pacificus* | dop-3, dop-6 — **no dop-5** |

First non-nematode member is still *Drosophila* Dop2R, entering at 16 tips.
**Note the two gaps**: *C. briggsae* contributes no dop-3 and *P. pacificus* no
dop-5 to this clade. Given W5 (RBH cannot resolve this radiation) and W4 (the
*P. pacificus* dop-5 "absence" was a Compara gap plus a truncated download),
**those two gaps fall squarely under charge B6 and should not be read as
absences until checked at the data level.** Five species is also a thin basis
for a phylum-level exclusivity claim, which is charge A2.
