# dop-6 / C24A8.6 — an inverted duplication, and whether DOP-6 is an octopamine receptor

Sub-project started 2026-09-16. Motivation: strains CHS1223 and CHS1256 (Pu et al. 2023,
*Nat Commun* 14:8410, PMID 38110404, Chen lab) both show an octopamine-related phenotype, and the
only genes they share are **dop-6** and **C24A8.6**.

---

## 1. The locus is a perfect inverted repeat, and C24A8.6 is half of dop-6

| | coordinates (WBcel235, chromosome X) | |
|---|---|---|
| arm A | 4,324,178–4,328,106 | 3,929 bp — dop-6 exons 10–17 |
| spacer | 4,328,107–4,329,302 | 1,196 bp — contains C24A8.5 |
| arm B | 4,329,303–4,333,231 | 3,929 bp — **reverse complement of arm A**; contains all of C24A8.6 |

**Identity 3,929/3,929 = 100.000 %, zero mismatches, sharp junctions.** The spacer is not itself
palindromic (26 % self-complementarity = background).

* **dop-6 = C24A8.1**, X:4,320,863–4,327,331 (+), 692 aa, 17 exons. Exons 1–9 (3,315 bp, codons
  1–405) are single-copy; exons 10–17 lie entirely in arm A.
* **C24A8.6** lies entirely in arm B: **0 bp of its 2,330 bp is unique**.
* **C24A8.6a (287 aa) is DOP-6 residues 406–692, exactly — 287/287 identical.** C24A8.6b adds a
  17 aa cryptic exon spliced from what is dop-6 intron 14. C24A8.6a's start codon *is* DOP-6's
  internal Met406, which makes the gene model look like an ORF-finder artefact on the inverted copy.
* That segment is end-of-ICL3 + TM6 + ECL3 + TM7 + the whole C-tail. The ligand pocket is **not**
  in it: the amine-binding TM3 aspartate **D3.32 = D103**, and separately the `DRL` motif
  (D120-R121-L122) at ~residue 120, are both in the unique half. *(An earlier version of this line
  conflated the two — D3.32 is D103, not D120 — which contradicted §5's correct identification of
  R121 as the DRY arginine. Both residues are in the N-half, so no conclusion changes.)* So C24A8.6 can encode at most a two-TM fragment with no binding pocket.
* The two mRNAs carry the **identical** 864 nt sequence (not reverse complements — verified), so no
  dsRNA/endo-siRNA mechanism is available, and reads from either are the same string.

**Consequences.** C24A8.6 cannot be measured apart from dop-6 by any short-read method, cannot be
mutated without also cutting dop-6, and has no diagnostic base anywhere. dop-6 *is* fully tractable
through exons 1–9. See `output/repeat_structure.txt`, `output/specific_tools.txt`,
`output/dop6_guides.tsv`.

## 2. It is real, and it is a common natural polymorphism

`survey_assemblies.py` over **56 assemblies** (`output/assembly_survey.tsv`):

| | 1 copy | 2 copies |
|---|---|---|
| N2 reference | | 1 |
| N2-derived lab strains | 0 | 7 (N2×3 pipelines, CGC1, UA44, BY250, yA, yA1244) |
| natural isolates | **19** | **29** |

Not an assembly artefact: present in three independent N2 assemblies (Sanger cosmid, U. Oregon
long-read, Wellcome Sanger). Not an assembler collapse in the non-carriers either — the U. Oregon
N2 and CB4856 assemblies are the same BioProject, same pipeline, same release date, and differ.
Replicates agree across pipelines: CB4856 1 copy in 3 assemblies, DL238 and ECA36 1 copy in 2 each,
EG4725 2 copies in 2. Absent from *C. briggsae*, *C. brenneri*, *C. remanei*, *C. japonica*.

**Arm-vs-arm identity is 100.000 % in the reference and 99.87–99.98 % in long-read carriers** — the
residual is most likely long-read basecalling error, not divergence.

### It sits on a ~2 Mb haplotype shared with N2 — not contamination, but not independent either

`haplotype_context.py` + `haplotype_extent.py`. Mean % identity to N2 in single-copy windows:

| window | carriers (n=29) | non-carriers (n=19) |
|---|---|---|
| dop-6 flank, left (X) | **99.977** | **99.477** |
| dop-6 flank, right (X) | **99.974** | **99.626** |
| X, 5.7 Mb away | 99.702 | 99.686 |
| chromosome II | 99.802 | 99.728 |
| chromosome III | 99.951 | 99.847 |

**Not whole-genome N2 contamination**: carriers show ordinary wild divergence from N2 on II and III.
**But carriage is identity-by-descent, not independent**: scanning outward along the X, the
carrier-specific N2-like block runs from about **−1 Mb to +1 Mb**, has decayed by ±3 Mb and is gone
by ±6 Mb. So the 29 carriers are essentially *one* haplotype observed 29 times, and whether the
duplication itself is under selection or merely hitchhiking on a swept chromosome-X haplotype
cannot be separated with these data. *C. elegans* selfs, so blocks persist far longer than in an
outcrosser and the block length must not be converted into a date.

## 3. Everything short-read at this locus is blind

CaeNDR 20250625, 684 isotypes, variant density per kb: unique flank **24.75**, spacer **34.28**,
arm A **3.82**, arm B **0.25**. The arms are a mapping-ambiguity desert (MAPQ 0), not a conserved
region. Therefore:

* Any πN/πS or polymorphism-based constraint estimate for codons 406–692 from short reads is
  meaningless — it will read as spuriously ultra-conserved.
* **Pu et al. validated their alleles by BWA-MEM against WS235 plus Pindel**, which cannot resolve
  these arms; they explicitly have a "disrupted larger CNV" class covering *"duplications,
  inversions, large deletions, or the combination of these events"*. The `c24a8.6` guide necessarily
  had two perfect targets 9 kb apart in inverted orientation. **Get Supplementary Data 1d (guides)
  and 1f (per-allele calls) and check `yum2349`, `yum2350`, `yum2603`, `yum2604` before treating
  dop-6 and c24a8.6 as independent lesions.**
* The fix is the one the assemblies already implement: read each isolate in its own coordinates.
  `polymorphism_from_assemblies.py` rebuilds the dop-6 CDS exon-by-exon from 56 assemblies and
  recovers what mapping cannot (`output/dop6_polymorphism.tsv`):

  | | segregating sites, N-half exons | segregating sites, C-half exons |
  |---|---|---|
  | CaeNDR short reads, 684 isotypes | 10 | **0** |
  | assemblies, 56 isolates | 3 | **7** |

  Short-read mapping reports **zero** variation across the whole duplicated half; the assemblies
  find seven segregating sites there from a sample an order of magnitude smaller. Biologically the
  result is also striking: **piN = 0 in both halves — all 10 segregating sites are synonymous**,
  across 56 natural isolates. Carriers are nearly invariant (1 site) versus non-carriers (8),
  which is the same shared-haplotype signal as above, not extra constraint.

## 3b. The strains: one guide, two cut sites, and a 75 kb duplication

Resolved 2026-09-16 from Pu et al. Supplementary Data 1 (`41467_2023_44177_MOESM4_ESM.xlsx`,
sheets d, e, f), supplied by MOD. This is the most consequential result in the sub-project.

**dop-6 and C24A8.6 were never targeted independently.** Sheet d gives *the same* guide, PAM and
ssODN for both:

```
guide  CAAGTTGGTTGCAACTGTAA   PAM TGG   ssODN ...CCAAGTTGGTT TAA AAGCTT TAA actgaaact...
```

Only the genotyping primers differ, because those sit in the unique flanks either side. That guide
has **exactly two perfect full-length matches in WBcel235** — X:4,325,604-4,325,623 (arm A, in
dop-6) and X:4,331,786-4,331,805 (arm B, in C24A8.6, minus strand) — as it must, since the arms are
identical. So one cutting reagent made two double-strand breaks 6.2 kb apart in inverted
orientation, and `yum2349`/`yum2350` (and `yum2603`/`yum2604`) are the two ends of one event, not
two alleles.

**CHS1223 is a complex rearrangement, not a double null.** Sheet f, both alleles, class
*"disrupted larger CNV"*:

| allele | gene | detected variation | location |
|---|---|---|---|
| yum2349 | dop-6 | **Inversion with insertion** (50 bp) | X:4,325,629-4,330,411 |
| yum2350 | C24A8.6 | — | — |

with the shared comment: ***"4kb deletion into dop-6 and C24A8.5, plus 75kb TD spanning multiple
genes"***. So the strain carries a 4 kb deletion taking out part of dop-6 *and* C24A8.5, an
inversion with a 50 bp insertion, and a **75 kb tandem duplication over an unspecified set of
genes**. The strain's other three alleles (F38E1.6, F09C6.16, K01A6.6, on V/V/IV) are clean.

**CHS1256's dop-6 and C24A8.6 alleles were never actually detected.** Sheet f records both as
*"not detected by Pindel"*, *"manually curated via IGV"*, with no mutation location and no detected
variation — in precisely the interval where short reads map at MAPQ 0 (§3). Two of that strain's
other alleles are also CNV-class (`ser-6`: "large deletion, part of ser-6 (2xINV)"; `dop-4`:
"2xINV, deletion of dop-4 exon 1"). It also carries `dop-1(yum2598)` at X:4,415,971, ~90 kb from
the dop-6 cut sites — three breaks inside 90 kb of chromosome X in one injection.

For scale, *"disrupted larger CNV"* is 82 of 1,725 alleles (4.8 %) across the whole resource;
CHS1223 has one and CHS1256 has three.

**Consequence for the project.** The shared-gene argument — "the only genes CHS1223 and CHS1256
have in common are dop-6 and C24A8.6" — does not survive this. The two are one target, so they
cannot be separated genetically in these strains; CHS1223 additionally carries a 75 kb duplication
over multiple genes and a deletion of C24A8.5; and CHS1256's dop-6/C24A8.6 lesions rest on manual
inspection of unmappable reads. **An octopamine phenotype in these strains cannot presently be
attributed to dop-6.** Pu et al. document all of this honestly in their supplement; the CGC genotype
strings do not.

This makes the clean experiment in §"Next" not merely preferable but necessary: a fresh dop-6 null
cut in the unique 5' half (a single-target guide, `output/dop6_guides.tsv`), outcrossed, plus cDNA
rescue.

## 4. DOP-6 is not a dopamine receptor — it is an orphan

* No pharmacology of any kind. Named by homology in Chase & Koelle 2004. **Koelle 2018 WormBook
  (NBK535514) Table 1 lists it as `C24A8.1` under "Poorly-studied putative small-molecule GPCRs",
  with every ligand, coupling and assay column empty** — and lists `C24A8.6` separately, propagating
  the reference artefact. The full orphan set is T02E9.3 (= dop-5), C24A8.1 (= dop-6), C24A8.6,
  F35H10.10, F59D12.1, T21B4.4.
* dN/dS vs *C. briggsae* (`output/dnds_vs_briggsae.tsv`): dop-6 **dN = 0.0252**. dS is saturated
  throughout, so dN is the usable column. **Corrected 2026-09-17 after referee audit — an earlier
  version of this section quoted 9 of the 15 receptors in that file and called dop-6 the lowest.
  It is 4th of 15:**

  | rank | gene | dN | | rank | gene | dN |
  |---|---|---|---|---|---|---|
  | 1 | tyra-2 | 0.0116 | | 9 | ser-6 | 0.0411 |
  | 2 | **octr-1** | **0.0236** | | 10 | tyra-3 | 0.0457 |
  | 3 | ser-5 | 0.0246 | | 11 | ser-3 | 0.0500 |
  | 4 | **dop-6** | **0.0252** | | 12 | dop-2 | 0.0508 |
  | 5 | dop-5 | 0.0261 | | 13 | ser-4 | 0.0609 |
  | 6 | ser-2 | 0.0290 | | 14 | dop-3 | 0.0618 |
  | 7 | gar-3 | 0.0319 | | 15 | dop-4 | 0.0733 |
  | 8 | dop-1 | 0.0322 | | | | |

  All three receptors that beat dop-6 were inside the omitted six, and no rule selects the nine that
  were shown. **This matters for the argument, not just the ranking: OCTR-1 — the best-evidenced
  octopamine receptor in the genome — is more constrained than dop-6.** So constraint cannot be used
  to distinguish dop-6 from an OA receptor, and to the extent it points anywhere it points the other
  way. There are no confidence intervals in that table; dop-6's would overlap ranks 1–8, so even
  "4th" should be read as a band.
* The palindrome does **not** bias the estimate — but on **one** argument, not three
  (`output/mirror_bias_check.txt`). What carries it: **the duplication postdates the *C. briggsae*
  split by a wide margin**, so the *C. briggsae* comparison cannot see it. The other two legs were
  withdrawn on referee audit: the Fisher test (p = 0.258) has **power 0.20** at the observed effect,
  whose point estimate is a **1.49× deflation of the C-half — the direction the bias hypothesis
  predicts** (Woolf 95 % CI on the rate ratio 0.77–2.95; 80 % power would need ~3×); and the
  four-isolate recomputation is near-tautological, since those isolates differ from N2 by 0–2 nt in
  2,079, which forces the NG86 counts to match arithmetically.

## 5. Could the duplicated C-terminus bind a G protein on its own?

Interface defined from coordinates (≤4.5 Å) in **11 receptor–G protein complexes spanning Gs, Gi,
Go, Gq and G11**, with receptor and Gα chains identified by sequence, not by chain label
(`output/gprotein_interface_multi.tsv`). Consensus contacts (≥6 of 11) = 20 positions: **13 lost,
7 retained → the fragment keeps 7 of 20 consensus positions**.

**Corrected 2026-09-17 after referee audit.** The class breakdown was previously given as "29–41 %,
no class rescues it". Two things were wrong with that:

* **6WHA was misclassified as Gs.** Classing a Gα chain by whole-chain identity silently
  misclassifies **mini-G chimeras**: a mini-Gs/q is a Gs scaffold carrying Gq's α5 helix, so it
  scores as Gs on whole-chain identity while the part that actually inserts into the receptor — the
  part this analysis measures — is Gq. 6WHA scored 88.6 % to the Gs reference **while its own PDB
  title reads "…IN COMPLEX WITH A MINI-GALPHA-Q PROTEIN…"**. The diagnostic was sitting in the
  table: 88.6 % was the lowest identity in the panel and every genuine chain scored ≥ 92 %. A
  chimera guard is now in `gprotein_interface_multi.py`.
* **The sampling is thin and was not stated.** n = 4 Gi, 3 Gs, 2 Gq, **1 each for Go and G11** — so
  "per class" means a single structure for two of the five classes.

Corrected figures, per structure and per class:

| class | n | mean % kept | range |
|---|---|---|---|
| Gs | 3 | 24.9 | 15.4–37.0 |
| Gi1 | 4 | 33.2 | 21.4–45.8 |
| Gq | 2 | 39.6 | 37.5–41.7 |
| Go | 1 | 40.9 | — |
| G11 | 1 | 40.7 | — |

Per-structure range across all 11 is **15.4–45.8 %**, not 29–41 %. The conclusion is unchanged and
the direction is if anything clearer — no class approaches an intact interface — but note the
retained fraction is **highest for Gq/Go/G11 and lowest for Gs**, driven by β2AR (15.4 %), the most
compact receptor in the set. **The "~35 %" figure is a count of positions with no geometry behind
it**; what actually carries the conclusion is the next paragraph, which is a structural argument.

What is lost is decisive: **R121, the DRY arginine, is a Gα contact in 11 of 11 structures**, plus
the whole ICL2 patch (A124, A125, P128, A129, Y131, K132, T133) and TM5/ICL3 (I206, L210, R213).
Retained: R467, R470, A471, L474, V475 (TM6 cytoplasmic end) and I534, N535. The Gα α5 helix docks
into a cavity opened by TM6 swinging away from TM3; with no TM3 there is no cavity. The live
hypothesis is therefore **not** autonomous signalling but complementation with / competition against
full-length DOP-6, or GRK/arrestin sequestration via the serine-rich tail.

## 6. Tree — where DOP-6 sits, and how strongly an octopamine placement is disfavoured

**Rewritten 2026-09-17. The previous version of this section described a superseded 50-taxon panel
(Q.PFAM+F+R6, 1,034 columns) in the present tense, cited two files that no longer exist, and
restated a claim withdrawn 200 lines below it. It has been deleted rather than patched.** The
current analysis is in §6b (how the alignment is built) and §6c (what the tree supports). The panel
is **b120**: 173 taxa, boundary at monoamine HMM ≥ 120, trees in `output/tree/b120/trees/`.

### Where DOP-6 falls

Nested clades containing `cel_dop-6`, UFboot at gt0.1, with the gappyout value in brackets as the
low end of the range (`output/tree/b120/confirm.txt`):

```
100 [100]   dop-6 orthologues (Cel, Cbr, Cbn, Cre, Ppa)
 95  [84]   + dop-5
 90  [60]   + dop-3
 90  [62]   + dop-2
 95  [65]   + dmDop2R, hDRD2, hDRD3, hDRD4
100  [99]   + the 5-HT2 group
```

Composition is identical at all four trim levels (8/8 clades), and six independent unconstrained ML
re-searches with different seeds all recover the same 5→9→13→16→19→24 nesting. **DOP-6 sits inside
the D2-like radiation.** Support is trim-dependent and the 13-tip clade falls to **60** at gappyout,
which must be quoted alongside the 90.

*Caveat recorded for reuse: per-run logL across those re-searches spanned −77505.7 to −77519.0
purely from independent free-rate parameter optimisation at 0.96 sites per parameter. **Raw logL
from separate IQ-TREE runs in this project are not comparable** and the spread is easy to mistake
for topological disagreement.*

### How strongly is an octopamine placement disfavoured?

Tested rather than asserted, with `topology_test.py` — which had been written but **never run**.
Four hard constraints force the five dop-6 orthologues monophyletic with each characterised OA
group in turn, excluding dop-5/dop-3/dop-2/dmDop2R/hDRD2-4:

| forced topology | ΔlogL | c-ELW | p-AU (gt0.1) | p-AU (gappyout) |
|---|---|---|---|---|
| + SER-3/SER-6/dmOAMB | 26.9 | 0.069 | 0.132 | 0.157 |
| + OCTR-1/α2-adrenergic | 25.7 | 0.021 | 0.063 | 0.103 |
| + dmOctβ1/2/3R | 27.9 | 0.028 | 0.096 | 0.133 |
| + SER-2/dmOctTyrR | 26.0 | 0.006 | 0.060 | 0.091 |

**No OA placement is rejected at α = 0.05 on either trim level**, but each costs 23–28 logL units
and all four together hold only ~12–15 % of the likelihood weight. The honest statement is
therefore: **an octopamine placement for DOP-6 is consistently and substantially disfavoured, but
not excluded.**

### Why the earlier hedge was wrong

An earlier version of §6c argued that the tree was *nearly uninformative against* an OA assignment,
because octopamine is invertebrate-specific and a phylum-specific ligand could sit in a
phylum-specific lineage that would not be expected to group with insect OA receptors.
**That premise is refuted by this panel's own data.** Every characterised OA receptor here sits in a
well-supported **cross-phylum** clade:

| OA group | tips | nematode / other | SH-aLRT/UFboot |
|---|---|---|---|
| SER-3/SER-6 + dmOAMB + *Limulus*, *Mizuhopecten*, *Lingula*, *Saccoglossus* | 11 | 6 / 5 | 96.4/98 |
| OCTR-1 + hADRA2A/B/C + *Danio* + dmOctα2R + amphioxus, oyster, sea cucumber | 14 | 4 / 10 | 98.1/100 |
| SER-2 + dmOctTyrR + bmOAR1 + oyster, priapulid | 9 | 4 / 5 | 92.9/99 |
| dmOctβ1/2/3R | 3 | 0 / 3 | 98/100 |

Nematode OA receptors demonstrably **do** group with their insect, mollusc and vertebrate relatives
— three times, at 92–100 support. These are four internal positive controls showing the test works,
in the same panel where the hedge claimed it could not. The hedge was also ligand-agnostic: taken
seriously it would equally void "D2-like by descent".

**The real reason the tree cannot settle the ligand is different, and simpler: ligand is not a
topological property.** The project's own Fitch count is 15 transmitter changes across 56 labelled
tips, so a ligand switch on the dop-6 stem is a common-rate event that no topology test can see.

### Can the tree predict coupling at all?

The precedent for reading a receptor's G protein off its nearest characterised relative is SER-5's
Gs call (Harris et al. 2009, *J Neurosci* 29:1446, PMID 19193891). Its validation is weaker than the
citation implies: the tree was **5-HT receptors only**, rooted on GAR-1; the method was ClustalW in
MegAlign/DNAStar, **hand-tuned**, then neighbour-joining in PAUP with DNAStar bootstrapping — no ML,
no model selection; **hypervariable regions were deleted, including ICL3**, which is where part of
the Gα interface lives (§5); and the validation is one sentence of eyeball assessment with no
statistical test. *(The full text is Cloudflare-gated here — this rests on the abstract, the
figure legend and the Methods as far as they are accessible. Flagged as unaudited.)*

So the intuition was tested rather than inherited. Fitch parsimony vs 10,000 label permutations on
b120, all four trim levels (`output/tree/b120/confirm.txt`):

| trait | tips with data | observed changes | random mean | p |
|---|---|---|---|---|
| ligand | 56 | 15 | 35.2 | 0.0001 |
| Gα coupling | 52 | 8 | 23.7 | 0.0001 |

Both are clustered at every trim level, so tree position does carry information about both. **But
clustering is not monophyly** — see §6c, "Coupling is clustered, not blocked".

**Tips are labelled by pharmacology, not gene name**: SER-2 is a *tyramine* receptor, SER-3 and
SER-6 are *octopamine* receptors, dop-5 and dop-6 are orphans.

## Scripts

| script | does |
|---|---|
| `fetch_locus.py` | locus sequence + gene models from Ensembl |
| `analyze_duplication.py` | repeat boundaries, exon mapping, protein relationship |
| `crossspecies_palindrome.py` | other *Caenorhabditis* loci |
| `ortholog_divergence.py` | divergence per half, protein and nucleotide |
| `selection_dnds.py` | NG86 dN/dS vs *C. briggsae* for 17 receptors |
| `mirror_bias_check.py` | does the palindrome bias the dN estimate |
| `wild_isolates.py` | CaeNDR variant density across the locus |
| `assembly_copy_number.py`, `survey_assemblies.py` | copy number across 56 assemblies |
| `wild_locus_structure.py` | per-assembly locus structure, gene-anchored |
| `haplotype_context.py`, `haplotype_extent.py` | polymorphism vs contamination vs sweep; haplotype block length |
| `polymorphism_from_assemblies.py` | diversity without mapping to a duplicated reference |
| `design_specific_tools.py` | dop-6-specific guides, primers, junction PCRs |
| `gprotein_interface.py`, `gprotein_interface_multi.py` | Gα interface mapping |
| `fetch_tree_seqs.py`, `curate_tree_panel.py`, `build_tree.py` | panel, curation, alignment, trees |
| `trait_signal.py` | can the tree predict ligand / coupling (Fitch + permutation) |
| `render_tree.py`, `plot_tree.py` | text tree; publication figure |

## Next

1. **BLOCKED, needs a human:** Pu et al. Supplementary Data 1d and 1f (guides and per-allele CNV
   calls) — needed to know what `yum2349/2350/2603/2604` actually are. Also the Harris et al. 2009
   full text (jneurosci.org is Cloudflare-blocked here; PMC has abstract only).
3. **Deorphanising DOP-6 in a heterologous system may not work.** The Komuniecki lab reported in
   the same paper that direct characterisation of SER-5 was unsuccessful, that similar attempts on
   other *C. elegans* GPCRs had failed (Rex & Komuniecki 2002; Hobson et al. 2003), and that this
   "may be problematic in the absence of a nematode cell line" because heterologous expression of
   many *C. elegans* GPCRs is poor and may require accessory proteins. Budget for that, or plan an
   in-vivo readout (e.g. epistasis with *goa-1* / *egl-30* / *gsa-1*) instead.
4. Decisive experiment, independent of all of the above: **a dop-6 null made in the unique 5′ half
   (leaving both arms untouched) plus cDNA rescue.** That assigns the phenotype without ever having
   to separate the two copies. A C24A8.6-specific null is not cleanly possible — arm B's only
   unique flanks are the spacer and a point inside cst-2's last exon.

## 6b. How the tree alignment is built, and what is deliberately thrown away

**This is the single most consequential methodological choice in the sub-project and it is stated
here in full, because it determines what the tree can and cannot answer.**

### The problem

Aminergic GPCRs share a ~250-residue seven-transmembrane core but differ enormously outside it. In
this 164-taxon panel the third intracellular loop (ICL3) ranges from **14 to 821 residues**; DOP-6's
is ~240. Aligning full-length sequences forces DOP-6's 240-residue ICL3 into the same columns as
beta2AR's 34, which has no correct answer. The result was an alignment that was **84.4 % gaps with
all 164 sequences above 50 %** — unusable.

### What is actually done

`extract_domains.py`. Each sequence is scanned with the monoamine HMM. A 7TM receptor does **not**
align to the model as one block: the alignment breaks at the insertions, and hmmsearch reports the
receptor in pieces. For DOP-6:

```
block 1: residues  27-164   = model positions   2-122
block 2: residues 181-224   = model positions 119-162
block 3: residues 466-546   = model positions 152-224
                  ↑ 241 residues skipped between blocks 2 and 3 — this is ICL3
```

The model coordinates run continuously (2→224); only the *receptor* coordinates jump. Three options
exist and only one works:

| | DOP-6 | beta2AR | dmOctB3R | |
|---|---|---|---|---|
| best block only | 138 aa | 140 aa | 142 aa | loses TM6, TM7 |
| span, first→last | 520 aa | 303 aa | 1093 aa | ICL3 back inside |
| **blocks concatenated** | **263 aa** | **254 aa** | **258 aa** | ← used |

Sequences of 413, 692 and 1256 aa converge on 254–263 aa. A **15-residue margin** is added each side
of every block (an HMM envelope is a detection threshold, not a domain boundary) and overlapping
blocks are merged. This is the same operation Harris et al. 2009 performed by hand; the HMM finds
the boundaries automatically.

### Effect

| | full-length | 7TM excerpts |
|---|---|---|
| taxa × columns (gt0.1) | 49 × 1034 | 163 × 362 |
| total gaps | 48.5 % | **11.4 %** |
| parsimony-informative | 783 / 1034 (76 %) | **356 / 362 (98 %)** |
| composition chi2 failures | 15 | 8 |

### What is thrown away, and three checks that it is safe

**1. The Galpha interface survives intact** (`output/gprotein_survival.txt`). ICL3 contributes to G
protein binding, so excising it could have deleted the residues the coupling argument rests on. It
does not: **all 20 consensus Galpha-contact positions** (contacting in ≥6 of 11 receptor–G protein
complexes) fall inside retained blocks in DOP-6 — including the TM5/ICL3 positions I206, L210, R213
and the TM6 positions R467, R470, A471, L474, V475. The reason is structural: **ICL3 contacts Galpha
at its ends**, where it leaves TM5 and enters TM6, and those ends are the conserved helix boundaries
the HMM recognises. The excised region is the middle of the loop (DOP-6 240–450), which makes no
consensus contact. For compact receptors the excision is negligible — beta2AR loses 4 aa.

**2. The ligand pocket survives.** D3.32, located independently from three ligand-bound structures,
is present in the excerpt of every receptor known to carry it.

**3. ICL3's own signal was measured, not assumed** (`output/icl3_information.txt`). Two results,
pointing opposite ways:

* ICL3 **length is phylogenetically conserved** — mean |length difference| between sister tips is
  52.7 aa versus 110.8 aa when lengths are shuffled across tips, p = 0.0001. So the loop is not
  noise; it carries real lineage information.
* ICL3 **sequence stops being comparable with depth** — identity runs 76.8 % between orthologues
  (7TM core: 86.5 %), 28.4 % between *C. elegans* paralogues (core 31.9 %), 26.4 % *C. elegans* vs
  human (core 37.1 %). *Caveat: these are global pairwise identities counted only over positions
  where both sequences have residues, which inflates the figure when the two loops differ greatly
  in length — the alignable fraction, not the identity, is the real limit.*

**Honest conclusion:** ICL3 does carry lineage signal, and excising it does lose something — but the
signal is usable only within orthologue groups, far shallower than the nodes this tree exists to
resolve, and retaining it destroys the deep alignment outright. The right use of ICL3 is as a
separate shallow-depth marker, not as tree columns.

### The limit this leaves

356 parsimony-informative sites for 163 taxa is **~1.1 sites per free parameter** (an unrooted
163-taxon tree has 323 branches). Domain-only alignments are routinely under-determined, so this
ratio is reported alongside any formal topology test rather than left implicit, and deep nodes in
this tree should be expected to carry modest support.

## 6c. What the tree actually supports about DOP-6 — and what it doesn't

Written after several claims made from this tree had to be withdrawn. The corrections are kept
visible rather than edited away, because each marks a way the tree can mislead.

### Clade size is not evidence size

About 40 of the 163 tips are Compara orthologues of *C. elegans* genes. They resolve topology but
inherit ligand and coupling by descent, so counting them as evidence is pseudo-replication. Recounted
in independent genes (`output/tree/independent_evidence.txt`):

| tips | UFboot | independent genes | characterised | which |
|---|---|---|---|---|
| 5 | 100 | 2 | **0** | — |
| 9 | 71 | 3 | **0** | — |
| 13 | 69 | 4 | **1** | dop-3 |
| 16 | 71 | 6 | 3 | + dop-2, dmDop2R |
| 19 | 88 | 9 | 6 | + hDRD2/3/4 |

DOP-6's nearest neighbours are all orphans. The first characterised gene is **DOP-3, four genes out**
— and on re-reading the primary data, its dopamine assignment is not open and shut.

Sugiura et al. 2005 (*J Neurochem* 94:1146, PMID 16001968) is **not a cAMP measurement**: it is a
CRE-luciferase transcriptional reporter in transiently transfected CHO cells, 3 h with 1 µM
forskolin. In DOP-3-expressing cells **5-HT beat dopamine on both axes** — pEC50 8.90 vs 7.57, and
~48 % vs ~40 % maximal inhibition. The authors discount 5-HT because untransfected CHO cells also
respond (pEC50 8.96, n.s.) through endogenous 5-HT1B, but that compares **potency only**: in their
own inset control cells bottom out near 90 % (~10 % reduction) against ~52 % (~48 %) in DOP-3 cells,
a ~5-fold efficacy difference in the direction of the transfected receptor. Dopamine's specificity
rests on the D87N (TM3 Asp3.32) and S167A mutants shifting its pEC50 — **and those mutants were
never tested against 5-HT**. So "DA > TA > OA" should not be cited as a selectivity order.
DOP-3's **Gα o coupling is separately weak**: IMP, inferred from a mutant phenotype (Chase, Pepper &
Koelle 2004, *Nat Neurosci* 7:1096, PMID 15378064), not a coupling assay.

Tree tips keep their assigned ligand colour, with an asterisk marking assignments weaker than their
evidence code implies (DOP-3, DOP-4 — same reporter format). Full breakdown in
`WormPooledCloning/data/processed/gpcr/cel_monoamine_evidence.md`. A Gi/o prediction for
DOP-6 only carries weight at the 19-tip level, where it really says "DOP-6 falls inside the D2-like
radiation" — not "DOP-6's immediate relatives are Gi/o." **An earlier version of this sentence
added "which is Gi/o throughout". That was withdrawn on referee audit as an untested clade-scale
blocking claim of exactly the kind retracted elsewhere in this file, and it depended on stopping at
19 tips: at 19 tips (UFboot 95) all 6 characterised members are Gi/o-family, but at 24 tips
(UFboot 100 — the better-supported clade) it is 6 Gi/o against 5 Gq. The best-supported clade
containing DOP-6 is coupling-mixed.**"

### The clade is Rhabditina-restricted, and the evidence gap is partly structural, partly sampling

**Rewritten 2026-09-17 after referee audit.** The deepest clade containing DOP-6 with no outside
member is the 13-tip dop-3/dop-5/dop-6 group (b120: UFboot **90** at gt0.1, **60** at gappyout); the
first outside member, *Drosophila* Dop2R, appears at 16 tips
(`output/tree/b120/nematode_expansion.txt`).

**The clade is not "nematode-specific" — all five species in it are Rhabditina** (four
*Caenorhabditis* plus *Pristionchus*). `nematode_expansion.txt` prints "Rhabditina" in its own column
and the earlier text then concluded "nematodes". The panel does contain four non-Rhabditina
nematodes — *Trichinella* (clade I), *Brugia* (III), *Strongyloides* (IV), *Necator* (V) — and none
falls inside the clade; but none returned a dop-3/5/6 orthologue either, and that is the RBH
artefact retired in the corrections log below, not evidence of absence.

**The sampling is worse than that.** Referee audit found real orthologues that exist and are simply
not in the panel: *C. briggsae* dop-3 = **Q61H86** (reviewed SwissProt, RBH 1092 bits — and already
used in this project's own `dnds_vs_briggsae.tsv`); *P. pacificus* dop-5 = **A0A2A6CA79**, found by
this project's own script 51 minutes *before* the panel was built and never added; and clean
reciprocal best hits for dop-6 in *Necator* (1110 bits) and *Strongyloides*, and dop-5 in
*Haemonchus*, all sitting unused in `~/data/proteomes/`. **The five-species table is Compara
coverage, not biology.** Whether the radiation is Rhabditina-restricted or nematode-wide is open
until those sequences are placed phylogenetically.

What survives: "D2-like by descent" is defensible; "a dopamine receptor" is not, without
pharmacology — consistent with the non-reciprocal BLAST (DOP-6's best human hit is DRD2, but DRD2's
best worm hit is OCTR-1 and dop-6 ranks 12th of 20). The thin characterised evidence near DOP-6 is
real regardless: its first characterised neighbour is DOP-3, four genes out, and DOP-3's own ligand
call is the one weakened above.

### Coupling is clustered, not blocked

Fitch parsimony vs 10,000 permutations, all trim levels: ligand 15–17 changes vs 35.2 random,
coupling 8 vs 23.7, both p = 0.0001. But clustering is not monophyly — the smallest clade containing
all Gs receptors holds **125 of 163 tips**, Gi/o 160, Gq all 163. Receptor *subtype* groups are
monophyletic (human α1 and β adrenergic perfectly so); *transmitter* classes are not (human 5-HT and
dopamine each scattered). Any Gi/o inference for DOP-6 is therefore a **local** claim about its own
clade, not an appeal to a genome-wide block. There is no such block.

### Long-branch orphans

`cel_F59D12.1` has the longest terminal branch in the tree — 5.98 subs/site, robust z = **17.4**
against a median of 0.348; `cel_T21B4.4` is second at 4.35 (z = 12.3). DOP-6's own branch is 0.020.
Both sit near DOP-6's clade. Beyond branch length, F59D12.1 has three independent strikes: monoamine
HMM score 52.8 (lowest of the *C. elegans* set), **Q** rather than D at 3.32, and only 318 aa. Leave-
one-out rebuilds (drop F59D12.1 / drop T21B4.4 / drop both, same alignment, model fixed) are in
`output/tree/rogue_test.txt`; a long branch alone is not grounds for removal, so the test is whether
their presence changes DOP-6's placement at all.

### Corrections logged

**Second round, 2026-09-17 — from the two-referee panel** (`reviews/2026-09-16_dop6_conclusions/`,
full disposition in `DISPOSITION.md`). Both referees independently found the first item.

* "dop-6 has the lowest dN of the aminergic receptors" — **withdrawn**. The source file holds 15
  receptors; the claim quoted 9. dop-6 is **4th**, and **OCTR-1, the best-evidenced octopamine
  receptor in the genome, is more constrained than it**. See §4. This removes the constraint
  argument entirely, rather than weakening it.
* "The palindrome does not bias dN, verified three ways" — **reduced to one way**. The Fisher test
  has power 0.20 and its point estimate is a 1.49× deflation *in the predicted direction*; the
  four-isolate recomputation is arithmetically forced. Only the timing argument survives, and it is
  sufficient.
* "The 12 test tips test the boundary" — **withdrawn**. `ador-1` (score 85.2, family `other`) was
  removed from the test set by a hardcoded drop-list and replaced with `nmur-3` (78.9), leaving
  12/12 annotated neuropeptide receptors — **the annotation boundary b120 existed to abolish,
  reimposed on the set meant to test it**, removing the one tip that would actually have tested it.
  The clade is real (98.3/100); the test is not. What does justify the cut: **120 sits in the single
  largest gap in the 142-gene ranking (27.0 points, gar-1 124.4 → npr-35 97.4) and the panel is
  invariant to any cut in (97.4, 124.4]**, so the value carries no information.
* "The tree is nearly uninformative against an OA assignment because a phylum-specific ligand could
  sit in a phylum-specific lineage" — **withdrawn**. Refuted by this panel's own data: all four
  characterised OA groups sit in cross-phylum clades at 92–100 support. Replaced by an actual AU
  test: disfavoured by 23–28 logL, ~12–15 % combined weight, **not rejected at α = 0.05**. See §6.
* "The D2-like radiation is Gi/o throughout" — **withdrawn** as an untested blocking claim that
  required stopping at the less-supported 19-tip level. See §6c.
* "The clade is nematode-specific" — **narrowed to Rhabditina**, and the sampling behind it shown to
  be Compara coverage rather than biology. See §6c.
* "I dropped the three sub-boundary genes on HMM score before any tree was built" — **inaccurate**.
  True for any *b120* tree, but `rogue_test.txt`, which shows that dropping exactly two of the three
  lifts the clade 69→91, predates the drop by 8m38s. The cut is still defensible on the gap argument
  above, which does not depend on the timeline.
* Independent-gene counts in "Clade size is not evidence size" were **inflated by one at every
  level** (an underscore/hyphen mismatch counted dop-6 twice). True counts 1, 2, 3, 5, 8, 13. The
  *characterised* counts were unaffected, so the conclusion was understated, not wrong.
* `README.md` §6 was a superseded 50-taxon analysis still in the present tense, and restated the
  withdrawn "nearest relatives are Gq 5-HT receptors" 200 lines above its own retraction. Deleted
  and rewritten. The same claim was also live in the auto-loaded memory file, which has been rewritten.

Outstanding from the panel and **not yet done**: regenerate `accession_audit.tsv` against the live
panel (it audits a deleted one and still lists three fabricated accessions); remove the mislabelled
`lsOA` (a *Lymnaea* receptor tyrosine kinase), `smTAR` (a *Schistosoma* tetraspanin) and the
duplicated `plADRA1`/`plADRA2`; propagate the dating retraction into `output/nematode_dating.*`,
which still writes the literal string `absent` into 22 cells; add `Cbr_dop-3` and `Ppa_dop-5` to the
panel and rebuild; fix the `gene_of()` bug in `independent_evidence.py`. Harris et al. 2009 full
text remains **BLOCKED** and §6's critique of it is unaudited.

* "DOP-6's nearest characterised relatives are Gq-coupled 5-HT receptors" — **withdrawn**. That rested
  on UFboot 37–62 in a 48 %-gap alignment with 2 of 13 human serotonin receptors present.
* "The canonical G-protein blocks hold" — **withdrawn**; asserted from a thumbnail, refuted by the
  monophyly test above.
* "dop-5 is absent from *P. pacificus*" — **withdrawn**. A direct search of the proteome finds a clean
  reciprocal best hit (`A0A2A6CA79`, bits 764), so Compara simply missed it and the dop-5 duplication
  predates the *Caenorhabditis*/*Pristionchus* split. *(The first run of that search used a truncated
  download — 15,999 of 26,122 proteins — so its coverage statement was wrong and any negative from it,
  including an apparent dop-2 non-reciprocation, is not usable. Being re-run complete.)*
