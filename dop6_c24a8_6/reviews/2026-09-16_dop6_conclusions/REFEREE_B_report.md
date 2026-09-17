# Referee B — audit report

dop-6 / C24A8.6 sub-project, claim set `CLAIMS.md` (2026-09-16, incl. the b120 ADDENDUM).
Independent recomputation from primary files. Nothing in `README.md` or `output/` was edited.

---

## BLOCKED

**BLOCKED: Harris et al. 2009, *J Neurosci* 29:1446 (PMID 19193891 / PMC3418693) — full text, Materials and Methods "Sequence analysis" paragraph and the Fig. 1 phylogeny legend.**
Three routes tried and all refused: `jneurosci.org` is Cloudflare-gated; `pmc.ncbi.nlm.nih.gov/articles/PMC3418693/` serves a reCAPTCHA to automated fetch; NCBI `efetch db=pmc id=3418693` returns `<!--The publisher of this article does not allow downloading of the full text in XML form.-->`.
*What to fetch:* the PDF of J Neurosci 29(5):1446–1456, saved to `dop6_c24a8_6/data/harris2009_jneurosci.pdf`. Needed to complete the B7 diff row for the Harris panel and to audit C35's five specific methodological assertions (ClustalW/MegAlign hand-tuned; NJ in PAUP 4.0b10; rooted on GAR-1; hypervariable regions incl. ICL3 deleted; one-sentence eyeball validation). **C35 is currently unauditable and I neither confirm nor refute it.**

**NOT blocked — two items the project lists as blocked are in fact retrievable, and I retrieved both:**
- Pu et al. 2023 Supplementary Data 1 (`41467_2023_44177_MOESM4_ESM.xlsx`, 812,657 bytes) downloads without authentication from `static-content.springer.com`. README §"Next" item 1 still says *"BLOCKED, needs a human"*. I parsed sheets d, e and f in full; C12–C14 are audited below.
- Koelle 2018 WormBook Table 1 (NBK535514) parses out of the Bookshelf HTML. C37's orphan list is audited below.

---

## Headline finding

**The tree work reproduces essentially perfectly; the prose around it does not.**

Every tree number I was asked to check reproduced to the digit — the 5-tip clade at UFboot 100 on all eight trees, the b120 nesting 95/90/90/95 at gt0.1 and 84/60/62/65 at gappyout, identical clade composition at all four trim levels on both panels, the 12-tip test clade at 98.3/100, 356 informative sites over 323 branches. The ADDENDUM's revised support table is correct in every cell I could check, and where it differs from my recomputation it is because it left cells blank, not because it filled them wrongly.

Against that, **three claims fail on the data and one exculpatory ordering claim fails on the timestamps**:

1. **C16 is false as worded.** `dop-6` is *not* "the lowest of the aminergic receptors tested". `selection_dnds.py` tested 17 receptors and returned 15 dN values; the claim lists 9 of them. Three receptors are more constrained than dop-6 — **tyra-2 0.0116, octr-1 0.0236, ser-5 0.0246** vs dop-6 0.0252 — and all three are among the six the claim omits. dop-6 is 4th of 15. One of the three is **OCTR-1, the octopamine receptor with the strongest evidence in the genome**, which makes the omission consequential rather than cosmetic. The error has propagated to `README.md` §4 and into the persistent memory file.
2. **C26's boundary test is not the test it says it is.** The two files C26 cites list *different* 12-gene test sets. `boundary_testset.txt` includes `ador-1` (HMM 85.2, WormBase family `other`) and excludes `nmur-3` (78.9); the tree actually built includes `nmur-3` and excludes `ador-1`. `ador-1` was removed by a hardcoded `DROP_WB` set. The surviving test set is 12/12 *annotated neuropeptide receptors* — so the annotation-label boundary that b120 was created to abolish was reimposed on the test set, and the one excluded gene that would have tested a non-neuropeptide annotation was the one taken out.
3. **W2 and W5 did not stay withdrawn.** W2 is stated as live fact in `README.md:249`, in two script docstrings, and — worst — in the persistent memory file `dop6_is_an_orphan_receptor.md`, which is auto-loaded into every future session and restates the entire withdrawn 50-taxon result including the deleted squid-rhodopsin panel's numbers. W5 has not propagated to `output/` at all: `nematode_dating.tsv` still writes the literal string `absent` into 22 cells.
4. **The ADDENDUM's ordering claim is inaccurate.** "The b120 drop was decided on HMM score, **before any tree was built**" is true only of *b120* trees. `output/tree/rogue_test.txt` — which demonstrates that dropping exactly `cel_F59D12.1` and `cel_T21B4.4` raises the 13-tip clade from 69 to 91 — was written at **22:21:34**, nine minutes *before* `rebuild_boundary120.py` ran at **22:30:12** and dropped exactly those two genes.

And on B6, the two absences the coordinator flagged are both **panel-sampling gaps, and I settled both against the sequence data**. *C. briggsae* has a reviewed SwissProt dop-3 (`Q61H86`, DOPR3_CAEBR) at 1,092 bits, 89.1 % identity, reciprocal, with a 6.7-fold gap to the next hit — and the project's own `dnds_vs_briggsae.tsv` already uses its Compara ortholog. *P. pacificus* dop-5 (`A0A2A6CA79`) was established by the project itself in W4 and never added to the panel. I also found **three further nematode dop-5/dop-6 orthologues already sitting in `~/data/proteomes/`** that the panel does not contain.

---

# MAJOR findings

### MAJOR-1 — C16: dop-6 is not the lowest-dN aminergic receptor; the claim lists a subset that excludes every lower value

`output/dnds_vs_briggsae.tsv` contains 15 whole-gene dN values (17 receptors requested; `ser-1` and `ser-7` have no *C. briggsae* ortholog). I recomputed NG86 dN and dS from the file's own `S`, `N`, `Sd`, `Nd` columns with Jukes–Cantor: every one of the 15 reproduces to ±0.0001. The ranking:

| rank | gene | dN | listed in C16? |
|---|---|---|---|
| 1 | **tyra-2** | **0.0116** | **no** |
| 2 | **octr-1** | **0.0236** | **no** |
| 3 | **ser-5** | **0.0246** | **no** |
| 4 | dop-6 | 0.0252 | yes ("the lowest") |
| 5 | dop-5 | 0.0261 | yes |
| 6 | ser-2 | 0.0290 | yes |
| 7 | **gar-3** | **0.0319** | **no** |
| 8 | dop-1 | 0.0322 | yes |
| 9 | **ser-6** | **0.0411** | **no** |
| 10 | **tyra-3** | **0.0457** | **no** |
| 11–15 | ser-3, dop-2, ser-4, dop-3, dop-4 | 0.0500–0.0733 | yes |

The nine listed values are individually correct and correctly ordered among themselves. The set is the problem: six of fifteen are omitted and three of those six are lower than dop-6. There is no rule that selects the nine — it is not "dopamine receptors" (ser-2/3/4 are in), not "genes with a computable dS" (dop-1 and dop-4 have `nan` dS and are in, ser-5 has dS and is out). `README.md` §4 repeats it verbatim; `memory/dop6_is_an_orphan_receptor.md` repeats it as *"the lowest of the aminergic receptors tested (dop-4 is 0.0733)"*.

This is also the **second instance of the B3 pattern**: a subset chosen by no stated rule whose effect is to remove exactly the counterexamples.

**Required action (compute-executable now):** regenerate the sentence from the TSV rather than by hand — e.g. "dop-6 dN = 0.0252, 4th lowest of the 15 aminergic receptors tested; tyra-2 (0.0116), octr-1 (0.0236) and ser-5 (0.0246) are lower" — in `README.md` §4, `CLAIMS.md` C16, and `memory/dop6_is_an_orphan_receptor.md`. No recomputation is needed; the file already holds the right numbers.

---

### MAJOR-2 — B4/C22: two mislabelled accessions and one duplicated accession in the current breadth set, and `accession_audit.tsv` no longer describes any panel that exists

I re-fetched all 63 distinct label/accession pairs in `output/panel_build/v2_anchor_panel/provenance.tsv` and `v3_breadth_panel/provenance.tsv` from `rest.uniprot.org` and compared each record's organism, recommended name, alternative names and **gene name** against the label.

**Good news first, because C22 asks the right question and most of the answer is clean.** All 38 human accessions are correct. All 14 *Drosophila* accessions are correct *by gene name* — including the four that look wrong by protein name. `dm5HT1A`=P28285 has recommended name "5-HT receptor **2A**" but gene name `5-HT1A`; `dm5HT1B`=P28286 has recommended name "5-HT receptor **2B**", alt name "5-HT receptor 1B", gene name `5-HT1B`. These are legacy UniProt names from Saudou 1992, not errors. `dm5HT7`=P20905 (gene `5-HT7`) and `dm5HT2A`=Q9VN38 (gene `5-HT2A`) are likewise correct. **The three fabricated accessions of C22 are gone from the live provenance** — they were replaced by logged UniProt queries.

**Three problems remain, all in the `v3_breadth_panel` invertebrate set:**

| label | accession | what the record actually is | verdict |
|---|---|---|---|
| `lsOA` | **O76997** | *Lymnaea stagnalis* "**Putative neurotrophin receptor LTRK 1**" | a receptor **tyrosine kinase**, not a GPCR, labelled as an octopamine receptor |
| `smTAR` | **C4QM85** | *Schistosoma mansoni* "**Tetraspanning orphan receptor**" (gene `TOR`) | a tetra/trispanning orphan, not a tyramine receptor |
| `plADRA1` / `plADRA2` | **A0A1J0F5X6** (both) | *P. dumerilii* "Alpha-**2** adrenergic receptor" | **one accession entered twice under two labels**; `plADRA1`'s label contradicts the record |
| `tcOA1` | A0AB74UEC7 | "Octopamine receptor **A2**", gene `OctA2` | subtype-label mismatch (minor) |

`lsOA` and `smTAR` are the same failure class as P31356/P20906 — a non-class-A protein sitting in an aminergic panel under a receptor label — one generation on. **Mitigating fact I want on the record: none of these four reached the final tree.** `basis.tsv` and `basis_b120.tsv` contain `plADRA1` only; `lsOA`, `smTAR`, `tcOA1`, `plADRA2`, `amTAR1`, `apOA`, `ap5HT`, `bmTAR`, `amDop1` are all absent from both panels. The damage is latent, not realised.

**Separately, the audit artefact itself is stale and misleading.** `output/accession_audit.tsv` was produced by `verify_accessions.py` reading `output/tree2/provenance.tsv` — a path that no longer exists (it is now `output/panel_build/v2_anchor_panel/`). The file still lists `dm5HT2A → P31356 (squid rhodopsin)`, `dm5HT7 → P20906 (Pseudomonas decarboxylase)` and a `dmHisCl` entry that is not in any current panel, and it is missing `dm5HT1B` and the entire breadth set. **So the one file whose job is to answer C22 answers it for a panel that was deleted.** Note also that `verify_accessions.py`'s keyword test (`'5HT' → 'hydroxytrypt'`) cannot distinguish 5-HT1A from 5-HT2A and would pass a subtype swap silently; it has no test at all for labels like `lsOA`/`smTAR` whose prefix is not in `EXPECT_ORG`/`EXPECT_KW`.

**Required action (compute-executable now):** (a) point `verify_accessions.py` at both live provenance files, add gene-name matching and a rule that every label prefix must resolve to a known organism (unknown prefix = FAIL, not pass), re-run, and let it overwrite `accession_audit.tsv`; (b) add a duplicate-accession check — `plADRA1`/`plADRA2` would have been caught by one line; (c) drop or correctly relabel `lsOA` and `smTAR` in `v3_breadth_panel/provenance.tsv` before that file is ever used to build a panel.

---

### MAJOR-3 — B3/C26: the boundary test set was filtered by annotation label, and the two cited files disagree about its membership

The score table `WormPooledCloning/.../cel_gpcr_monoamine_scores.tsv` ranks the excluded genes 20th onward:

```
20 npr-35 97.4 nprec   26 ador-1 85.2  OTHER    31 tkr-1  79.0 nprec
21 ckr-1  94.9 nprec   27 npr-22 84.9  nprec    32 nmur-3 78.9 nprec
22 trhr-1 94.8 nprec   28 npr-34 82.6  nprec    33 npr-9  77.1 nprec
23 ckr-2  91.9 nprec   29 npr-24 82.0  nprec
24 nmur-2 87.9 nprec   30 npr-20 81.4  nprec
25 tkr-2  87.7 nprec
```

- `output/tree/boundary_testset.txt` (22:20) lists ranks 20–31 — **includes `ador-1`, excludes `nmur-3`** — and writes `panel_boundary.faa` (176 sequences).
- `output/tree/boundary_rebuild.txt` (22:30) lists ranks 20–25 and 27–32 — **excludes `ador-1`, includes `nmur-3`** — and writes `panel_b120.faa` (173 sequences), which is the panel the trees were built from.

C26 cites **both files** as its source. The tree contains the second set (I confirmed the 12 `TEST_` tips directly in all four b120 treefiles). So **C26's "the 12 highest-scoring excluded genes" is false**: `ador-1` at 85.2 outscores six of the twelve that are in.

The mechanism is explicit in `rebuild_boundary120.py`: `DROP_WB={'WBGene00011878', ...}` removes `ador-1` from the test-set candidate pool, commented *"matched by WBGene so ador-1 cannot duplicate"*. The rationale is that `ador-1` **is** the panel gene `cel_T21B4.4` (I confirmed the identity: ador-1 = WBGene00011878, 362 aa; `cel_T21B4.4` in `basis.tsv` = 362 aa). But `T21B4.4` was *dropped* from the panel in the same script, so there was nothing left to duplicate — the exclusion removes it from the panel *and* from the test of the panel's boundary.

**Why this matters and is not bookkeeping.** `ador-1`/`T21B4.4` is (i) the highest-scoring excluded gene whose WormBase family is `other` rather than `neuropeptide_receptor`, (ii) one of Koelle's six orphans, and (iii) the second-longest branch in the v1 tree. It is therefore the single most informative tip for the question C26 exists to answer — *is there an unannotated non-neuropeptide receptor sitting just below the cut?* With it removed, the test set is homogeneous by annotation, and "all 12 form one clade with 0 aminergic receptors inside it" becomes close to a foregone conclusion: neuropeptide GPCRs are a well-known coherent family. I confirmed the clade is real and strong (exact 12-tip clade, SH-aLRT/UFboot 98.1/100 gappyout, **98.3/100 gt0.1**, 98.1/100 gt0.2, 97.8/100 gt0.3) — the result is right, it just tests less than claimed.

**Required action (compute-executable now):** rebuild with `N_TEST=13` including `TEST_ador-1` (its protein is already in `output/tree/panel.faa` as `cel_T21B4.4`), re-run the four trim levels, and report whether it joins the neuropeptide clade, sits alone, or falls inside the aminergic ingroup. That is a ~2 h IQ-TREE job and it converts C26 from an assertion into a test. Until then C26 must be worded *"the 12 highest-scoring excluded genes that are annotated neuropeptide receptors"*, and the stale `boundary_testset.txt`/`panel_boundary.faa` pair should be deleted or marked superseded so it cannot be cited again.

---

### MAJOR-4 — B5/W2: the withdrawn headline claim is live in four places, including the persistent memory

W2 = *"DOP-6's nearest characterised relatives are Gq-coupled 5-HT receptors."* Grep of the whole sub-project plus the memory directory:

| location | status |
|---|---|
| `README.md:459` | correct withdrawal ✓ |
| **`README.md:249`** | **live, stated as fact**: *"Its nearest characterised relatives are two **5-HT receptors, both Gq**."* — 210 lines *above* its own withdrawal |
| **`plot_tree.py:7`** | **live** in the docstring of a figure-generating script |
| **`topology_test.py:10–11`** | **live**, and used as the *premise* of the topology test |
| **`memory/dop6_is_an_orphan_receptor.md`** | **live**, in full |

The memory entry is the serious one. It is auto-loaded into every future session and reads:

> **TREE RESULT (2026-09-16): DOP-6 does NOT group with octopamine receptors.** 50-taxon aminergic tree, IQ-TREE Q.PFAM+F+R6 … + F59D12.1, T21B4.4, ser-1 and hHTR2A (77.8/62). So four of Koelle's six orphans form a real lineage, and **DOP-6's nearest characterised relatives are two 5-HT receptors, both Gq**.

Every element of that paragraph is superseded: the 50-taxon panel is the one `output/MANIFEST.md` records as **deleted** for containing "a squid rhodopsin, a bacterial decarboxylase and a ubiquitin ligase"; `Q.PFAM+F+R6` is superseded by LG+F+R7 / LG+F+I+R7; the support values `77.8/62` are the 48 %-gap values W2 was withdrawn on; and `F59D12.1`/`T21B4.4` are now below the panel boundary. The same file also carries the MAJOR-1 dN error and a superseded Fitch test (*"34 tips … 7 changes vs 16.2"*; current: 52 tips, 8 vs 23.7).

`README.md` §6 as a whole is stale, not just line 249: it describes "50 sequences", model Q.PFAM+F+R6, embeds `![tree](output/tree/dop6_tree_gt0.1.png)` — **which does not exist** — and cites `output/tree/trait_signal.txt`, **which does not exist**. §6b and §6c beside it carry the current 163-taxon numbers, so the document contradicts itself.

W1 and W3 propagated correctly — no surviving statement of either anywhere. W4's *text* propagated (see MINOR-12 and MAJOR-6 for what did not).

**Required action (PI-only for the memory file; compute-executable for the rest):** the memory file is user-owned and I did not touch it — the PI must rewrite the TREE RESULT paragraph and the dN sentence in `memory/dop6_is_an_orphan_receptor.md`, since a stale memory re-injects a withdrawn claim into every future session and is the highest-leverage single fix in this report. For the repo: rewrite `README.md` §6 against the 163/173-taxon results or delete it in favour of §6b/§6c, fix the two dead file references, and strip the W2 sentence from `plot_tree.py` and `topology_test.py` docstrings.

---

### MAJOR-5 — B5/W5: the withdrawn dating matrix is still asserted in `output/`, in machine-readable form

W5 withdraws the nematode dating matrix as an RBH artefact. The withdrawal reached `CLAIMS.md` and nothing else:

- **`output/nematode_dating.tsv`** writes the literal string **`absent`** into **22 of 30** gene cells. A TSV is the form most likely to be re-read by a script or a collaborator, and it carries no withdrawal marker of any kind.
- **`output/nematode_dating.txt`** ends with *"READING: dop-2 present in Dorylaimia and shared with Drosophila = the ancestral lineage. The clade at which dop-3/dop-5/dop-6 first appear brackets the radiation."* — the inference W5 retracts, stated as the file's conclusion.
- **`nematode_dating2.py`**'s docstring (mtime 23:26, the newest version) still argues the method: *"Every call is a reciprocal best hit in both directions… Trichinella is the decisive taxon."*

I confirmed W5's own stated rationale at the data level: in *Haemonchus*, `cel_dop-5 → A0A7I4YXV5` at **976 bits** (W5 says 969; the difference is BLAST settings) and the reverse hit lands on a dop-5 isoform, while `cel_dop-6`'s best *Haemonchus* hit is `A0A7I4Z6U1`, which is also `cel_dop-4`'s best hit at 692 bits. RBH does pile onto the wrong paralogue here. **W5 is right to be withdrawn; it just has not been withdrawn where it is written down.**

The matrix is also now *stale as well as wrong*: it prints `Ppac: not in the store yet` and `Name: not in the store yet`, but `~/data/proteomes/MANIFEST.tsv` lists both as complete (`Ppac` 26,122, `Name` 19,125) with file mtimes 23:43 and 23:31 — after the dating run at 23:26. Those are exactly the two species where I found clean dop-5/dop-6 orthologues (MAJOR-6).

**Required action (compute-executable now):** prepend a withdrawal banner to `output/nematode_dating.txt`, and either replace every `absent` in `output/nematode_dating.tsv` with `not_resolved_by_RBH` or delete the TSV. Do not re-run `nematode_dating2.py` with `Ppac`/`Name` added — that would produce a better-populated version of a method W5 has already ruled inadequate. The stated replacement (phylogenetic placement) is the right route.

---

### MAJOR-6 — B6: both flagged "absences" are panel-sampling gaps, and at least three more nematode dop-5/dop-6 orthologues are already on disk and not in the panel

Checked at the data level with `blastp` RBH against complete proteomes, both directions, not from the tree.

**(a) *C. briggsae* dop-3 — present, unambiguous.** Downloaded UP000008549 (21,744 proteins) and ran RBH for `cel_dop-1…dop-6`:

| query | best *C. briggsae* hit | bits | %id | len | reverse best hit | reciprocal? |
|---|---|---|---|---|---|---|
| cel_dop-3 | **sp\|Q61H86\|DOPR3_CAEBR** | **1092** | 89.05 | 612 | sp\|Q6RYS9\|DOPR3_CAEEL (1070) | **YES** |
| cel_dop-1 | A8XK01 | 914 | 96.5 | 460 | DOPR1_CAEEL | YES |
| cel_dop-2 | A8XVC3 | 1399 | 89.2 | 854 | DOPR2_CAEEL | YES |
| cel_dop-4 | A8X5P4 | 909 | 88.4 | 517 | DOPR4_CAEEL | YES |
| cel_dop-5 | A8Y4I4 | 1508 | 95.9 | 758 | O45732 (dop-5) | YES |
| cel_dop-6 | A8XJW6 | 1278 | 90.9 | 694 | P91096 (dop-6) | YES |

`Q61H86` is a **reviewed SwissProt entry explicitly named "Dopamine receptor 3"**, and the bit-score gap to the next *C. briggsae* hit is 1092 → 163, i.e. 6.7-fold. This is not a marginal RBH in a radiation; it is as clean as the dop-1 and dop-4 calls. The project already relies on this orthologue: `output/dnds_vs_briggsae.tsv` records `dop-3 → WBGene00032108`, 604 codons, 90.73 % aa identity. **The "no dop-3" cell in the ADDENDUM's species table is a Compara/panel gap; *C. briggsae* is not missing dop-3.** The panel's *C. briggsae* set is 7 genes (dop-4, dop-5, dop_6, octr-1, ser-3, ser-4, ser-6) against *C. brenneri*'s 11 and *C. remanei*'s 11 — the asymmetry is in the sampling, not the genomes.

**(b) *P. pacificus* dop-5 — present; the project established this itself in W4 and never applied it.** I re-ran the RBH against the **complete** 26,122-protein `Ppac_UP000005239.faa` and reproduced `ppa_dop5_search.txt` exactly: `cel_dop-5 → A0A2A6CA79` at **764 bits**, reverse `A0A2A6CA79 → O45732 (dop-5)` at **775 bits**. `A0A2A6CA79` is annotated in UniProt as `Dop-5 OS=Pristionchus pacificus GN=WBGene00091462`. `basis_b120.tsv` contains `Ppa_dop-3`, `Ppa_dop_6`, `Ppa_ser-4`, `Ppa_ser-6` — **no `Ppa_dop-5`**. So the b120 tree reproduces the exact absence W4 withdrew, and the ADDENDUM re-presents it as a feature of the clade.

**(c) New: three further nematode orthologues already in `~/data/proteomes/` and absent from the panel.** RBH, both directions, complete proteomes:

| species | gene | accession | fwd bits / %id | reverse best | in panel? |
|---|---|---|---|---|---|
| *Necator americanus* | **dop-6** | W2T1T3 | **1110 / 81.7 %** (667 aa) | P91096 = dop-6, 1121 | **no** |
| *Necator americanus* | **dop-5** | W2TY13 | 958 / 68.1 % | O45732 = dop-5, 961 | **no** |
| *Strongyloides ratti* | **dop-6** | A0A090LQI8 | 687 / 51.6 % | P91096 = dop-6, 719 | **no** |
| *Haemonchus contortus* | **dop-5** | A0A7I4YXV5 | 976 / 69.1 % | A0A486WU00 = dop-5, 970 | **no** (species absent) |

The panel *does* contain *Necator* (`Nam_ser-3`, `Nam_ser-7`) and *Strongyloides* (`Sra_ser-2`, `Sra_ser-6`) — so these species were sampled, but for other genes. **The 13-tip clade is 5 species because Compara supplied 5 species for these genes, not because the radiation is confined to them.** Every entry in the ADDENDUM's species-composition table is therefore a statement about Compara coverage, and the table should not be read as biology.

For completeness on the other B6 items: the same six-proteome RBH found **no** dop-3 reciprocal partner outside *Caenorhabditis*/*Pristionchus* (in *Haemonchus*, cel_dop-3's best hit reciprocates to `ser-2`; in *Necator*, to dop-6) — so C29's dop-3 restriction survives this test, subject to the W5 caveat that RBH cannot prove absence inside a radiation. And nothing in the tree places a non-nematode inside the 13-tip clade at any trim level on either panel; the first non-nematode is `dmDop2R` at 16 tips, exactly as C29 states. **C29's nematode-exclusivity claim is intact; only its species table is not.**

**Required action (compute-executable now):** add `A0A2A6CA79` (Ppa dop-5), `Q61H86` (Cbr dop-3), `W2T1T3` + `W2TY13` (Nam dop-6, dop-5), `A0A090LQI8` (Sra dop-6) and `A0A7I4YXV5` (Hco dop-5) to the nematode orthologue set and rebuild b120 at all four trim levels. This raises the nematode-exclusive clade from 5 species to 7 and tests C29 far harder than the current panel does; my prediction is that exclusivity strengthens and the species table stops being quotable. Until then, delete the species table from the ADDENDUM or relabel it "Compara coverage in this panel".

---

### MAJOR-7 — the ADDENDUM's ordering claim: no *b120* tree preceded the drop, but the *v1 rogue test* did, by nine minutes

The ADDENDUM asks referees to check, independently, that the b120 taxon drop was decided on HMM score and not on the support it produced. Filesystem mtimes, all 2026-09-16:

```
20:23–21:53  v1 163-taxon trees built (gt0.1, gappyout, gt0.2, gt0.3)
21:53        output/tree/build_report.txt
22:20:16     add_boundary_testset.py → boundary_testset.txt, panel_boundary.faa
22:21:34     output/tree/rogue_test.txt      ← "drop F59D12.1: 9-tip 71→91, 13-tip 69→91"
22:30:12     rebuild_boundary120.py → boundary_rebuild.txt, panel_b120.faa, basis_b120.tsv
22:49:31     b120 gt0.1.treefile   ← EARLIEST b120 tree, 19 min AFTER the panel was fixed
23:13:52     b120 gappyout.treefile
23:32:17     b120 gt0.2.treefile
23:54:16     b120 gt0.3.treefile
```

**The narrow claim holds:** no b120 tree file predates `panel_b120.faa`, and no b120 IQ-TREE checkpoint exists outside `b120/trees/`. **The claim as worded does not:** `rogue_test.txt`, which is precisely a demonstration that removing `cel_F59D12.1` and `cel_T21B4.4` lifts the clade of interest by ~20 points, was written nine minutes earlier. The information "dropping these two improves support" was in hand when the drop was made. I am flagging this because the ADDENDUM invites the check and the sentence "before any tree was built" reads as broader than the evidence supports.

**Three facts that weigh the other way, and I think they are decisive on the substance:**
1. The rule is **external and pre-existing**. The ≥120 cut and the 27.0-point gap are properties of `cel_gpcr_monoamine_scores.tsv` in a *different repository* (`WormPooledCloning`), computed before this sub-project. I verified independently: the largest gap in the top 60 is **27.0**, between `gar-1` (124.4) and `npr-35` (97.4), and the 19 genes above it are exactly 16 `known_monoamine` + `gar-1/2/3`. Nothing about the tree influences that number.
2. **The rule reproduces the drop exactly.** `F59D12.1`/`pcdr-1` = 52.8, `T21B4.4`/`ador-1` = 85.2, `F35H10.10` = `no_hit` on the monoamine HMM. All three are below 120 by the external table, with no discretion applied.
3. **A third gene was dropped that was not a rogue.** `F35H10.10` has an ordinary branch length and does not appear in the rogue test at all — indeed it had already been dropped from the v1 *tree* by the domain-excerpting step (`domain_report.txt`: "1 still have no envelope and are dropped: cel_F35H10.10"), so `basis.tsv` has 164 rows and the tree has 163. The drop set is not "the rogues".

**One thing that weakens the code's own testimony:** `rebuild_boundary120.py` implements the exclusions as a hardcoded name list, `DROP_CEL={'T21B4.4','F59D12.1','F35H10.10'}`, checked *before* any score lookup and `continue`d unconditionally. So the script does not *execute* the ≥120 rule for the three contested genes; it asserts it, and (see MINOR-1) prints `nan` where their scores should be. An auditor has to reconstruct the rule from an external file to confirm it — which I did, and it holds — but the code should demonstrate it.

**Required action (compute-executable now):** replace `DROP_CEL` with an actual score filter (`s = score_by_wbgene[wb]; if s is None or s < CUT: drop`), keyed on WBGene so it cannot miss, and let the script print the real score next to each dropped gene. Re-running it must reproduce `panel_b120.faa` byte-identically; if it does, the ordering question is closed by construction rather than by timestamp. **Also (PI-only):** commit the b120 work — `git status` shows `output/tree/b120/`, `basis_b120.tsv`, `boundary_rebuild.txt`, `bundles/` and the new scripts as **untracked**, so mtimes are currently the *only* evidence of ordering and mtimes are mutable.

---

### MAJOR-8 — C32: every "independent genes" count is inflated by exactly one (a hyphen/underscore bug)

`independent_evidence.py` collapses orthologues with:

```python
def gene_of(t):
    if t.startswith(NEMPRE):
        return 'cel_'+t.split('_',1)[1]        # orthologue -> its C. elegans gene
    return t
```

Panel tips use `cel_dop-6` (hyphen) but `Cbn_dop_6`, `Cbr_dop_6`, `Cre_dop_6`, `Ppa_dop_6` (underscore). So `gene_of('Cbn_dop_6')` returns `'cel_dop_6'`, which never equals `'cel_dop-6'`, and **dop-6 is counted as two independent genes at every level**. The symptom is already visible in `output/tree/b120/nematode_expansion.txt`, which prints *"genes in it: dop-3, dop-5, dop-6, dop_6"* — four names for three genes. `dop-5` and `dop-3` are unaffected because their orthologue tips use hyphens.

Recomputed on `b120/trees/gt0.1.treefile`:

| clade | as reported | correct |
|---|---|---|
| 5 tips | 2 | **1** |
| 9 tips | 3 | **2** |
| 13 tips | 4 | **3** |
| 16 tips | 6 | **5** |
| 19 tips | 9 | **8** |
| 24 tips | 14 | **13** |

The **characterised** counts (0, 0, 1, 3, 6, 10) are unaffected — dop-6 is not characterised — so the ADDENDUM's headline that higher support bought no extra evidence is **correct, and in fact understated**: the neighbourhood is one gene thinner than reported at every level. At 5 tips the honest statement is "one gene, observed in five species", not two.

**Required action (compute-executable now):** normalise separators in `gene_of` (`.replace('_','-')` after the prefix split) and re-run `independent_evidence.py` on both panels; update the C32 tables in `CLAIMS.md` and `README.md` §6c. One-line fix, ~1 s runtime.

---

# MINOR findings

### MINOR-1 — B2/C27: confirmed, cause identified, and the PI's exculpation is correct
`boundary_rebuild.txt` prints `HMM nan` for all three dropped genes although `boundary_testset.txt` gives `F59D12.1 = 52.8`. **Confirmed.** Cause: `score` is keyed on the **public** gene name (`pcdr-1`, `ador-1`) taken from columns 2–3 of `cel_gpcr_monoamine_scores.tsv`, while `DROP_CEL` and the panel keys use **sequence** names (`F59D12.1`, `T21B4.4`), so `score.get('F59D12.1')` misses and falls through to `float('nan')`.

**Did the `nan` change which genes were dropped? No.** `DROP_CEL` is tested first and `continue`s unconditionally, before any score is consulted; the score is used only in the printed line. No gene anywhere in the run was dropped for `no HMM score` (the report shows exactly 3 dropped, all from `DROP_CEL`), and all 19 ingroup genes resolved their scores because their panel keys are public names. C27 is correct as stated.

**One thing C27 does not say and should.** The three `nan`s do not mean the same thing. `F59D12.1` (52.8) and `T21B4.4` (85.2) have real scores that the lookup missed; `F35H10.10` = WBGene00018073 has **`no_hit`** in the score table — it does not hit the monoamine HMM at all, at 1,404 aa, consistent with C37's class-C assignment. Printing one `nan` for both situations hides the fact that the strongest exclusion in the set is the best-justified one.
**Required action (compute-executable now):** covered by the MAJOR-7 fix — key the lookup on WBGene and print `no_hit` distinctly from a lookup miss.

### MINOR-2 — C34: the "robust z" is not a robust z, and has no producing code
Recomputed from `trees/gt0.1.treefile`, 163 terminal branches: median **0.3478** ✓ (claim: 0.348), `cel_F59D12.1` **5.9814** ✓ (5.98), `cel_T21B4.4` **4.3448** ✓ (4.35), `cel_dop-6` **0.0201** ✓ (0.020). Ranks confirmed: longest and second-longest terminal branches in the tree.

The z values do not reproduce. Standard robust z = (x − median)/(1.4826·MAD) gives **11.7** and **8.3**, not 17.4 and 12.3. The only formula that returns the claimed values is (x − median)/MAD with the raw MAD (1.4826·MAD = 0.4805; raw MAD = 0.3241; (5.9814−0.3478)/0.3241 = **17.38**; (4.3448−0.3478)/0.3241 = **12.33**). Omitting the 1.4826 consistency constant inflates the figure by ~48 % relative to the normal-calibrated scale the label "z" implies. Grep across the whole repo: `17.4` appears **only** in the `rogue_test.py` docstring (both copies) and `README.md:449` — **no script computes it and no output file contains it.** It is an unsourced number.
**Required action (compute-executable now):** either add the computation to `rogue_test.py` so it emits into `rogue_test.txt`, using (x−median)/(1.4826·MAD) and reporting 11.7/8.3, or drop the z and keep the defensible rank statement ("longest and second-longest terminal branch; 5.98 and 4.35 vs a median of 0.348"), which is what the argument actually uses. The branch-length distribution is nowhere near normal, so no z of any flavour carries much meaning here.

### MINOR-3 — C34: "roughly 20 points" holds for one of the three leave-one-out runs
From `rogue_test.txt`, 9-tip clade: full 71 → **drop_F59 91** (+20) → **drop_T21 80** (+9) → **drop_both 81** (+10). The 13-tip clade: 69 → 91 / 77 / 79. So the ~20-point figure is the `drop_F59` result specifically, and — counter-intuitively — dropping **both** rogues improves support *less* than dropping `F59D12.1` alone. That is worth a sentence, since it is mild evidence that the two long branches are not simply additive attractors. The "**8/8 clades identical**" claim is confirmed in all three runs.
**Required action (compute-executable now):** reword C34 and `README.md` §6c to give the range (+9 to +20) and note the non-additivity.

### MINOR-4 — C7/C8: 29/19 counts assemblies, not isotypes
`assembly_survey.tsv` has exactly 56 rows, 1 `ref` + 7 `lab` + 48 `wild`, and the wild split is **29 two-copy / 19 one-copy** ✓ — C7 reproduces at the assembly level, and the four replicate pairs agree exactly as claimed (CB4856 1 copy ×3, DL238 1 ×2, ECA36 1 ×2, EG4725 2 ×2; zero discordant isotypes). But those replicates mean the 48 wild assemblies are **43 unique isotypes**, and the isotype-level split is **28 carriers / 15 non-carriers**. C7 says "among natural isolates" and C8 says "the 29 carriers are essentially one haplotype observed 29 times" — both read as isotype counts. It is 28.
**Required action (compute-executable now):** add an isotype-collapsed column to `assembly_survey.tsv` and quote both numbers ("29 of 48 assemblies, 28 of 43 isotypes").

### MINOR-5 — C11: the "56 isolates" denominator contains N2 twice; all substantive numbers survive
`dop6_cds_isolates.fna` holds **57** CDSs: the 56 assemblies plus a separate `N2_WBcel235_reference` entry. `BY250` is excluded for an indel, leaving 56 "usable". I verified that `N2_WBcel235_reference` is **byte-identical** to `N2_reference_(Sanger_cosmids)` (WBcel235 *is* GCA_000002985.3), so N2 is counted twice and the "all isolates n=56" group is 55 assemblies + a duplicate. It is also why the groups do not add up: carriers 36 + non-carriers 19 = **55**, because `N2_WBcel235_reference` has no row in the survey and falls into neither group.

**None of it changes a number.** I recomputed segregating sites four ways:

| set | n | N-half S | C-half S | nonsyn |
|---|---|---|---|---|
| as the script does | 56 | 3 | 7 | **0** |
| wild assemblies only | 48 | 3 | 7 | **0** |
| wild carriers | 29 | 0 | 1 | 0 |
| wild non-carriers | 19 | 3 | 5 | 0 |
| unique wild isotypes | 43 | 3 | 7 | **0** |

All 10 segregating sites are at **codon position 3** and all are synonymous, so **πN = 0 in both halves is robust** to every reasonable choice of denominator, and "carriers 1 site vs non-carriers 8" reproduces exactly. One caveat worth stating: the carrier group as coded includes 8 laboratory N2 derivatives, which trivially lowers carrier diversity — but restricting to wild carriers gives the same 0 and 1, so the C8 interpretation survives.
**Required action (compute-executable now):** drop the duplicate reference sequence from `dop6_cds_isolates.fna` or exclude it from the "all isolates" group, and relabel the denominator "55 assemblies (43 wild isotypes)".

### MINOR-6 — README/memory: exons 1–9 do not encode codons 1–405
`README.md` §1 and `memory/dop6_c24a8_6_duplication.md` both say *"Exons 1–9 (3,315 bp, codons 1–405) are single-copy"*. From `repeat_structure.txt`, exons 1–9 sum to **1,200 nt = codons 1–400**; 3,315 bp is the *genomic* unique span (4,320,863–4,324,177). Codons **401–405 sit in exon 10**, which is inside arm A and therefore duplicated. The protein-level statement ("residues 1–405 are single-copy in the protein, because C24A8.6a starts at Met406") is correct; attributing them to exons 1–9 is not.
**Required action (compute-executable now):** reword to "exons 1–9 (1,200 nt, codons 1–400) lie in the 3,315 bp unique genomic span; codons 401–405 fall in exon 10, inside arm A."

### MINOR-7 — C18: one structure's Gα class contradicts its own PDB title
`gprotein_interface_multi.tsv` assigns **6WHA** to **Gs** at 88.6 % identity, while the title it stores in the same row reads *"HTR2A BOUND TO 25-CN-NBOH IN COMPLEX WITH A **MINI-GALPHA-Q** PROTEIN…"*. The 88.6 % match to Gs is consistent with 6WHA's mini-Gs/q chimera (Gs backbone, Gq α5), so the sequence-based call is defensible — but the row is self-contradictory as written, and it changes the class tally to **Gs n=4, Gi1 n=4, Go n=1, G11 n=1, Gq n=1**. C18's *"11 complexes spanning Gs, Gi, Go, Gq, G11"* and *"per class 29–41 %, no class rescues it"* therefore rest on **n = 1 for three of the five classes**.
**Required action (compute-executable now):** add a `chimera` column and record 6WHA as `mini-Gs/q`; state the per-class n alongside the per-class percentage so the reader sees that Go, Gq and G11 are single structures.

### MINOR-8 — C20: the positive control is 35 of ~44 human class-A aminergic receptors
Recounted from `pocket_332_v2.tsv`: **38 human tips, 35 return D**, the 3 non-D being `hGPR50` (T), `hMTNR1A` (M), `hMTNR1B` (M) — correctly, all melatonin-family. All 80 tips map to alignment column **694** and 74/80 carry D. **C20 reproduces exactly.** But the human set omits **hCHRM4 (P08173)** entirely — the panel has M1, M2, M3, M5 — and omits hTAAR2–9 and hHTR5B. A positive control is the instrument that caught the fabricated accessions, so its completeness matters: "35/35 human aminergic receptors" should read "35/35 of the 35 human aminergic receptors in the panel (hCHRM4 and hTAAR2–9 are not included)".
**Required action (compute-executable now):** add `P08173` (hCHRM4) to `v2_anchor_panel/provenance.tsv` and re-run `pocket_332_v2.py`; it will take ~10 min and makes the control 36/36.

### MINOR-9 — README §3b: the "disrupted larger CNV" tallies for both strains are wrong
Verified against Pu et al. sheet f, which I parsed in full (1,725 alleles; class distribution **1,641 `disrupted` / 82 `disrupted larger CNV` / 2 `possibly disrupted`**, so README's "82 of 1,725 (4.8 %)" ✓ exact). But:

| strain | README says | sheet f says |
|---|---|---|
| CHS1223 | "has one" | **two** — `yum2349` (dop-6) and `yum2350` (C24A8.6) are **both** `disrupted larger CNV` |
| CHS1256 | "has three" | **two** — `yum2601` (dop-4) and `yum2605` (ser-6); its `yum2603`/`yum2604` are class **`disrupted`**, not CNV |

Reading the two dop-6/C24A8.6 alleles as "one event" is the project's own thesis, but then CHS1256's count cannot be three on the same logic. The literal sheet counts are 2 and 2 — which, against a 4.8 % base rate, is 2/5 alleles for CHS1223 and 2/8 for CHS1256, i.e. **a stronger version of the point the paragraph is making**.
**Required action (compute-executable now):** correct both counts from the sheet.

### MINOR-10 — C14 omits a fourth chromosome-X lesion in CHS1256
Sheet f gives CHS1256 **eight** alleles. Besides `dop-1(yum2598)` at X:4,415,971, it carries **`dop-3(yum2600)` at CHROMOSOME_X:6,569,413** — a second X break beyond the dop-6 pair, 2.24 Mb away. C14's "three breaks inside 90 kb of X" is literally true but the strain's X burden is four lesions, not three. This again strengthens rather than weakens C15.
**Required action (compute-executable now):** add the dop-3 X locus to C14 and README §3b.

### MINOR-11 — `ppa_dop5_search.txt` self-certifies a coverage it did not have
The file states *"proteome UP000005239: **13,985 proteins** ← coverage check: a negative against a proteome this size is a fair test, not a database gap"*. `~/data/proteomes/MANIFEST.tsv` records `Ppac UP000005239` as **26,122 / 26,122, complete**, and I counted 26,122 `>` records in the FASTA. So the run used **54 %** of the proteome, having been written (21:39) before the complete download landed (23:43). W4 said the first run used 15,999/26,122 and was "being re-run complete"; the file in `output/` is a *third*, different, still-truncated run.

**The conclusion is unaffected** — I re-ran the full six-query RBH against the complete 26,122-protein file and reproduced the output *exactly* (same five accessions, bits 329/173/513/501/764/848, same reciprocation pattern including dop-2 mapping back to dop-3). A positive found in a subset is still a positive. But the sentence certifying coverage is false, and W4's stated remedy was never applied to the file.
**Required action (compute-executable now):** re-run `ppa_dop5_search.py` against the complete store file so the printed protein count matches the manifest.

### MINOR-12 — the boundary is scored on different protein isoforms than the panel contains
The HMM scores that define the ≥120 rule come from a table whose `protein_len` often differs from the panel sequence: `gar-1` 682 aa scored vs **713 aa** in `basis_b120.tsv`; `dop-2` 705 vs **849**; `ser-2` 432 vs **455**; `dop-3` 605 vs **607**. `gar-1` at **124.4** is the gene nearest the cut, 4.4 points above it, scored on a shorter isoform than the one in the tree. The rule is defensible but is not applied to the objects it admits.
**Required action (deferred-with-owner — MOD, alongside the `WormPooledCloning` HMM work):** rescore the panel sequences themselves with the monoamine HMM and confirm `gar-1` stays above 120; if it does not, the boundary needs restating rather than the gene removing.

### MINOR-13 — the CHAR evidence tiers contradict the table they cite
`independent_evidence.py` records `cel_tyra-2` and `cel_tyra-3` as **"IBA only"**, citing *"Koelle 2018 WormBook Table 1 and the lab's curated pharmacology table"*. Koelle Table 1, which I retrieved verbatim, records for **TYRA-2**: *"Binding studies, heterologous cell signaling"*, Gαi/o, tyramine EC50 50 nM (Rex et al. 2005); for **TYRA-3**: *"Binding studies, mutant resistant to tyramine"*, Gαq, tyramine EC50 70.8 nM (Hapiak et al. 2013). Both are direct-assay entries, not IBA. Neither gene is near dop-6, so the 5/9/13-tip counts are unaffected, but the tiers feed the trait-signal labelling.
**Required action (compute-executable now):** correct the two CHAR entries against Koelle Table 1 and re-run `confirm_tree.py` / `independent_evidence.py`.

### MINOR-14 — the Pu supplement is not deposited, and the README still calls it blocked
C12–C15 are "the most consequential result in the sub-project" and rest entirely on a file that is not in the repository (`dop6_c24a8_6/data/` is empty; no `.xlsx` anywhere in the repo relates to Pu et al.). `README.md` §"Next" item 1 still lists it as *"BLOCKED, needs a human"*. It downloads without authentication (I fetched all 812 KB and parsed every relevant row).
**Required action (compute-executable now):** `curl` the supplement to `dop6_c24a8_6/data/41467_2023_44177_MOESM4_ESM.xlsx`, commit it (or a checksum plus a fetch script if licensing forbids redistribution), and delete the stale BLOCKED entry.

---

## What I verified positively, and one thing the project under-claims

Stating this explicitly, because a report that only lists problems is not an audit.

**Locus structure (C1–C5) is exact.** I recomputed the palindrome from the raw FASTA `locus_X_4318000-4334500.fa` rather than from `repeat_structure.txt`: arm A 4,324,178–4,328,106 reverse-complemented against arm B 4,329,303–4,333,231 gives **3,929/3,929, zero mismatches, 100.000000 %**. Spacer 1,196 bp, containing C24A8.5 (4,328,400–4,328,916) ✓. `C24A8.6a.1` is an **exact substring** of DOP-6 at positions 406–692, **287/287**, and DOP-6[406] is **M** ✓. CDS lengths 2,079 / 864 confirm the 864-nt shared mRNA of C4 ✓. Exons 10–17 lie wholly inside arm A ✓.

**C12 is confirmed at the data level, and the project under-claims it.** Pu et al. sheet d gives `dop-6` and `c24a8.6` the identical guide `CAAGTTGGTTGCAACTGTAA`, identical PAM `TGG`, identical ssODN (`CGATTCGTCGAGAGACTCGTGCAACCAAGTTGGTT TAA AAGCTT TAA actgaaact…`) and identical restriction enzyme (HindIII), with only the genotyping primers differing — exactly as claimed. I searched the locus and found **exactly two** perfect full-length matches, at **X:4,325,604–4,325,623 (plus strand, PAM TGG)** and **X:4,331,786–4,331,805 (minus strand, PAM TGG)** — the stated coordinates to the base, cut sites 6,169 bp apart ✓ ("6.2 kb"). **The under-claim:** across all 1,672 distinct guides in sheet d, only **19** serve more than one gene, and every other one of those 19 is an adjacent paralogue pair (`srh-227/228`, `str-5/6`, `sri-13/14`, `srd-71/72`…). dop-6/c24a8.6 is one of 19 documented one-guide-two-genes cases in Pu's own design sheet — that is a stronger, independent confirmation than the argument currently makes.

**C13 is confirmed verbatim.** Both CHS1223 alleles class `disrupted larger CNV`; `yum2349` = "Inversion with insertion", insertion size **50**, `CHROMOSOME_X:4325629-4330411`; shared comment *"4kb deletion into dop-6 and C24A8.5, plus 75kb TD spanning multiple genes"* — exact. The other three alleles (F38E1.6 on V, F09C6.16 on V, K01A6.6 on IV) are all plain `disrupted` with perfect ssODN matches ✓.

**C14's core is confirmed verbatim.** `yum2603` and `yum2604` both read `not detected by Pindel`, comment `manually curated via IGV (see methods)`, mutation location **empty**, detected variation **empty** ✓. `dop-1(yum2598)` at `CHROMOSOME_X:4415971-4415987` ✓, 90,351 bp from the arm-A cut site ✓. `ser-6(yum2605)` = `disrupted larger CNV`, *"large deletion, part of ser-6 (2xINV)"*; `dop-4(yum2601)` = `disrupted larger CNV`, *"2xINV, deletion of dop-4 exon 1"* ✓. **A point stronger than the one made:** despite having no location, no called variation and no Pindel detection, Pu et al. class `yum2603`/`yum2604` as **`disrupted`** — the *confident* class, not the CNV class. The genotype is asserted at full confidence on manual inspection of unmappable reads. C15 is if anything conservative.

**C37 is confirmed verbatim from the primary source.** Koelle 2018 Table 1, group *"Poorly-studied putative small-molecule GPCRs"*, contains exactly six rows with **every** ligand/coupling/assay/reference column empty: `T02E9.3`, `C24A8.1`, `C24A8.6`, `F35H10.10`, `F59D12.1`, `T21B4.4`. `C24A8.6` **is** listed separately ✓. The identifications behind C34/C37 also check out by protein length: `pcdr-1` = WBGene00010329, 52.8, **318 aa** = `cel_F59D12.1` (318 aa in `basis.tsv`); `ador-1` = WBGene00011878, 85.2, **362 aa** = `cel_T21B4.4` (362 aa); `F35H10.10` = WBGene00018073, **`no_hit`**, 1,404 aa = `cel_F35H10.10` (1,404 aa).

**The b120 rule's own justification checks out.** Largest gap in the top 60 of the ranking = **27.0 points** between `gar-1` (124.4) and `npr-35` (97.4) ✓; the 19 genes above it are exactly the 16 `known_monoamine` plus `gar-1/2/3` ✓; the `cel_` ingroup in `basis_b120.tsv` is exactly those 19 ✓.

**C23/C25 alignment figures reproduce.** Full-length L-INS-i = **164 × 3,270 at 84.4 % gaps** ✓; trimmed gt0.1 = **163 × 362 at 11.4 % gaps** ✓; **356 parsimony-informative sites** from `gt0.1.iqtree` ✓; unrooted 163-taxon tree has 2(163)−3 = **323** branches, 356/323 = **1.10 sites per free parameter** ✓.

**C18's arithmetic is exact.** The 20 consensus positions in `gprotein_survival.txt` split at residue 405 into **13 lost** (57, 120, 121, 124, 125, 128, 129, 131, 132, 133, 206, 210, 213) and **7 retained** (467, 470, 471, 474, 475, 534, 535) = **35.0 %** ✓. Per-structure `lost + kept = n_contacts` in all 11 rows ✓. Per-class means 29–41 % ✓. (README's prose lists the lost ICL2/TM5 residues but silently omits D120 and T57 from the 13.)

---

## Reproduction table

| claim | asserted | my recomputation | verdict |
|---|---|---|---|
| **C1** | arm A vs rc(arm B) 3,929/3,929 = 100.000 %, 0 mismatches | 3,929/3,929, 0 mismatches, 100.000000 % (from raw FASTA) | **CONFIRMED** |
| **C1** | spacer 1,196 bp, contains C24A8.5 | 1,196 bp; C24A8.5 at 4,328,400–4,328,916 inside it | **CONFIRMED** |
| **C2** | C24A8.6a = DOP-6 406–692, 287/287, start = Met406 | exact substring, 287/287, DOP-6[406] = M | **CONFIRMED** |
| **C4** | identical 864 nt in both mRNAs | C24A8.6a CDS = 864 nt = 287×3+3 | **CONFIRMED** |
| **C7** | 56 assemblies; 29 carriers / 19 non-carriers (natural isolates) | 56 rows; wild 29/19 ✓ — but **28/15 across 43 unique isotypes** | CONFIRMED at assembly level; **MINOR-4** at isotype level |
| **C7** | replicates agree (CB4856 ×3, DL238 ×2, ECA36 ×2, EG4725 ×2) | all four agree; zero discordant isotypes | **CONFIRMED** |
| **C11** | C-half segregating sites 0 (short reads) vs 7 (assemblies) | assemblies 7 ✓ (robust to wild-only and isotype-collapsed sets) | **CONFIRMED** |
| **C11** | πN = 0 both halves; all 10 sites synonymous | 3 + 7 = 10, all at codon position 3, all synonymous, nonsyn = 0 | **CONFIRMED** |
| **C11** | carriers 1 site vs non-carriers 8 | 0+1 = 1 and 3+5 = 8 (same in wild-only) | **CONFIRMED** |
| **C11** | "56 isolates" | 55 assemblies + a duplicate N2 reference; groups sum to 55 | **DISCREPANT (MINOR-5)** |
| **C16** | dop-6 dN 0.0252, **the lowest** of the aminergic receptors tested | 0.0252 ✓ but **4th of 15**; tyra-2 0.0116, octr-1 0.0236, ser-5 0.0246 lower | **REFUTED (MAJOR-1)** |
| **C16** | the 8 other listed values | all 8 reproduce to ±0.0001 from S/N/Sd/Nd + JC | **CONFIRMED** |
| **C18** | 20 consensus contacts; 7 retained, 13 lost, ~35 % | 7/20 = 35.0 %; 13 lost; per-structure sums consistent | **CONFIRMED** |
| **C18** | 11 complexes spanning Gs/Gi/Go/Gq/G11 | 11 ✓, but 6WHA labelled Gs vs its own "mini-Gαq" title; n=1 for Go, Gq, G11 | CONFIRMED with **MINOR-7** |
| **C20** | 3 structures → one column; 35/35 human return D | all 80 tips at column **694**; 38 human tips, **35 D**, 3 non-D = GPR50/MTNR1A/MTNR1B | **CONFIRMED** (set incomplete — MINOR-8) |
| **C21** | dop-6 D ×5, dop-5 G ×5, F59D12.1 Q, T21B4.4 L, melatonin M/M/T | all as stated in `pocket_332_v2.tsv` | **CONFIRMED** |
| **C23** | full-length 84.4 % gaps, 164 seqs | 164 × 3,270, **84.4 %** | **CONFIRMED** |
| **C23** | excerpts: 11.4 % gaps, 98 % informative | 163 × 362, **11.4 %**; 356/362 = 98.3 % | **CONFIRMED** |
| **C24(a)** | all 20 Gα contacts inside retained blocks | 20/20 RETAINED in `gprotein_survival.txt`; excised region 240–450 | **CONFIRMED** |
| **C25** | 356 sites / 323 branches ≈ 1.1 | 2(163)−3 = 323; 356/323 = 1.10 | **CONFIRMED** |
| **C26** | 12 test tips form one clade, 98/100, 0 aminergic inside | exact 12-tip clade: **98.3/100** gt0.1, 98.1/100 gappyout & gt0.2, 97.8/100 gt0.3 | **CONFIRMED** (98.3, not 98) |
| **C26** | "the 12 **highest-scoring** excluded genes" | `ador-1` (85.2) excluded, `nmur-3` (78.9) in; two cited files disagree | **REFUTED (MAJOR-3)** |
| **C27** | `HMM nan` defect; did not change which genes were dropped | confirmed; cause = public-name vs sequence-name key; `DROP_CEL` short-circuits before the lookup | **CONFIRMED (B2 upheld)** |
| **C28** | 5-tip clade UFboot 100 at every trim level, both panels | 100 on all 8 trees (SH-aLRT 99.9–100) | **CONFIRMED** |
| **C28** | b120 gt0.1: +dop-5 95/95, +dop-3 85.7/90, +dop-2 88/90 | 95/95, 85.7/90, 88/90 | **CONFIRMED** |
| **C28** | b120 gappyout: 96.1/84, 83.3/60, 89/62 | 96.1/84, 83.3/60, 89/62 | **CONFIRMED** |
| **C28** | identical topology and tip counts 5→9→13→16 | clade **composition identical on all 8 trees**, also at 19 tips | **CONFIRMED** |
| **ADDENDUM** | b120 UFboot table (95/90/90/95/100; gappyout 84/60/62/65/99) | gt0.1 95/90/90/95/100; gappyout 84/60/62/65/**99**; gt0.2 87/79/79/**87**/98; gt0.3 94/82/82/**89**/100 | **CONFIRMED** (two "—" cells are computable: 87 and 89) |
| **ADDENDUM** | 8/8 composition agreement at all three other trim levels | 8/8 at gt0.2, gt0.3, gappyout | **CONFIRMED** |
| **C29** | 13-tip clade exclusively nematode; first outsider dmDop2R at 16 | 13-tip = 5 spp. dop-3/5/6 only; 16-tip adds Tsp_dop-2, cel_dop-2, **dmDop2R** | **CONFIRMED** |
| **C29 spp. table** | Cbr no dop-3; Ppa no dop-5 | true of the panel; **false of the genomes** (Q61H86; A0A2A6CA79) | **REFUTED as biology (MAJOR-6)** |
| **C32** | indep genes 2 / 3 / 4 at 5 / 9 / 13 tips | **1 / 2 / 3** — hyphen/underscore bug splits dop-6 in two | **DISCREPANT (MAJOR-8)** |
| **C32** | characterised 0 / 0 / 1; first is dop-3, four genes out | 0 / 0 / 1; dop-3 | **CONFIRMED** |
| **C34** | F59D12.1 branch 5.98, longest; T21B4.4 4.35; dop-6 0.020; median 0.348 | 5.9814 (longest), 4.3448 (2nd), 0.0201, median 0.3478 | **CONFIRMED** |
| **C34** | robust z = 17.4 and 12.3 | **11.7** and **8.3** with (x−med)/(1.4826·MAD); 17.4/12.3 uses raw MAD; no producing code | **DISCREPANT (MINOR-2)** |
| **C34** | leave-one-out: no topology change, 8/8 clades identical | 8/8 in all three runs | **CONFIRMED** |
| **C34** | "lift support by roughly 20 points (71 → 91)" | drop_F59 +20; drop_T21 +9; drop_both +10 | CONFIRMED for one of three (**MINOR-3**) |
| **C12** | one guide, two perfect matches, stated coordinates | exactly 2 matches, both PAM TGG, at the stated bases; identical ssODN in sheet d | **CONFIRMED** |
| **C13** | both alleles "disrupted larger CNV"; inversion + 50 bp at X:4,325,629–4,330,411 | verbatim match in sheet f | **CONFIRMED** |
| **C14** | both CHS1256 alleles "not detected by Pindel" / IGV, no location, no variation | verbatim match | **CONFIRMED** |
| README §3b | CNV class = 82/1,725 (4.8 %); CHS1223 one, CHS1256 three | 82/1,725 ✓; **CHS1223 two, CHS1256 two** | **DISCREPANT (MINOR-9)** |
| **C37** | Koelle Table 1 orphan group, six entries, all columns empty, C24A8.6 separate | exact: T02E9.3, C24A8.1, C24A8.6, F35H10.10, F59D12.1, T21B4.4 | **CONFIRMED** |
| **W4** | Ppa dop-5 = A0A2A6CA79, bits 764, clean RBH | reproduced on the **complete** 26,122-protein proteome: 764 fwd / 775 rev | **CONFIRMED** |
| **C35** | Harris 2009 methods (5 specific assertions) | — | **BLOCKED** |

**Tally: 38 asserted values confirmed, 8 discrepant** (C16 ordering; C26 test-set membership; C32 ×6 counts, scored as one; C34 z ×2, scored as one; C34 "20 points"; C11 denominator; C7 isotype denominator; README CNV tallies), **1 blocked** (C35).

---

## B7 — dataset diff

### (i) Pu et al. 2023, PMID 38110404 — Supplementary Data 1 (obtained and parsed: sheets a–f)

| dimension | Pu et al. | this project | direction |
|---|---|---|---|
| guide designs | **1,766 rows / 1,672 distinct guides**, whole GPCR family | 2 rows read (`dop-6`, `c24a8.6`) | theirs ⊃ ours |
| guides serving >1 gene | **19** (18 adjacent-paralogue pairs + dop-6/c24a8.6) | 1 discussed | **theirs ⊃ ours — the other 18 are not cited and are the natural control** |
| per-allele calls | **1,725 alleles**, 287 strains | 13 rows read (CHS1223 ×5, CHS1256 ×8) | theirs ⊃ ours |
| class labels | `disrupted` 1,641 / `disrupted larger CNV` 82 / `possibly disrupted` 2 | 82/1,725 quoted ✓ | agree |
| CHS1223 alleles | 5: yum2346 (F38E1.6, V), yum2347 (F09C6.16, V), yum2348 (K01A6.6, IV), **yum2349 (dop-6, X)**, **yum2350 (C24A8.6)** | all 5 | **agree** |
| CHS1256 alleles | 8: yum2598 (dop-1, X), yum2599 (dop-2, V), **yum2600 (dop-3, X:6,569,413)**, yum2601 (dop-4, CNV), yum2602 (dop-5, V), yum2603 (C24A8.6), yum2604 (dop-6), yum2605 (ser-6, IV, CNV) | 5 discussed | **theirs ⊃ ours — dop-3 on X is not mentioned (MINOR-10)** |
| guide sequence | `CAAGTTGGTTGCAACTGTAA`, PAM `TGG` | identical | agree |
| ssODN | identical for both genes, HindIII | identical | agree |
| genotyping primers | differ between the two genes | noted | agree |
| deposition | public, unauthenticated download | **not in the repo**; README calls it BLOCKED | **ours ⊂ theirs (MINOR-14)** |
| ours not in theirs | — | the 2-perfect-match genome search for the guide; the palindrome structure; the MAPQ-0 argument | **ours ⊄ theirs — this is the project's genuine contribution** |

### (ii) Koelle 2018 WormBook Table 1 (NBK535514) — obtained

| dimension | Koelle Table 1 | this project | direction |
|---|---|---|---|
| characterised small-molecule GPCRs | 22 rows: SER-1/4/5/7, DOP-1/2/3/4, OCTR-1, SER-3, SER-6, SER-2, TYRA-2, TYRA-3, GAR-1/2/3, GBB-1/2, MGL-1/2/3 | 19 in the b120 ingroup (16 monoamine + gar-1/2/3) | **theirs ⊃ ours — GBB-1/2 and MGL-1/2/3 correctly excluded (not class A)** |
| orphan group | **exactly 6, all columns empty**: T02E9.3, C24A8.1, **C24A8.6**, F35H10.10, F59D12.1, T21B4.4 | same 6 quoted | **agree exactly** |
| TYRA-2 evidence | "Binding studies, heterologous cell signaling"; Gαi/o; tyramine EC50 50 nM (Rex 2005) | `independent_evidence.py`: **"IBA only"** | **disagree (MINOR-13)** |
| TYRA-3 evidence | "Binding studies, mutant resistant to tyramine"; Gαq; EC50 70.8 nM (Hapiak 2013) | **"IBA only"** | **disagree (MINOR-13)** |
| SER-5 coupling | "Gαs **predicted from sequence analysis**, Harris et al. 2009" | C35's premise | **agree — Koelle's own hedge corroborates C35** |
| DOP-3 | "Mutant resistant to dopamine, heterologous cell signaling"; Gαo in ventral cord, Gαi/o in cells; **DA 27 nM, TA 500 nM, OA 250 µM** (Sugiura 2005) | C36 re-reads Sugiura; CHAR = ('DA','Go','IMP from mutant phenotype') | agree; C36's critique is of Sugiura, not of Koelle |
| DOP-1 coupling | Gαq (Chase 2004, Kindt 2007) | CHAR = Gq ✓ | agree |
| ours not in theirs | — | HMM scores for all 1,537 genes; the ≥120 boundary; the observation that C24A8.6 is not a gene | **ours ⊄ theirs — C37's substantive point** |

### (iii) Harris et al. 2009 panel, PMID 19193891 — **BLOCKED**

| dimension | Harris et al. | this project | direction |
|---|---|---|---|
| receptors in tree | **unverified** (C35: "5-HT receptors only") | 173 taxa, 5 ligand classes + outgroup | unverifiable |
| alignment | **unverified** (C35: ClustalW/MegAlign, hand-tuned) | MAFFT L-INS-i, no hand edits | unverifiable |
| tree method | **unverified** (C35: NJ in PAUP 4.0b10, no model selection) | IQ-TREE 3.1.3, ModelFinder, 1000 UFBoot + 1000 SH-aLRT, 4 trim levels | unverifiable |
| rooting | **unverified** (C35: GAR-1) | muscarinic outgroup (hCHRM1/2/3/5 + gar-1/2/3) | unverifiable |
| ICL3 | **unverified** (C35: hypervariable regions deleted) | ICL3 middle excised by HMM blocks; the excision is itself tested (C24) | unverifiable |
| validation | **unverified** (C35: one eyeball sentence) | Fitch vs 10,000 permutations, p = 0.0001, all 4 trims | unverifiable |
| corroborated from elsewhere | Koelle Table 1 records SER-5's Gαs as *"predicted from sequence analysis, Harris et al., 2009"* | C35's premise that the Gs call is sequence-based, not assayed | **one element independently corroborated** |

**Everything in C35 except the sequence-prediction premise rests on a reading I could not repeat.** See the BLOCKED block at the top.

---

## Summary of required actions by class

**Compute-executable now (18):** MAJOR-1 (rewrite C16 from the TSV), MAJOR-2a/b/c (fix + re-run `verify_accessions.py`; duplicate check; drop `lsOA`/`smTAR`), MAJOR-3 (rebuild with `TEST_ador-1`; delete the stale boundary files), MAJOR-4 repo half (README §6, two docstrings, two dead file refs), MAJOR-5 (withdrawal banner + relabel the `absent` cells), MAJOR-6 (add 6 nematode orthologues, rebuild b120), MAJOR-7 first half (replace `DROP_CEL` with a score filter), MAJOR-8 (one-line `gene_of` fix), MINOR-1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -13, -14.

**PI-only (2):** MAJOR-4 memory half — `memory/dop6_is_an_orphan_receptor.md` is user-owned and carries a fully withdrawn tree result plus the C16 error into every future session; this is the single highest-leverage fix in the report. MAJOR-7 second half — commit the b120 work, which is entirely untracked, so that ordering rests on git rather than on mutable mtimes.

**Deferred-with-owner (2):** MINOR-12, rescoring panel isoforms against the monoamine HMM — owner MOD, alongside the existing `WormPooledCloning` HMM work. The BLOCKED Harris 2009 retrieval — owner: whoever has institutional access to *J Neurosci*; C35 stays unauditable until then and should carry that caveat.

---

*Referee B. All recomputation scripts are in the session scratchpad; none wrote to `README.md` or `output/`.*
