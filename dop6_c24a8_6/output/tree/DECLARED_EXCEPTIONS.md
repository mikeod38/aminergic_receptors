# Tips kept below the monoamine HMM ≥ 120 boundary, and why

The panel boundary is **monoamine HMM ≥ 120**. Fifteen tips in the current panels score below it
and are kept deliberately. They are marked **§** in the figures and listed in
`declared_exceptions.tsv` so no sub-cutoff tip is silently present — which is precisely how
`F59D12.1` survived in the superseded v1 panel.

| tip | HMM | why it is kept |
|---|---|---|
| hMTNR1A | 102.5 | melatonin reference. DOP-5 carries **G at 3.32**, the melatonin-receptor state rather than the amine-binding aspartate; without these three tips that observation has no reference group in the tree. |
| hGPR50 | 98.3 | melatonin-related orphan, same comparison |
| hMTNR1B | 95.3 | melatonin reference, same comparison |
| 12 × `TEST_*` | 78.9–97.5 | the **boundary test set** — the 12 highest-scoring excluded genes, included to test the boundary from outside. They form one 12-tip clade at 98/100 with zero aminergic receptors inside it, which is the evidence that the boundary is sound. Removing them would remove that evidence. |

## Why F59D12.1 is NOT among them

`F59D12.1` scores **52.8** — far below the boundary, and below every tip above — and it failed on
three further independent grounds:

- the **longest terminal branch** in the v1 tree, 5.98 subs/site, robust z ≈ 11.7 against a median of 0.348
- **Q rather than D at 3.32**, so it lacks the amine-binding aspartate
- only **318 aa**, short for the family

Its exclusion does not rest on the HMM score alone, which is what separates it from the fifteen
above. `T21B4.4` (L at 3.32, second-longest branch) and `F35H10.10` (class C — 7tm_3 141.6 plus
ANF/VFT 57.9, no class-A signal) were dropped on the same combined grounds.

## Sensitivity check

Strict panels containing **only** tips ≥ 120 — b120 158 tips, placement 253 — are built in
`output/tree/strict_b120/` and `output/tree/strict_placement/`. They exist to answer one question:
does keeping the fifteen change DOP-6's placement? Compare the nested clade composition and support
against the primary trees. If they agree, the exceptions are cosmetic; if they do not, the primary
tree is the one that needs defending.
