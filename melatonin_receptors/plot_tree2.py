"""plot_tree2.py — the 178-taxon tree as built, nothing collapsed.

Round 2, gt0.1 trim (the trim level where the melatonin clade is well supported).
Every tip shown. Colour is the only summarising done.
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio import Phylo

HERE = Path(__file__).parent
NP = {"hNPY2R", "hNPFFR2", "hTACR3", "hQRFPR", "hGPR52", "hGPR21"}
HMEL = {"hMTNR1A", "hMTNR1B", "hGPR50"}
COL = {"mel": "#c0392b", "melhit": "#e08030", "pcdr1": "#1a5fb4",
       "nbr": "#7aa6d8", "worm": "#2e7d32", "np": "#999999", "amine": "#555555"}
LEG = [("mel", "melatonin receptors (human + zebrafish)"),
       ("pcdr1", "PCDR-1 / F59D12.1"),
       ("worm", "C. elegans genes outscoring PCDR-1 (71)"),
       ("nbr", "PCDR-1 nearest neighbours (best BLAST hit per species)"),
       ("melhit", "other reciprocal-BLAST melatonin hits"),
       ("amine", "aminergic reference receptors"),
       ("np", "neuropeptide / orphan outgroup")]


def cls(n):
    if n in HMEL or (n.startswith("MEL_Drer")):
        return "mel"
    if "pcdr-1" in n:
        return "pcdr1"
    if n.startswith("CEL_"):
        return "worm"
    if n.startswith("PCDRfam_"):
        return "nbr"
    if n in NP:
        return "np"
    if n.startswith("MEL_"):
        return "melhit"
    return "amine"


def pretty(n):
    if "pcdr-1" in n:
        return "PCDR-1 / F59D12.1  (C. elegans)"
    if n.startswith("CEL_"):
        return f"C. elegans {n[4:]}"
    for p in ("MEL_", "PCDRfam_"):
        if n.startswith(p):
            return n[len(p):].replace("_", " ")
    return n


tree = Phylo.read(HERE / "output/tree2/gt0.1.treefile", "newick")
tree.root_with_outgroup(*[{"name": n} for n in NP])
terms = tree.get_terminals()

fig, ax = plt.subplots(figsize=(9.5, 20))
Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False,
           label_func=lambda c: pretty(c.name) if c in terms else "")
for t in ax.texts:
    raw = t.get_text().strip()
    m = next((x.name for x in terms if pretty(x.name) == raw), None)
    k = cls(m) if m else "amine"
    t.set_color(COL[k])
    t.set_fontsize(4.6)
    if k in ("pcdr1", "mel"):
        t.set_fontweight("bold")
        t.set_fontsize(6.0 if k == "pcdr1" else 5.4)

mel = tree.common_ancestor(list(HMEL))
pc = next(x for x in terms if "pcdr-1" in x.name)
for clade, lab, col in ((mel, "melatonin clade  100/100", COL["mel"]),):
    ys = [i + 1 for i, t in enumerate(terms) if t in clade.get_terminals()]
    ax.annotate(lab, (tree.distance(clade) - 0.1, sum(ys) / len(ys)), fontsize=7,
                color=col, ha="right", va="center", fontweight="bold")
yp = terms.index(pc) + 1
ax.annotate("PCDR-1 clade  99.8/100", (tree.distance(pc.name) * 0.45, yp + 2.5),
            fontsize=7, color=COL["pcdr1"], ha="right", va="center", fontweight="bold")

ax.legend(handles=[Line2D([], [], color=COL[k], lw=3, label=v) for k, v in LEG],
          loc="lower left", fontsize=7.5, frameon=True, framealpha=0.95, edgecolor="#ccc")
ax.get_legend().set_zorder(20)
ax.set_xlabel("substitutions per site", fontsize=8)
ax.set_ylabel("")
ax.set_yticks([])
ax.tick_params(labelsize=7)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.set_title("Melatonin receptors and PCDR-1 / F59D12.1 — 178 taxa, nothing collapsed\n"
             "gt0.1 trim, MAFFT L-INS-i, LG+F+I+R7, rooted on neuropeptide receptors",
             fontsize=10)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(HERE / f"output/tree2/tree_full.{ext}", dpi=170, bbox_inches="tight")
print("wrote", HERE / "output/tree2/tree_full.png")
