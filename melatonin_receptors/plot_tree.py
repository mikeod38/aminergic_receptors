"""plot_tree.py — small figure of where PCDR-1 falls relative to the melatonin receptors.

The 120-taxon gt0.1 tree with the aminergic backbone collapsed to one node, so the
comparison that matters is legible: the melatonin clade, every sequence reciprocal
BLAST had called a melatonin hit, the PCDR-1 family, and the neuropeptide outgroup.

Style follows ../spatial_filter/blend_render.py.
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
import matplotlib.pyplot as plt
from Bio import Phylo

HERE = Path(__file__).parent
NP = {"hNPY2R", "hNPFFR2", "hTACR3", "hQRFPR", "hGPR52", "hGPR21"}
HMEL = {"hMTNR1A", "hMTNR1B", "hGPR50"}
SPECIES = {"Drer": "zebrafish", "Dmel": "fly", "Pcau": "priapulid", "Cvir": "oyster",
           "Lana": "lingula", "Myes": "scallop", "Nvec": "anemone", "Spur": "urchin",
           "Bflo": "amphioxus", "Cint": "ciona"}
COL = {"mel": "#c0392b", "rbh": "#e08030", "pcdr1": "#1a5fb4",
       "fam": "#7aa6d8", "np": "#888888", "amine": "#555555"}


def cls(n):
    if n in HMEL:
        return "mel"
    if n.startswith("PCDR-1"):
        return "pcdr1"
    if n.startswith("PCDRfam_"):
        return "fam"
    if n in NP:
        return "np"
    if n.startswith("MEL_"):
        return "rbh"
    return "amine"


def pretty(n):
    if n.startswith("PCDR-1"):
        return "C. elegans PCDR-1 / F59D12.1"
    if n.startswith("aminergic"):
        return n
    for pre in ("MEL_", "PCDRfam_"):
        if n.startswith(pre):
            sp, acc = n[len(pre):].split("_", 1)
            tag = "PCDR-1 family" if pre == "PCDRfam_" else "BLAST melatonin hit"
            return f"{SPECIES.get(sp, sp)} {acc[-4:]}  ({tag})"
    return n


tree = Phylo.read(HERE / "output/tree/gt0.1.treefile", "newick")
tree.root_with_outgroup(*[{"name": n} for n in NP])
for t in tree.get_terminals():
    if "pcdr-1" in t.name:
        t.name = "PCDR-1"

# collapse the largest clade that is entirely aminergic reference receptors
cand = [c for c in tree.get_nonterminals()
        if all(cls(x.name) == "amine" for x in c.get_terminals())]
big = max(cand, key=lambda c: len(c.get_terminals()))
n_amine = len(big.get_terminals())
big.clades = []
big.name = f"aminergic reference receptors  ({n_amine} tips, collapsed)"

terms = tree.get_terminals()
fig, ax = plt.subplots(figsize=(7.0, 6.4))
Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False,
           label_func=lambda c: pretty(c.name) if c in terms else "")

for t in ax.texts:
    raw = t.get_text().strip()
    match = next((x.name for x in terms if pretty(x.name).strip() == raw), None)
    k = cls(match) if match else "amine"
    t.set_color(COL[k])
    t.set_fontsize(5.2)
    if k in ("pcdr1", "mel"):
        t.set_fontweight("bold")
        t.set_fontsize(6.6)
    if match and match.startswith("aminergic"):
        t.set_fontsize(6.6)
        t.set_style("italic")

mel = tree.common_ancestor(list(HMEL))
fam = tree.common_ancestor(["PCDR-1", "PCDRfam_Dmel_Q9VZD1", "PCDRfam_Pcau_A0ABM1EQA9",
                            "PCDRfam_Cvir_A0A8B8D717", "PCDRfam_Lana_A0A1S3KFU7"])
for clade, label, col in ((mel, "100/100", COL["mel"]),
                          (fam, "99.7/100", COL["pcdr1"])):
    ys = [i + 1 for i, t in enumerate(terms) if t in clade.get_terminals()]
    x = min(tree.distance(t) for t in clade.get_terminals())
    ax.annotate(label, (x - 0.08, sum(ys) / len(ys)), fontsize=7, color=col,
                ha="right", va="center", fontweight="bold")

ax.set_xlabel("substitutions per site", fontsize=8)
ax.set_ylabel("")
ax.set_yticks([])
ax.tick_params(labelsize=7)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.set_title("PCDR-1 / F59D12.1 does not group with the melatonin receptors\n"
             "120 taxa, gt0.1 trim, LG+F+I+R7, rooted on neuropeptide receptors; "
             "same topology at gappyout",
             fontsize=8.5)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(HERE / f"output/tree/pcdr1_placement.{ext}", dpi=200, bbox_inches="tight")
print("wrote", HERE / "output/tree/pcdr1_placement.png", "| collapsed", n_amine, "aminergic tips")
