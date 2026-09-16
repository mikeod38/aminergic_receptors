#!/usr/bin/env python3
"""plot_tree.py — figure of the aminergic GPCR tree, with DOP-6 in its phylogenetic context.

Draws the ML phylogram with every tip labelled by its experimentally determined transmitter and
Galpha coupling, so DOP-6's placement can be read against both at a glance. The point the figure has
to make is that the octopamine receptors are scattered across three separate parts of the tree and
DOP-6 is in none of them; its nearest characterised relatives are two Gq-coupled 5-HT receptors.

House style (matching sos/sos_style.py): Arial, TrueType-embedded PDF, gene names italic.

Usage: plot_tree.py [trim_level]   (default gt0.1)
Writes output/tree/dop6_tree_<trim>.pdf and .png
"""
import os, re, sys
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__)); T = os.path.join(HERE, 'output', 'tree')
trim = sys.argv[1] if len(sys.argv) > 1 else 'gt0.1'

LIG = {
 'cel_dop-1':'DA','cel_dop-2':'DA','cel_dop-3':'DA','cel_dop-4':'DA','hDRD1':'DA','hDRD2':'DA',
 'dmDop1R1':'DA','dmDop2R':'DA',
 'cel_ser-1':'5-HT','cel_ser-4':'5-HT','cel_ser-7':'5-HT','hHTR1A':'5-HT','hHTR2A':'5-HT',
 'cel_ser-5':'5-HT?',
 'cel_ser-3':'OA','cel_ser-6':'OA','cel_octr-1':'OA','dmOAMB':'OA','dmOctA2R':'OA',
 'dmOctB1R':'OA','dmOctB2R':'OA','dmOctB3R':'OA','bmOAR1':'OA',
 'cel_ser-2':'TA','cel_tyra-2':'TA','dmTyrR_Q9VEG1':'TA','dmTyrRII_Q9VEG2':'TA',
 'cel_tyra-3':'TA/OA','dmOctTyrR_P22270':'TA/OA',
 'hADRA1A':'NE','hADRA1B':'NE','hADRA1D':'NE','hADRA2A':'NE','hADRA2C':'NE','hADRB1':'NE',
 'hADRB2':'NE','plADRA1':'NE',
 'hCHRM1':'ACh','hCHRM3':'ACh','cel_gar-3':'ACh',
}
LIGCOL = {'OA':'#C15C3E', 'TA':'#E0A458', 'TA/OA':'#D4813F', 'NE':'#8E3E29',
          'DA':'#1A2C7A', '5-HT':'#3C8DAD', '5-HT?':'#9EC5D4', 'ACh':'#7A7A7A',
          'orphan':'#FFFFFF'}
COUP = {}
for ln in open(f'{T}/coupling.tsv').read().strip().split('\n')[1:]:
    p = ln.split('\t')
    if len(p) > 1 and p[1].strip() and p[1] != '?': COUP[p[0]] = p[1]
# pretty tip names; gene names get italics, accessions and species tags do not
PRETTY = {'dmOctTyrR_P22270':'Oct-TyrR', 'dmTyrR_Q9VEG1':'TyrR', 'dmTyrRII_Q9VEG2':'TyrRII'}
# explicit species tag + display name; avoids regex-stripping accessions out of the ortholog labels
SPECIES = {
 'Cbn-dop-6_CBN19707.1':    ('Cbn', 'dop-6'),
 'Cbr-dop-6_CBG14393b.1':   ('Cbr', 'dop-6'),
 'Cre-dop-6_XM_053735786.1':('Cre', 'dop-6'),
 'Ppa-dop-6_PPA00437.1':    ('Ppa', 'dop-6'),
}
NOITAL = {'N2'}

def parse(nw):
    nw = nw.strip().rstrip(';'); pos=[0]
    def node():
        if nw[pos[0]] == '(':
            pos[0]+=1; kids=[]
            while True:
                kids.append(node())
                if nw[pos[0]]==',': pos[0]+=1; continue
                if nw[pos[0]]==')': pos[0]+=1; break
            m=re.match(r'[^,\)\(:]*',nw[pos[0]:]); lab=m.group(0); pos[0]+=len(lab)
            bl=0.0
            if pos[0]<len(nw) and nw[pos[0]]==':':
                m2=re.match(r':([-\d.eE+]*)',nw[pos[0]:]); pos[0]+=len(m2.group(0)); bl=float(m2.group(1) or 0)
            return {'tip':False,'lab':lab,'kids':kids,'bl':bl}
        m=re.match(r'[^,\)\(:]*',nw[pos[0]:]); lab=m.group(0); pos[0]+=len(lab)
        bl=0.0
        if pos[0]<len(nw) and nw[pos[0]]==':':
            m2=re.match(r':([-\d.eE+]*)',nw[pos[0]:]); pos[0]+=len(m2.group(0)); bl=float(m2.group(1) or 0)
        return {'tip':True,'lab':lab,'kids':[],'bl':bl}
    return node()

tree = parse(open(f'{T}/trees/{trim}.treefile').read())

# ── layout: y by tip order, x by cumulative branch length ────────────────────
ys = {}; order = []
def assign(n, depth):
    n['x'] = depth + n['bl']
    if n['tip']:
        order.append(n); n['y'] = len(order) - 1.0
    else:
        for k in n['kids']: assign(k, n['x'])
        n['y'] = sum(k['y'] for k in n['kids']) / len(n['kids'])
assign(tree, 0.0)
ntip = len(order)
xmax = max(n['x'] for n in order)

fig_h = 0.155 * ntip + 0.9
fig_w = 7.2
fig = plt.figure(figsize=(fig_w, fig_h))
ax = fig.add_axes([0.012, 0.055/fig_h*6, 0.50, 1 - 0.30/fig_h])
ax.set_axis_off()

def draw(n):
    if not n['tip']:
        y0 = min(k['y'] for k in n['kids']); y1 = max(k['y'] for k in n['kids'])
        ax.plot([n['x'], n['x']], [y0, y1], color='#333333', lw=0.7, solid_capstyle='round')
        for k in n['kids']:
            ax.plot([n['x'], k['x']], [k['y'], k['y']], color='#333333', lw=0.7, solid_capstyle='round')
            draw(k)
draw(tree)

# support labels: SH-aLRT/UFboot, shown when UFboot >= 50
def sup(n):
    if n['tip'] or not n['lab']: return
    m = re.match(r'([\d.]+)/([\d.]+)', n['lab'])
    if m and float(m.group(2)) >= 50:
        ax.text(n['x'] - 0.012*xmax, n['y'] + 0.28, f"{float(m.group(1)):.0f}/{float(m.group(2)):.0f}",
                fontsize=3.6, ha='right', va='center', color='#777777')
    for k in n['kids']: sup(k)
sup(tree)

# ── tips ─────────────────────────────────────────────────────────────────────
xlab = xmax * 1.035
for n in order:
    lab = n['lab']
    lig = LIG.get(lab, 'orphan')
    col = LIGCOL.get(lig, '#FFFFFF')
    ax.plot([n['x'], xlab - xmax*0.012], [n['y'], n['y']], color='#CCCCCC', lw=0.35, zorder=0)
    ax.add_patch(Rectangle((xlab, n['y'] - 0.30), xmax*0.022, 0.60, facecolor=col,
                           edgecolor='#555555', lw=0.4, clip_on=False, zorder=3))
    if lab in SPECIES:
        species, disp = SPECIES[lab]
    else:
        disp = PRETTY.get(lab, lab)
        species = ''
        if lab.startswith('cel_'):  species, disp = 'Cel', disp[4:]
        elif lab.startswith('dm'):  species, disp = 'Dme', disp[2:] if disp.startswith('dm') else disp
        elif lab.startswith('bm'):  species, disp = 'Bmo', disp[2:]
        elif lab.startswith('pl'):  species, disp = 'Pdu', disp[2:]
        elif lab.startswith('h'):   species, disp = 'Hsa', disp[1:]
    bold = (lab == 'cel_dop-6')
    if species:
        ax.text(xlab + xmax*0.032, n['y'], species, fontsize=4.4, va='center', ha='left',
                color='#A0A0A0', style='italic')
    ax.text(xlab + xmax*0.085, n['y'], disp, fontsize=5.4, va='center', ha='left',
            style='normal' if disp in NOITAL else 'italic',
            weight='bold' if bold else 'normal',
            color='#000000' if not bold else '#B03030')
    c = COUP.get(lab, '')
    if c:
        ax.text(xlab + xmax*0.255, n['y'], c.replace('*', ''), fontsize=4.6, va='center',
                ha='left', color='#555555', style='italic' if c.endswith('*') else 'normal')
    if bold:
        ax.annotate('', xy=(xlab - xmax*0.030, n['y']), xytext=(xlab - xmax*0.095, n['y']),
                    arrowprops=dict(arrowstyle='-|>', color='#B03030', lw=1.0))

ax.set_ylim(ntip - 0.4, -0.9)
ax.set_xlim(-xmax*0.01, xmax*1.40)

# scale bar
sb = round(xmax/6, 1) or 0.1
ax.plot([0, sb], [ntip - 0.1, ntip - 0.1], color='#333333', lw=0.9)
ax.text(sb/2, ntip + 0.35, f'{sb} substitutions/site', fontsize=4.8, ha='center', va='top',
        color='#333333')

# legend
handles = [Line2D([], [], marker='s', linestyle='none', markersize=4.2,
                  markerfacecolor=LIGCOL[k], markeredgecolor='#555555', markeredgewidth=0.4,
                  label=lab)
           for k, lab in [('OA','octopamine'), ('TA','tyramine'), ('TA/OA','tyramine/OA'),
                          ('NE','norepinephrine'), ('DA','dopamine'), ('5-HT','serotonin'),
                          ('ACh','acetylcholine'), ('orphan','orphan — no ligand data')]]
lg = fig.legend(handles=handles, loc='lower right', bbox_to_anchor=(0.995, 0.012),
                frameon=False, fontsize=4.9, handletextpad=0.45, labelspacing=0.34,
                title='experimentally determined ligand', title_fontsize=5.2)
lg.get_title().set_color('#333333')

fig.text(0.512, 0.985, f'Aminergic GPCR phylogeny  ·  {trim}, IQ-TREE Q.PFAM+F+R6',
         fontsize=6.2, va='top', ha='left', weight='bold')
fig.text(0.512, 0.963,
         'node labels SH-aLRT/UFboot (shown when UFboot ≥ 50); right column = Gα coupling\n'
         '(Koelle 2018 WormBook Table 1; * = sequence prediction, not measured)',
         fontsize=4.9, va='top', ha='left', color='#555555')
fig.text(0.512, 0.915,
         'DOP-6 sits in a clade of orphans (DOP-5, F59D12.1, T21B4.4); its nearest\n'
         'characterised relatives are two Gq-coupled 5-HT receptors. The octopamine\n'
         'receptors are scattered across three other parts of the tree.',
         fontsize=5.0, va='top', ha='left', color='#B03030')

for ext in ('pdf', 'png'):
    fig.savefig(f'{T}/dop6_tree_{trim}.{ext}', dpi=400, bbox_inches='tight',
                facecolor='white')
print(f'-> output/tree/dop6_tree_{trim}.pdf / .png   ({ntip} tips)')
