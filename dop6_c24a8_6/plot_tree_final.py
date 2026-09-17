#!/usr/bin/env python3
"""plot_tree_final.py — figure of the rebuilt 163-taxon aminergic tree.

Tips are coloured by experimentally determined ligand where one exists, and DOP-6's clade is
highlighted. Proteome-scan tips carry only an accession, so they are shown as species + accession
and deliberately left uncoloured - they are breadth, not evidence about ligand.

Usage: plot_tree_final.py [trim]   (default gt0.1)
"""
import os, re, sys
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype']=42; matplotlib.rcParams['ps.fonttype']=42
matplotlib.rcParams['font.family']='sans-serif'
matplotlib.rcParams['font.sans-serif']=['Arial','Helvetica','DejaVu Sans']
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

HERE=os.path.dirname(os.path.abspath(__file__)); T=os.path.join(HERE,'output','tree')
trim=sys.argv[1] if len(sys.argv)>1 else 'gt0.1'
LIG={'cel_dop-1':'DA','cel_dop-2':'DA','cel_dop-3':'DA','cel_dop-4':'DA',
 'hDRD1':'DA','hDRD2':'DA','hDRD3':'DA','hDRD4':'DA','hDRD5':'DA','dmDop1R1':'DA','dmDop2R':'DA',
 'cel_ser-1':'5-HT','cel_ser-4':'5-HT','cel_ser-7':'5-HT','dm5HT1A':'5-HT','dm5HT1B':'5-HT',
 'dm5HT2A':'5-HT','dm5HT7':'5-HT','cel_ser-5':'5-HT?',
 'cel_ser-3':'OA','cel_ser-6':'OA','cel_octr-1':'OA','dmOAMB':'OA','dmOctA2R':'OA',
 'dmOctB1R':'OA','dmOctB2R':'OA','dmOctB3R':'OA','bmOAR1':'OA','amOA1':'OA','tcOA1':'OA','apOA':'OA',
 'cel_ser-2':'TA','cel_tyra-2':'TA','dmTyrR':'TA','dmTyrRII':'TA','amTAR1':'TA','bmTAR':'TA',
 'cel_tyra-3':'TA/OA','dmOctTyrR':'TA/OA',
 'hCHRM1':'ACh','hCHRM2':'ACh','hCHRM3':'ACh','hCHRM5':'ACh','cel_gar-1':'ACh','cel_gar-2':'ACh',
 'cel_gar-3':'ACh','hMTNR1A':'MEL','hMTNR1B':'MEL','hGPR50':'MEL',
 'hHRH1':'HA','hHRH2':'HA','hHRH3':'HA','hHRH4':'HA','hTAAR1':'TA'}
for k in ['hHTR1A','hHTR1B','hHTR1D','hHTR1E','hHTR1F','hHTR2A','hHTR2B','hHTR2C','hHTR4',
          'hHTR5A','hHTR6','hHTR7']: LIG[k]='5-HT'
for k in ['hADRA1A','hADRA1B','hADRA1D','hADRA2A','hADRA2B','hADRA2C','hADRB1','hADRB2','hADRB3',
          'plADRA1','plADRA2','am5HT1A','amDop1']: LIG.setdefault(k,'NE')
LIG['am5HT1A']='5-HT'; LIG['amDop1']='DA'
COL={'OA':'#C15C3E','TA':'#E0A458','TA/OA':'#D4813F','NE':'#8E3E29','DA':'#1A2C7A',
     '5-HT':'#3C8DAD','5-HT?':'#9EC5D4','ACh':'#7A7A7A','MEL':'#6B4E8E','HA':'#4F8A5B',None:'#FFFFFF'}
SUBJECT={'cel_dop-6','cel_dop-5','cel_F59D12.1','cel_T21B4.4'}
# Galpha coupling, verbatim from Koelle 2018 WormBook Table 1 for C. elegans, textbook for the
# vertebrate/insect references. Blank where no coupling has been determined.
COUP={'cel_ser-1':'Gq','cel_ser-4':'Gi/o','cel_ser-7':'Gs','cel_dop-1':'Gq','cel_dop-2':'Gi/o',
 'cel_dop-3':'Go','cel_dop-4':'Gs','cel_octr-1':'Go','cel_ser-3':'Gq','cel_ser-6':'Gs',
 'cel_ser-2':'Gi/o','cel_tyra-2':'Gi/o','cel_tyra-3':'Gq',
 'hDRD1':'Gs','hDRD5':'Gs','hDRD2':'Gi/o','hDRD3':'Gi/o','hDRD4':'Gi/o',
 'hADRB1':'Gs','hADRB2':'Gs','hADRB3':'Gs','hADRA1A':'Gq','hADRA1B':'Gq','hADRA1D':'Gq',
 'hADRA2A':'Gi/o','hADRA2B':'Gi/o','hADRA2C':'Gi/o',
 'hHTR1A':'Gi/o','hHTR1B':'Gi/o','hHTR1D':'Gi/o','hHTR1E':'Gi/o','hHTR1F':'Gi/o',
 'hHTR2A':'Gq','hHTR2B':'Gq','hHTR2C':'Gq','hHTR4':'Gs','hHTR6':'Gs','hHTR7':'Gs',
 'hCHRM1':'Gq','hCHRM3':'Gq','hCHRM5':'Gq','hCHRM2':'Gi/o',
 'dmDop1R1':'Gs','dmDop2R':'Gi/o','dmOctB1R':'Gs','dmOctB2R':'Gs','dmOctB3R':'Gs',
 'dmOAMB':'Gq','dmOctA2R':'Gi/o','dmOctTyrR':'Gi/o','dmTyrR':'Gi/o','dmTyrRII':'Gi/o'}
POCK={}
_pt=os.path.join(T,'pocket.tsv')
if os.path.exists(_pt):
    for _ln in open(_pt).read().strip().split('\n')[1:]:
        _p=_ln.split('\t'); POCK[_p[0]]=_p[1]
def parse(nw):
    nw=nw.strip().rstrip(';'); pos=[0]
    def node():
        if nw[pos[0]]=='(':
            pos[0]+=1; kids=[]
            while True:
                kids.append(node())
                if nw[pos[0]]==',': pos[0]+=1; continue
                if nw[pos[0]]==')': pos[0]+=1; break
            m=re.match(r'[^,\)\(:]*',nw[pos[0]:]); lab=m.group(0); pos[0]+=len(lab); bl=0.0
            if pos[0]<len(nw) and nw[pos[0]]==':':
                m2=re.match(r':([-\d.eE+]*)',nw[pos[0]:]); pos[0]+=len(m2.group(0)); bl=float(m2.group(1) or 0)
            return {'tip':False,'lab':lab,'kids':kids,'bl':bl}
        m=re.match(r'[^,\)\(:]*',nw[pos[0]:]); lab=m.group(0); pos[0]+=len(lab); bl=0.0
        if pos[0]<len(nw) and nw[pos[0]]==':':
            m2=re.match(r':([-\d.eE+]*)',nw[pos[0]:]); pos[0]+=len(m2.group(0)); bl=float(m2.group(1) or 0)
        return {'tip':True,'lab':lab,'kids':[],'bl':bl}
    return node()
t=parse(open(f'{T}/trees/{trim}.treefile').read())
# optional second argument: focus on the smallest clade containing cel_dop-6 with >= N tips,
# because 163 tips in one image is not readable
FOCUS=int(sys.argv[2]) if len(sys.argv)>2 else 0
if FOCUS:
    def tips_(n): return [n['lab']] if n['tip'] else [x for k in n['kids'] for x in tips_(k)]
    best=[None]
    def walk_(n):
        if n['tip']: return
        tl=tips_(n)
        if 'cel_dop-6' in tl and len(tl)>=FOCUS:
            if best[0] is None or len(tl)<len(tips_(best[0])): best[0]=n
        for k in n['kids']: walk_(k)
    walk_(t)
    if best[0] is not None:
        t=best[0]; t['bl']=0.0
        print(f'focused on the {len(tips_(t))}-tip clade containing cel_dop-6')
order=[]
def assign(n,d):
    n['x']=d+n['bl']
    if n['tip']: order.append(n); n['y']=len(order)-1.0
    else:
        for k in n['kids']: assign(k,n['x'])
        n['y']=sum(k['y'] for k in n['kids'])/len(n['kids'])
assign(t,0.0)
N=len(order); xmax=max(n['x'] for n in order)
_tall = 0.132 if N>80 else 0.30
fig=plt.figure(figsize=(13.5 if N>80 else 12.0, _tall*N+1.3))
ax=fig.add_axes([0.006,0.035,0.74,0.90]); ax.set_axis_off()
def draw(n):
    if n['tip']: return
    y0=min(k['y'] for k in n['kids']); y1=max(k['y'] for k in n['kids'])
    ax.plot([n['x'],n['x']],[y0,y1],color='#333',lw=0.55,solid_capstyle='round')
    for k in n['kids']:
        ax.plot([n['x'],k['x']],[k['y'],k['y']],color='#333',lw=0.55,solid_capstyle='round')
        draw(k)
draw(t)
def sup(n):
    if n['tip'] or not n['lab']: return
    m=re.match(r'([\d.]+)/([\d.]+)',n['lab'])
    if m and float(m.group(2))>=70:
        ax.text(n['x']-0.006*xmax,n['y']+0.42,f"{float(m.group(2)):.0f}",fontsize=2.9 if N>80 else 5.4,
                ha='right',va='center',color='#888')
    for k in n['kids']: sup(k)
sup(t)
xlab=xmax*1.02
# column positions, as multiples of xmax. Generous and fixed, so long accessions cannot collide
# with the species or coupling columns - the previous version placed them relatively and they did.
# Columns are multiples of xmax measured from xlab. Species names run to ~15 characters
# ("D. melanogaster", "S. kowalevskii") and accessions to ~11, so the gene column has to start well
# clear of the species column rather than just after it - previous spacings still collided.
# Column order, as multiples of xmax from xlab:
#   ligand colour box (at xlab) | Galpha coupling | residue 3.32 | gene name | species
# The three functional annotations sit together next to the colour box; the species name goes last
# because it is the longest and least often needed.
CCOUP, C332, CGENE, CSP = (0.05, 0.26, 0.36, 0.95)
SPCODE={'Ajap':'A. japonica','Apla':'A. planci','Bflo':'B. floridae','Cint':'C. intestinalis',
 'Cvir':'C. virginica','Drer':'D. rerio','Ebur':'E. diaphana','Hvul':'H. vulgaris',
 'Lana':'L. anatina','Lpol':'L. polyphemus','Mlei':'M. leidyi','Myes':'M. yessoensis',
 'Nvec':'N. vectensis','Pcau':'P. caudatus','Pmar':'P. marinus','Pmin':'P. miniata',
 'Skow':'S. kowalevskii','Spur':'S. purpuratus','Tadh':'T. adhaerens','Aque':'A. queenslandica',
 'Cbr':'C. briggsae','Cbn':'C. brenneri','Cre':'C. remanei','Ppa':'P. pacificus',
 'Asu':'A. suum','Bma':'B. malayi','Sra':'S. ratti','Hco':'H. contortus','Tsp':'T. spiralis',
 'Nam':'N. americanus'}
_SCAN=re.compile(r'^([A-Z][a-z]{3})_(?:tr|sp)_([A-Za-z0-9]+)_')
_NEM=re.compile(r'^([A-Z][a-z]{2})_(.+)$')
_MLEI=re.compile(r'^([A-Z][a-z]{3})_(\S+)$')
def pretty(lab):
    """-> (species label, short name, italicise the name?)"""
    if lab.startswith('cel_'): return 'C. elegans', lab[4:], True
    m=_SCAN.match(lab)                      # proteome-scan hit: Xxxx_tr_ACC_ACC_ORG
    if m:
        acc=m.group(2)
        if acc.startswith('A0A'): acc=acc[3:]   # universal TrEMBL prefix, dropped for width
        return SPCODE.get(m.group(1),m.group(1)), acc, False
    m=_NEM.match(lab)                       # nematode ortholog: Xxx_gene
    if m and m.group(1) in SPCODE: return SPCODE[m.group(1)], m.group(2).replace('_','-'), True
    m=_MLEI.match(lab)                      # Mlei_KAL... style
    if m and m.group(1) in SPCODE: return SPCODE[m.group(1)], m.group(2), False
    for pre,sp in (('hGPR','H. sapiens'),('hMTNR','H. sapiens'),('h','H. sapiens'),
                   ('dm','D. melanogaster'),('am','A. mellifera'),('bm','B. mori'),
                   ('tc','T. castaneum'),('ap','A. californica'),('pl','P. dumerilii')):
        if lab.startswith(pre) and len(lab)>len(pre): return sp, lab[len(pre):], True
    return '', lab, True

for n in order:
    lab=n['lab']; lg=LIG.get(lab); c=COL.get(lg,'#FFFFFF')
    ax.plot([n['x'],xlab-xmax*0.006],[n['y'],n['y']],color='#DDD',lw=0.28,zorder=0)
    ax.add_patch(Rectangle((xlab,n['y']-0.34),xmax*0.014,0.68,facecolor=c,edgecolor='#666',
                           lw=0.3,clip_on=False,zorder=3))
    sp,gene,ital=pretty(lab)
    bold = lab=='cel_dop-6'
    # residue at 3.32 - the amine-binding aspartate. Red where it is absent.
    r332=POCK.get(lab,'')
    cp=COUP.get(lab,'')
    if cp:
        ax.text(xlab+xmax*CCOUP,n['y'],cp,fontsize=3.4 if N>80 else 6.4,va='center',ha='left',
                color='#444')
    if r332:
        ax.text(xlab+xmax*C332,n['y'],r332,fontsize=3.8 if N>80 else 7.0,va='center',ha='center',
                family='monospace',
                color='#000' if r332=='D' else '#C0392B',
                weight='normal' if r332=='D' else 'bold')
    ax.text(xlab+xmax*CGENE,n['y'],gene,fontsize=3.9 if N>80 else 7.2,va='center',ha='left',
            style='italic' if ital else 'normal',
            weight='bold' if bold else 'normal',
            color='#B03030' if bold else ('#7A3030' if lab in SUBJECT else '#000'))
    ax.text(xlab+xmax*CSP,n['y'],sp,fontsize=3.1 if N>80 else 5.6,va='center',ha='left',
            color='#AAA',style='italic')
    if bold:
        ax.annotate('',xy=(xlab-xmax*0.016,n['y']),xytext=(xlab-xmax*0.062,n['y']),
                    arrowprops=dict(arrowstyle='-|>',color='#B03030',lw=0.9))
ax.set_ylim(N-0.4,-1.2); ax.set_xlim(-xmax*0.01,xmax*2.35)
sb=round(xmax/6,1) or 0.1
ax.plot([0,sb],[N+0.2,N+0.2],color='#333',lw=0.8)
ax.text(sb/2,N+0.95,f'{sb} subs/site',fontsize=3.6,ha='center',va='top',color='#333')
h=[Line2D([],[],marker='s',linestyle='none',markersize=3.6,markerfacecolor=COL[k],
          markeredgecolor='#666',markeredgewidth=0.3,label=v)
   for k,v in [('OA','octopamine'),('TA','tyramine'),('TA/OA','tyramine/OA'),('NE','adrenergic'),
               ('DA','dopamine'),('5-HT','serotonin'),('HA','histamine'),('MEL','melatonin'),
               ('ACh','acetylcholine'),(None,'no ligand data')]]
lg=fig.legend(handles=h,loc='upper right',bbox_to_anchor=(0.998,0.985),frameon=False,fontsize=4.4,
              handletextpad=0.4,labelspacing=0.3,title='experimentally determined ligand',
              title_fontsize=4.8)
fig.text(0.60,0.995,f'Aminergic GPCR phylogeny, rebuilt  ·  163 taxa, 7TM domain excerpts  ·  {trim}',
         fontsize=5.6,va='top',ha='left',weight='bold')
fig.text(0.60,0.978,'node labels = UFboot (shown when ≥ 70). 22 C. elegans, 40 other nematode,\n'
         '38 human, 14 Drosophila, 50 other metazoan.',fontsize=4.2,va='top',ha='left',color='#555')
_sfx = f'_clade{N}' if FOCUS else ''
for ext in ('pdf','png'):
    fig.savefig(f'{T}/tree_final_{trim}{_sfx}.{ext}',dpi=420,bbox_inches='tight',facecolor='white')
print(f'-> output/tree/tree_final_{trim}{_sfx}.pdf / .png   ({N} tips)')
