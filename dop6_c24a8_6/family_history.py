#!/usr/bin/env python3
"""family_history.py — the duplication history of the C. elegans dopamine-receptor-like family.

Reads the order and the timing of the duplications off the 268-taxon placement tree, where every
nematode monoamine-HMM candidate >= 120 from 8 proteomes spanning clades I, III, IV and V was
included and NO orthology call was made - the topology assigns paralogy.

A duplication is dated by the deepest lineage present on BOTH sides of it: if clade IV sequences
sit inside both daughter clades, the duplication predates the clade IV/V split. This is the
inference reciprocal best hits cannot make, and is why the earlier RBH dating matrix was retracted.

Nematode phylogeny used for bounding (Blaxter clades):
    I   Dorylaimia    (Trichinella, Trichuris)          earliest-diverging
    III Spirurina     (Brugia)
    IV  Tylenchina    (Strongyloides)
    V   Rhabditina    (Caenorhabditis, Pristionchus, Necator, Haemonchus)
Order of divergence: I first, then III, then IV/V.

Writes output/tree/placement/family_history.txt
"""
import csv, os
from collections import Counter
from Bio import Phylo
HERE=os.path.dirname(os.path.abspath(__file__)); T=f'{HERE}/output/tree/placement'
TRIMS=('gt0.1','gappyout','gt0.2','gt0.3')
meta={r['tip']:r for r in csv.DictReader(open(f'{T}/candidates.tsv'),delimiter='\t')}
LEG={'cel':'V','Cbr':'V','Cbn':'V','Cre':'V','Ppa':'V','Bma':'III','Sra':'IV','Tsp':'I',
     'Nam':'V','Hco':'V'}
DEPTH={'I':1,'III':2,'IV':3,'V':4}
def clade_of(t):
    m=meta.get(t)
    if m: return m['clade'].split()[0]
    return LEG.get(t.split('_',1)[0])
GENES=['dop-1','dop-2','dop-3','dop-4','dop-5','dop-6']
rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)

say('='*100); say('DUPLICATION HISTORY OF THE C. ELEGANS DOPAMINE-RECEPTOR-LIKE FAMILY'); say('='*100)
for trim in TRIMS:
    if not os.path.exists(f'{T}/trees/{trim}.iqtree'):
        say(f'\n--- {trim}: not finished, skipped'); continue
    t=Phylo.read(f'{T}/trees/{trim}.treefile','newick')
    og=[x for x in t.get_terminals() if x.name in
        ('hCHRM1','hCHRM2','hCHRM3','hCHRM5','cel_gar-1','cel_gar-2','cel_gar-3')]
    t.root_with_outgroup(*og)
    tips=[x.name for x in t.get_terminals()]
    say(f'\n{"="*100}'); say(f'--- {trim}'); say('='*100)
    # 1. each gene's own orthologue clade and how deep it reaches
    say('\n1. HOW FAR BACK DOES EACH GENE GO? (deepest lineage inside its orthologue clade)')
    say(f'   {"gene":<8} {"clade":>6} {"support":>11}  {"lineages inside":<18} earliest lineage')
    anc={}
    for g in GENES:
        seed=[n for n in tips if n.replace('_','-').endswith('-'+g)]
        if not seed: say(f'   {g:<8} not on tree'); continue
        # walk outward until a second C. elegans dop gene appears -> that is the orthologue clade
        best=None
        for c in t.find_clades():
            tl=[x.name for x in c.get_terminals()]
            if not set(seed)<=set(tl): continue
            others={o for o in GENES if o!=g and any(x.replace('_','-').endswith('-'+o) for x in tl)}
            if others: continue
            if best is None or len(tl)>len(best[1]): best=(c,tl)
        if not best: say(f'   {g:<8} no exclusive clade'); continue
        c,tl=best
        lin=sorted({clade_of(x) for x in tl if clade_of(x)},key=lambda z:DEPTH.get(z,9))
        anc[g]=set(lin)
        say(f'   {g:<8} {len(tl):>6} {str(c.name):>11}  {", ".join(lin):<18} {lin[0] if lin else "-"}')
    # 2. pairwise duplication dating
    say('\n2. WHEN DID EACH PAIR SPLIT?')
    say('   RECONCILIATION BOUND = deepest lineage anywhere in the duplication node\'s subtree.')
    say('   If a clade-III sequence sits below the node, both daughter lineages existed before')
    say('   clade III diverged, so the duplication predates that split.')
    say('   CONSERVATIVE BOUND = deepest lineage present on BOTH sides independently, which')
    say('   survives one misplaced sequence but under-dates whenever sampling is uneven.')
    say('')
    say(f'   {"pair":<16} {"MRCA":>6} {"support":>11}  {"reconciliation":<22} conservative')
    for i,a in enumerate(GENES):
        for b in GENES[i+1:]:
            sa=[n for n in tips if n.replace('_','-').endswith('-'+a)]
            sb=[n for n in tips if n.replace('_','-').endswith('-'+b)]
            if not sa or not sb: continue
            m=t.common_ancestor(sa+sb)
            sub={clade_of(x.name) for x in m.get_terminals()}
            sub={x for x in sub if x}
            rec=sorted(sub,key=lambda z:DEPTH.get(z,9))
            recv=f'predates {rec[0]} split' if rec else '-'
            shared=sorted(anc.get(a,set())&anc.get(b,set()),key=lambda z:DEPTH.get(z,9))
            cons=f'predates {shared[0]} split' if shared else 'unbounded'
            say(f'   {a+"/"+b:<16} {len(m.get_terminals()):>6} {str(m.name):>11}  {recv:<22} {cons}')
say('')
say('='*100); say('READING'); say('='*100)
say('"Predates the X split" means sequences from lineage X sit inside BOTH daughter clades, so the')
say('duplication must be older than the divergence of X. It is a LOWER bound: absence of a lineage')
say('inside a clade is not evidence the gene was absent there, because these proteomes are of')
say('uneven quality and a real orthologue can simply not have been called.')
open(f'{T}/family_history.txt','w').write('\n'.join(rep)+'\n')
say(f'\n-> {T}/family_history.txt')
