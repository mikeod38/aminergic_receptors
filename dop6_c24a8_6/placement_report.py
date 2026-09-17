#!/usr/bin/env python3
"""placement_report.py — read the dating off the 268-taxon placement tree.

Replaces the retracted RBH matrix (see nematode_dating2.py, retracted). No orthology call is made
anywhere: every nematode monoamine-HMM candidate >= 120 from 8 proteomes spanning clades I, III, IV
and V went into the panel, and the TOPOLOGY assigns paralogy. A lineage is "in" a clade if a
sequence from it falls inside, which is what a best-hit table cannot establish.

Writes output/tree/placement/dating.txt
"""
import csv, os, sys
from Bio import Phylo
HERE=os.path.dirname(os.path.abspath(__file__)); T=f'{HERE}/output/tree/placement'
TRIMS=('gt0.1','gappyout','gt0.2','gt0.3')
meta={r['tip']:r for r in csv.DictReader(open(f'{T}/candidates.tsv'),delimiter='\t')}
LEG={'cel':'V','Cbr':'V','Cbn':'V','Cre':'V','Ppa':'V','Bma':'III','Sra':'IV','Tsp':'I',
     'Nam':'V','Hco':'V'}
NEMPRE=tuple(LEG)+('Tspi','Tmur','Bmal','Srat','Hcon','Name','Ppac','Cbri')
def clade_of(t):
    m=meta.get(t)
    if m: return m['clade'].split()[0]
    return LEG.get(t.split('_',1)[0])
def is_nem(t): return clade_of(t) is not None
rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)

say('='*100); say('DATING THE dop-3/dop-5/dop-6 RADIATION BY PLACEMENT (268-taxon panel)'); say('='*100)
for trim in TRIMS:
    tf=f'{T}/trees/{trim}.treefile'; iq=f'{T}/trees/{trim}.iqtree'
    if not os.path.exists(iq):
        say(f'\n--- {trim}: not finished (no .iqtree) - SKIPPED, its treefile is not final'); continue
    t=Phylo.read(tf,'newick')
    og=[x for x in t.get_terminals() if x.name in
        ('hCHRM1','hCHRM2','hCHRM3','hCHRM5','cel_gar-1','cel_gar-2','cel_gar-3')]
    t.root_with_outgroup(*og)
    tips=[x.name for x in t.get_terminals()]
    d6=[n for n in tips if n.replace('_','-').endswith('-dop-6')]
    sizes={}
    for c in t.find_clades():
        tl=[x.name for x in c.get_terminals()]
        if set(d6)<=set(tl): sizes[len(tl)]=(c,tl)
    say(f'\n--- {trim}')
    say(f'  {"tips":>5} {"SH/UFboot":>11}  {"nematode clades inside":<26} first outside member')
    deepest=None
    for n in sorted(sizes):
        if n>60: break
        c,tl=sizes[n]
        cl=sorted({clade_of(x) for x in tl if clade_of(x)})
        out=[x for x in tl if not is_nem(x)]
        if not out: deepest=(n,c,tl,cl)
        say(f'  {n:>5} {str(c.name):>11}  {", ".join(cl):<26} {sorted(out)[0] if out else "-"}')
    if deepest:
        n,c,tl,cl = deepest
        say(f'\n  deepest clade with no non-nematode member: {n} tips, support {c.name}')
        say(f'  lineages represented: {", ".join(cl)}')
    # the dop-6 orthologue clade specifically
    for n in sorted(sizes):
        c,tl=sizes[n]
        cl={clade_of(x) for x in tl if clade_of(x)}
        if len(tl)>=6 and 'IV' in cl:
            say(f'\n  SHALLOWEST clade containing dop-6 AND a clade-IV sequence: {n} tips, support {c.name}')
            for x in sorted(tl):
                m=meta.get(x)
                say(f'     {x:<24} {m["species"] if m else "Compara/panel tip"}')
            break
say('')
say('='*100)
say('READING')
say('='*100)
say('A sequence from a lineage falling INSIDE a clade dates that clade to before the lineage split.')
say('This is the claim a reciprocal-best-hit table cannot make, and the reason the previous dating')
say('matrix was retracted: RBH sends a real orthologue to whichever paralogue scores highest, so a')
say('gene that is present is recorded as absent.')
open(f'{T}/dating.txt','w').write('\n'.join(rep)+'\n')
say(f'\n-> {T}/dating.txt')
