#!/usr/bin/env python3
"""assemble_final_panel.py — the panel the tree is finally built from, with an auditable basis
recorded for every taxon.

Membership is not one HMM cutoff (the score bands of wanted outgroups and unwanted families
interleave). Each sequence carries one of three explicit bases instead:

  anchor     verified identity and known ligand/coupling - the human, Drosophila and C. elegans
             reference receptors, each accession checked against its own UniProt record
  ortholog   Compara ortholog of a characterised C. elegans receptor in another nematode, filtered
             on length and alignment coverage
  scan       recovered by scanning a downloaded reference proteome with the monoamine HMM, i.e.
             found by sequence rather than by anyone's annotation text
  subject    the uncharacterised genes the analysis exists to place, included by definition

Deep/outgroup lineages are sampled more sparsely than the ingroup, and the model organisms already
covered by the anchor set are not re-sampled from their own proteomes.

Writes output/tree/panel.faa and basis.tsv
"""
import os, re
from collections import defaultdict

HERE=os.path.dirname(os.path.abspath(__file__))
OUT=os.path.join(HERE,'output','tree'); os.makedirs(OUT,exist_ok=True)
T2=f'{HERE}/output/tree2'; T3=f'{HERE}/output/tree3'; SCAN=f'{HERE}/output/proteome_scan'
SUBJECT={'cel_dop-6','cel_dop-5','cel_F35H10.10','cel_F59D12.1','cel_T21B4.4'}
SKIP_SCAN={'Hsap','Dmel','Cele','Mmus'}          # already covered by the verified anchor set
PER={'Drer':2,'Ajap':2,'Pmar':2,'Bflo':2,'Cint':2,'Skow':2,'Spur':2,'Apla':2,'Pmin':2,
     'Cvir':3,'Myes':3,'Lana':3,'Lpol':3,'Pcau':3,
     'Nvec':3,'Hvul':3,'Ebur':3,'Tadh':3,'Mlei':3,'Aque':3}
def rd(p):
    d={}
    if not os.path.exists(p): return d
    for r in open(p).read().split('>')[1:]:
        h,s=r.split('\n',1); d[h.strip()]=s.replace('\n','')
    return d
seqs={}; basis={}
# anchors + subjects from the corrected v2 panel
v2=rd(f'{T2}/panel.faa')
for k,v in v2.items():
    seqs[k]=v; basis[k]='subject' if k in SUBJECT else 'anchor'
# nematode orthologs from the breadth round
v3=rd(f'{T3}/panel.faa')
NEM=('Cbr_','Cbn_','Cre_','Ppa_','Asu_','Bma_','Sra_','Hco_','Tsp_','Nam_')
for k,v in v3.items():
    if k in seqs: continue
    if k.startswith(NEM): seqs[k]=v; basis[k]='ortholog'
# proteome scan breadth
sc=rd(f'{SCAN}/candidates.faa')
bysp=defaultdict(list)
score={}
for ln in open(f'{SCAN}/hits.tsv').read().strip().split('\n')[1:]:
    p=ln.split('\t'); score[(p[0],p[1])]=float(p[2])
for k,v in sc.items():
    sp=k.split('_')[0]
    if sp in SKIP_SCAN: continue
    bysp[sp].append((k,v))
for sp,rows in bysp.items():
    rows.sort(key=lambda kv: -score.get((sp,kv[0][len(sp)+1:].replace('_','|')), 0))
    for k,v in rows[:PER.get(sp,2)]:
        if k in seqs: continue
        seqs[k]=v; basis[k]='scan'
with open(f'{OUT}/panel.faa','w') as f:
    for k in sorted(seqs):
        f.write(f'>{k}\n')
        for i in range(0,len(seqs[k]),60): f.write(seqs[k][i:i+60]+'\n')
with open(f'{OUT}/basis.tsv','w') as f:
    f.write('taxon\tinclusion_basis\tlength\n')
    for k in sorted(seqs): f.write(f'{k}\t{basis[k]}\t{len(seqs[k])}\n')
from collections import Counter
c=Counter(basis.values())
print(f'{len(seqs)} sequences -> output/tree/panel.faa')
for k,v in sorted(c.items()): print(f'   {k:<10} {v}')
grp=Counter('nematode' if k.startswith(NEM) else 'human' if k.startswith('h') else
            'C. elegans' if k.startswith('cel_') else 'Drosophila' if k.startswith('dm') else
            'other metazoan' for k in seqs)
print()
for k,v in sorted(grp.items()): print(f'   {k:<16} {v}')
