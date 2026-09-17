#!/usr/bin/env python3
"""rebuild_boundary120.py — rebuild the panel on a stated, uniform inclusion rule.

The previous panel's boundary was provenance, not evidence: the C. elegans members were whatever
Koelle 2018 listed plus the named aminergic receptors, which let in ador-1 (monoamine HMM 85.2,
rank 26 - below six neuropeptide receptors) and pcdr-1 (52.8, rank 56) while excluding 36 genes that
score higher. Both then behaved as rogues, which is a symptom of the boundary rather than an
argument for deleting them individually.

New rule, applied uniformly:

  IN   monoamine HMM score >= 120. The largest gap in the whole ranking is 27.0 points, between
       gar-1 (124.4) and npr-35 (97.4); above it are exactly the 19 aminergic + muscarinic genes,
       below it the neuropeptide receptors begin and descend without further structure.
  IN   orthologues of those genes in other nematodes - for topology, not as evidence
  IN   reference receptors of known ligand/coupling (human, Drosophila, other metazoan)
  IN   a TEST SET of the 12 highest-scoring excluded genes (79-97), labelled as such and included
       once, to verify that the 120 boundary is real rather than assumed
  OUT  ador-1 / T21B4.4 (85.2, below boundary; its own name says "homolog")
  OUT  pcdr-1 / F59D12.1 (52.8, below boundary; named for a phenotype)
  OUT  F35H10.10 - class C receptor (7tm_3 141.6 + ANF/VFT 57.9, no class-A signal at all)

Writes output/tree/panel_b120.faa, basis_b120.tsv
"""
import os, re
HERE=os.path.dirname(os.path.abspath(__file__)); T=os.path.join(HERE,'output','tree')
SC='/Users/mikeodonnell/git/WormPooledCloning/data/processed/gpcr/cel_gpcr_monoamine_scores.tsv'
PEP=('/private/tmp/claude-501/-Users-mikeodonnell-git-Malaiwong-Oglu-2026/'
     '93015e40-e331-495a-8362-20ac7ea07c75/scratchpad/cel_pep.fa')
CUT=120.0; N_TEST=12
DROP_CEL={'T21B4.4','F59D12.1','F35H10.10'}          # by sequence name, as they appear in the panel
DROP_WB={'WBGene00011878','WBGene00010329','WBGene00018073'}
rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)

score={}; fam={}
for l in open(SC).read().strip().split('\n')[1:]:
    p=l.split('\t')
    try: s=float(p[4])
    except (ValueError,IndexError): continue
    score[p[2]]=s; fam[p[2]]=p[3]; score[p[1]]=s; fam[p[1]]=p[3]   # by gene name AND WBGene

old={}
for r in open(f'{T}/panel.faa').read().split('>')[1:]:
    h,s=r.split('\n',1); old[h.strip()]=s.replace('\n','')
NEM=('Cbr_','Cbn_','Cre_','Ppa_','Asu_','Bma_','Sra_','Hco_','Tsp_','Nam_')
keep={}; basis={}
say('='*94); say('PANEL REBUILD, BOUNDARY AT MONOAMINE HMM >= 120'); say('='*94)
dropped=[]
for k,v in old.items():
    if k.startswith('cel_'):
        g=k[4:]
        if g in DROP_CEL:
            dropped.append((k,score.get(g,float('nan')),'below boundary / wrong class')); continue
        s=score.get(g)
        if s is None:
            dropped.append((k,float('nan'),'no HMM score')); continue
        if s<CUT:
            dropped.append((k,s,'below boundary')); continue
        keep[k]=v; basis[k]=f'ingroup (HMM {s:.1f})'
    elif k.startswith(NEM):
        g=k.split('_',1)[1].replace('_','-')
        s=score.get(g) or score.get(g.replace('-','_'))
        if g in DROP_CEL or (s is not None and s<CUT):
            dropped.append((k,s if s else float('nan'),'orthologue of an excluded gene')); continue
        keep[k]=v; basis[k]='nematode orthologue (topology only)'
    else:
        keep[k]=v; basis[k]='reference receptor'
say(f'\nDROPPED ({len(dropped)}):')
for k,s,why in sorted(dropped,key=lambda x:-(x[1] if x[1]==x[1] else -1)):
    say(f'   {k:<24} HMM {s if s==s else float("nan"):>6.1f}  {why}')

# test set: the 12 highest-scoring genes still excluded, matched by WBGene so ador-1 cannot duplicate
have_wb=set()
for l in open(SC).read().strip().split('\n')[1:]:
    p=l.split('\t')
    if f'cel_{p[2]}' in keep: have_wb.add(p[1])
rows=[]
for l in open(SC).read().strip().split('\n')[1:]:
    p=l.split('\t')
    try: s=float(p[4])
    except (ValueError,IndexError): continue
    if p[1] in have_wb or p[1] in DROP_WB: continue
    rows.append((p[1],p[2],p[3],s))
rows.sort(key=lambda x:-x[3]); test=rows[:N_TEST]
byg={}
for r in open(PEP).read().split('>')[1:]:
    h,s=r.split('\n',1); seq=s.replace('\n','')
    m=re.search(r'gene:(\S+)',h)
    if not m: continue
    g=m.group(1)
    if g not in byg or len(seq)>len(byg[g]): byg[g]=seq
say(f'\nTEST SET ({N_TEST} highest-scoring excluded genes, to verify the boundary):')
for wb,g,f_,s in test:
    if wb not in byg: say(f'   {g}: no protein'); continue
    keep[f'TEST_{g}']=byg[wb]; basis[f'TEST_{g}']=f'boundary test (HMM {s:.1f}, labelled {f_})'
    say(f'   TEST_{g:<12} HMM {s:>6.1f}  {f_}')
with open(f'{T}/panel_b120.faa','w') as f:
    for k in sorted(keep):
        f.write(f'>{k}\n')
        for i in range(0,len(keep[k]),60): f.write(keep[k][i:i+60]+'\n')
with open(f'{T}/basis_b120.tsv','w') as f:
    f.write('taxon\tbasis\tlength\n')
    for k in sorted(keep): f.write(f'{k}\t{basis[k]}\t{len(keep[k])}\n')
from collections import Counter
c=Counter(b.split(' (')[0] for b in basis.values())
say(f'\n{len(keep)} sequences -> output/tree/panel_b120.faa')
for k,v in sorted(c.items()): say(f'   {k:<32} {v}')
open(f'{T}/boundary_rebuild.txt','w').write('\n'.join(rep)+'\n')
