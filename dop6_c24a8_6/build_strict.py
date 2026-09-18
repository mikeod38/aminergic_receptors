#!/usr/bin/env python3
"""build_strict.py — panels containing ONLY sequences that pass the monoamine HMM >= 120 boundary.

MOD: "The final tree should include only those that passed our cutoffs."

Scoring every tip of the existing panels found 15 below the boundary:
  12 TEST_* neuropeptide receptors (78.9-97.5) - these are the DECLARED boundary-test set and are
     supposed to be below it; they demonstrate the boundary works, and that result is reported
     separately. They do not belong in a panel presented as the aminergic family.
  hMTNR1A 102.5, hGPR50 98.3, hMTNR1B 95.3 - the melatonin receptors, carried as ordinary
     reference members. These were added because DOP-5 has G at 3.32, the melatonin-receptor
     state, but that rationale was never declared as an exception - so they sat below the same
     cutoff that excluded F59D12.1 at 52.8. That is the curation asymmetry the referee panel
     flagged under charge B3, and this script removes it.

Everything else in both panels passes. Builds strict versions of b120 (173 -> 158) and placement
(268 -> 253) and re-runs IQ-TREE at the four trim levels.
"""
import os, re, subprocess, sys
import pyhmmer
HERE=os.path.dirname(os.path.abspath(__file__)); T=f'{HERE}/output/tree'
MONO='/Users/mikeodonnell/git/WormPooledCloning/data/raw/MonoamineHMM_Model.hmm'
IQ='/Users/mikeodonnell/git/receptor_phylogeny/scripts/bin/iqtree3-mt'
CUT=120.0; THREADS='3'
OUTG=['hCHRM1','hCHRM2','hCHRM3','hCHRM5','cel_gar-1','cel_gar-2','cel_gar-3']
def readfa(p):
    d={};k=None
    for l in open(p):
        if l.startswith('>'): k=l[1:].strip(); d[k]=''
        elif k: d[k]+=l.strip()
    return d
with pyhmmer.plan7.HMMFile(MONO) as fh: hmm=next(fh)
alpha=hmm.alphabet
def score_all(seqs):
    digs=[pyhmmer.easel.TextSequence(name=k.encode(),
           sequence=re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]','X',v)).digitize(alpha)
          for k,v in seqs.items()]
    best={}
    for hits in pyhmmer.hmmsearch([hmm],digs,E=10.0,cpus=6):
        for h in hits:
            n=h.name.decode() if isinstance(h.name,bytes) else h.name
            best[n]=max(best.get(n,-99.0),h.score)
    return best
for tag,src,dst in [('b120',f'{T}/panel_b120.faa',f'{T}/strict_b120'),
                    ('placement',f'{T}/placement/panel.faa',f'{T}/strict_placement')]:
    os.makedirs(f'{dst}/trees',exist_ok=True)
    seqs=readfa(src); sc=score_all(seqs)
    keep={k:v for k,v in seqs.items() if sc.get(k,-99)>=CUT}
    drop=sorted(((sc.get(k,-99),k) for k in seqs if sc.get(k,-99)<CUT))
    print(f'\n=== {tag}: {len(seqs)} -> {len(keep)} tips ({len(drop)} dropped)')
    for s,k in drop: print(f'    dropped {s:7.1f}  {k}')
    with open(f'{dst}/panel.faa','w') as f:
        for k,v in keep.items(): f.write(f'>{k}\n{v}\n')
    open(f'{dst}/dropped.tsv','w').write('tip\thmm_score\n'+
        ''.join(f'{k}\t{s:.1f}\n' for s,k in drop))
print('\npanels written; run build_strict_trees.sh to build')
