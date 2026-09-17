#!/usr/bin/env python3
"""build_placement.py — panel = b120 + every nematode monoamine-HMM candidate >= 120, deduplicated.

The point is to let the TOPOLOGY assign paralogy across the phylum, instead of asking a
best-hit heuristic that cannot resolve a recent radiation. Compara orthologue tips are kept as
named anchors (they carry gene names the scan hits do not), but any scan hit identical to one of
them is dropped so no protein enters twice.
"""
import os, re, subprocess, sys
HERE=os.path.dirname(os.path.abspath(__file__)); T=f'{HERE}/output/tree'
P=f'{T}/placement'; os.makedirs(f'{P}/trees',exist_ok=True)

def readfa(p):
    d={}; k=None
    for l in open(p):
        if l.startswith('>'): k=l[1:].strip(); d[k]=''
        elif k: d[k]+=l.strip()
    return d

base=readfa(f'{T}/panel_b120.faa'); cand=readfa(f'{P}/candidates.faa')
byseq={v:k for k,v in base.items()}
added=0; dropped=[]
for k,v in cand.items():
    if v in byseq: dropped.append((k,byseq[v])); continue
    base[k]=v; added+=1
print(f'b120 {len(base)-added} tips + {added} new scan candidates = {len(base)}')
print(f'{len(dropped)} scan hits dropped as exact duplicates of an existing tip')
for a,b in dropped[:12]: print(f'   {a:<22} == {b}')
with open(f'{P}/panel.faa','w') as f:
    for k,v in base.items(): f.write(f'>{k}\n{v}\n')
open(f'{P}/dedup.txt','w').write('\n'.join(f'{a}\t{b}' for a,b in dropped)+'\n')

# 7TM excerpts, then align/trim/tree with the same settings as b120
sys.argv=['extract_domains.py']
os.environ['DOP6_PANEL']=f'{P}/panel.faa'; os.environ['DOP6_OUT']=f'{P}/panel_domains.faa'
print('\nextracting 7TM excerpts ...')
r=subprocess.run([sys.executable,f'{HERE}/extract_domains.py'],capture_output=True,text=True,
                 env={**os.environ})
print(r.stdout[-1500:] or r.stderr[-1500:])
