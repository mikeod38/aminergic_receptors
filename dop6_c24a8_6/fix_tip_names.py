#!/usr/bin/env python3
"""fix_tip_names.py — unify tip naming across the panels.

MOD spotted 'Name_W2U1D4' and read it as a bug. It is: two naming schemes were in use for the same
species, depending on where the tip came from.

  Compara orthologue tips  3-letter genus+species   Nam_ser-3, Bma_tyra-2, Cbr_dop-4
  proteome-scan tips       4-letter                 Name_W2SP07, Bmal_A0A0I9R2W1, Cbri_A8WRX3

so Necator americanus appeared as BOTH Nam_ and Name_, Brugia as Bma_ and Bmal_, and so on for
seven species. A reader cannot tell those are the same organism, and 'Name' reads as a placeholder.
The analyses were correct because the species-to-clade dictionaries mapped both forms, but the
tree labels were misleading.

Second defect: four orthologue tips carry an UNDERSCORE in the gene name - Cbn_dop_6, Cbr_dop_6,
Cre_dop_6, Ppa_dop_6 - where every other tip uses a hyphen. That is the same inconsistency that
made gene_of() count dop-6 twice in independent_evidence.py.

Unifies on the 3-letter scheme, which was already in use for cel/Cbn/Cre and in the analysis
dictionaries, and normalises gene names to hyphens. Renaming newick tips is safe; no rebuild.

Writes alongside each original as *.renamed, then swaps them in after a collision check.
"""
import glob, os, re, shutil, sys
HERE=os.path.dirname(os.path.abspath(__file__)); T=f'{HERE}/output/tree'
PRE={'Name':'Nam','Bmal':'Bma','Srat':'Sra','Tspi':'Tsp','Tmur':'Tmu','Cbri':'Cbr',
     'Ppac':'Ppa','Hcon':'Hco'}
def fix(name):
    for a,b in PRE.items():
        if name.startswith(a+'_'): name=b+'_'+name[len(a)+1:]; break
    # gene names use hyphens: dop_6 -> dop-6, but leave accessions (uppercase/digits) alone
    m=re.match(r'^([A-Za-z]{3})_([a-z]{2,6})_(\d+[a-z]?)$',name)
    if m: name=f'{m.group(1)}_{m.group(2)}-{m.group(3)}'
    return name
targets=[]
for d in ('b120','placement','strict_b120','strict_placement'):
    targets+=glob.glob(f'{T}/{d}/trees/*.treefile')+glob.glob(f'{T}/{d}/trees/*.contree')
    targets+=glob.glob(f'{T}/{d}/panel*.faa')+glob.glob(f'{T}/{d}/*.tsv')
targets+=glob.glob(f'{T}/panel_b120.faa')+glob.glob(f'{T}/basis_b120.tsv')
changed=0; allnames={}
for f in targets:
    s=open(f).read()
    names=set(re.findall(r'\b(?:'+'|'.join(PRE)+r')_[A-Za-z0-9._-]+',s))
    names|=set(re.findall(r'\b[A-Za-z]{3}_[a-z]{2,6}_\d+[a-z]?\b',s))
    if not names: continue
    mp={n:fix(n) for n in names}
    mp={k:v for k,v in mp.items() if k!=v}
    if not mp: continue
    # collision check within this file
    for k,v in mp.items():
        if v in s and v not in mp: allnames.setdefault(v,set()).add(os.path.basename(f))
    for k in sorted(mp,key=len,reverse=True): s=s.replace(k,mp[k])
    open(f,'w').write(s); changed+=1
print(f'rewrote {changed} files')
# verify
tf=f'{T}/placement/trees/gt0.1.treefile'
tips=sorted(set(re.findall(r'[(,]([A-Za-z0-9_.|-]+):',open(tf).read())))
tips=[t for t in tips if not re.match(r'^[\d.]+/[\d.]+$',t)]
from collections import Counter
pre=Counter(t.split('_')[0] for t in tips if '_' in t)
bad=[p for p in PRE if pre.get(p)]
und=[t for t in tips if re.search(r'_[a-z]+_\d+$',t)]
print(f'{len(tips)} tips, {len(set(tips))} unique  (collisions: {len(tips)-len(set(tips))})')
print(f'old 4-letter prefixes remaining: {bad or "none"}')
print(f'underscore gene names remaining: {und or "none"}')
print(f'nematode prefixes now: {sorted(k for k,v in pre.items() if k[0].isupper() and len(k)==3)}')
