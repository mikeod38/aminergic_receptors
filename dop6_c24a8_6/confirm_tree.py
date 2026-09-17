#!/usr/bin/env python3
"""confirm_tree.py — the checks that decide whether the rebuilt tree's DOP-6 placement is usable.

Three things, run over every trim level so no conclusion rests on one trimming choice:

  1. PLACEMENT STABILITY - the composition of the nested clades containing cel_dop-6, with support.
     A placement that changes between trim levels is not a placement.
  2. TRAIT SIGNAL - Fitch parsimony vs 10,000 label permutations for ligand and for Galpha coupling,
     on the receptors where each is experimentally known. This says whether tree position carries
     information about either trait at all, which has to be established before any inference is
     drawn from DOP-6's position.
  3. LOCAL COUPLING HOMOGENEITY - the Gi/o prediction for DOP-6 rests on its immediate clade being
     coupling-uniform, NOT on a genome-wide Gi/o block (there is no such block: the smallest clade
     containing all Gs receptors holds 125 of 163 tips). So the claim is stated and tested locally:
     of the characterised receptors within N nodes of DOP-6, what fraction share one coupling?

Writes output/tree/confirm.txt
"""
import os, random, re, sys
from collections import Counter

HERE=os.path.dirname(os.path.abspath(__file__)); T=os.environ.get('DOP6_TREEDIR') or os.path.join(HERE,'output','tree')
TRIMS=['gt0.1','gt0.2','gt0.3','gappyout']; NPERM=10000
random.seed(0)
rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)

LIG={'cel_dop-1':'DA','cel_dop-2':'DA','cel_dop-3':'DA','cel_dop-4':'DA',
 'hDRD1':'DA','hDRD2':'DA','hDRD3':'DA','hDRD4':'DA','hDRD5':'DA','dmDop1R1':'DA','dmDop2R':'DA',
 'cel_ser-1':'5HT','cel_ser-4':'5HT','cel_ser-7':'5HT',
 'cel_ser-3':'OA','cel_ser-6':'OA','cel_octr-1':'OA','dmOAMB':'OA','dmOctA2R':'OA',
 'dmOctB1R':'OA','dmOctB2R':'OA','dmOctB3R':'OA','bmOAR1':'OA','amOA1':'OA','tcOA1':'OA','apOA':'OA',
 'cel_ser-2':'TA','cel_tyra-2':'TA','dmTyrR':'TA','dmTyrRII':'TA','amTAR1':'TA','bmTAR':'TA',
 'hCHRM1':'ACh','hCHRM2':'ACh','hCHRM3':'ACh','hCHRM5':'ACh',
 'hHRH1':'HA','hHRH2':'HA','hHRH3':'HA','hHRH4':'HA'}
for k in ['hHTR1A','hHTR1B','hHTR1D','hHTR1E','hHTR1F','hHTR2A','hHTR2B','hHTR2C','hHTR4',
          'hHTR5A','hHTR6','hHTR7']: LIG[k]='5HT'
for k in ['hADRA1A','hADRA1B','hADRA1D','hADRA2A','hADRA2B','hADRA2C','hADRB1','hADRB2','hADRB3']:
    LIG[k]='NE'
COUP={'cel_ser-1':'Gq','cel_ser-4':'Gi/o','cel_ser-7':'Gs','cel_dop-1':'Gq','cel_dop-2':'Gi/o',
 'cel_dop-3':'Gi/o','cel_dop-4':'Gs','cel_octr-1':'Gi/o','cel_ser-3':'Gq','cel_ser-6':'Gs',
 'cel_ser-2':'Gi/o','cel_tyra-2':'Gi/o','cel_tyra-3':'Gq',
 'hDRD1':'Gs','hDRD5':'Gs','hDRD2':'Gi/o','hDRD3':'Gi/o','hDRD4':'Gi/o',
 'hADRB1':'Gs','hADRB2':'Gs','hADRB3':'Gs','hADRA1A':'Gq','hADRA1B':'Gq','hADRA1D':'Gq',
 'hADRA2A':'Gi/o','hADRA2B':'Gi/o','hADRA2C':'Gi/o',
 'hHTR1A':'Gi/o','hHTR1B':'Gi/o','hHTR1D':'Gi/o','hHTR1E':'Gi/o','hHTR1F':'Gi/o',
 'hHTR2A':'Gq','hHTR2B':'Gq','hHTR2C':'Gq','hHTR4':'Gs','hHTR6':'Gs','hHTR7':'Gs',
 'hCHRM1':'Gq','hCHRM3':'Gq','hCHRM5':'Gq','hCHRM2':'Gi/o',
 'dmDop1R1':'Gs','dmDop2R':'Gi/o','dmOctB1R':'Gs','dmOctB2R':'Gs','dmOctB3R':'Gs',
 'dmOAMB':'Gq','dmOctA2R':'Gi/o','dmOctTyrR':'Gi/o','dmTyrR':'Gi/o','dmTyrRII':'Gi/o'}

def parse(nw):
    nw=nw.strip().rstrip(';'); pos=[0]
    def node():
        if nw[pos[0]]=='(':
            pos[0]+=1; k=[]
            while True:
                k.append(node())
                if nw[pos[0]]==',': pos[0]+=1; continue
                if nw[pos[0]]==')': pos[0]+=1; break
            m=re.match(r'[^,\)\(:]*',nw[pos[0]:]); lab=m.group(0); pos[0]+=len(lab)
            if pos[0]<len(nw) and nw[pos[0]]==':':
                m2=re.match(r':[-\d.eE+]*',nw[pos[0]:]); pos[0]+=len(m2.group(0))
            return {'tip':False,'lab':lab,'kids':k}
        m=re.match(r'[^,\)\(:]*',nw[pos[0]:]); lab=m.group(0); pos[0]+=len(lab)
        if pos[0]<len(nw) and nw[pos[0]]==':':
            m2=re.match(r':[-\d.eE+]*',nw[pos[0]:]); pos[0]+=len(m2.group(0))
        return {'tip':True,'lab':lab,'kids':[]}
    return node()
def tips(n): return [n['lab']] if n['tip'] else [x for k in n['kids'] for x in tips(k)]
def fitch(n,st):
    if n['tip']:
        v=st.get(n['lab']); return ({v} if v else None),0
    sets=[];ch=0
    for k in n['kids']:
        s,c=fitch(k,st); ch+=c
        if s: sets.append(s)
    if not sets: return None,ch
    i=set.intersection(*sets)
    return (i,ch) if i else (set.union(*sets),ch+1)

say('='*94); say('CONFIRMATION CHECKS ON THE REBUILT TREE'); say('='*94)

# ── 1. placement stability ───────────────────────────────────────────────────
say('\n1. IS DOP-6\'s PLACEMENT STABLE ACROSS TRIM LEVELS?')
comp={}
for trim in TRIMS:
    tf=f'{T}/trees/{trim}.treefile'
    if not os.path.exists(tf): say(f'   {trim}: missing'); continue
    t=parse(open(tf).read()); path=[]
    def find(n,anc):
        if n['tip']:
            if n['lab']=='cel_dop-6': path.extend(anc)
            return
        for k in n['kids']: find(k,anc+[n])
    find(t,[])
    rows=[]
    for n in reversed(path):
        tl=tips(n)
        if len(tl)>30: break
        m=re.match(r'([\d.]+)/([\d.]+)',n['lab'] or '')
        rows.append((len(tl),f'{float(m.group(2)):.0f}' if m else '-',frozenset(tl)))
    comp[trim]=rows
    cel=lambda S: sorted(x[4:] for x in S if x.startswith('cel_'))
    say(f'   {trim}:')
    for sz,sup,S in rows[:6]:
        say(f'      {sz:>3} tips  UFboot {sup:>4}   C. elegans: {", ".join(cel(S)) or "-"}')
ref=comp.get('gt0.1')
if ref and len(comp)>1:
    say('\n   agreement with gt0.1 on clade composition:')
    for trim,rows in comp.items():
        if trim=='gt0.1': continue
        shared=sum(1 for _,_,S in rows if any(S==S2 for _,_,S2 in ref))
        say(f'      {trim:<10} {shared}/{len(rows)} of its dop-6 clades are identical to a gt0.1 clade')

# ── 2. trait signal ──────────────────────────────────────────────────────────
say('\n2. DOES TREE POSITION CARRY INFORMATION ABOUT LIGAND / COUPLING?')
for name,tbl in (('ligand',LIG),('Galpha coupling',COUP)):
    say(f'   --- {name}')
    for trim in TRIMS:
        tf=f'{T}/trees/{trim}.treefile'
        if not os.path.exists(tf): continue
        t=parse(open(tf).read()); tp=set(tips(t))
        st={k:v for k,v in tbl.items() if k in tp}
        _,obs=fitch(t,st)
        keys=list(st); vals=list(st.values()); null=[]
        for _ in range(NPERM):
            random.shuffle(vals); null.append(fitch(t,dict(zip(keys,vals)))[1])
        p=(sum(1 for x in null if x<=obs)+1)/(NPERM+1)
        say(f'      {trim:<10} {len(st):>3} labelled tips  observed {obs:>3} changes  '
            f'random {sum(null)/len(null):>5.1f}  p={p:.4f}  '
            f'{"clustered" if p<0.05 else "NO SIGNAL"}')

# ── 3. local coupling homogeneity around DOP-6 ───────────────────────────────
say('\n3. IS DOP-6\'s IMMEDIATE CLADE COUPLING-UNIFORM?')
say('   (the Gi/o prediction is a LOCAL claim - there is no genome-wide Gi/o block)')
for trim in TRIMS:
    tf=f'{T}/trees/{trim}.treefile'
    if not os.path.exists(tf): continue
    t=parse(open(tf).read()); path=[]
    def find2(n,anc):
        if n['tip']:
            if n['lab']=='cel_dop-6': path.extend(anc)
            return
        for k in n['kids']: find2(k,anc+[n])
    find2(t,[])
    for n in reversed(path):
        tl=tips(n)
        if len(tl)>30: break
        known=[COUP[x] for x in tl if x in COUP]
        if len(known)>=3:
            c=Counter(known); top,cnt=c.most_common(1)[0]
            say(f'   {trim:<10} smallest clade with >=3 characterised members: {len(tl)} tips, '
                f'{len(known)} characterised, {cnt}/{len(known)} are {top}')
            break
open(f'{T}/confirm.txt','w').write('\n'.join(rep)+'\n')
say('\n-> output/tree/confirm.txt')
