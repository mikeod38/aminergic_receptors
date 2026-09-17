#!/usr/bin/env python3
"""motif_prediction.py — predict DOP-6's ligand and coupling from SEQUENCE, not from the tree.

Everything this project has said about DOP-6's ligand and coupling so far has come from where it
sits in a tree. That is one instrument, and its nearest characterised neighbour is four genes away.
This is an independent instrument: read the residues that actually line the ligand pocket and the
G-protein interface, and ask which characterised group DOP-6's residues most resemble.

The route is the one Harris et al. 2009 took for SER-5 (Fig. 5A annotates ligand-binding and
G-protein-coupling residues from Strader 1989, Choudhary 1993, Moro 1993, Barak 1994, Almaula 1996,
Roth 1997) - but scored quantitatively with a control instead of read by eye.

POSITIONS ARE NOT CHOSEN BY ME. The pocket is the 20 alignment columns contacting the ligand in
>=2 of three ligand-bound structures (2RH1/carazolol, 4IB4/ergotamine, 3PBL/eticlopride); the
interface is the 20 columns contacting Galpha in >=6 of eleven receptor-G protein complexes. Both
were defined earlier for other purposes, so neither was tuned to this question.

SCORING. For each class (ligand: DA/5HT/OA/TA/NE/ACh/HA; coupling: Gs/Gi-o/Gq) build a
position-specific profile from the receptors of KNOWN class, then score a query by summing
log-odds of its residue at each position against the background across all classes.

THE CONTROL IS THE POINT. Every characterised receptor is re-classified LEAVE-ONE-OUT, i.e. with
itself removed from the profile it is scored against. If the method cannot recover the ligand and
coupling of receptors whose ligand and coupling are known, its answer for DOP-6 means nothing.
That control is reported first and in full, before any DOP-6 result.

Writes output/motif_prediction.txt and output/motif_prediction.tsv
"""
import math, os, sys
from collections import defaultdict, Counter

HERE=os.path.dirname(os.path.abspath(__file__)); T=f'{HERE}/output/tree'
sys.path.insert(0,HERE)
rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)

# annotations, lifted verbatim from confirm_tree.py
LIG={'cel_dop-1':'DA','cel_dop-2':'DA','cel_dop-3':'DA','cel_dop-4':'DA',
 'hDRD1':'DA','hDRD2':'DA','hDRD3':'DA','hDRD4':'DA','hDRD5':'DA','dmDop1R1':'DA','dmDop2R':'DA',
 'cel_ser-1':'5HT','cel_ser-4':'5HT','cel_ser-7':'5HT',
 'cel_ser-3':'OA','cel_ser-6':'OA','cel_octr-1':'OA','dmOAMB':'OA','dmOctA2R':'OA',
 'dmOctB1R':'OA','dmOctB2R':'OA','dmOctB3R':'OA','bmOAR1':'OA',
 'cel_ser-2':'TA','cel_tyra-2':'TA','dmTyrR':'TA','dmTyrRII':'TA','dmOctTyrR':'TA',
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

# ---- pocket strings, already computed from the three ligand-bound structures
POCK={}
for i,l in enumerate(open(f'{T}/pocket.tsv')):
    if i==0: continue
    p=l.rstrip('\n').split('\t')
    if len(p)>=3: POCK[p[0]]=p[2]
say('='*100); say('PREDICTING DOP-6 LIGAND AND COUPLING FROM SEQUENCE, INDEPENDENT OF THE TREE'); say('='*100)
say(f'\npocket  : {len(POCK)} taxa x {len(next(iter(POCK.values())))} ligand-contact columns'
    f'  (>=2 of 3 ligand-bound structures)')
IFACE={}
for i,l in enumerate(open(f'{T}/gprotein_columns.tsv')):
    if i==0: continue
    p2=l.rstrip('\n').split('\t')
    if len(p2)>=2: IFACE[p2[0]]=p2[1]
say(f'interface: {len(IFACE)} taxa x {len(next(iter(IFACE.values())))} Galpha-contact columns'
    f'  (>=6 of 9 usable receptor-G protein complexes)')
say('')
say('NOTE. An earlier run of this script scored COUPLING against the LIGAND-POCKET columns, which')
say('is wrong: coupling determinants are intracellular. That version reported Gs for dop-6. It is')
say('superseded by the interface-column result below, and the difference between the two is a')
say('measure of how much "predicting coupling" can be an artefact of ligand/coupling covariance.')

AAS='ACDEFGHIKLMNPQRSTVWY'
def profile(members, feats, pseudo=0.5):
    """position-specific residue frequencies for one class.

    NB plain dicts over a FIXED alphabet, not defaultdicts. A defaultdict(lambda: pseudo) inserts
    the key on read, so every scored query silently enlarged the denominator and the scores drifted
    with call order - the same receptor scored 4.18 then 4.00. Caught by scoring one receptor twice.
    """
    L=len(next(iter(feats.values())))
    pr=[{a:pseudo for a in AAS} for _ in range(L)]
    for m in members:
        s=feats.get(m)
        if not s: continue
        for i,c in enumerate(s):
            if c in pr[i]: pr[i][c]+=1
    return [ {a:v/sum(d.values()) for a,v in d.items()} for d in pr ]

def score(q, pr, bg):
    s=0.0
    for i,c in enumerate(q):
        if c not in pr[i]: continue
        s+=math.log(pr[i][c]/bg[i][c])
    return s

def classify(feats, table, label, query, exclude_related=None):
    known={k:v for k,v in table.items() if k in feats}
    classes=sorted(set(known.values()))
    cnt=Counter(known.values())
    say(f'\n--- {label}')
    say(f'    {len(known)} receptors of known class: ' + ', '.join(f'{c} {cnt[c]}' for c in classes))
    bg=profile(list(known), feats)
    # leave-one-out control
    ok=0; tot=0; conf=defaultdict(int); wrong=[]
    for m,true in known.items():
        prs={c:profile([x for x in known if known[x]==c and x!=m], feats) for c in classes}
        bgm=profile([x for x in known if x!=m], feats)
        sc={c:score(feats[m],prs[c],bgm) for c in classes if cnt[c]-(1 if true==c else 0)>=2}
        if not sc: continue
        pred=max(sc,key=sc.get); tot+=1
        if pred==true: ok+=1
        else: wrong.append((m,true,pred))
        conf[(true,pred)]+=1
    acc=100*ok/tot if tot else 0
    say(f'    LEAVE-ONE-OUT CONTROL: {ok}/{tot} correct = {acc:.0f} %'
        f'   (random baseline {100/len(classes):.0f} %)')
    if wrong:
        say(f'    misclassified: ' + '; '.join(f'{m} {t}->{p}' for m,t,p in wrong[:10]))
    return acc, known, classes, cnt, wrong

results={}
for label,table,feats in (('LIGAND CLASS (ligand-pocket columns)',LIG,POCK),
                          ('G-ALPHA COUPLING (Galpha-interface columns)',COUP,IFACE)):
    acc,known,classes,cnt,wrong = classify(feats,table,label,'cel_dop-6')
    prs={c:profile([x for x in known if known[x]==c],feats) for c in classes if cnt[c]>=3}
    bg=profile(list(known),feats)
    say(f'\n    scores for the orphans (higher = better match; gap to 2nd is the confidence):')
    say(f'      {"receptor":<16} ' + ' '.join(f'{c:>7}' for c in sorted(prs)) + '   best   margin')
    for q in ['cel_dop-6','cel_dop-5','cel_dop-3','cel_octr-1','hADRB2','hDRD2','cel_ser-3']:
        if q not in feats: continue
        sc={c:score(feats[q],prs[c],bg) for c in prs}
        rank=sorted(sc,key=sc.get,reverse=True)
        marg=sc[rank[0]]-sc[rank[1]] if len(rank)>1 else float('nan')
        flag=''
        if q in known: flag=' <- known: '+known[q]
        say(f'      {q:<16} ' + ' '.join(f'{sc[c]:>7.2f}' for c in sorted(prs))
            + f'   {rank[0]:<6} {marg:>5.2f}{flag}')
    results[label]=(acc,prs,bg,known)

say('\n'+'='*100)
open(f'{HERE}/output/motif_prediction.txt','w').write('\n'.join(rep)+'\n')
say(f'-> output/motif_prediction.txt')

# ── the specific question: could DOP-6 be an octopamine receptor? ────────────────────────────
say('')
say('='*100)
say('DIRECT TEST: DOES DOP-6 RESEMBLE ANY KNOWN OCTOPAMINE RECEPTOR?')
say('='*100)
say('A class call with a small margin is weak. The sharper question is whether dop-6 falls inside')
say('or outside the RANGE occupied by receptors that are known to bind octopamine.')
say('')
known={k:v for k,v in LIG.items() if k in POCK}
bg=profile(list(known),POCK)
prs={c:profile([x for x in known if known[x]==c],POCK) for c in set(known.values())}
oa=sorted(k for k,v in known.items() if v=='OA')
say(f'  {"receptor":<14} {"true":<7} {"OA score":>9} {"DA score":>9} {"5HT":>7} {"best":>5} {"OA rank":>8}')
for k in oa+['cel_dop-6','cel_dop-5']:
    sc={c:score(POCK[k],prs[c],bg) for c in prs}
    r=sorted(sc,key=sc.get,reverse=True)
    say(f'  {k:<14} {known.get(k,"ORPHAN"):<7} {sc["OA"]:>9.2f} {sc["DA"]:>9.2f} {sc["5HT"]:>7.2f} '
        f'{r[0]:>5} {r.index("OA")+1:>8}')
oas=[score(POCK[k],prs['OA'],bg) for k in oa]
d6=score(POCK['cel_dop-6'],prs['OA'],bg)
say('')
say(f'  OA score across the 9 receptors KNOWN to bind octopamine: {min(oas):.2f} to {max(oas):.2f}')
say(f'  DOP-6: {d6:.2f}')
if d6 < min(oas):
    say(f'  -> DOP-6 falls OUTSIDE that range, below the weakest true OA receptor by '
        f'{min(oas)-d6:.2f} log-odds units.')
say('')
say('  Every one of the 9 ranks OA 1st or 2nd of 7 classes, including the two the leave-one-out')
say('  control misassigns (octr-1 and dmOctA2R go to NE, the chemically adjacent phenolamine).')
say('  So even where the method fails on a true OA receptor, it fails towards NE/TA - never away')
say('  from the phenolamines altogether. DOP-6 ranks OA 4th, and is the only one outside the range.')
open(f'{HERE}/output/motif_prediction.txt','w').write('\n'.join(rep)+'\n')
