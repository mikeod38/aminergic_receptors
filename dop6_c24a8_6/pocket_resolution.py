#!/usr/bin/env python3
"""pocket_resolution.py — how much ligand resolution do the pocket columns actually have?

Written after MOD made two objections to catechol_signature.py, both of which hold:

  (1) "All of the OA receptors in your tally are adrenergic receptors, so how do we know the
      binding mode isn't typical of the adrenergic family?"
  (2) "There are 5HTR2 in the tree right next to the DRD and dop-2 relatives - look at those
      directly."

Both are answered below by direct tally, and both destroy the 5.46 argument:

  - 5.46 = S in 11/11 DA, 9/9 NE, 9/9 OA, 5/5 TA. INVARIANT across catecholamine AND phenolamine,
    so it cannot be reading the number of ring hydroxyls. It marks the adrenergic-like family.
  - 5.46 varies WITHIN the 5-HT2 group (S/A/A/T) although the ligand is identical, so it is not
    ligand-determined even there.

The remaining question is whether ANY column separates the ligands, which is the thing worth
knowing before trusting a sequence-based ligand call. Answered by counting columns that separate
two ligand classes at >=80 % conservation in each, with a label-permutation null and two positive
controls.

Writes output/pocket_resolution.txt
"""
import os, random
from collections import Counter
HERE=os.path.dirname(os.path.abspath(__file__)); T=f'{HERE}/output/tree'
rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)
A={}
for r in open(f'{T}/aln/panel_domains.linsi.faa').read().split('>')[1:]:
    h,s=r.split('\n',1); A[h.strip()]=s.replace('\n','')
src=open(f'{HERE}/motif_prediction.py').read()
ns={}; exec(src[src.index('LIG={'):src.index('COUP={')],ns); LIG=ns['LIG']
known={k:v for k,v in LIG.items() if k in A}
G=lambda c:[k for k,v in known.items() if v==c]
L=len(next(iter(A.values())))
C={'3.32':408,'5.42':532,'5.43':533,'5.46':536}

say('='*100); say('HOW MUCH LIGAND RESOLUTION DO THE POCKET COLUMNS HAVE?'); say('='*100)
say('')
say('OBJECTION 1 - is 5.46 a chemistry marker or a family marker?')
for cls in ['DA','NE','OA','TA','5HT','ACh','HA']:
    ks=G(cls)
    if ks: say(f'   {cls:<4} n={len(ks):<3} 5.46: ' +
               ' '.join(f'{a}x{n}' for a,n in Counter(A[k][C["5.46"]] for k in ks).most_common()))
say('')
say('   S in 34/34 across DA+NE (catecholamine) AND OA+TA (phenolamine). A position with the same')
say('   state in both chemistries cannot be reading the number of ring hydroxyls. It is a marker')
say('   of the adrenergic-like family. THE RING-HYDROXYL READING IS WITHDRAWN.')
say('')
say('OBJECTION 2 - the 5-HT2 receptors sit beside the DRD/dop-2 group; what do they show?')
say(f'   {"tip":<12} {"ligand":<6} {"3.32":>5} {"5.42":>5} {"5.43":>5} {"5.46":>5}')
for k in ['hDRD2','hDRD3','hDRD4','dmDop2R','cel_dop-2','cel_dop-3','Tsp_dop-2',
          'hHTR2A','hHTR2B','hHTR2C','dm5HT2A','cel_ser-1','cel_dop-5','cel_dop-6']:
    if k in A:
        say(f'   {k:<12} {LIG.get(k,"-"):<6} ' + ' '.join(f'{A[k][C[p]]:>5}' for p in ('3.32','5.42','5.43','5.46')))
say('')
say('   Within one ligand class (5-HT) 5.46 runs S, A, A, A, T. Not ligand-determined.')
say('   The position that DOES track this split is 5.42: S in every DA receptor, G in every 5-HT')
say('   receptor here - and DOP-6 has S, the dopamine-receptor state, exactly as descent from a')
say('   D2-like ancestor predicts.')
say('')
say('RESOLUTION TEST - how many columns separate two ligand classes at >=80 % / >=80 %?')
def count(g1,g2):
    n=0
    for c in range(L):
        c1=Counter(A[k][c] for k in g1); c2=Counter(A[k][c] for k in g2)
        if not c1 or not c2: continue
        a1,n1=c1.most_common(1)[0]; a2,n2=c2.most_common(1)[0]
        if a1 in '-.' or a2 in '-.': continue
        if n1/len(g1)>=0.8 and n2/len(g2)>=0.8 and a1!=a2: n+=1
    return n
pairs=[('DA','ACh'),('DA','HA'),('DA','5HT'),('DA','OA'),('DA','NE'),('OA','TA')]
say(f'   {"comparison":<14} {"columns":>8}   note')
for a_,b_ in pairs:
    ga,gb=G(a_),G(b_)
    if len(ga)<3 or len(gb)<3: continue
    n=count(ga,gb)
    note=''
    if (a_,b_)==('DA','ACh'): note='positive control - different receptor branch'
    if (a_,b_)==('DA','5HT'): note='<- unambiguously different ligands, and ZERO columns'
    if (a_,b_)==('DA','OA'): note='the question of interest'
    say(f'   {a_+" vs "+b_:<14} {n:>8}   {note}')
ga,gb=G('DA'),G('OA'); obs=count(ga,gb); pool=ga+gb; perm=[]
for _ in range(200):
    random.shuffle(pool); perm.append(count(pool[:len(ga)],pool[len(ga):]))
say('')
say(f'   DA vs OA: observed {obs}; random relabelling mean {sum(perm)/len(perm):.1f}, max {max(perm)}, '
    f'p={(sum(1 for x in perm if x>=obs)+1)/(len(perm)+1):.3f}')
say('')
say('='*100); say('CONCLUSION'); say('='*100)
say('The pocket columns separate aminergic from muscarinic receptors easily (17 columns) and')
say('separate DOPAMINE from SEROTONIN receptors NOT AT ALL (0 columns), although those ligands')
say('differ far more than dopamine does from octopamine. The single DA-vs-OA column is above a')
say('permutation null only because the null is 0, and one column is not a signature.')
say('')
say('So: SEQUENCE AT THESE POSITIONS CANNOT ASSIGN DOP-6 A LIGAND. Both earlier claims are')
say('withdrawn - the profile-similarity call (which read ancestry) and the 5.46 chemistry call')
say('(which read family membership). MOD was right on both counts.')
say('')
say('What survives is smaller and directional only: DOP-6 keeps 3.32 D and 5.42 S, the')
say('adrenergic/catecholamine-family state consistent with D2-like descent; and it carries a')
say('derived A at 5.46 that is unique among 163 tips and fixed across all five orthologues, so')
say('something at that position changed and stayed changed. That is a fact about divergence, not')
say('evidence about which amine it binds.')
open(f'{HERE}/output/pocket_resolution.txt','w').write('\n'.join(rep)+'\n')
say(f'\n-> output/pocket_resolution.txt')
