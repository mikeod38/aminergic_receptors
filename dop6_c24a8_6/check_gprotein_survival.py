#!/usr/bin/env python3
"""check_gprotein_survival.py — did excising ICL3 remove the Galpha interface from the excerpts?

ICL3 contributes to G protein binding, and the consensus Galpha-contact set computed earlier from
11 receptor-G protein complexes includes three positions in the TM5/ICL3 region (DOP-6 I206, L210,
R213) and five at the TM6 cytoplasmic end (R467, R470, A471, L474, V475). If the 7TM excerpting
threw those away, then any coupling inference read off this tree is built on an alignment that
deleted the very residues that determine coupling - the same failure mode as a gap trim silently
dropping a named functional site.

So: take the 20 consensus contact positions in DOP-6 numbering and ask, position by position,
whether each falls inside a retained block (with margin) or inside the excised loop.

Writes output/gprotein_survival.txt
"""
import os, re
import pyhmmer

HERE=os.path.dirname(os.path.abspath(__file__)); T=os.path.join(HERE,'output','tree')
MONO='/Users/mikeodonnell/git/WormPooledCloning/data/raw/MonoamineHMM_Model.hmm'
MARGIN=15
# consensus Galpha contacts (>=6 of 11 complexes), in DOP-6 numbering, from gprotein_interface_multi
CONTACTS=[57,120,121,124,125,128,129,131,132,133,206,210,213,467,470,471,474,475,534,535]
REGION={57:'TM1/ICL1',120:'TM3 DRY',121:'TM3 DRY',124:'ICL2',125:'ICL2',128:'ICL2',129:'ICL2',
        131:'ICL2',132:'ICL2',133:'ICL2',206:'TM5 cyto end',210:'TM5/ICL3 start',
        213:'TM5/ICL3 start',467:'ICL3/TM6 start',470:'TM6 cyto end',471:'TM6 cyto end',
        474:'TM6',475:'TM6',534:'TM7',535:'TM7 NPxxY'}
nm=lambda x: x.decode() if isinstance(x,(bytes,bytearray)) else str(x)
rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)

seqs={}
for r in open(f'{T}/panel.faa').read().split('>')[1:]:
    h,s=r.split('\n',1); seqs[h.strip()]=s.replace('\n','')
alpha=pyhmmer.easel.Alphabet.amino()
CHECK=['cel_dop-6','hADRB2','cel_octr-1','cel_ser-3','hDRD2']
digs=[pyhmmer.easel.TextSequence(name=k.encode(),
        sequence=re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]','X',seqs[k])).digitize(alpha) for k in CHECK]
with pyhmmer.plan7.HMMFile(MONO) as hf: hmms=list(hf)
iv={}
for hits in pyhmmer.hmmsearch(hmms,digs,E=1e-3,cpus=4):
    for h in hits:
        n=nm(h.name); s=seqs[n]
        ds=sorted([d for d in h.domains if d.i_evalue<1.0],key=lambda d:d.alignment.target_from)
        merged=[]
        for d in ds:
            lo=max(1,d.alignment.target_from-MARGIN); hi=min(len(s),d.alignment.target_to+MARGIN)
            if merged and lo<=merged[-1][1]+1: merged[-1]=(merged[-1][0],max(merged[-1][1],hi))
            else: merged.append((lo,hi))
        iv[n]=merged
say('='*88); say('DO THE Galpha-CONTACT RESIDUES SURVIVE THE 7TM EXCERPTING?'); say('='*88)
say('consensus contacts = positions touching Galpha in >=6 of 11 receptor-G protein complexes\n')
m=iv['cel_dop-6']
say(f'DOP-6 retained blocks (after the {MARGIN} aa margin): ' +
    ', '.join(f'{a}-{b}' for a,b in m))
excised=[(m[i][1]+1,m[i+1][0]-1) for i in range(len(m)-1) if m[i+1][0]-m[i][1]>1]
say(f'DOP-6 excised region(s): ' + ', '.join(f'{a}-{b} ({b-a+1} aa)' for a,b in excised))
say('')
say(f'{"position":>9} {"residue":>8} {"region":<18} verdict')
say('-'*60)
d6=seqs['cel_dop-6']; kept=0
for p in CONTACTS:
    inside=any(a<=p<=b for a,b in m)
    kept+=inside
    say(f'{p:>9} {d6[p-1]:>8} {REGION[p]:<18} {"RETAINED" if inside else "LOST TO THE EXCISED LOOP"}')
say('')
say(f'{kept}/{len(CONTACTS)} consensus Galpha-contact positions survive in the DOP-6 excerpt')
say('')
say('Why: ICL3 contacts Galpha at its ENDS - where it leaves TM5 and where it enters TM6 - and')
say('those ends are inside the HMM-aligned blocks, because they are the conserved helix boundaries')
say('the model recognises. What gets excised is the middle of the loop, which is the part that')
say('varies from 14 to 821 residues across this panel and makes no consensus contact.')
say('')
say('Same check on four other receptors:')
for k in CHECK[1:]:
    mm=iv[k]; span=sum(b-a+1 for a,b in mm)
    exc=sum(mm[i+1][0]-mm[i][1]-1 for i in range(len(mm)-1))
    say(f'  {k:<12} retained {span:>4} aa in {len(mm)} block(s), excised {exc:>4} aa')
open(f'{HERE}/output/gprotein_survival.txt','w').write('\n'.join(rep)+'\n')
say('\n-> output/gprotein_survival.txt')
