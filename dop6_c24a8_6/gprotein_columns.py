#!/usr/bin/env python3
"""gprotein_columns.py — map the Galpha interface into alignment columns, for every taxon.

The ligand-pocket equivalent of this already exists (pocket_positions_final.py). This does the
same operation on the INTRACELLULAR side, which is where coupling determinants actually live:
predicting a receptor's G protein from its ligand pocket only works to the extent ligand and
coupling covary, which is exactly the assumption under test.

For each of the 11 receptor-Galpha complexes: identify the receptor chain by sequence, find every
receptor residue within 4.5 A of the Galpha chain, align that receptor's SEQUENCE AS OBSERVED IN
THE STRUCTURE to the panel alignment, and carry the contact residues into alignment columns.
Columns supported by >= 6 of 11 structures are the consensus interface.

Writes output/tree/gprotein_columns.tsv
"""
import os, re, sys, gzip, urllib.request
from collections import defaultdict
HERE=os.path.dirname(os.path.abspath(__file__)); T=f'{HERE}/output/tree'
WORK=('/private/tmp/claude-501/-Users-mikeodonnell-git-Malaiwong-Oglu-2026/'
      '93015e40-e331-495a-8362-20ac7ea07c75/scratchpad/pdb')
ALN=f'{T}/aln/panel_domains.linsi.faa'
AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
     'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
     'TYR':'Y','VAL':'V','MSE':'M'}
# pdb -> (receptor chain, galpha chain) taken from output/gprotein_interface_multi.tsv
import csv
MULTI={r['pdb']:r for r in csv.DictReader(open(f'{HERE}/output/gprotein_interface_multi.tsv'),delimiter='\t')}

GA_REF={'Gs':'P63092','Gi1':'P63096','Go':'P09471','Gq':'P50148','G11':'P29992','G12':'Q03113'}
def up(acc):
    f=f'{WORK}/{acc}.fa'
    if not os.path.exists(f):
        with urllib.request.urlopen(f'https://rest.uniprot.org/uniprotkb/{acc}.fasta',timeout=60) as r:
            open(f,'w').write(r.read().decode())
    return ''.join(l.strip() for l in open(f) if not l.startswith('>'))
from Bio import Align
al=Align.PairwiseAligner(); al.open_gap_score=-11; al.extend_gap_score=-1
from Bio.Align import substitution_matrices
al.substitution_matrix=substitution_matrices.load('BLOSUM62')
def pid_local(a,b):
    """identical to gprotein_interface_multi.py's - local alignment, identity over aligned pairs"""
    al.mode='local'
    x=al.align(a,b)[0]; A_,B_=str(x[0]),str(x[1])
    both=[(p_,q_) for p_,q_ in zip(A_,B_) if p_!='-' and q_!='-']
    if not both: return 0,0
    return 100*sum(1 for p_,q_ in both if p_==q_)/len(both), len(both)
GAS={k:up(v) for k,v in GA_REF.items()}

rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)
A={}
for r in open(ALN).read().split('>')[1:]:
    h,s=r.split('\n',1); A[h.strip()]=s.replace('\n','')
W=len(next(iter(A.values())))
say(f'alignment {len(A)} x {W}')

def getpdb(p):
    f=f'{WORK}/{p}.pdb'
    if not os.path.exists(f):
        os.makedirs(WORK,exist_ok=True)
        with urllib.request.urlopen(f'https://files.rcsb.org/download/{p}.pdb',timeout=60) as r:
            open(f,'wb').write(r.read())
    return f

def parse(f):
    ch=defaultdict(dict); at=defaultdict(list)
    for l in open(f):
        if l[:6] not in ('ATOM  ','HETATM'): continue
        rn=l[17:20].strip()
        if rn not in AA3: continue
        c=l[21]; ri=l[22:27].strip()
        ch[c][ri]=AA3[rn]
        at[(c,ri)].append((float(l[30:38]),float(l[38:46]),float(l[46:54])))
    return ch,at

def align_to(seq, target):
    """ungapped structure sequence -> column index per residue, by simple global alignment
       against the target's ungapped alignment row"""
    tg=[(i,c) for i,c in enumerate(target) if c not in '-.']
    t=''.join(c for _,c in tg)
    # Needleman-Wunsch, linear gap
    n,m=len(seq),len(t)
    if n==0 or m==0: return {}
    import numpy as np
    S=np.zeros((n+1,m+1)); P=np.zeros((n+1,m+1),dtype=int)
    g=-4
    S[:,0]=np.arange(n+1)*g; S[0,:]=np.arange(m+1)*g
    P[:,0]=1; P[0,:]=2; P[0,0]=0
    for i in range(1,n+1):
        for j in range(1,m+1):
            d=S[i-1,j-1]+(4 if seq[i-1]==t[j-1] else -1)
            u=S[i-1,j]+g; l=S[i,j-1]+g
            b=max(d,u,l); S[i,j]=b
            P[i,j]=0 if b==d else (1 if b==u else 2)
    i,j=n,m; out={}
    while i>0 or j>0:
        p=P[i,j]
        if p==0: out[i-1]=tg[j-1][0]; i-=1; j-=1
        elif p==1: i-=1
        else: j-=1
    return out

col_hits=defaultdict(set); res_at={}
for pdb,row in MULTI.items():
    rc=row['receptor_chain']; gc=None
    f=getpdb(pdb); ch,at=parse(f)
    # galpha chain = the one whose sequence best matches the recorded class; fall back: largest
    # non-receptor chain. We only need contacts, and the tsv already vetted the pairing.
    # Galpha chain identified BY SEQUENCE against the UniProt references, as the original
    # interface script does. Picking "largest non-receptor chain" selects Gbeta about half the
    # time - Gbeta is 340 aa against Galpha's ~350 - and silently produces a near-empty contact
    # set, which is how the first version of this script returned 1 consensus column.
    gc=None; gbest=0; gcls=None
    for c in ch:
        if c==rc: continue
        order_c=sorted(ch[c],key=lambda x:(int(re.sub(r'\D','',x) or 0),x))
        cs=''.join(ch[c][r] for r in order_c)
        for cls,ref in GAS.items():
            pid_,n_=pid_local(cs,ref)
            if n_>=150 and pid_>=70 and pid_>gbest: gbest=pid_; gc=c; gcls=cls
    if gc is None: say(f'  {pdb}: no Galpha chain at >=70% over >=150 aa - skipped'); continue
    gpts=[p for ri in ch[gc] for p in at[(gc,ri)]]
    import numpy as np
    G=np.array(gpts)
    contact=[]
    for ri in ch[rc]:
        P2=np.array(at[(rc,ri)])
        d=np.sqrt(((P2[:,None,:]-G[None,:,:])**2).sum(-1)).min()
        if d<=4.5: contact.append(ri)
    # structure sequence in order
    order=sorted(ch[rc],key=lambda x:(int(re.sub(r'\D','',x) or 0),x))
    sseq=''.join(ch[rc][r] for r in order)
    # which panel tip is this receptor? use the tsv title/heuristic: best identity match
    best=None;bs=-1
    for k,v in A.items():
        u=v.replace('-','').replace('.','')
        if not u: continue
        # crude: count shared 6-mers
        s6={u[i:i+6] for i in range(len(u)-5)}
        q6={sseq[i:i+6] for i in range(len(sseq)-5)}
        sc=len(s6&q6)
        if sc>bs: bs=sc; best=k
    if bs<5: say(f'  {pdb}: receptor not in panel (best {best}, {bs} shared 6-mers) - skipped'); continue
    mp=align_to(sseq,A[best])
    n=0
    for ri in contact:
        idx=order.index(ri)
        if idx in mp: col_hits[mp[idx]].add(pdb); n+=1
    say(f'  {pdb}  receptor={best:<10} class={gcls:<5} '
        f'{len(contact)} contacts, {n} mapped')

cons=sorted(c for c,s in col_hits.items() if len(s)>=6)
say(f'\nconsensus Galpha-interface columns (contact in >=6 of {len(MULTI)}): {len(cons)}')
with open(f'{T}/gprotein_columns.tsv','w') as f:
    f.write('taxon\tinterface\n')
    for k in sorted(A):
        f.write(k+'\t'+''.join(A[k][c] for c in cons)+'\n')
say(f'-> {T}/gprotein_columns.tsv')
open(f'{T}/gprotein_columns.txt','w').write('\n'.join(rep)+'\n')
