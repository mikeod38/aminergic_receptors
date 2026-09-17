#!/usr/bin/env python3
"""pocket_positions_final.py — ligand-pocket residues for every tip of the final tree.

Defines the pocket from atomic contacts in three ligand-bound aminergic structures spanning three
ligand classes, maps each contact into the final 7TM alignment, keeps the positions that at least
two of the three structures agree on, and reads every taxon's residue there. The amine-binding
aspartate (3.32) is reported separately because it is the single diagnostic position; the rest are
reported as a pocket string so the figure can show what each receptor's pocket actually looks like.

Positive control: every human aminergic receptor must return D at 3.32.

Writes output/tree/pocket.tsv
"""
import gzip, os, sys, urllib.request
from collections import defaultdict

HERE=os.path.dirname(os.path.abspath(__file__)); T=os.path.join(HERE,'output','tree')
WORK=('/private/tmp/claude-501/-Users-mikeodonnell-git-Malaiwong-Oglu-2026/'
      '93015e40-e331-495a-8362-20ac7ea07c75/scratchpad/pdb')
ALN=f'{T}/aln/panel_domains.linsi.faa'
REFS=[('2RH1','CAU','hADRB2'),('4IB4','ERM','hHTR2B'),('3PBL','ETQ','hDRD3')]
AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
     'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
     'TYR':'Y','VAL':'V','MSE':'M'}
rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)
A={}
for r in open(ALN).read().split('>')[1:]:
    h,s=r.split('\n',1); A[h.strip()]=s.replace('\n','')
W=len(next(iter(A.values())))
say(f'alignment {len(A)} x {W}')

def contacts(pdbid,het):
    p=f'{WORK}/{pdbid}.pdb'
    if not os.path.exists(p):
        with urllib.request.urlopen(f'https://files.rcsb.org/download/{pdbid}.pdb.gz',timeout=120) as r:
            open(p,'w').write(gzip.decompress(r.read()).decode('utf-8','replace'))
    prot=defaultdict(list); lig=[]
    for ln in open(p):
        if ln.startswith(('ATOM','HETATM')):
            rn=ln[17:20].strip(); ch=ln[21]
            try: ri=int(ln[22:26])
            except ValueError: continue
            xyz=(float(ln[30:38]),float(ln[38:46]),float(ln[46:54]))
            if rn==het: lig.append(xyz)
            elif rn in AA3: prot[ch].append((ri,AA3[rn],xyz))
    best=None
    for ch,ats in prot.items():
        hits={}
        for ri,aa,(x,y,z) in ats:
            for (lx,ly,lz) in lig:
                if (x-lx)**2+(y-ly)**2+(z-lz)**2<=4.5**2: hits[ri]=aa; break
        if hits and (best is None or len(hits)>len(best)): best=hits
    return best or {}

def col_of(label,resnum):
    """alignment column of UniProt residue `resnum` of `label` - the excerpt is a CONCATENATION of
    blocks, so the mapping must go through the excerpt's own residue order, not raw numbering."""
    ex=A[label].replace('-','')
    full=FULL[label]
    i=full.find(ex[:25])
    if i<0: return None
    # walk the excerpt, tracking position in the full sequence across concatenation joins
    return None
FULL={}
for r in open(f'{T}/panel.faa').read().split('>')[1:]:
    h,s=r.split('\n',1); FULL[h.strip()]=s.replace('\n','')

# map via pairwise alignment of the full reference sequence to its own excerpt
from Bio.Align import PairwiseAligner, substitution_matrices
pa=PairwiseAligner(); pa.substitution_matrix=substitution_matrices.load('BLOSUM62')
pa.open_gap_score,pa.extend_gap_score,pa.mode=-11,-1,'global'
def resnum_to_col(label,resnum):
    ex=A[label]; exs=ex.replace('-','')
    aln=pa.align(FULL[label],exs)[0]; X,Y=str(aln[0]),str(aln[1])
    fi=ei=0; target=None
    for a,b in zip(X,Y):
        if a!='-': fi+=1
        if b!='-': ei+=1
        if a!='-' and b!='-' and fi==resnum: target=ei; break
    if target is None: return None
    c=0
    for idx,ch in enumerate(ex):
        if ch!='-':
            c+=1
            if c==target: return idx
    return None

votes=defaultdict(list)
say('')
for pdbid,het,lab in REFS:
    hits=contacts(pdbid,het)
    if not hits or lab not in A: say(f'{pdbid}: skipped'); continue
    n=0
    for ri,aa in sorted(hits.items()):
        c=resnum_to_col(lab,ri)
        if c is not None and A[lab][c]==aa: votes[c].append((pdbid,ri,aa)); n+=1
    say(f'{pdbid} ({lab}): {len(hits)} ligand contacts, {n} mapped into the alignment')
cons=sorted(c for c,v in votes.items() if len({x[0] for x in v})>=2)
say(f'\nconsensus pocket columns (contact in >=2 of 3 structures): {len(cons)}')
# which one is 3.32? the column where hADRB2 has D113
c332=resnum_to_col('hADRB2',113)
say(f'3.32 = alignment column {c332}; hADRB2 residue there = {A["hADRB2"][c332]}')
HUM=[k for k in A if k.startswith(('hHTR','hDRD','hADR','hCHRM','hHRH','hTAAR'))]
bad=[k for k in HUM if A[k][c332]!='D']
say(f'POSITIVE CONTROL: {len(HUM)-len(bad)}/{len(HUM)} human aminergic receptors have D at that '
    f'column.  {"PASS" if not bad else "FAIL "+str(bad)}')
if bad: say('control failed; not writing pocket calls'); open(f'{T}/pocket.txt','w').write('\n'.join(rep)); sys.exit(1)
with open(f'{T}/pocket.tsv','w') as f:
    f.write('taxon\tres_3_32\tpocket\n')
    for k in sorted(A):
        pk=''.join(A[k][c] if A[k][c]!='-' else '.' for c in cons)
        f.write(f'{k}\t{A[k][c332]}\t{pk}\n')
nonD=[k for k in sorted(A) if A[k][c332]!='D']
say(f'\n{len(A)-len(nonD)}/{len(A)} tips carry D at 3.32; without it:')
for k in nonD: say(f'   {k:<28} {A[k][c332]}')
open(f'{T}/pocket.txt','w').write('\n'.join(rep)+'\n')
say(f'\n-> output/tree/pocket.tsv')
