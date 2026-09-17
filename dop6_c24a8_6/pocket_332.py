#!/usr/bin/env python3
"""pocket_332.py — what sits at position 3.32 in DOP-6 and the other orphans?

Position 3.32 (Ballesteros-Weinstein) is the aspartate in transmembrane helix 3 that forms the salt
bridge with the protonated amine of every biogenic-amine ligand. It is the single most diagnostic
residue for "is this an amine receptor at all". Melatonin receptors, whose ligand is an N-acetylated
indoleamine with no free amine, carry a small residue there instead — so the identity of this one
position separates amine receptors from the melatonin branch.

Rather than borrowing a residue number from a paper (which breaks across species and across
structures with unmodelled loops), 3.32 is located here from a real ligand-bound structure: the
beta2-adrenergic receptor bound to its ligand, where the TM3 aspartate contacting the ligand is
identified from atomic contacts, then carried to every panel sequence through a single alignment.

Writes output/pocket_332.tsv and .txt
"""
import gzip, os, re, subprocess, urllib.request
from collections import defaultdict
from Bio.Align import PairwiseAligner, substitution_matrices

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'output')
WORK = ('/private/tmp/claude-501/-Users-mikeodonnell-git-Malaiwong-Oglu-2026/'
        '93015e40-e331-495a-8362-20ac7ea07c75/scratchpad/pdb')
os.makedirs(WORK, exist_ok=True)
PDB = '2RH1'          # beta2AR with carazolol bound
LIGCODE = 'CAU'       # carazolol
AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V','MSE':'M'}
rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)

path=f'{WORK}/{PDB}.pdb'
if not os.path.exists(path):
    with urllib.request.urlopen(f'https://files.rcsb.org/download/{PDB}.pdb.gz',timeout=120) as r:
        open(path,'w').write(gzip.decompress(r.read()).decode('utf-8','replace'))
prot=defaultdict(list); lig=[]
for ln in open(path):
    if ln.startswith(('ATOM','HETATM')):
        rn=ln[17:20].strip(); ch=ln[21]
        try: ri=int(ln[22:26])
        except ValueError: continue
        xyz=(float(ln[30:38]),float(ln[38:46]),float(ln[46:54]))
        if rn==LIGCODE: lig.append(xyz)
        elif rn in AA3: prot[ch].append((ri,AA3[rn],xyz))
say('='*86); say('POSITION 3.32 — THE AMINE-BINDING ASPARTATE'); say('='*86)
say(f'located from {PDB} (beta2AR + carazolol): receptor residues within 4.5 A of the ligand\n')
best=None
for ch,ats in prot.items():
    hits={}
    for ri,aa,(x,y,z) in ats:
        for (lx,ly,lz) in lig:
            if (x-lx)**2+(y-ly)**2+(z-lz)**2 <= 4.5**2: hits[ri]=aa; break
    if hits and (best is None or len(hits)>len(best[1])): best=(ch,hits)
ch,hits=best
asps=[ri for ri,aa in sorted(hits.items()) if aa=='D']
say(f'chain {ch}: {len(hits)} residues contact the ligand; aspartates among them: {asps}')
anchor=asps[0]
say(f'-> taking Asp{anchor} as 3.32 (the TM3 aspartate contacting the ligand amine)\n')

b2=''.join(hits[k] for k in sorted(hits))  # not used for alignment; fetch full seq instead
with urllib.request.urlopen('https://rest.uniprot.org/uniprotkb/P07550.fasta',timeout=60) as r:
    B2=''.join(r.read().decode().split('\n')[1:])
# structure numbering for 2RH1 chain A is UniProt numbering for beta2AR
al=PairwiseAligner(); al.substitution_matrix=substitution_matrices.load('BLOSUM62')
al.open_gap_score,al.extend_gap_score,al.mode=-11,-1,'global'
seqs={}
for f in (f'{OUT}/tree2/panel.faa',):
    for r_ in open(f).read().split('>')[1:]:
        h,s=r_.split('\n',1); seqs[h.strip()]=s.replace('\n','')
say(f'{"receptor":<18} {"3.32":>5}  {"ligand class":<14} note')
say('-'*86)
rows=[]
GROUP={'h5HT':'5-HT','hHTR':'5-HT','hDRD':'dopamine','hADR':'adrenergic','hCHRM':'ACh',
       'hHRH':'histamine','hTAAR':'trace amine','hMTNR':'MELATONIN','hGPR50':'melatonin-related'}
for k in sorted(seqs):
    aln=al.align(B2,seqs[k])[0]; A,Bs=str(aln[0]),str(aln[1])
    i=j=0; res='-'
    for x,y in zip(A,Bs):
        if x!='-': i+=1
        if y!='-': j+=1
        if i==anchor and x!='-':
            res=y if y!='-' else '-'; break
    cls=''
    for pre,c in GROUP.items():
        if k.startswith(pre): cls=c; break
    rows.append((k,res,cls))
for k,res,cls in rows:
    flag=''
    if res=='D': flag='amine-binding D present'
    elif res in 'GAST': flag='small residue — melatonin-receptor-like'
    elif res=='-': flag='not alignable at this column'
    else: flag=f'neither D nor small ({res})'
    mark=' <<<' if k in ('cel_dop-6','cel_dop-5') else ''
    say(f'{k:<18} {res:>5}  {cls:<14} {flag}{mark}')
with open(f'{OUT}/pocket_332.tsv','w') as f:
    f.write('receptor\tres_3_32\tligand_class\n')
    for r_ in rows: f.write('\t'.join(r_)+'\n')
open(f'{OUT}/pocket_332.txt','w').write('\n'.join(rep)+'\n')
say('\n-> output/pocket_332.tsv')
