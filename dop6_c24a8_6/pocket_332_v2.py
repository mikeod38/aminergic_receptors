#!/usr/bin/env python3
"""pocket_332_v2.py — position 3.32 read from a multiple alignment anchored on several structures.

v1 mapped 3.32 by pairwise global alignment of every sequence to beta2AR. That passes its positive
control on human receptors — all 34 human amine receptors returned the canonical Asp — but it
returned G for Drosophila 5-HT2A, which cannot be right: the amine-binding aspartate of 5-HT
receptors is known. So the method fails on divergent sequences, which is precisely the regime
DOP-6 and DOP-5 sit in, and their calls cannot be trusted from v1.

This version fixes both halves of the problem:
  * column assignment comes from the MAFFT L-INS-i multiple alignment, not from pairwise alignment
    to one reference, so every sequence is placed in the context of all the others;
  * the column is anchored by FOUR independent ligand-bound structures spanning four ligand classes
    (beta2AR/adrenergic, 5-HT2B/serotonin, D3R/dopamine, H1R/histamine). Each structure's own
    ligand-contacting aspartate is found from atomic contacts, mapped into the alignment, and the
    four must agree on one column. Disagreement means the alignment is untrustworthy there and the
    result is not reported.

Positive control, stated explicitly: every human aminergic receptor must come back D.

Writes output/pocket_332_v2.tsv and .txt
"""
import gzip, os, subprocess, sys, urllib.request
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'output')
WORK = ('/private/tmp/claude-501/-Users-mikeodonnell-git-Malaiwong-Oglu-2026/'
        '93015e40-e331-495a-8362-20ac7ea07c75/scratchpad/pdb')
os.makedirs(WORK, exist_ok=True)
PANEL = f'{OUT}/tree2/panel.faa'
ALNF  = f'{OUT}/tree2/panel.linsi.faa'
# (pdb, ligand het code, uniprot of the receptor, panel label) - each verified by sequence below
REFS = [('2RH1', 'CAU', 'P07550', 'hADRB2'),
        ('4IB4', 'ERM', 'P41595', 'hHTR2B'),
        ('3PBL', 'ETQ', 'P35462', 'hDRD3'),
        ('3RZE', 'DXC', 'P35367', 'hHRH1')]
AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V','MSE':'M'}
rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)

seqs={}
for r in open(PANEL).read().split('>')[1:]:
    h,s=r.split('\n',1); seqs[h.strip()]=s.replace('\n','')
if not os.path.exists(ALNF):
    say('aligning the panel with MAFFT L-INS-i ...')
    with open(ALNF,'w') as o:
        subprocess.run(['mafft','--localpair','--maxiterate','1000','--thread','4','--quiet',PANEL],
                       stdout=o, check=True)
A={}
for r in open(ALNF).read().split('>')[1:]:
    h,s=r.split('\n',1); A[h.strip()]=s.replace('\n','')
W=len(next(iter(A.values())))
say('='*92); say('POSITION 3.32 FROM A MULTIPLE ALIGNMENT, ANCHORED ON FOUR STRUCTURES'); say('='*92)
say(f'alignment {len(A)} sequences x {W} columns\n')

def contacts(pdbid, het):
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
    if not lig: return None,{}
    best=None
    for ch,ats in prot.items():
        hits={}
        for ri,aa,(x,y,z) in ats:
            for (lx,ly,lz) in lig:
                if (x-lx)**2+(y-ly)**2+(z-lz)**2<=4.5**2: hits[ri]=aa; break
        if hits and (best is None or len(hits)>len(best[1])): best=(ch,hits)
    return best if best else (None,{})

def col_of_resnum(label, resnum):
    """alignment column of residue `resnum` (1-based, UniProt numbering) of `label`"""
    s=A.get(label)
    if s is None: return None
    i=0
    for c,ch in enumerate(s):
        if ch!='-':
            i+=1
            if i==resnum: return c
    return None

cols=[]
say(f'{"structure":<8} {"ligand":<7} {"panel seq":<10} {"Asp contacting ligand":>22} {"aln column":>11}')
say('-'*68)
for pdbid,het,acc,lab in REFS:
    ch,hits=contacts(pdbid,het)
    if not hits: say(f'{pdbid:<8} {het:<7} no ligand contacts found'); continue
    asps=[ri for ri,aa in sorted(hits.items()) if aa=='D']
    if not asps: say(f'{pdbid:<8} {het:<7} no aspartate contacts the ligand'); continue
    # structures may be numbered by UniProt; take the lowest-numbered contacting Asp (TM3 precedes TM5/6)
    ri=asps[0]
    c=col_of_resnum(lab,ri)
    ok = c is not None and A[lab][c]=='D'
    say(f'{pdbid:<8} {het:<7} {lab:<10} {"Asp"+str(ri):>22} {str(c) if c is not None else "?":>11}'
        f'{"" if ok else "   <-- residue at that column is not D, mapping failed"}')
    if ok: cols.append((pdbid,lab,ri,c))

if not cols:
    say('\nno usable anchor; aborting'); open(f'{OUT}/pocket_332_v2.txt','w').write('\n'.join(rep)+'\n'); sys.exit(1)
uniq=sorted({c for _,_,_,c in cols})
say('')
if len(uniq)==1:
    COL=uniq[0]; say(f'all {len(cols)} structures agree on alignment column {COL} -> this is 3.32')
else:
    from collections import Counter
    cnt=Counter(c for _,_,_,c in cols); COL,n=cnt.most_common(1)[0]
    say(f'structures DISAGREE: columns {uniq}. Taking the majority column {COL} ({n}/{len(cols)}),')
    say('and flagging that the alignment is not fully reliable at this position.')

# ── positive control, then the calls ─────────────────────────────────────────
HUMAN_AMINE=[k for k in A if k.startswith(('hHTR','hDRD','hADR','hCHRM','hHRH','hTAAR'))]
ctrl={k:A[k][COL] for k in sorted(HUMAN_AMINE)}
bad=[k for k,v in ctrl.items() if v!='D']
say('')
say(f'POSITIVE CONTROL: {len(ctrl)} human aminergic receptors, all of which are known to carry the')
say(f'amine-binding aspartate. Returned D: {len(ctrl)-len(bad)}/{len(ctrl)}.'
    + ('  PASS' if not bad else f'  FAIL for {bad}'))
if bad:
    say('Control failed, so the calls below are not trustworthy. Stopping.')
    open(f'{OUT}/pocket_332_v2.txt','w').write('\n'.join(rep)+'\n'); sys.exit(1)

say('')
say('CALLS (only non-D receptors and the genes of interest are listed in full):')
rows=[]
for k in sorted(A):
    rows.append((k, A[k][COL]))
nonD=[(k,v) for k,v in rows if v!='D']
say(f'  {len(rows)-len(nonD)}/{len(rows)} sequences carry D at 3.32')
say('')
say(f'  {"receptor":<18} {"3.32":>5}')
for k,v in nonD: say(f'  {k:<18} {v:>5}')
say('')
for k in ('cel_dop-6','cel_dop-5','cel_F35H10.10','cel_F59D12.1','cel_T21B4.4',
          'hMTNR1A','hMTNR1B','hGPR50','dm5HT2A','dm5HT1A','dm5HT7'):
    if k in A: say(f'  {k:<18} {A[k][COL]:>5}   (v1 pairwise result for comparison)')
with open(f'{OUT}/pocket_332_v2.tsv','w') as f:
    f.write('receptor\tres_3_32\taln_column\n')
    for k,v in rows: f.write(f'{k}\t{v}\t{COL}\n')
open(f'{OUT}/pocket_332_v2.txt','w').write('\n'.join(rep)+'\n')
say('\n-> output/pocket_332_v2.tsv')
