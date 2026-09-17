#!/usr/bin/env python3
"""build_placement_tree.py — same pipeline as build_b120, run on the phylum-wide placement panel."""
import os, re, subprocess, sys
import pyhmmer
HERE=os.path.dirname(os.path.abspath(__file__)); T=os.path.join(HERE,'output','tree')
B=f'{T}/placement'; os.makedirs(B,exist_ok=True); os.makedirs(f'{B}/trees',exist_ok=True)
MONO='/Users/mikeodonnell/git/WormPooledCloning/data/raw/MonoamineHMM_Model.hmm'
TM7='/Users/mikeodonnell/git/WormPooledCloning/data/raw/PF00001.hmm'
IQ='/Users/mikeodonnell/git/receptor_phylogeny/scripts/bin/iqtree3-mt'
MARGIN=15; THREADS='6'
OUTG=['hCHRM1','hCHRM2','hCHRM3','hCHRM5','cel_gar-1','cel_gar-2','cel_gar-3']
nm=lambda x: x.decode() if isinstance(x,(bytes,bytearray)) else str(x)
rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)
seqs={}
for r in open(f'{T}/placement/panel.faa').read().split('>')[1:]:
    h,s=r.split('\n',1); seqs[h.strip()]=s.replace('\n','')
say(f'{len(seqs)} sequences; 7TM excerpts, {MARGIN} aa margin, blocks concatenated (ICL3 excised)')
alpha=pyhmmer.easel.Alphabet.amino()
digs=[pyhmmer.easel.TextSequence(name=k.encode(),
        sequence=re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]','X',v)).digitize(alpha)
      for k,v in sorted(seqs.items())]
def env(path):
    with pyhmmer.plan7.HMMFile(path) as hf: hmms=list(hf)
    out={}
    for hits in pyhmmer.hmmsearch(hmms,digs,E=1e-3,cpus=6):
        for h in hits:
            n=nm(h.name)
            ds=[d for d in h.domains if d.i_evalue<1.0]
            if not ds: continue
            segs=sorted((d.alignment.target_from,d.alignment.target_to) for d in ds)
            tot=sum(b-a+1 for a,b in segs)
            if n not in out or tot>out[n][1]: out[n]=(segs,tot)
    return out
E=env(MONO)
miss=[k for k in seqs if k not in E]
if miss:
    say(f'{len(miss)} without a monoamine envelope; falling back to 7tm_1: {", ".join(miss[:6])}')
    E2=env(TM7)
    for k in miss:
        if k in E2: E[k]=E2[k]
out={}
for k,s in seqs.items():
    if k not in E: say(f'  dropped {k}: no envelope'); continue
    iv=[]
    for a,b in E[k][0]:
        lo=max(1,a-MARGIN); hi=min(len(s),b+MARGIN)
        if iv and lo<=iv[-1][1]+1: iv[-1]=(iv[-1][0],max(iv[-1][1],hi))
        else: iv.append((lo,hi))
    ex=''.join(s[lo-1:hi] for lo,hi in iv)
    if len(ex)<150: say(f'  dropped {k}: excerpt {len(ex)} aa'); continue
    out[k]=ex
L=sorted(len(v) for v in out.values())
say(f'{len(out)} excerpts; min {L[0]}, median {L[len(L)//2]}, max {L[-1]}')
fa=f'{B}/panel_domains.faa'
with open(fa,'w') as f:
    for k in sorted(out):
        f.write(f'>{k}\n')
        for i in range(0,len(out[k]),60): f.write(out[k][i:i+60]+'\n')
aln=f'{B}/panel.linsi.faa'
if not os.path.exists(aln):
    say('MAFFT L-INS-i ...')
    with open(aln,'w') as o:
        subprocess.run(['mafft','--localpair','--maxiterate','1000','--thread',THREADS,'--quiet',fa],
                       stdout=o,check=True)
A={}
for r in open(aln).read().split('>')[1:]:
    h,s=r.split('\n',1); A[h.strip()]=s.replace('\n','')
W=len(next(iter(A.values())))
gap=sum(s.count('-') for s in A.values())/(len(A)*W)
say(f'alignment {len(A)} x {W}, {100*gap:.1f} % gaps')
trims=[]
for name,args in [('gappyout',['-gappyout']),('gt0.1',['-gt','0.1']),
                  ('gt0.2',['-gt','0.2']),('gt0.3',['-gt','0.3'])]:
    o=f'{B}/panel.{name}.faa'
    r=subprocess.run(['trimal','-in',aln,'-out',o]+args,capture_output=True,text=True)
    if r.returncode: continue
    B2={}
    for rr in open(o).read().split('>')[1:]:
        h,s=rr.split('\n',1); B2[h.strip()]=s.replace('\n','')
    say(f'  {name:<10} {len(next(iter(B2.values()))):>5} columns')
    trims.append((name,o))
say('\nIQ-TREE ...')
for name,path in [t for t in trims if t[0]=='gt0.1']+[t for t in trims if t[0]!='gt0.1']:
    pre=f'{B}/trees/{name}'
    if os.path.exists(pre+'.iqtree'): say(f'  {name}: done already'); continue
    og=[o for o in OUTG if o in A]
    r=subprocess.run([IQ,'-s',path,'-m','MFP','-B','1000','-alrt','1000','-T',THREADS,
                      '--prefix',pre,'-o',','.join(og),'--quiet','-redo'],capture_output=True,text=True)
    if r.returncode: say(f'  {name}: FAILED {(r.stderr or r.stdout)[-200:]}')
    else:
        m=''
        for ln in open(pre+'.iqtree'):
            if ln.startswith('Best-fit model'): m=ln.split(':')[1].strip(); break
        say(f'  {name}: done, {m}')
open(f'{B}/build.txt','w').write('\n'.join(rep)+'\n')
say(f'\n-> {B}/trees/')
