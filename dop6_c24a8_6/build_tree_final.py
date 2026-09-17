#!/usr/bin/env python3
"""build_tree_final.py — the tree, built on the audited 164-taxon panel.

Differences from the earlier rounds, each one forced by something that went wrong:
  * every accession re-verified against its own record (three were a squid rhodopsin, a bacterial
    decarboxylase and a ubiquitin ligase)
  * breadth recovered by scanning downloaded proteomes with the monoamine HMM rather than by
    UniProt text query (which returned a kinase and a tetraspanin)
  * nematode depth from Compara orthologs, so DOP-6 is no longer the only worm gene in its region
    of the tree
  * the multithreaded IQ-TREE binary, not the single-core one on PATH
  * four trim levels, and every conclusion has to hold across all of them

Writes output/tree/aln/*, output/tree/trees/*
"""
import os, subprocess, sys
HERE=os.path.dirname(os.path.abspath(__file__)); T=os.path.join(HERE,'output','tree')
ALN=f'{T}/aln'; TRE=f'{T}/trees'
os.makedirs(ALN,exist_ok=True); os.makedirs(TRE,exist_ok=True)
IQ='/Users/mikeodonnell/git/receptor_phylogeny/scripts/bin/iqtree3-mt'
THREADS='6'
OUTGROUP=['hCHRM1','hCHRM2','hCHRM3','hCHRM5','cel_gar-1','cel_gar-2','cel_gar-3']
rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)

v=subprocess.run([IQ,'--version'],capture_output=True,text=True).stdout.split('\n')[0]
say(f'IQ-TREE binary: {v}')
if 'single-core' in v.lower(): say('refusing to run multithreaded on a sequential build'); sys.exit(1)

aln=f'{ALN}/panel_domains.linsi.faa'
if not os.path.exists(aln):
    say('MAFFT L-INS-i on the 7TM domain excerpts ...')
    with open(aln,'w') as o:
        r=subprocess.run(['mafft','--localpair','--maxiterate','1000','--thread',THREADS,'--quiet',
                          f'{T}/panel_domains.faa'],stdout=o,stderr=subprocess.PIPE)
    if r.returncode: say('MAFFT failed'); sys.exit(1)
A={}
for r_ in open(aln).read().split('>')[1:]:
    h,s=r_.split('\n',1); A[h.strip()]=s.replace('\n','')
W=len(next(iter(A.values())))
gap=sum(s.count('-') for s in A.values())/(len(A)*W)
say(f'alignment {len(A)} x {W}, {100*gap:.1f} % gaps')
say(f'sequences >50 % gaps: {sum(1 for s in A.values() if s.count("-")/W>0.5)}')

trims=[]
for name,args in [('gappyout',['-gappyout']),('gt0.1',['-gt','0.1']),
                  ('gt0.2',['-gt','0.2']),('gt0.3',['-gt','0.3'])]:
    out=f'{ALN}/panel.{name}.faa'
    r=subprocess.run(['trimal','-in',aln,'-out',out]+args,capture_output=True,text=True)
    if r.returncode: say(f'  {name}: trimal failed'); continue
    B={}
    for rr in open(out).read().split('>')[1:]:
        h,s=rr.split('\n',1); B[h.strip()]=s.replace('\n','')
    w=len(next(iter(B.values())))
    say(f'  {name:<10} {w:>5} columns')
    trims.append((name,out,w))

say('\nIQ-TREE: ModelFinder + 1000 UFBoot + 1000 SH-aLRT, rooted on the muscarinic outgroup')
for name,path,w in [t for t in trims if t[0]=='gt0.1']+[t for t in trims if t[0]!='gt0.1']:
    pre=f'{TRE}/{name}'
    if os.path.exists(pre+'.iqtree'): say(f'  {name}: already done'); continue
    r=subprocess.run([IQ,'-s',path,'-m','MFP','-B','1000','-alrt','1000','-T',THREADS,
                      '--prefix',pre,'-o',','.join(OUTGROUP),'--quiet','-redo'],
                     capture_output=True,text=True)
    if r.returncode: say(f'  {name}: FAILED {(r.stderr or r.stdout)[-250:]}')
    else:
        m=''
        for ln in open(pre+'.iqtree'):
            if ln.startswith('Best-fit model'): m=ln.split(':')[1].strip(); break
        say(f'  {name}: done, {m}')
open(f'{T}/build_report.txt','w').write('\n'.join(rep)+'\n')
say('\n-> output/tree/trees/')
