#!/usr/bin/env python3
"""rogue_test.py — do the two long-branch orphans distort DOP-6's placement?

cel_F59D12.1 has the longest terminal branch in the tree (5.98 subs/site, robust z = 17.4) and
cel_T21B4.4 the second (4.35, z = 12.3); DOP-6's own branch is 0.020. Long branches attract each
other and can drag unrelated taxa together, and both of these sit near DOP-6's clade.

A long branch is not by itself grounds for exclusion, so the question here is not whether to drop
them but whether their presence CHANGES anything. Three rebuilds, one variable each, against the
same alignment and model:

   drop_F59   F59D12.1 removed
   drop_T21   T21B4.4 removed
   drop_both  both removed

If DOP-6's clade composition is identical in all four trees, the rogues are harmless and stay in
(their instability is then itself a finding). If it moves, the placement was partly an artefact.

Writes output/tree/rogue_test.txt
"""
import os, re, subprocess, sys

HERE=os.path.dirname(os.path.abspath(__file__)); T=os.environ.get('DOP6_TREEDIR') or os.path.join(HERE,'output','tree')
IQ='/Users/mikeodonnell/git/receptor_phylogeny/scripts/bin/iqtree3-mt'
ALN=f'{T}/aln/panel.gt0.1.faa'
MODEL='LG+F+R7'          # the model ModelFinder chose for this alignment; fixed so only the taxon
                         # set differs between runs
W=f'{T}/rogue'; os.makedirs(W,exist_ok=True)
OUTG=['hCHRM1','hCHRM2','hCHRM3','hCHRM5','cel_gar-1','cel_gar-2','cel_gar-3']
RUNS=[('drop_F59',{'cel_F59D12.1'}),('drop_T21',{'cel_T21B4.4'}),
      ('drop_both',{'cel_F59D12.1','cel_T21B4.4'})]
rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)
seqs={}
for r in open(ALN).read().split('>')[1:]:
    h,s=r.split('\n',1); seqs[h.strip()]=s.replace('\n','')
say('='*90); say('ROGUE-TAXON CONTROL: does removing the long branches move DOP-6?'); say('='*90)
say(f'alignment {ALN.split("/")[-1]}, {len(seqs)} taxa, model {MODEL} fixed across runs\n')
for name,drop in RUNS:
    fa=f'{W}/{name}.faa'
    with open(fa,'w') as f:
        for k,v in seqs.items():
            if k in drop: continue
            f.write(f'>{k}\n{v}\n')
    if os.path.exists(f'{W}/{name}.treefile'):
        say(f'{name}: already built'); continue
    og=[o for o in OUTG if o not in drop]
    say(f'{name}: building ({len(seqs)-len(drop)} taxa) ...')
    r=subprocess.run([IQ,'-s',fa,'-m',MODEL,'-B','1000','-alrt','1000','-T','6',
                      '--prefix',f'{W}/{name}','-o',','.join(og),'--quiet','-redo'],
                     capture_output=True,text=True)
    if r.returncode: say(f'  FAILED: {(r.stderr or r.stdout)[-250:]}')
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
def clades(tf):
    t=parse(open(tf).read()); path=[]
    def find(n,anc):
        if n['tip']:
            if n['lab']=='cel_dop-6': path.extend(anc)
            return
        for k in n['kids']: find(k,anc+[n])
    find(t,[])
    out=[]
    for n in reversed(path):
        tl=tips(n)
        if len(tl)>26: break
        m=re.match(r'([\d.]+)/([\d.]+)',n['lab'] or '')
        out.append((len(tl),f'{float(m.group(2)):.0f}' if m else '-',frozenset(tl)))
    return out
say('')
base=clades(f'{T}/trees/gt0.1.treefile')
cel=lambda S: ', '.join(sorted(x[4:] for x in S if x.startswith('cel_'))) or '-'
say('FULL TREE (163 taxa):')
for sz,sup,S in base: say(f'   {sz:>3} tips  UFboot {sup:>4}   C. elegans: {cel(S)}')
for name,drop in RUNS:
    tf=f'{W}/{name}.treefile'
    if not os.path.exists(tf): say(f'\n{name}: no tree'); continue
    cl=clades(tf)
    say(f'\n{name.upper()} ({163-len(drop)} taxa):')
    for sz,sup,S in cl: say(f'   {sz:>3} tips  UFboot {sup:>4}   C. elegans: {cel(S)}')
    # compare composition ignoring the dropped taxa
    b=[frozenset(S-drop) for _,_,S in base]
    same=sum(1 for _,_,S in cl if S in b)
    say(f'   -> {same}/{len(cl)} clades identical to the full tree (dropped taxa ignored)')
open(f'{T}/rogue_test.txt','w').write('\n'.join(rep)+'\n')
say('\n-> output/tree/rogue_test.txt')
