#!/usr/bin/env python3
"""assess_resolution.py — does this matrix actually resolve the tree, and where does it fail?

356 parsimony-informative sites for 163 taxa is ~1.1 sites per free parameter, which is the regime
where domain-only alignments are routinely under-determined. So before reading anything off the
topology, measure the resolution directly: the distribution of node support, split by how deep the
node is. A matrix can resolve shallow structure (orthologue groups, species pairs) perfectly well
while leaving every deep node a coin flip - and it is the deep nodes that DOP-6's placement depends
on.

Reports support by clade size, and specifically the support of the nodes on the path from DOP-6 to
the root, since those are the ones any placement claim rests on.

Writes output/tree/resolution.txt
"""
import os, re, sys
from collections import defaultdict

HERE=os.path.dirname(os.path.abspath(__file__)); T=os.path.join(HERE,'output','tree')
rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)

def parse(nw):
    nw=nw.strip().rstrip(';'); pos=[0]
    def node():
        if nw[pos[0]]=='(':
            pos[0]+=1; kids=[]
            while True:
                kids.append(node())
                if nw[pos[0]]==',': pos[0]+=1; continue
                if nw[pos[0]]==')': pos[0]+=1; break
            m=re.match(r'[^,\)\(:]*',nw[pos[0]:]); lab=m.group(0); pos[0]+=len(lab)
            if pos[0]<len(nw) and nw[pos[0]]==':':
                m2=re.match(r':[-\d.eE+]*',nw[pos[0]:]); pos[0]+=len(m2.group(0))
            return {'tip':False,'lab':lab,'kids':kids}
        m=re.match(r'[^,\)\(:]*',nw[pos[0]:]); lab=m.group(0); pos[0]+=len(lab)
        if pos[0]<len(nw) and nw[pos[0]]==':':
            m2=re.match(r':[-\d.eE+]*',nw[pos[0]:]); pos[0]+=len(m2.group(0))
        return {'tip':True,'lab':lab,'kids':[]}
    return node()
def tips(n): return [n['lab']] if n['tip'] else [x for k in n['kids'] for x in tips(k)]

for trim in ['gt0.1','gappyout','gt0.2','gt0.3']:
    tf=f'{T}/trees/{trim}.treefile'
    if not os.path.exists(tf): continue
    t=parse(open(tf).read()); N=len(tips(t))
    nodes=[]
    def walk(n):
        if n['tip']: return
        m=re.match(r'([\d.]+)/([\d.]+)',n['lab'] or '')
        if m: nodes.append((len(tips(n)), float(m.group(1)), float(m.group(2))))
        for k in n['kids']: walk(k)
    walk(t)
    if not nodes: say(f'{trim}: no support labels'); continue
    say('='*86)
    say(f'{trim}: {N} taxa, {len(nodes)} internal nodes with support')
    say('='*86)
    say(f'{"clade size":<16} {"n nodes":>8} {"UFboot>=95":>11} {"UFboot>=70":>11} {"median UFboot":>14}')
    bands=[('2 (tip pairs)',2,2),('3-5 (shallow)',3,5),('6-15',6,15),
           ('16-40 (deep)',16,40),('>40 (backbone)',41,10**6)]
    for lab,lo,hi in bands:
        sel=[n for n in nodes if lo<=n[0]<=hi]
        if not sel: continue
        b=[n[2] for n in sel]; b.sort()
        say(f'{lab:<16} {len(sel):>8} {100*sum(1 for x in b if x>=95)/len(b):>10.0f}% '
            f'{100*sum(1 for x in b if x>=70)/len(b):>10.0f}% {b[len(b)//2]:>14.0f}')
    allb=sorted(n[2] for n in nodes)
    say(f'{"ALL":<16} {len(allb):>8} {100*sum(1 for x in allb if x>=95)/len(allb):>10.0f}% '
        f'{100*sum(1 for x in allb if x>=70)/len(allb):>10.0f}% {allb[len(allb)//2]:>14.0f}')
    # path from dop-6 to root
    path=[]
    def find(n,anc):
        if n['tip']:
            if n['lab']=='cel_dop-6': path.extend(anc)
            return
        for k in n['kids']: find(k,anc+[k] if False else anc+[n])
    find(t,[])
    say('')
    say(f'nodes on the path from cel_dop-6 outward (the ones a placement claim rests on):')
    shown=0
    for n in reversed(path):
        sz=len(tips(n))
        if sz>60 or shown>8: break
        m=re.match(r'([\d.]+)/([\d.]+)',n['lab'] or '')
        sup=f'{float(m.group(1)):.0f}/{float(m.group(2)):.0f}' if m else 'n/a'
        tl=sorted(tips(n))
        cel=[x for x in tl if x.startswith('cel_')]
        say(f'   {sz:>3} tips  support {sup:<10}  C. elegans members: '
            f'{", ".join(c[4:] for c in cel[:8])}{" ..." if len(cel)>8 else ""}')
        shown+=1
    say('')
open(f'{T}/resolution.txt','w').write('\n'.join(rep)+'\n')
say('-> output/tree/resolution.txt')
