#!/usr/bin/env python3
"""topology_test.py — is HTR2A in the right place, and does it matter for DOP-6?

In the ML tree HTR2A pairs with C. elegans SER-1 (88/99) and that pair sits inside the block that
also holds the D2-like dopamine receptors — so HTR2A ends up loosely grouped with DRD2. The
textbook vertebrate amine-receptor topology instead splits the family largely by G protein:
ADRA1 + 5-HT2 together (Gq), ADRA2 + DRD2 + 5-HT1 together (Gi/o), ADRB + DRD1 together (Gs).
On that expectation HTR2A should sit with ADRA1A/B/D, NOT near DRD2.

This matters beyond bookkeeping: SER-1 and HTR2A are DOP-6's nearest characterised relatives, so if
their placement is wrong, the statement "DOP-6's nearest characterised relatives are Gq-coupled
5-HT receptors" is wrong too.

Test: constrain the canonical Gq cluster (hADRA1A, hADRA1B, hADRA1D, hHTR2A, cel_ser-1) to be
monophyletic, re-search under the same model, and compare against the unconstrained tree with the
AU and SH tests (10,000 RELL replicates). Reported as verdicts plus dLogL, never as a comparison of
p-value magnitudes.

Writes output/tree/topology_test.txt
"""
import os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__)); T = os.path.join(HERE, 'output', 'tree')
IQ = '/Users/mikeodonnell/git/receptor_phylogeny/scripts/bin/iqtree3-mt'
ALN = f'{T}/aln/panel.gt0.1.faa'
WORK = f'{T}/topo'; os.makedirs(WORK, exist_ok=True)
MODEL = 'Q.PFAM+F+R6'
GQ = ['hADRA1A', 'hADRA1B', 'hADRA1D', 'hHTR2A', 'cel_ser-1']
rep = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); rep.append(s)

names = [l[1:].strip() for l in open(ALN) if l.startswith('>')]
say('=' * 88)
say('IS HTR2A MISPLACED, AND DOES IT MOVE DOP-6?')
say('=' * 88)
say(f'alignment {ALN.split("/")[-1]}, {len(names)} taxa, model {MODEL} (fixed, from ModelFinder)')
say(f'constraint: {" + ".join(GQ)} monophyletic  (the canonical Gq cluster)')

# A partial constraint file must still CONTAIN A SPLIT. Writing just "(A,B,C,D,E);" is a star of
# five taxa with no internal edge, so it encodes no constraint at all and IQ-TREE silently returns
# the unconstrained tree - which is how the first version of this test ended up comparing a tree
# with itself (RF = 0, dLogL = 5e-06). The group must be nested against at least one other taxon.
OTHERS = ['hCHRM1', 'hCHRM3']
open(f'{WORK}/gq.constr', 'w').write('((' + ','.join(GQ) + '),' + ','.join(OTHERS) + ');\n')
say(f'constraint file: (({" ,".join(GQ)}), {", ".join(OTHERS)}) -> {len(GQ)+len(OTHERS)} of {len(names)} taxa.')
say('PARTIAL constraint (silent about the other taxa), but it does contain a real split, which a')
say('bare "(A,B,C,D,E);" does not.')

if True:
    say('\nrunning constrained ML search ...')
    r = subprocess.run([IQ, '-s', ALN, '-m', MODEL, '-g', f'{WORK}/gq.constr', '-T', '4',
                        '--prefix', f'{WORK}/gq', '--quiet', '-redo'], capture_output=True, text=True)
    if r.returncode: say('constrained search FAILED: ' + (r.stderr or r.stdout)[-300:]); sys.exit(1)

# SANITY CHECK: did the constraint actually change the topology, and is the group now a clade?
def _splits(path):
    nw=open(path).read().strip().rstrip(';'); pos=[0]
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
            return ('int',lab,kids)
        m=re.match(r'[^,\)\(:]*',nw[pos[0]:]); lab=m.group(0); pos[0]+=len(lab)
        if pos[0]<len(nw) and nw[pos[0]]==':':
            m2=re.match(r':[-\d.eE+]*',nw[pos[0]:]); pos[0]+=len(m2.group(0))
        return ('tip',lab,[])
    t=node()
    def tp(n): return [n[1]] if n[0]=='tip' else [x for k in n[2] for x in tp(k)]
    allt=frozenset(tp(t)); S=set()
    def walk(n):
        if n[0]=='tip': return
        s=frozenset(tp(n))
        if 1<len(s)<len(allt):
            o=frozenset(allt-s); S.add(s if sorted(s)[0]<sorted(o)[0] else o)
        for k in n[2]: walk(k)
    walk(t); return S
import re
sa=_splits(f'{T}/trees/gt0.1.treefile'); sb=_splits(f'{WORK}/gq.treefile')
rf=len(sa-sb)+len(sb-sa)
took = frozenset(GQ) in sb
say(f'\nsanity check: Robinson-Foulds(ML, constrained) = {rf}; Gq group is a clade in the '
    f'constrained tree: {took}')
if rf == 0 or not took:
    say('THE CONSTRAINT DID NOT TAKE - the two trees are the same, so any AU test on them is')
    say('meaningless. Aborting rather than reporting a null result that is really a bug.')
    open(f'{T}/topology_test.txt','w').write('\n'.join(rep)+'\n'); sys.exit(1)

# collect the two topologies and run AU/SH
with open(f'{WORK}/both.trees', 'w') as f:
    f.write(open(f'{T}/trees/gt0.1.treefile').read().strip() + '\n')
    f.write(open(f'{WORK}/gq.treefile').read().strip() + '\n')
say('running AU / SH topology tests (10,000 RELL replicates) ...')
r = subprocess.run([IQ, '-s', ALN, '-m', MODEL, '-z', f'{WORK}/both.trees', '-n', '0', '-zb', '10000',
                    '-au', '-T', '4', '--prefix', f'{WORK}/au', '--quiet', '-redo'],
                   capture_output=True, text=True)
if r.returncode: say('AU test FAILED: ' + (r.stderr or r.stdout)[-400:]); sys.exit(1)
out = open(f'{WORK}/au.iqtree').read()
blk = out[out.find('USER TREES'):]
say('')
for ln in blk.split('\n'):
    if ln.strip() and ('logL' in ln or ln.strip().startswith(('1', '2', '-')) or 'deltaL' in ln):
        say('  ' + ln.rstrip())

# where does dop-6 sit in the constrained tree?
import re
def tips_of_clade_containing(nw, target, maxn=12):
    nw = nw.strip().rstrip(';'); pos=[0]
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
            return ('int',lab,kids)
        m=re.match(r'[^,\)\(:]*',nw[pos[0]:]); lab=m.group(0); pos[0]+=len(lab)
        if pos[0]<len(nw) and nw[pos[0]]==':':
            m2=re.match(r':[-\d.eE+]*',nw[pos[0]:]); pos[0]+=len(m2.group(0))
        return ('tip',lab,[])
    t=node()
    def tips(n): return [n[1]] if n[0]=='tip' else [x for k in n[2] for x in tips(k)]
    best=None
    def walk(n):
        nonlocal best
        if n[0]=='tip': return
        tl=tips(n)
        if target in tl and len(tl)<=maxn and (best is None or len(tl)>len(best)): best=tl
        for k in n[2]: walk(k)
    walk(t)
    return best
say('')
say('DOP-6\'s neighbourhood in each tree (largest clade containing it with <= 12 tips):')
for lab, f in (('unconstrained ML', f'{T}/trees/gt0.1.treefile'), ('Gq-constrained', f'{WORK}/gq.treefile')):
    tl = tips_of_clade_containing(open(f).read(), 'cel_dop-6')
    say(f'  {lab:<18}: {", ".join(sorted(tl)) if tl else "?"}')
open(f'{T}/topology_test.txt', 'w').write('\n'.join(rep) + '\n')
say('\n-> output/tree/topology_test.txt')
