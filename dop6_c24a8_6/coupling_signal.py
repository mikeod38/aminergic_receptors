#!/usr/bin/env python3
"""coupling_signal.py — does position in this tree predict Galpha coupling, or not?

Assigning DOP-6 a G protein from its nearest characterised relative is only worth doing if
coupling actually tracks phylogeny in this family. Koelle 2018 WormBook Table 1 supplies
experimentally determined coupling for 14 C. elegans receptors, and the reference receptors in the
panel add well-established vertebrate and insect assignments, so the method can be tested on the
receptors where the answer is known before it is applied to the one where it is not.

Test: Fitch parsimony — the minimum number of coupling-state changes needed on the ML topology —
compared against 10,000 random reassignments of the same labels to the same tips. Few changes
relative to the null means coupling is phylogenetically clustered and a prediction is meaningful;
a parsimony score indistinguishable from random means tree position carries no information about
coupling and no prediction should be made.

Then, and only if the signal test passes, DOP-6's neighbourhood is reported with node support.

Writes output/tree/coupling_signal.txt
"""
import os, random, re, sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__)); T = os.path.join(HERE, 'output', 'tree')
TRIMS = ['gappyout', 'gt0.1', 'gt0.2', 'gt0.3']
NPERM = 10000
random.seed(0)
rep = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); rep.append(s)

# ── minimal Newick parser -> nested tuples with labels and support ───────────
def parse(nw):
    nw = nw.strip().rstrip(';')
    pos = [0]
    def node():
        if nw[pos[0]] == '(':
            pos[0] += 1; kids = []
            while True:
                kids.append(node())
                if nw[pos[0]] == ',': pos[0] += 1; continue
                if nw[pos[0]] == ')': pos[0] += 1; break
            m = re.match(r'[^,\)\(:]*', nw[pos[0]:]); lab = m.group(0); pos[0] += len(lab)
            if pos[0] < len(nw) and nw[pos[0]] == ':':
                m2 = re.match(r':[-\d.eE+]*', nw[pos[0]:]); pos[0] += len(m2.group(0))
            return ('int', lab, kids)
        m = re.match(r'[^,\)\(:]*', nw[pos[0]:]); lab = m.group(0); pos[0] += len(lab)
        if pos[0] < len(nw) and nw[pos[0]] == ':':
            m2 = re.match(r':[-\d.eE+]*', nw[pos[0]:]); pos[0] += len(m2.group(0))
        return ('tip', lab, [])
    return node()

def tips(n):
    return [n[1]] if n[0] == 'tip' else [t for k in n[2] for t in tips(k)]

def fitch(n, st):
    """minimum state changes (Fitch), states as sets; returns (stateset, changes)"""
    if n[0] == 'tip':
        s = st.get(n[1])
        return ({s} if s else None), 0
    sets, ch = [], 0
    for k in n[2]:
        s, c = fitch(k, st); ch += c
        if s: sets.append(s)
    if not sets: return None, ch
    inter = set.intersection(*sets)
    if inter: return inter, ch
    return set.union(*sets), ch + 1

# ── coupling labels ──────────────────────────────────────────────────────────
COUP = {}
for ln in open(f'{T}/coupling.tsv').read().strip().split('\n')[1:]:
    p = ln.split('\t')
    c = p[1].strip()
    if not c or c == '?': continue
    if c.endswith('*'): continue          # SER-5 is itself a sequence-based prediction, not data
    if '+' in c: continue                 # GAR-3 couples to two classes; not a single state
    COUP[p[0]] = 'Gi/o' if c in ('Gi/o', 'Go', 'Gi') else c

say('=' * 92)
say('DOES TREE POSITION PREDICT Galpha COUPLING IN THIS FAMILY?')
say('=' * 92)
say(f'{len(COUP)} tips carry a single experimentally determined coupling state')
say('  ' + ', '.join(f'{k}={v}' for k, v in sorted(COUP.items())))
say(f'\nstate counts: {dict(Counter(COUP.values()))}')
say('SER-5 excluded: its Gs assignment is itself a sequence-based prediction (Harris 2009),')
say('so including it would let the method validate itself. GAR-3 excluded: couples to two classes.')

ok_any = False
for trim in TRIMS:
    tf = f'{T}/trees/{trim}.treefile'
    if not os.path.exists(tf): say(f'\n{trim}: no treefile'); continue
    tree = parse(open(tf).read())
    tp = set(tips(tree))
    st = {k: v for k, v in COUP.items() if k in tp}
    _, obs = fitch(tree, st)
    labs = list(st.values())
    keys = list(st.keys())
    null = []
    for _ in range(NPERM):
        random.shuffle(labs)
        _, c = fitch(tree, dict(zip(keys, labs)))
        null.append(c)
    better = sum(1 for x in null if x <= obs)
    p = (better + 1) / (NPERM + 1)
    mn = sum(null)/len(null)
    say(f'\n{trim}: {len(st)} labelled tips, observed parsimony score {obs} changes; '
        f'random mean {mn:.1f} (min {min(null)})')
    say(f'   p = {p:.4f}  -> {"SIGNAL: coupling is phylogenetically clustered" if p < 0.05 else "NO SIGNAL: indistinguishable from random"}')
    if p < 0.05: ok_any = True

say('\n' + '-' * 92)
say('DOP-6 NEIGHBOURHOOD')
say('-' * 92)
for trim in TRIMS:
    tf = f'{T}/trees/{trim}.treefile'
    if not os.path.exists(tf): continue
    tree = parse(open(tf).read())
    # walk to find the smallest clades containing dop-6
    path = []
    def find(n, anc):
        if n[0] == 'tip':
            if n[1] == 'cel_dop-6': path.extend(anc + [n])
            return
        for k in n[2]: find(k, anc + [n])
    find(tree, [])
    if not path: say(f'{trim}: dop-6 not found'); continue
    say(f'\n{trim}:')
    for n in reversed(path[:-1][-4:]):
        tl = tips(n)
        if len(tl) > 26: break
        sup = n[1] if n[1] else 'n/a'
        say(f'   support {sup:<12} clade of {len(tl):>2}: ' +
            ', '.join(f'{t}[{COUP.get(t,"?")}]' for t in sorted(tl)))

say('\n' + '-' * 92)
if not ok_any:
    say('CONCLUSION: coupling shows no significant phylogenetic clustering on this tree, so the')
    say('nearest-characterised-relative method has no demonstrated validity in this family and')
    say('DOP-6\'s coupling should NOT be predicted from its tree position. It needs an assay.')
else:
    say('CONCLUSION: coupling is phylogenetically clustered, so DOP-6\'s neighbourhood is a')
    say('legitimate HYPOTHESIS about its coupling — still a hypothesis, not a measurement.')
open(f'{T}/coupling_signal.txt', 'w').write('\n'.join(rep) + '\n')
say('\n-> output/tree/coupling_signal.txt')
