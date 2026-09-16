#!/usr/bin/env python3
"""trait_signal.py — can this tree predict a receptor's LIGAND, its Galpha COUPLING, both, or neither?

The project's question (a) is whether DOP-6 is an octopamine receptor, and the tempting way to
answer it is to look at which clade DOP-6 lands in. That only works if ligand identity actually
tracks phylogeny in this family — which is a separate question from whether coupling does, and has
to be tested before either is used. Aminergic receptors are known to group by receptor subtype
across species rather than by transmitter, so the two traits can behave very differently.

For each trait: Fitch parsimony (minimum state changes on the ML topology) for the tips where the
trait is experimentally known, against 10,000 random reassignments of the same labels to the same
tips. Run on every trim level, so the answer cannot be a trimming artefact.

Writes output/tree/trait_signal.txt
"""
import os, random, re, sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__)); T = os.path.join(HERE, 'output', 'tree')
TRIMS = ['gappyout', 'gt0.1', 'gt0.2', 'gt0.3']
NPERM = 10000
random.seed(0)
rep = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); rep.append(s)

LIGAND = {  # experimentally determined transmitter; orphans and inferences excluded
 'cel_dop-1': 'DA', 'cel_dop-2': 'DA', 'cel_dop-3': 'DA', 'cel_dop-4': 'DA',
 'hDRD1': 'DA', 'hDRD2': 'DA', 'dmDop1R1': 'DA', 'dmDop2R': 'DA',
 'cel_ser-1': '5HT', 'cel_ser-4': '5HT', 'cel_ser-7': '5HT', 'hHTR1A': '5HT', 'hHTR2A': '5HT',
 'cel_ser-3': 'OA', 'cel_ser-6': 'OA', 'cel_octr-1': 'OA',
 'dmOAMB': 'OA', 'dmOctA2R': 'OA', 'dmOctB1R': 'OA', 'dmOctB2R': 'OA', 'dmOctB3R': 'OA',
 'bmOAR1': 'OA',
 'cel_ser-2': 'TA', 'cel_tyra-2': 'TA', 'dmTyrR_Q9VEG1': 'TA', 'dmTyrRII_Q9VEG2': 'TA',
 'hADRA1A': 'NE', 'hADRA1B': 'NE', 'hADRA1D': 'NE', 'hADRA2A': 'NE', 'hADRA2C': 'NE',
 'hADRB1': 'NE', 'hADRB2': 'NE',
 'hCHRM1': 'ACh', 'hCHRM3': 'ACh', 'cel_gar-3': 'ACh',
}
COUP = {}
for ln in open(f'{T}/coupling.tsv').read().strip().split('\n')[1:]:
    p = ln.split('\t'); c = p[1].strip()
    if not c or c == '?' or c.endswith('*') or '+' in c: continue
    COUP[p[0]] = 'Gi/o' if c in ('Gi/o', 'Go', 'Gi') else c

def parse(nw):
    nw = nw.strip().rstrip(';'); pos = [0]
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
def tips(n): return [n[1]] if n[0] == 'tip' else [t for k in n[2] for t in tips(k)]
def fitch(n, st):
    if n[0] == 'tip':
        s = st.get(n[1]); return ({s} if s else None), 0
    sets, ch = [], 0
    for k in n[2]:
        s, c = fitch(k, st); ch += c
        if s: sets.append(s)
    if not sets: return None, ch
    inter = set.intersection(*sets)
    return (inter, ch) if inter else (set.union(*sets), ch + 1)

say('=' * 96)
say('WHAT CAN THIS TREE ACTUALLY PREDICT?')
say('=' * 96)
say('Fitch parsimony on the ML topology vs 10,000 random label permutations, per trait, per trim.')
say('Fewer changes than random = the trait is phylogenetically clustered and a prediction for an')
say('orphan is meaningful. Indistinguishable from random = the tree cannot answer that question.\n')

results = {}
for trait, tbl in (('LIGAND', LIGAND), ('Galpha COUPLING', COUP)):
    say(f'--- {trait}: {len(tbl)} tips with experimental data, states {dict(Counter(tbl.values()))}')
    for trim in TRIMS:
        tf = f'{T}/trees/{trim}.treefile'
        if not os.path.exists(tf): say(f'    {trim:<10} no treefile'); continue
        tree = parse(open(tf).read()); tp = set(tips(tree))
        st = {k: v for k, v in tbl.items() if k in tp}
        _, obs = fitch(tree, st)
        keys, labs = list(st.keys()), list(st.values())
        null = []
        for _ in range(NPERM):
            random.shuffle(labs)
            null.append(fitch(tree, dict(zip(keys, labs)))[1])
        p = (sum(1 for x in null if x <= obs) + 1) / (NPERM + 1)
        mn = sum(null)/len(null)
        verdict = 'CLUSTERED' if p < 0.05 else 'no signal'
        say(f'    {trim:<10} {len(st):>2} tips   observed {obs:>2} changes   random mean {mn:>5.1f}   '
            f'p = {p:.4f}   {verdict}')
        results.setdefault(trait, []).append((trim, obs, mn, p))
    say('')

say('-' * 96)
say('READING')
say('-' * 96)
for trait, rs in results.items():
    ps = [r[3] for r in rs]
    if all(p < 0.05 for p in ps):
        say(f'{trait}: clustered at every trim level -> tree position is informative about it.')
    elif all(p >= 0.05 for p in ps):
        say(f'{trait}: no signal at any trim level -> tree position says NOTHING about it.')
    else:
        say(f'{trait}: signal is trim-dependent (p = {", ".join(f"{p:.3f}" for p in ps)}) -> weak, '
            f'do not rely on it.')
say('')
say('If ligand is unpredictable but coupling is, then this tree can suggest what G protein DOP-6')
say('uses but cannot tell you whether its ligand is octopamine. That question would need a binding')
say('or activation assay, not a tree.')
open(f'{T}/trait_signal.txt', 'w').write('\n'.join(rep) + '\n')
say('-> output/tree/trait_signal.txt')
