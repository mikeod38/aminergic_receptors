#!/usr/bin/env python3
"""orphan_clade.py — read the orphan-family tree: who is in, how deep, and is the
alignment trustworthy at this depth?

Runs after build_orphan_tree.py. Four questions, in the order that decides whether
the later ones mean anything:

  0. POSITIVE CONTROL (methodology 7.2). The melatonin family is monophyletic on
     evidence independent of this panel - 100/100 in output/tree2, and both lamprey
     members are reciprocal best hits to human MTNR1A. If it fragments here, the
     alignment cannot resolve this depth and nothing below is interpretable. This is
     asserted: the script says FAIL loudly rather than carrying on quietly.
  1. Does the seed family come out monophyletic, and with what support?
  2. Which of the 93 above-floor hits are actually inside it? The floor was set by
     human ADRB2, so other species' aminergic receptors could clear it; the ones
     that fall among the aminergic backbone instead are bleed-through, and that is
     a number worth having rather than a nuisance.
  3. How deep does the clade reach - which is just "what is the oldest lineage with
     a member inside it", read off the tips rather than asserted.

Writes output/orphan/tree/clade.txt
"""
import os, re, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'output', 'orphan', 'tree')

SEED = ['ORPH_Hsap_Q9Y2T5', 'ORPH_Hsap_Q99679', 'ORPH_Spur_A0A7M7GHW5',
        'ORPH_Dmel_Q9VZD1', 'ORPH_Pcau_A0ABM1EQA9', 'ORPH_Cvir_A0A8B8D717',
        'ORPH_Lana_A0A1S3KFU7']
# Matched by accession suffix rather than by full label: the first run of this script
# "failed" its own positive control because the panel had MEL_MEL_Drer_* and this list
# said MEL_Drer_*, so 5 of 10 members silently went missing and the 5 real zebrafish
# melatonin receptors were reported as intruders. Never anchor a control on a label
# that another script composes.
MEL_KEYS = ['hMTNR1A', 'hMTNR1B', 'hGPR50', 'A0A8M1NRI5', 'A0A8M9PN98', 'A2BGT9',
            'P51046', 'P51049', 'A0AAJ7T8G6', 'A0ACM8ASR7']
SEED_KEYS = ['Q9Y2T5', 'Q99679', 'A0A7M7GHW5', 'Q9VZD1', 'A0ABM1EQA9', 'A0A8B8D717',
             'A0A1S3KFU7']


def resolve(keys, allt):
    """key -> the one tip whose label ends with it (or equals it)."""
    out, miss = [], []
    for k in keys:
        m = [t for t in allt if t == k or t.endswith('_' + k)]
        if len(m) == 1:
            out.append(m[0])
        else:
            miss.append(f'{k}({len(m)} matches)')
    return out, miss
NP = ['hNPY2R', 'hNPFFR2', 'hTACR3']

# oldest-first, so the deepest member inside the clade can be read off directly
DEPTH = [('Aque', 'sponge'), ('Tadh', 'Trichoplax'), ('Mlei', 'ctenophore'),
         ('Nvec', 'cnidarian (Nematostella)'), ('Hvul', 'cnidarian (Hydra)'),
         ('Cele', 'nematode'), ('Dmel', 'insect'), ('Lpol', 'chelicerate'),
         ('Pcau', 'priapulid'), ('Cvir', 'mollusc (oyster)'),
         ('Myes', 'mollusc (scallop)'), ('Lana', 'brachiopod'),
         ('Ebur', 'hagfish'), ('Pmin', 'echinoderm (starfish)'),
         ('Apla', 'echinoderm (crown-of-thorns)'), ('Ajap', 'echinoderm (sea cucumber)'),
         ('Spur', 'echinoderm (sea urchin)'), ('Skow', 'hemichordate'),
         ('Bflo', 'amphioxus'), ('Cint', 'tunicate'), ('Pmar', 'lamprey'),
         ('Drer', 'teleost'), ('Mmus', 'mouse'), ('Hsap', 'human')]


class N:
    __slots__ = ('par', 'ch', 'name', 'sup', 'bl')

    def __init__(self, par=None):
        self.par, self.ch, self.name, self.sup, self.bl = par, [], None, None, None


def parse(t):
    t = t.strip().rstrip(';')
    root = N(); cur = root; buf = ''

    def flush():
        nonlocal buf
        if not buf:
            return
        lab, _, bl = buf.partition(':')
        if cur.ch:
            cur.sup = lab or None
        else:
            cur.name = lab or None
        try:
            cur.bl = float(bl)
        except ValueError:
            pass
        buf = ''

    for c in t:
        if c == '(':
            n = N(cur); cur.ch.append(n); cur = n
        elif c == ',':
            flush(); n = N(cur.par); cur.par.ch.append(n); cur = n
        elif c == ')':
            flush(); cur = cur.par
        else:
            buf += c
    flush()
    return root


def tips(n):
    if not n.ch:
        return {n.name}
    s = set()
    for c in n.ch:
        s |= tips(c)
    return s


def walk(n, out=None):
    out = [] if out is None else out
    out.append(n)
    for c in n.ch:
        walk(c, out)
    return out


def sup(n):
    if n.sup:
        return n.sup
    if n.par is None:
        return 'root'
    sib = [c for c in n.par.ch if c is not n]
    if n.par.par is None and len(sib) == 1:
        return sib[0].sup or 'unlabelled'
    return 'unlabelled'


def mrca(root, target):
    best = None
    for n in walk(root):
        tp = tips(n)
        if target <= tp and (best is None or len(tp) < len(tips(best))):
            best = n
    return best


def species_of(tip):
    m = re.match(r'^(?:ORPH|MEL|AM)_([A-Z][a-z]{2,3})_', tip)
    if m:
        return m.group(1)
    if tip.startswith('h'):
        return 'Hsap'
    if tip.startswith('cel_'):
        return 'Cele'
    if tip.startswith('dm'):
        return 'Dmel'
    m = re.match(r'^(?:AM_)?([A-Z][a-z]{2,3})_', tip)
    return m.group(1) if m else '?'


L = []
for tag in ('gt0.1', 'gappyout'):
    path = f'{OUT}/{tag}.treefile'
    if not os.path.exists(path):
        L.append(f'=== {tag}: not built ==='); continue
    root = parse(open(path).read())
    allt = tips(root)
    L.append('=' * 78)
    L.append(f'=== {tag} — {len(allt)} tips ===')
    L.append('')

    # 0. positive control
    mel, melmiss = resolve(MEL_KEYS, allt)
    if melmiss:
        L.append(f'  LABEL PROBLEM, control cannot be judged: {", ".join(melmiss)}')
    node = mrca(root, set(mel))
    inside = tips(node)
    extra = sorted(inside - set(mel))
    ok = not extra
    L.append(f'POSITIVE CONTROL — melatonin family monophyly: '
             f'{"PASS" if ok else "FAIL"}')
    L.append(f'  {len(mel)} members, smallest containing clade has {len(inside)} tips, '
             f'support {sup(node)}')
    if extra:
        L.append(f'  intruders ({len(extra)}): {", ".join(extra[:12])}'
                 + (' ...' if len(extra) > 12 else ''))
        L.append('  The alignment does not cleanly resolve a group that independent '
                 'evidence says is monophyletic. Read everything below with that in mind.')
    L.append('')

    # 1. the seed family
    seed, seedmiss = resolve(SEED_KEYS, allt)
    if seedmiss:
        L.append(f'  NOTE: seed members not uniquely resolvable: {", ".join(seedmiss)}')
    fam = mrca(root, set(seed))
    fam_tips = tips(fam)
    L.append(f'ORPHAN FAMILY — smallest clade containing all {len(seed)} seed members:')
    L.append(f'  {len(fam_tips)} tips, support {sup(fam)}')
    L.append(f'  exactly monophyletic: {len(fam_tips) == len(seed)}')
    L.append('')

    # 2. bleed-through
    orph = {t for t in allt if t.startswith('ORPH_')}
    inside_o = sorted(orph & fam_tips)
    outside_o = sorted(orph - fam_tips)
    L.append(f'ABOVE-FLOOR HITS — {len(orph)} in the panel, '
             f'{len(inside_o)} inside the family clade, {len(outside_o)} outside '
             f'(aminergic bleed-through past the ADRB2 floor):')
    byspec = defaultdict(lambda: [0, 0])
    for t in inside_o:
        byspec[species_of(t)][0] += 1
    for t in outside_o:
        byspec[species_of(t)][1] += 1
    L.append(f'  {"sp":8} {"in":>4} {"out":>4}')
    for sp in sorted(byspec, key=lambda s: -byspec[s][0]):
        i, o = byspec[sp]
        L.append(f'  {sp:8} {i:4} {o:4}')
    L.append('')
    nonorph = sorted(fam_tips - orph)
    if nonorph:
        L.append(f'  Also inside the family clade, from other parts of the panel '
                 f'({len(nonorph)}): {", ".join(nonorph[:14])}'
                 + (' ...' if len(nonorph) > 14 else ''))
        L.append('')

    # 3. PCDR-1
    pc = [t for t in allt if 'pcdr-1' in t or t == 'ORPH_Cele_O17899']
    for t in pc:
        L.append(f'PCDR-1 ({t}): inside the family clade = {t in fam_tips}')
        if t in fam_tips:
            n = next(x for x in walk(root) if not x.ch and x.name == t)
            cur = n.par
            lvl = 0
            while cur is not None and lvl < 5:
                tp = sorted(tips(cur))
                L.append(f'  + {len(tp):3} tips @ {sup(cur):12}  {", ".join(tp)}'
                         if len(tp) <= 12 else
                         f'  + {len(tp):3} tips @ {sup(cur):12}  '
                         f'{", ".join(tp[:12])} ...')
                cur = cur.par
                lvl += 1
    L.append('')

    # 4. depth
    present = {species_of(t) for t in fam_tips}
    L.append('DEPTH — oldest lineages with a member inside the family clade:')
    hit = False
    for code, label in DEPTH:
        if code in present:
            mark = '<-- deepest' if not hit else ''
            L.append(f'  {code:6} {label:30} in clade  {mark}')
            hit = True
    absent = [f'{c} ({l})' for c, l in DEPTH if c not in present]
    L.append('')
    L.append(f'  Not in the clade: {", ".join(absent)}')
    L.append('')

txt = '\n'.join(L)
print(txt)
open(f'{OUT}/clade.txt', 'w').write(txt + '\n')
