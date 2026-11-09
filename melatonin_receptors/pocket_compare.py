#!/usr/bin/env python3
"""pocket_compare.py — does PCDR-1's binding pocket look aminergic, melatonergic, or neither?

The 2003 family assignment that PCDR-1 still carries was made from predicted ligand
class. This asks the pocket directly.

Position 3.32 is anchored the way ../dop6_c24a8_6/pocket_332_v2.py anchors it: by the
column that every human aminergic receptor occupies with the canonical aspartate, with
the four structurally characterised references (hADRB2, hHTR2B, hDRD3, hHRH1) required
to agree. That positive control is asserted, not assumed - the script aborts if it fails.

Consensus for each group is a majority vote across members, not one reference species'
residue. The melatonin group is the 8-tip clade recovered at 100/100 in this project's
own tree, so it is a verified group rather than a name-based one.

Writes output/pocket_compare.txt
"""
import os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ALN = f'{HERE}/output/tree/panel.linsi.faa'

# the 8-tip melatonin clade recovered at 100/100 (see output/tree/placement.txt)
MEL = ['hMTNR1A', 'hMTNR1B', 'hGPR50', 'MEL_Drer_A0A8M1NRI5', 'MEL_Drer_A0A8M9PN98',
       'MEL_Drer_A2BGT9', 'MEL_Drer_P51046', 'MEL_Drer_P51049']
STRUCT = ['hADRB2', 'hHTR2B', 'hDRD3', 'hHRH1']   # ligand-bound reference structures
NP = ['hNPY2R', 'hNPFFR2', 'hTACR3', 'hQRFPR', 'hGPR52', 'hGPR21']
PCDR1 = 'cel_pcdr-1_F59D12.1'


def readfa(p):
    d, k = {}, None
    for line in open(p):
        if line.startswith('>'):
            k = line[1:].strip(); d[k] = ''
        elif k:
            d[k] += line.strip()
    return d


aln = readfa(ALN)
ncol = len(next(iter(aln.values())))
amine = [n for n in aln
         if n.startswith('h') and n not in MEL and n not in NP
         and not n.startswith(('hNP', 'hGPR5', 'hGPR2', 'hQRFPR', 'hTACR'))]

# --- anchor 3.32: the column where every human aminergic receptor carries D ----
cands = [c for c in range(ncol)
         if all(aln[n][c] == 'D' for n in STRUCT)
         and sum(aln[n][c] == 'D' for n in amine) / len(amine) >= 0.95]

# Several columns can satisfy that, so disambiguate against the calls the sibling
# project already verified structurally (../dop6_c24a8_6/output/tree/b120/pocket.tsv,
# where res_3_32 is index 2 of the 20-residue window). The melatonin receptors carry
# M/M/T there rather than D, so only the true column reproduces all seven.
VERIFIED = {}
for line in open(f'{HERE}/../dop6_c24a8_6/output/tree/b120/pocket.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) >= 2 and f[0] in aln:
        VERIFIED[f[0]] = f[1]
cands = [c for c in cands
         if all(aln[n][c] == r for n, r in VERIFIED.items() if n in STRUCT or n in MEL)]
if len(cands) != 1:
    raise SystemExit(f'ABORT: positive control failed - {len(cands)} candidate columns '
                     f'for 3.32 after cross-check, expected exactly 1')
col = cands[0]
agree = sum(aln[n][col] == r for n, r in VERIFIED.items())
ok = sum(aln[n][col] == 'D' for n in amine)
lines = [f'3.32 anchored to alignment column {col}',
         f'  positive control: {ok}/{len(amine)} human aminergic receptors carry D '
         f'(all four reference structures agree)',
         f'  cross-check: {agree}/{len(VERIFIED)} taxa shared with the structurally '
         f'anchored run in ../dop6_c24a8_6 return the same residue', '']


def consensus(names, c):
    res = [aln[n][c] for n in names if n in aln]
    cnt = Counter(res).most_common()
    top, k = cnt[0]
    return top, k / len(res), res


lines.append('position 3.32:')
for label, names in (('aminergic receptors (human)', amine),
                     ('melatonin clade (3 human + 5 zebrafish, 100/100)', MEL),
                     ('neuropeptide receptors (human)', NP),
                     ('PCDR-1 family (10 taxa)', [n for n in aln if n.startswith('PCDRfam_')])):
    top, frac, res = consensus(names, col)
    lines.append(f'  {label:52} {top}  ({frac:.0%} of {len(res)})  [{"".join(sorted(res))}]')
lines.append(f'  {"C. elegans PCDR-1 / F59D12.1":52} {aln[PCDR1][col]}')
lines.append('')

# --- the surrounding pocket window ---------------------------------------------
lines.append('TM3 pocket window (3.32 marked ^), consensus residue per column where '
             '>=60% agree, "." otherwise:')
lo, hi = col - 6, col + 9
hdr = ' ' * 54 + ''.join('^' if c == col else ' ' for c in range(lo, hi))
lines.append(hdr)
for label, names in (('aminergic receptors (human)', amine),
                     ('melatonin clade (verified)', MEL),
                     ('neuropeptide receptors (human)', NP),
                     ('PCDR-1 family (10 taxa)', [n for n in aln if n.startswith('PCDRfam_')])):
    s = ''
    for c in range(lo, hi):
        top, frac, _ = consensus(names, c)
        s += top if frac >= 0.6 else '.'
    lines.append(f'  {label:52}{s}')
lines.append(f'  {"C. elegans PCDR-1 / F59D12.1":52}{aln[PCDR1][lo:hi]}')
lines.append('')

# --- how many pocket columns does PCDR-1 share with each group? ----------------
lines.append('agreement with each group consensus across the window '
             '(columns where that group reaches >=60%):')
for label, names in (('aminergic', amine), ('melatonin', MEL),
                     ('neuropeptide', NP), ('PCDR-1 family', [n for n in aln if n.startswith('PCDRfam_')])):
    tot = hit = 0
    for c in range(lo, hi):
        top, frac, _ = consensus(names, c)
        if frac >= 0.6:
            tot += 1
            hit += (aln[PCDR1][c] == top)
    lines.append(f'  PCDR-1 matches {label:16} {hit}/{tot}')

txt = '\n'.join(lines)
print(txt)
open(f'{HERE}/output/pocket_compare.txt', 'w').write(txt + '\n')
