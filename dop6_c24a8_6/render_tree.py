#!/usr/bin/env python3
"""render_tree.py — draw the aminergic GPCR tree as annotated text.

Each tip carries its experimentally determined ligand and Galpha coupling, so DOP-6's placement can
be read against both. Nodes show SH-aLRT/UFboot support. Usage: render_tree.py [trim_level]
"""
import os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__)); T = os.path.join(HERE, 'output', 'tree')
trim = sys.argv[1] if len(sys.argv) > 1 else 'gt0.1'
tf = f'{T}/trees/{trim}.treefile'
LIG = {  # experimentally determined ligand; '' = orphan / no pharmacology
 'cel_dop-1': 'DA', 'cel_dop-2': 'DA', 'cel_dop-3': 'DA', 'cel_dop-4': 'DA',
 'cel_dop-5': 'ORPHAN', 'cel_dop-6': 'ORPHAN', 'cel_F35H10.10': 'ORPHAN',
 'cel_F59D12.1': 'ORPHAN', 'cel_T21B4.4': 'ORPHAN',
 'cel_ser-1': '5HT', 'cel_ser-4': '5HT', 'cel_ser-5': '5HT?', 'cel_ser-7': '5HT',
 'cel_ser-2': 'TA', 'cel_tyra-2': 'TA', 'cel_tyra-3': 'TA/OA',
 'cel_ser-3': 'OA', 'cel_ser-6': 'OA', 'cel_octr-1': 'OA', 'cel_gar-3': 'ACh',
 'hADRA1A': 'NE', 'hADRA1B': 'NE', 'hADRA1D': 'NE', 'hADRA2A': 'NE', 'hADRA2C': 'NE',
 'hADRB1': 'NE', 'hADRB2': 'NE', 'hDRD1': 'DA', 'hDRD2': 'DA',
 'hHTR1A': '5HT', 'hHTR2A': '5HT', 'hCHRM1': 'ACh', 'hCHRM3': 'ACh',
 'dmOAMB': 'OA', 'dmOctA2R': 'OA', 'dmOctB1R': 'OA', 'dmOctB2R': 'OA', 'dmOctB3R': 'OA',
 'dmOctTyrR_P22270': 'TA/OA', 'dmTyrR_Q9VEG1': 'TA', 'dmTyrRII_Q9VEG2': 'TA',
 'dmDop1R1': 'DA', 'dmDop2R': 'DA', 'bmOAR1': 'OA', 'plADRA1': 'NE',
}
COUP = {}
if os.path.exists(f'{T}/coupling.tsv'):
    for ln in open(f'{T}/coupling.tsv').read().strip().split('\n')[1:]:
        p = ln.split('\t')
        if len(p) > 1 and p[1].strip() and p[1] != '?': COUP[p[0]] = p[1]

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
            bl = 0.0
            if pos[0] < len(nw) and nw[pos[0]] == ':':
                m2 = re.match(r':([-\d.eE+]*)', nw[pos[0]:]); pos[0] += len(m2.group(0))
                bl = float(m2.group(1) or 0)
            return {'tip': False, 'lab': lab, 'kids': kids, 'bl': bl}
        m = re.match(r'[^,\)\(:]*', nw[pos[0]:]); lab = m.group(0); pos[0] += len(lab)
        bl = 0.0
        if pos[0] < len(nw) and nw[pos[0]] == ':':
            m2 = re.match(r':([-\d.eE+]*)', nw[pos[0]:]); pos[0] += len(m2.group(0))
            bl = float(m2.group(1) or 0)
        return {'tip': True, 'lab': lab, 'kids': [], 'bl': bl}
    return node()

if not os.path.exists(tf): print(f'no treefile for {trim}'); sys.exit(1)
tree = parse(open(tf).read())
out = []
def draw(n, pre='', last=True, root=False):
    conn = '' if root else ('└─ ' if last else '├─ ')
    if n['tip']:
        lab = n['lab']
        ann = []
        if LIG.get(lab): ann.append(LIG[lab])
        if COUP.get(lab): ann.append(COUP[lab])
        tag = f"  [{' | '.join(ann)}]" if ann else ''
        mark = '  <<<<<' if lab == 'cel_dop-6' else ''
        out.append(f'{pre}{conn}{lab}{tag}{mark}')
    else:
        sup = n['lab']
        out.append(f'{pre}{conn}' + (f'({sup})' if sup else '·'))
    if not n['tip']:
        np_ = pre + ('' if root else ('   ' if last else '│  '))
        for i, k in enumerate(n['kids']):
            draw(k, np_, i == len(n['kids'])-1)
draw(tree, root=True)
print(f'=== {trim}  ({len(open(tf).read())} bytes)   node labels = SH-aLRT / UFboot ===')
print('tips: [ligand | Galpha coupling]; coupling from Koelle 2018 WormBook Table 1\n')
print('\n'.join(out))
