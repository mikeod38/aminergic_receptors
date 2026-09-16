#!/usr/bin/env python3
"""gprotein_interface.py — could the duplicated C-terminal fragment (DOP-6 406-692, = C24A8.6)
engage a G protein on its own?

The fragment contains the end of ICL3, TM6, ECL3, TM7 and the whole cytoplasmic tail. Some of the
Galpha interface lives there; the rest does not. Rather than reasoning from textbook cartoons or
borrowing a paper's residue numbers, this defines the interface from a real receptor-G protein
complex directly:

  1. download the beta2-adrenergic receptor / Gs heterotrimer complex (PDB 3SN6) from RCSB
  2. compute every beta2AR residue with any atom within 4.5 A of any Galpha-s atom, and separately
     within 4.5 A of Gbeta, straight from the coordinates
  3. align beta2AR to DOP-6 and carry those contact positions across
  4. count how many land in DOP-6 1-405 (lost from the fragment) vs 406-692 (retained)

Also locates the class-A motifs (D/ERY, CWxP, NPxxY) and helix 8 in DOP-6 and says which side of
the 405/406 boundary each falls on.

Writes output/gprotein_interface.txt
"""
import gzip, io, os, urllib.request
from collections import defaultdict
from Bio.Align import PairwiseAligner, substitution_matrices

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'output')
WORK = ('/private/tmp/claude-501/-Users-mikeodonnell-git-Malaiwong-Oglu-2026/'
        '93015e40-e331-495a-8362-20ac7ea07c75/scratchpad')
SPLIT = 405
PDB = '3SN6'          # beta2AR-Gs; receptor chain R, Galpha-s chain A, Gbeta chain B
CUT = 4.5
AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V','MSE':'M'}
rep = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); rep.append(s)

# ── 1. structure ─────────────────────────────────────────────────────────────
path = f'{WORK}/{PDB}.pdb'
if not os.path.exists(path):
    url = f'https://files.rcsb.org/download/{PDB}.pdb.gz'
    with urllib.request.urlopen(url, timeout=120) as r: raw = r.read()
    open(path, 'w').write(gzip.decompress(raw).decode('utf-8', 'replace'))
atoms = defaultdict(list); seqres = defaultdict(list)
for ln in open(path):
    if ln.startswith(('ATOM', 'HETATM')):
        ch = ln[21]; rn = ln[17:20].strip()
        if rn not in AA3: continue
        try: ri = int(ln[22:26])
        except ValueError: continue
        xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        atoms[ch].append((ri, AA3[rn], xyz))
say('=' * 96)
say('CAN DOP-6 406-692 (= C24A8.6) BIND A G PROTEIN?')
say('=' * 96)
say(f'interface defined from {PDB} (beta2AR-Gs complex), contacts at <= {CUT} A, from coordinates')
say(f'chains present: {", ".join(f"{c}:{len(set(r for r,_,_ in v))} res" for c, v in sorted(atoms.items()))}\n')

REC, GA, GB = 'R', 'A', 'B'
if REC not in atoms or GA not in atoms:
    say('expected chains not found; aborting'); raise SystemExit

def contacts(chain_r, chain_p):
    P = [x for _, _, x in atoms[chain_p]]
    out = {}
    c2 = CUT * CUT
    # bucket the partner atoms on a 4.5 A grid for speed
    grid = defaultdict(list)
    for x in P: grid[(int(x[0]//CUT), int(x[1]//CUT), int(x[2]//CUT))].append(x)
    for ri, aa, (x, y, z) in atoms[chain_r]:
        if ri in out: continue
        gx, gy, gz = int(x//CUT), int(y//CUT), int(z//CUT)
        hit = False
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for (px, py, pz) in grid.get((gx+dx, gy+dy, gz+dz), ()):
                        if (x-px)**2 + (y-py)**2 + (z-pz)**2 <= c2: hit = True; break
                    if hit: break
                if hit: break
            if hit: break
        if hit: out[ri] = aa
    return out
ca = contacts(REC, GA); cb = contacts(REC, GB) if GB in atoms else {}
say(f'beta2AR residues contacting Galpha-s: {len(ca)}')
say(f'beta2AR residues contacting Gbeta   : {len(cb)}')

# receptor sequence as observed in the structure, with its numbering
robs = {}
for ri, aa, _ in atoms[REC]: robs[ri] = aa
rmin, rmax = min(robs), max(robs)
say(f'receptor chain observed residues {rmin}-{rmax} ({len(robs)} modelled; gaps are unmodelled loops)\n')

# ── 2. align beta2AR (full UniProt sequence) to DOP-6 ────────────────────────
import json, urllib.request as U
def up(acc):
    with U.urlopen(f'https://rest.uniprot.org/uniprotkb/{acc}.fasta', timeout=60) as r:
        return ''.join(r.read().decode().split('\n')[1:])
b2 = up('P07550')
dop6 = ''.join(open(f'{OUT}/dop6_c24a8_6_proteins.faa').read().split('>')[1].split('\n')[1:])
al = PairwiseAligner(); al.substitution_matrix = substitution_matrices.load('BLOSUM62')
al.open_gap_score, al.extend_gap_score, al.mode = -11, -1, 'global'
aln = al.align(b2, dop6)[0]; A, B = str(aln[0]), str(aln[1])
b2i = d6i = 0; m = {}
for x, y in zip(A, B):
    if x != '-': b2i += 1
    if y != '-': d6i += 1
    if x != '-' and y != '-': m[b2i] = (d6i, x, y)
say(f'beta2AR {len(b2)} aa vs DOP-6 {len(dop6)} aa: {len(m)} aligned positions')
say('NOTE: 3SN6 numbering is the beta2AR UniProt numbering for the receptor portion, but the')
say('      construct has an unmodelled ICL3; only modelled residues can contribute contacts.\n')

# ── 3. carry contacts across ─────────────────────────────────────────────────
say('-' * 96)
say('WHERE THE G-PROTEIN INTERFACE LANDS ON DOP-6')
say('-' * 96)
for lab, cs in (('Galpha-s', ca), ('Gbeta', cb)):
    if not cs: continue
    kept, lost, unmapped = [], [], []
    for ri in sorted(cs):
        if ri not in m: unmapped.append(ri); continue
        d6, x, y = m[ri]
        (kept if d6 > SPLIT else lost).append((ri, x, d6, y))
    tot = len(kept) + len(lost)
    say(f'\n{lab}: {len(cs)} contact residues, {tot} of them alignable to DOP-6 '
        f'({len(unmapped)} fall in unalignable stretches)')
    say(f'  in DOP-6 1-405  (LOST from the fragment): {len(lost)}')
    say(f'  in DOP-6 406-692 (RETAINED in C24A8.6)  : {len(kept)}')
    if tot: say(f'  -> the fragment retains {100*len(kept)/tot:.0f} % of the alignable interface')
    def fmt(v): return ', '.join(f'{x}{ri}->{y}{d6}' for ri, x, d6, y in v[:22])
    say(f'  lost    : {fmt(lost)}{" ..." if len(lost) > 22 else ""}')
    say(f'  retained: {fmt(kept)}{" ..." if len(kept) > 22 else ""}')

# ── 4. class-A motifs and helix 8 in DOP-6 ───────────────────────────────────
say('\n' + '-' * 96)
say('CLASS-A MOTIFS IN DOP-6 AND WHICH SIDE OF RESIDUE 405 THEY SIT ON')
say('-' * 96)
import re
MOT = [('D(E)RY  — TM3 cytoplasmic end, R3.50 is the central Galpha contact', r'D[RK][YLFH]'),
       ('CWxP    — TM6 rotamer toggle', r'CW.P'),
       ('NPxxY   — TM7 cytoplasmic end', r'N[PR]..Y'),
       ('helix 8 — amphipathic FxxxFxxxL-type run just after TM7', r'F..[KR]..[FL]')]
for name, pat in MOT:
    for mm in re.finditer(pat, dop6):
        p = mm.start() + 1
        say(f'  {name:<62} {mm.group()} at {p:>3}  ->  {"RETAINED (406-692)" if p > SPLIT else "LOST (1-405)"}')

say('\n' + '-' * 96)
say('READING')
say('-' * 96)
say('The alpha5 helix of Galpha inserts into a cavity whose walls are TM3, TM5, TM6 and ICL2/ICL3.')
say('A fragment starting at 406 brings TM6, TM7 and the tail, but not TM3 (with its D/ERY arginine),')
say('not ICL2, and not TM5 — and without TM1-5 there is no cavity for alpha5 to enter, because the')
say('cavity is formed by the outward swing of TM6 AWAY from TM3. Contact residues are necessary,')
say('not sufficient: they only work assembled into the seven-helix bundle.')
open(f'{OUT}/gprotein_interface.txt', 'w').write('\n'.join(rep) + '\n')
say('\n-> output/gprotein_interface.txt')
