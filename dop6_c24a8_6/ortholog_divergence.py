#!/usr/bin/env python3
"""ortholog_divergence.py — where along DOP-6 does cross-species divergence actually sit?

The palindrome duplicates DOP-6 residues 406-692. This script asks whether that segment
diverges from other Caenorhabditis at the same rate as the unique N-terminal half
(residues 1-405), at protein AND nucleotide level. Synonymous sites cannot be conserved
across 20-40 Myr, so a nucleotide-identical segment is diagnostic of contamination or a
mis-assembly rather than of selection.

Writes output/ortholog_divergence.txt
"""
import json, os, time, urllib.request
from Bio.Align import PairwiseAligner, substitution_matrices

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'output')
ENS = 'https://rest.ensembl.org'
BOUND = 405                      # last residue of the unique (non-duplicated) half
ORTH = [('Cbr-dop-6', 'caenorhabditis_briggsae', 'CBG14393b.1'),
        ('Cbn-dop-6', 'caenorhabditis_brenneri', 'CBN19707.1'),
        ('Cre-dop-6', 'caenorhabditis_remanei_gca010183535v1rs', 'XM_053735786.1'),
        ('Ppa-dop-6', 'pristionchus_pacificus', 'PPA00437.1')]

def get(u, c='application/json', tries=6):
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(u, headers={'Content-Type': c}), timeout=90) as r:
                t = r.read().decode()
            return json.loads(t) if c == 'application/json' else t.strip()
        except Exception:
            if k == tries - 1: raise
            time.sleep(2 * (k + 1))

rep = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); rep.append(s)

pa = PairwiseAligner(); pa.substitution_matrix = substitution_matrices.load('BLOSUM62')
pa.open_gap_score, pa.extend_gap_score, pa.mode = -11, -1, 'global'
na = PairwiseAligner(); na.match_score, na.mismatch_score = 2, -3
na.open_gap_score, na.extend_gap_score, na.mode = -10, -0.5, 'global'

seqs = {}
for r in open(f'{OUT}/dop6_orthologs.faa').read().split('>')[1:]:
    h, s = r.split('\n', 1); seqs[h.strip().split('|')[0]] = s.replace('\n', '')
cel = seqs['Cel-dop-6']
cel_cds = get(f'{ENS}/sequence/id/C24A8.1.1?type=cds;content-type=text/plain', 'text/plain')
say(f'Cel DOP-6 {len(cel)} aa; CDS {len(cel_cds)} nt (3x{len(cel)}+3 = {3*len(cel)+3})')
say(f'\nsplit at residue {BOUND}: N-half 1-{BOUND} = unique; C-half {BOUND+1}-692 = the segment '
    f'duplicated as C24A8.6\n')
say(f"{'ortholog':<12} {'':<6} {'protein':<26} {'nucleotide (CDS)':<26}")
say(f"{'':<12} {'half':<6} {'ident/aln':>12} {'%id':>6}   {'ident/aln':>12} {'%id':>6}")
say('-' * 78)

for label, sp, tid in ORTH:
    if label not in seqs: say(f'{label}: protein missing'); continue
    other = seqs[label]
    try:
        ocds = get(f'{ENS}/sequence/id/{tid}?type=cds;content-type=text/plain', 'text/plain')
    except Exception as e:
        ocds = None; say(f'{label}: CDS fetch failed ({e})')
    # protein alignment -> map Cel residue index to alignment column
    x = pa.align(cel, other)[0]; A, B = str(x[0]), str(x[1])
    ci = 0; halves = {'N': [0, 0], 'C': [0, 0]}
    for p, q in zip(A, B):
        if p != '-': ci += 1
        if p == '-' or q == '-': continue
        h = 'N' if ci <= BOUND else 'C'
        halves[h][1] += 1
        if p == q: halves[h][0] += 1
    # nucleotide alignment on the CDS, split at the same codon boundary
    ncell = {}
    if ocds:
        celN, celC = cel_cds[:3*BOUND], cel_cds[3*BOUND:]
        # align full CDS then split by Cel coordinate
        y = na.align(cel_cds, ocds)[0]; C, D = str(y[0]), str(y[1])
        ni = 0; nh = {'N': [0, 0], 'C': [0, 0]}
        for p, q in zip(C, D):
            if p != '-': ni += 1
            if p == '-' or q == '-': continue
            h = 'N' if ni <= 3*BOUND else 'C'
            nh[h][1] += 1
            if p == q: nh[h][0] += 1
        ncell = nh
    for h in ('N', 'C'):
        i, n = halves[h]
        if ncell:
            ni_, nn = ncell[h]
            say(f'{label if h=="N" else "":<12} {h:<6} {i}/{n:<8} {100*i/n if n else 0:>6.1f}   '
                f'{ni_}/{nn:<8} {100*ni_/nn if nn else 0:>6.1f}')
        else:
            say(f'{label if h=="N" else "":<12} {h:<6} {i}/{n:<8} {100*i/n if n else 0:>6.1f}   '
                f'{"n/a":>12} {"":>6}')
    say('')

say('Interpretation: synonymous sites saturate over Caenorhabditis divergence times, so a')
say('C-half nucleotide identity near 100 % cannot be produced by purifying selection.')
open(f'{OUT}/ortholog_divergence.txt', 'w').write('\n'.join(rep) + '\n')
say('-> output/ortholog_divergence.txt')
