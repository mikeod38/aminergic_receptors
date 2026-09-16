#!/usr/bin/env python3
"""analyze_duplication.py — what exactly is duplicated between dop-6 (C24A8.1) and C24A8.6?

Reads output/locus_X_4318000-4334500.fa and output/gene_models.json (fetch_locus.py),
self-aligns the locus, and reports the repeat's boundaries in genomic, exonic and protein
coordinates. Writes output/repeat_structure.txt.
"""
import json, os, subprocess, sys
from Bio.Align import PairwiseAligner, substitution_matrices

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'output')
L0 = 4318000                                              # locus fasta starts here (1-based genomic)
seq = ''.join(l.strip() for l in open(f'{OUT}/locus_X_4318000-4334500.fa') if not l.startswith('>'))
G = lambda a, b: seq[a - L0:b - L0 + 1]                   # genomic 1-based inclusive slice
RC = lambda s: s[::-1].translate(str.maketrans('ACGTacgt', 'TGCAtgca'))
models = json.load(open(f'{OUT}/gene_models.json'))
rep = []

def say(*a):
    line = ' '.join(str(x) for x in a); print(line); rep.append(line)

# ── 1. repeat arms from the self-blast ───────────────────────────────────────
ARM_A = (4324178, 4328106)      # locus 6179-10107   (plus strand)
ARM_B = (4329303, 4333231)      # locus 11304-15232  (matches RC of arm A)
a, b = G(*ARM_A), RC(G(*ARM_B))
say('=' * 78); say('1. THE REPEAT')
say(f'arm A  X:{ARM_A[0]}-{ARM_A[1]}  {len(a)} bp')
say(f'arm B  X:{ARM_B[0]}-{ARM_B[1]}  {len(b)} bp   (reverse complement of arm A)')
mm = sum(1 for x, y in zip(a, b) if x != y)
say(f'identity {len(a) - mm}/{len(a)} = {100 * (len(a) - mm) / len(a):.4f} %   mismatches {mm}')
spacer = (ARM_A[1] + 1, ARM_B[0] - 1)
say(f'spacer X:{spacer[0]}-{spacer[1]}  {spacer[1] - spacer[0] + 1} bp   '
    f'(between the arms; contains C24A8.5)')
say(f'whole palindrome X:{ARM_A[0]}-{ARM_B[1]} = {ARM_B[1] - ARM_A[0] + 1} bp')
say('\njunction sequences (20 bp each side; | = repeat boundary)')
say(f'  arm A start  ...{G(ARM_A[0]-20, ARM_A[0]-1)}|{G(ARM_A[0], ARM_A[0]+19)}...')
say(f'  arm A end    ...{G(ARM_A[1]-19, ARM_A[1])}|{G(ARM_A[1]+1, ARM_A[1]+20)}...')
say(f'  arm B start  ...{G(ARM_B[0]-20, ARM_B[0]-1)}|{G(ARM_B[0], ARM_B[0]+19)}...')
say(f'  arm B end    ...{G(ARM_B[1]-19, ARM_B[1])}|{G(ARM_B[1]+1, ARM_B[1]+20)}...')

# is the spacer itself palindromic (i.e. does the hairpin extend with mismatches)?
sp = G(*spacer); half = len(sp) // 2
mmsp = sum(1 for x, y in zip(sp[:half], RC(sp)[:half]) if x != y)
say(f'spacer self-complementarity: {half - mmsp}/{half} bp match between its two halves '
    f'({100*(half-mmsp)/half:.1f} %) — {"palindromic too" if mmsp < half*0.2 else "not palindromic"}')

# ── 2. genes vs. repeat ──────────────────────────────────────────────────────
say('\n' + '=' * 78); say('2. WHICH PARTS OF EACH GENE SIT IN THE REPEAT')
def overlap(s, e, arm): return max(0, min(e, arm[1]) - max(s, arm[0]) + 1)
for name, g in models.items():
    inA, inB = overlap(g['start'], g['end'], ARM_A), overlap(g['start'], g['end'], ARM_B)
    tot = g['end'] - g['start'] + 1
    say(f"\n{name:<9} X:{g['start']}-{g['end']} ({tot} bp) strand {g['strand']:+d}  "
        f"in armA {inA} bp, in armB {inB} bp, unique {tot - inA - inB} bp")
    for t in g['transcripts']:
        if t['biotype'] != 'protein_coding': continue
        exs = t['exons'] if g['strand'] > 0 else t['exons'][::-1]
        say(f"  {t['id']} ({t['protein_length']} aa)")
        cum = 0
        for i, (es, ee) in enumerate(exs, 1):
            ln = ee - es + 1; oA, oB = overlap(es, ee, ARM_A), overlap(es, ee, ARM_B)
            tag = 'in repeat' if (oA == ln or oB == ln) else ('PARTLY in repeat' if (oA or oB) else 'unique')
            say(f"     exon {i:>2} X:{es}-{ee} {ln:>5} bp  {tag}")
            cum += ln

# ── 3. proteins ──────────────────────────────────────────────────────────────
say('\n' + '=' * 78); say('3. PROTEIN-LEVEL RELATIONSHIP')
prots = {}
for rec in open(f'{OUT}/dop6_c24a8_6_proteins.faa').read().split('>')[1:]:
    h, s = rec.split('\n', 1); prots[h.strip()] = s.replace('\n', '')
dop6 = [v for k, v in prots.items() if k.startswith('dop-6')][0]
say(f'dop-6 / C24A8.1  {len(dop6)} aa')
al = PairwiseAligner(); al.substitution_matrix = substitution_matrices.load('BLOSUM62')
al.open_gap_score, al.extend_gap_score, al.mode = -11, -1, 'local'
for k, v in sorted(prots.items()):
    if k.startswith('dop-6'): continue
    if not k.startswith('C24A8.6'): continue
    aln = al.align(dop6, v)[0]
    A, B = str(aln[0]), str(aln[1])
    idn = sum(1 for x, y in zip(A, B) if x == y and x != '-')
    ln = len(A)
    qs, qe = aln.aligned[0][0][0] + 1, aln.aligned[0][-1][-1]
    ss, se = aln.aligned[1][0][0] + 1, aln.aligned[1][-1][-1]
    say(f'\n{k}  {len(v)} aa')
    say(f'  local aln {ln} cols, {idn} identical = {100*idn/ln:.1f} %')
    say(f'  dop-6 {qs}-{qe}  vs  {k.split("|")[1]} {ss}-{se}')
    say(f'  -> {k.split("|")[1]} is {"an EXACT substring" if v[ss-1:se] == dop6[qs-1:qe] else "NOT an exact substring"} '
        f'of dop-6 over the aligned region')
    exact = v in dop6
    say(f'  -> full {k.split("|")[1]} protein is {"" if exact else "NOT "}an exact substring of the dop-6 protein')

open(f'{OUT}/repeat_structure.txt', 'w').write('\n'.join(rep) + '\n')
print(f'\n-> output/repeat_structure.txt')
