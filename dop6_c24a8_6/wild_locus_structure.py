#!/usr/bin/env python3
"""wild_locus_structure.py — the dop-6 locus in each assembly, located from the gene itself.

Rather than trusting anchor arithmetic, this locates dop-6 exon 4 (single copy, outside the
palindrome) in each assembly, extracts a generous window around it in the gene's own orientation,
self-aligns that window to detect an inverted repeat, and reports how many copies of the 3,929 bp
repeat unit are present and at what identity. No identity threshold is applied, so a diverged
second copy would still be found.

Writes output/wild_locus_structure.txt
"""
import os, subprocess

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'output')
WORK = ('/private/tmp/claude-501/-Users-mikeodonnell-git-Malaiwong-Oglu-2026/'
        '93015e40-e331-495a-8362-20ac7ea07c75/scratchpad/asm')
L0 = 4318000
ARM = (4324178, 4328106)
n2 = ''.join(l.strip() for l in open(f'{OUT}/locus_X_4318000-4334500.fa') if not l.startswith('>')).upper()
G = lambda a, b: n2[a-L0:b-L0+1]
RC = lambda s: s[::-1].translate(str.maketrans('ACGTacgt', 'TGCAtgca'))
ASM = [('GCA_028201515_1', 'N2 (U. Oregon long-read)'), ('GCA_052911875_1', 'CGC1'),
       ('GCA_028201415_1', 'CB4856 (same pipeline as the N2 above)'), ('GCA_016989125_1', 'XZ1516'),
       ('GCA_016989505_1', 'DL238'), ('GCA_016989235_1', 'MY2147')]
UP, DOWN = 6000, 16000          # window upstream / downstream of dop-6 exon 4, in gene orientation
EX4 = G(4321367, 4321664)       # dop-6 exon 4: single copy, 5' of the palindrome
REPEAT = G(*ARM)

rep = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); rep.append(s)

open(f'{WORK}/ex4.fa', 'w').write(f'>ex4\n{EX4}\n')
say('=' * 98)
say('THE dop-6 LOCUS IN EACH ASSEMBLY, LOCATED FROM dop-6 EXON 4')
say('=' * 98)
say(f'repeat unit = N2 X:{ARM[0]}-{ARM[1]} ({len(REPEAT)} bp) = dop-6 exons 10-17.')
say('In N2 it is present twice, in inverted orientation (the second copy is annotated C24A8.6).')
say('Question for each assembly: how many copies, in what orientation, at what identity?\n')

for tag, desc in ASM:
    db = f'{WORK}/{tag}'
    if not os.path.exists(db + '.nsq'): say(f'{desc}: no db'); continue
    r = subprocess.run(['blastn', '-query', f'{WORK}/ex4.fa', '-db', db, '-dust', 'no',
                        '-evalue', '1e-30', '-outfmt', '6 sseqid sstart send pident length'],
                       capture_output=True, text=True)
    hs = [h.split('\t') for h in r.stdout.strip().split('\n') if h]
    hs = [h for h in hs if int(h[4]) > 250]
    if len(hs) != 1:
        say(f'{desc}: exon 4 found {len(hs)} times — skipping'); continue
    contig, s0, e0 = hs[0][0], int(hs[0][1]), int(hs[0][2])
    fwd = e0 > s0
    lo = (s0 - UP) if fwd else (s0 - DOWN)
    hi = (e0 + DOWN) if fwd else (e0 + UP)
    lo = max(1, lo)
    rr = subprocess.run(['seqkit', 'subseq', '--chr', contig, '-r', f'{lo}:{hi}', f'{WORK}/{tag}.fna'],
                        capture_output=True, text=True)
    win = ''.join(l for l in rr.stdout.split('\n')[1:] if not l.startswith('>')).upper()
    if not win: say(f'{desc}: extraction failed'); continue
    if not fwd: win = RC(win)                       # put dop-6 on the plus strand of the window
    fa = f'{WORK}/{tag}_win.fa'
    open(fa, 'w').write(f'>{tag}\n{win}\n')
    subprocess.run(['makeblastdb', '-in', fa, '-dbtype', 'nucl', '-out', f'{WORK}/{tag}_win'],
                   capture_output=True)
    # copies of the repeat unit
    open(f'{WORK}/rep.fa', 'w').write(f'>repeat\n{REPEAT}\n')
    r2 = subprocess.run(['blastn', '-query', f'{WORK}/rep.fa', '-db', f'{WORK}/{tag}_win', '-dust', 'no',
                         '-evalue', '1e-20', '-outfmt', '6 qstart qend sstart send pident length sstrand'],
                        capture_output=True, text=True)
    cps = [h.split('\t') for h in r2.stdout.strip().split('\n') if h and int(h.split('\t')[5]) > 500]
    # self-alignment: inverted repeat present?
    r3 = subprocess.run(['blastn', '-query', fa, '-db', f'{WORK}/{tag}_win', '-dust', 'no',
                         '-evalue', '1e-20', '-outfmt', '6 qstart qend sstart send pident length sstrand'],
                        capture_output=True, text=True)
    inv = [h.split('\t') for h in r3.stdout.strip().split('\n')
           if h and h.split('\t')[6] == 'minus' and int(h.split('\t')[5]) > 500]
    say(f'{desc}')
    say(f'  contig {contig}  window {lo}-{hi} ({len(win)} bp), dop-6 on the plus strand')
    say(f'  copies of the 3,929 bp repeat unit: {len(cps)}')
    for h in sorted(cps, key=lambda x: int(x[2])):
        say(f'     window {h[2]}-{h[3]}  {h[6]:<5}  {h[4]} % identity over {h[5]} bp')
    if inv:
        b = max(inv, key=lambda x: int(x[5]))
        say(f'  self-alignment: INVERTED REPEAT {b[0]}-{b[1]} vs {b[2]}-{b[3]}, '
            f'{b[4]} % over {b[5]} bp')
    else:
        say('  self-alignment: no inverted repeat > 500 bp')
    say('')

open(f'{OUT}/wild_locus_structure.txt', 'w').write('\n'.join(rep) + '\n')
say('-> output/wild_locus_structure.txt')
