#!/usr/bin/env python3
"""haplotype_extent.py — how far does the N2-like haplotype carrying the dop-6 duplication extend?

haplotype_context.py showed carriers are far more N2-like than non-carriers immediately around
dop-6 (99.98 % vs 99.48 % over 16 kb) but not at a window 5.7 Mb away, and not on chromosomes II
and III. That rules out whole-genome N2 contamination and says carriage is identity-by-descent on a
shared haplotype. What is left is the length of that haplotype, which separates the two remaining
readings:

  * a RECENT laboratory introgression of N2 sequence into strain collections would leave a long,
    barely recombined block — hundreds of kb to megabases
  * an OLDER segregating haplotype that rose in frequency would be eroded by recombination to a
    shorter block

So: scan windows outward from dop-6 along chromosome X, measure identity to N2 in carriers and in
non-carriers separately, and find where the two converge.

Important caveat carried into the output: C. elegans is a selfer, so effective recombination is
very low and haplotype blocks persist far longer than in an outcrosser. Block length therefore
bounds the event loosely and must not be turned into a confident date.

Writes output/haplotype_extent.tsv and .txt
"""
import os, subprocess

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'output')
WORK = ('/private/tmp/claude-501/-Users-mikeodonnell-git-Malaiwong-Oglu-2026/'
        '93015e40-e331-495a-8362-20ac7ea07c75/scratchpad/asm')
GENOME = ('/private/tmp/claude-501/-Users-mikeodonnell-git-Malaiwong-Oglu-2026/'
          '93015e40-e331-495a-8362-20ac7ea07c75/scratchpad/genome/ce.fa')
DOP6 = 4324000           # the palindrome sits at X:4,324,178-4,333,231
WIN = 10000
# offsets from the locus, negative = leftward, in bp
OFFSETS = [-3000000, -1000000, -500000, -200000, -100000, -50000, -20000,
           20000, 50000, 100000, 200000, 500000, 1000000, 3000000, 6000000]
rep = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); rep.append(s)

qf = f'{WORK}/_extent.fa'
wins = []
with open(qf, 'w') as f:
    for off in OFFSETS:
        a = DOP6 + off
        if a < 1 or a > 17_600_000 - WIN: continue
        r = subprocess.run(['seqkit', 'subseq', '--chr', 'X', '-r', f'{a}:{a+WIN}', GENOME],
                           capture_output=True, text=True)
        s = ''.join(l for l in r.stdout.split('\n')[1:] if not l.startswith('>')).upper()
        if len(s) < WIN * 0.9: continue
        name = f'off{off:+d}'
        f.write(f'>{name}\n{s}\n'); wins.append((name, off, a))
say(f'{len(wins)} windows of {WIN} bp along chromosome X, offsets relative to X:{DOP6}\n')

carrier = {}
for ln in open(f'{OUT}/assembly_survey.tsv').read().strip().split('\n')[1:]:
    p = ln.split('\t')
    if len(p) >= 4 and p[3].isdigit(): carrier[p[2]] = (p[0], p[1], int(p[3]))

vals = {n: {1: [], 2: []} for n, _, _ in wins}
for acc, (strain, kind, cp) in carrier.items():
    if kind != 'wild': continue
    tag = acc.replace('.', '_')
    if not os.path.exists(f'{WORK}/{tag}.nsq'): continue
    r = subprocess.run(['blastn', '-query', qf, '-db', f'{WORK}/{tag}', '-dust', 'no',
                        '-evalue', '1e-50', '-outfmt', '6 qseqid pident length'],
                       capture_output=True, text=True)
    hs = [x.split('\t') for x in r.stdout.strip().split('\n') if x]
    for name, off, a in wins:
        sel = [h for h in hs if h[0] == name and int(h[2]) > 2000]
        if not sel: continue
        tot = sum(int(h[2]) for h in sel[:6])
        pid = sum(float(h[1])*int(h[2]) for h in sel[:6]) / tot
        vals[name][cp].append(pid)

say(f'{"offset from dop-6":>20} {"carriers":>18} {"non-carriers":>18} {"difference":>12}')
say('-' * 74)
rows = []
for name, off, a in sorted(wins, key=lambda w: w[1]):
    c = vals[name][2]; nc = vals[name][1]
    if not c or not nc: continue
    mc, mn = sum(c)/len(c), sum(nc)/len(nc)
    lab = f'{off/1000:+.0f} kb' if abs(off) < 1_000_000 else f'{off/1_000_000:+.1f} Mb'
    star = '  <-- still N2-like' if mc - mn > 0.15 else ''
    say(f'{lab:>20} {mc:>17.3f}% {mn:>17.3f}% {mc-mn:>11.3f}{star}')
    rows.append((off, a, len(c), mc, len(nc), mn, mc-mn))
say('')
say('at the locus itself (from haplotype_context.py): carriers 99.977 %, non-carriers 99.477 %,')
say('difference 0.500 over the 16 kb window immediately left of the palindrome')
say('')
big = [r for r in rows if r[6] > 0.15]
if big:
    lo = min(r[0] for r in big); hi = max(r[0] for r in big)
    say(f'the carrier-specific N2-like haplotype is still detectable from {lo/1000:+.0f} kb to '
        f'{hi/1000:+.0f} kb around the locus')
say('')
say('CAVEAT: C. elegans reproduces mainly by selfing, so effective recombination is very low and')
say('haplotype blocks persist far longer than in an outcrossing species. Block length bounds the')
say('event only loosely and should not be converted into a date.')
with open(f'{OUT}/haplotype_extent.tsv', 'w') as f:
    f.write('offset\tX_start\tn_carriers\tcarrier_pid\tn_noncarriers\tnoncarrier_pid\tdifference\n')
    for r in rows: f.write('\t'.join(str(x) for x in r) + '\n')
open(f'{OUT}/haplotype_extent.txt', 'w').write('\n'.join(rep) + '\n')
say('-> output/haplotype_extent.tsv')
