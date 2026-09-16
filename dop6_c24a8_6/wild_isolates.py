#!/usr/bin/env python3
"""wild_isolates.py — do wild isolates carry the dop-6 palindrome?

Short reads cannot be placed inside a perfect inverted repeat, so a variant caller sees the two
arms as one ambiguous locus. That gives a readable signature: if the palindrome is present in the
wild isolates, both arms should be a variant desert relative to the unique sequence around them,
because the reads there map with MAPQ 0 and get filtered. This script profiles CaeNDR variant
density in 500 bp windows across X:4,314,000-4,338,000 and compares the arms with the unique
flanks and the spacer.

Source: CaeNDR release 20250625 hard-filter isotype VCF (S3, queried remotely with bcftools).
Writes output/wild_isolate_variants.tsv and output/wild_isolate_variants.txt
"""
import os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'output')
VCF = ('https://caendr-open-access-data-bucket.s3.us-east-2.amazonaws.com/dataset_release/'
       'c_elegans/20250625/variation/WI.20250625.hard-filter.isotype.vcf.gz')
REG = ('X', 4314000, 4338000); WIN = 500
ARM_A = (4324178, 4328106); ARM_B = (4329303, 4333231); SPACER = (4328107, 4329302)

rep = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); rep.append(s)

chrom, a, b = REG
cmd = ['bcftools', 'query', '-r', f'{chrom}:{a}-{b}', '-f', '%POS\t%REF\t%ALT\t%AC\t%AN\n', VCF]
say(f'querying CaeNDR 20250625 hard-filter isotype VCF over {chrom}:{a}-{b} ...')
try:
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
except subprocess.TimeoutExpired:
    say('bcftools timed out'); sys.exit(1)
if r.returncode != 0:
    say('bcftools failed:'); say(r.stderr[-800:]); sys.exit(1)
rows = [l.split('\t') for l in r.stdout.strip().split('\n') if l]
say(f'{len(rows)} variant sites returned\n')

def label(p):
    if ARM_A[0] <= p <= ARM_A[1]: return 'armA'
    if ARM_B[0] <= p <= ARM_B[1]: return 'armB'
    if SPACER[0] <= p <= SPACER[1]: return 'spacer'
    return 'unique_flank'
from collections import Counter
cnt = Counter(label(int(x[0])) for x in rows)
span = {'armA': ARM_A[1]-ARM_A[0]+1, 'armB': ARM_B[1]-ARM_B[0]+1,
        'spacer': SPACER[1]-SPACER[0]+1,
        'unique_flank': (b-a+1) - (ARM_A[1]-ARM_A[0]+1) - (ARM_B[1]-ARM_B[0]+1) - (SPACER[1]-SPACER[0]+1)}
say(f'{"region":<14} {"bp":>7} {"variants":>9} {"per kb":>8}')
say('-' * 42)
for k in ['unique_flank', 'armA', 'spacer', 'armB']:
    n = cnt.get(k, 0); say(f'{k:<14} {span[k]:>7} {n:>9} {1000*n/span[k]:>8.2f}')

say('\n500 bp windows (| = one variant, capped at 40):')
with open(f'{OUT}/wild_isolate_variants.tsv', 'w') as f:
    f.write('win_start\twin_end\tregion\tn_variants\tper_kb\n')
    for w in range(a, b, WIN):
        n = sum(1 for x in rows if w <= int(x[0]) < w+WIN)
        lab = label(w + WIN//2)
        f.write(f'{w}\t{w+WIN-1}\t{lab}\t{n}\t{2*n:.1f}\n')
        say(f'  X:{w}-{w+WIN-1} {lab:<13} {n:>3} {"|"*min(n,40)}')
open(f'{OUT}/wild_isolate_variants.txt', 'w').write('\n'.join(rep) + '\n')
say('\n-> output/wild_isolate_variants.tsv')
