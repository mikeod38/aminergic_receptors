#!/usr/bin/env python3
"""haplotype_context.py — is duplication carriage in "wild" isolates real, N2 contamination,
or the chromosome-X selective sweep?

Carriers among the natural isolates could mean three different things, and they make different
predictions about divergence from N2 measured OUTSIDE the repeat:

  (1) genuine segregating polymorphism  -> carriers diverge from N2 normally both at the dop-6
                                           flanks and at unlinked control loci
  (2) N2 contamination / N2 escape      -> carrier is N2-like EVERYWHERE, including control loci
                                           on chromosomes the sweep does not affect
  (3) the global selective sweep        -> carriers are N2-like around dop-6 (chromosome X is a
                                           swept chromosome) but diverge normally on II and III,
                                           which the sweep largely spares

So: measure percent identity to the N2 reference in unique windows flanking dop-6, and at control
windows on II, III (largely unswept) and I, V (swept), for every assembly, and tabulate it against
carrier status. The repeat arms themselves are excluded — only single-copy sequence is compared.

Writes output/haplotype_context.tsv and output/haplotype_context.txt
"""
import json, os, subprocess, time, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'output')
WORK = ('/private/tmp/claude-501/-Users-mikeodonnell-git-Malaiwong-Oglu-2026/'
        '93015e40-e331-495a-8362-20ac7ea07c75/scratchpad/asm')
GENOME = ('/private/tmp/claude-501/-Users-mikeodonnell-git-Malaiwong-Oglu-2026/'
          '93015e40-e331-495a-8362-20ac7ea07c75/scratchpad/genome/ce.fa')
# windows of N2 reference sequence used as queries. 'swept' follows Andersen et al. 2012:
# chromosomes I, IV, V and X carry the global sweep; II and III largely do not.
WINDOWS = [
    ('dop6_left',   'X',  4304000, 4320000, 'X, 4-20 kb left of the palindrome (swept chr)'),
    ('dop6_right',  'X',  4334000, 4350000, 'X, just right of the palindrome (swept chr)'),
    ('X_far',       'X', 10000000, 10016000, 'X, 5.7 Mb away (swept chr)'),
    ('II_ctrl',    'II',  8000000,  8016000, 'chromosome II (sweep largely absent)'),
    ('III_ctrl',  'III',  8000000,  8016000, 'chromosome III (sweep largely absent)'),
    ('I_ctrl',      'I',  8000000,  8016000, 'chromosome I (swept chr)'),
    ('V_ctrl',      'V', 10000000, 10016000, 'chromosome V (swept chr)'),
]
rep = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); rep.append(s)

# pull the query windows out of the local WBcel235 fasta
qf = f'{WORK}/_hapwin.fa'
if not os.path.exists(qf):
    with open(qf, 'w') as f:
        for name, chrom, a, b, _ in WINDOWS:
            r = subprocess.run(['seqkit', 'subseq', '--chr', chrom, '-r', f'{a}:{b}', GENOME],
                               capture_output=True, text=True)
            s = ''.join(l for l in r.stdout.split('\n')[1:] if not l.startswith('>')).upper()
            if not s: say(f'WARNING: could not extract {name}'); continue
            f.write(f'>{name}\n{s}\n')
say(f'query windows: {", ".join(w[0] for w in WINDOWS)}\n')

carrier = {}
order = []
for ln in open(f'{OUT}/assembly_survey.tsv').read().strip().split('\n')[1:]:
    p = ln.split('\t')
    if len(p) >= 4 and p[3].isdigit():
        carrier[p[2]] = (p[0], p[1], int(p[3])); order.append(p[2])

say('=' * 118)
say('DIVERGENCE FROM N2 OUTSIDE THE REPEAT, BY CARRIER STATUS')
say('=' * 118)
for n, c, a, b, d in WINDOWS: say(f'  {n:<12} {c}:{a}-{b}  {d}')
say('')
hdr = f'{"strain":<20} {"kind":<5} {"cp":>3} ' + ' '.join(f'{w[0]:>11}' for w in WINDOWS)
say(hdr); say('-' * len(hdr))
rows = []
for acc in order:
    strain, kind, cp = carrier[acc]
    tag = acc.replace('.', '_')
    if not os.path.exists(f'{WORK}/{tag}.nsq'): continue
    r = subprocess.run(['blastn', '-query', qf, '-db', f'{WORK}/{tag}', '-dust', 'no',
                        '-evalue', '1e-50', '-max_target_seqs', '5',
                        '-outfmt', '6 qseqid pident length'], capture_output=True, text=True)
    hs = [x.split('\t') for x in r.stdout.strip().split('\n') if x]
    vals = []
    for name, chrom, a, b, _ in WINDOWS:
        sel = [h for h in hs if h[0] == name and int(h[2]) > 3000]
        if not sel: vals.append(None); continue
        # coverage-weighted identity over the best HSPs
        tot = sum(int(h[2]) for h in sel[:6])
        pid = sum(float(h[1])*int(h[2]) for h in sel[:6]) / tot
        vals.append(pid)
    say(f'{strain[:20]:<20} {kind:<5} {cp:>3} ' +
        ' '.join(f'{v:>11.3f}' if v is not None else f'{"-":>11}' for v in vals))
    rows.append((strain, kind, cp, vals))

# ── summary by carrier status ────────────────────────────────────────────────
say('')
say('-' * 118)
say('MEAN % IDENTITY TO N2, BY CARRIER STATUS (natural isolates only; lab strains excluded)')
say('-' * 118)
say(f'{"group":<26} {"n":>3} ' + ' '.join(f'{w[0]:>11}' for w in WINDOWS))
for lab, sel in [('carriers (2 copies)', [r for r in rows if r[1] == 'wild' and r[2] == 2]),
                 ('non-carriers (1 copy)', [r for r in rows if r[1] == 'wild' and r[2] == 1])]:
    if not sel: continue
    ms = []
    for i in range(len(WINDOWS)):
        v = [r[3][i] for r in sel if r[3][i] is not None]
        ms.append(sum(v)/len(v) if v else None)
    say(f'{lab:<26} {len(sel):>3} ' +
        ' '.join(f'{m:>11.3f}' if m is not None else f'{"-":>11}' for m in ms))
say('')
say('How to read it:')
say('  carriers N2-like at II/III control windows too -> N2 contamination or escape')
say('  carriers N2-like only on X (and other swept chromosomes) -> the duplication sits on the')
say('     swept haplotype and rose in frequency with it; carriage is real but not independent')
say('  carriers diverging normally everywhere -> an ordinary segregating structural polymorphism')

with open(f'{OUT}/haplotype_context.tsv', 'w') as f:
    f.write('strain\tkind\trepeat_copies\t' + '\t'.join(w[0] for w in WINDOWS) + '\n')
    for s, k, c, v in rows:
        f.write(f'{s}\t{k}\t{c}\t' + '\t'.join('' if x is None else f'{x:.4f}' for x in v) + '\n')
open(f'{OUT}/haplotype_context.txt', 'w').write('\n'.join(rep) + '\n')
say('\n-> output/haplotype_context.tsv')
