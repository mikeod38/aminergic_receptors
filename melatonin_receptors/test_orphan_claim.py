#!/usr/bin/env python3
"""test_orphan_claim.py — formally test the claim this project put in writing.

The 178-taxon round (output/tree2) reported PCDR-1 grouping with hGPR52, hGPR21 and
five invertebrate receptors at 99.8/100, and that went into the README and the
collaborator brief. The 258-taxon round, which samples the family properly instead of
one member per species, does not reproduce it: PCDR-1's well-supported clade is
protostome-only (99.7/100), and the node that adds the human orphans is 36.2/91 at
gt0.1 and 3.2/52 at gappyout.

A weak support value is not the same as a rejected hypothesis, so this tests the
written claim directly: force the 8 taxa to be exclusively monophyletic, and ask
whether the data reject that tree against the unconstrained ML tree.

Construction follows the methodology's section 7:
  - 7.1 EXCLUSIVE, not partial. The constraint tree names every taxon in the panel,
    with the 8 as one group, so "monophyletic" means exclusive of the other 250. A
    partial constraint listing only the 8 constrains nothing (7.4a), and an earlier
    round of this project wasted a run learning that.
  - 7.2 CALIBRATION. The same constraint is applied to a control group whose answer
    is known: the 10 melatonin receptors, which are monophyletic at 100/100 in this
    very tree. Forcing a true group should cost ~0 log units. Without that, a large
    dLogL cannot be read as "rejected" rather than "this test rejects everything".
  - A second calibration at the other end: the 8-tip protostome clade that IS well
    supported here. Forcing it should also be near-free.

Writes output/orphan/constraint/
"""
import os, subprocess, re

HERE = os.path.dirname(os.path.abspath(__file__))
TREE = os.path.join(HERE, 'output', 'orphan', 'tree')
OUT = os.path.join(HERE, 'output', 'orphan', 'constraint')
IQ = '/Users/mikeodonnell/git/receptor_phylogeny/scripts/bin/iqtree3-mt'
os.makedirs(OUT, exist_ok=True)
MODEL = 'LG+F+I+R7'
THREADS = max(1, os.sysconf('SC_NPROCESSORS_ONLN') // 2)

# the claim as written: the 7 corroborated seed members plus PCDR-1
CLAIM = ['ORPH_Hsap_Q9Y2T5', 'ORPH_Hsap_Q99679', 'ORPH_Spur_A0A7M7GHW5',
         'ORPH_Dmel_Q9VZD1', 'ORPH_Pcau_A0ABM1EQA9', 'ORPH_Cvir_A0A8B8D717',
         'ORPH_Lana_A0A1S3KFU7', 'ORPH_Cele_O17899']
# calibration 1: known TRUE here (100/100 in this tree) - should be ~free to force
MEL = ['hMTNR1A', 'hMTNR1B', 'hGPR50', 'MEL_MEL_Drer_A0A8M1NRI5',
       'MEL_MEL_Drer_A0A8M9PN98', 'MEL_MEL_Drer_A2BGT9', 'MEL_MEL_Drer_P51046',
       'MEL_MEL_Drer_P51049', 'MEL_Pmar_A0AAJ7T8G6', 'MEL_Pmar_A0ACM8ASR7']
# calibration 2: the protostome clade that IS well supported (99.7/100)
PROTO = ['ORPH_Cele_O17899', 'ORPH_Dmel_Q9VZD1', 'ORPH_Pcau_A0ABM1EQA9',
         'ORPH_Cvir_A0A8B8D717', 'ORPH_Lana_A0A1S3KFU7', 'ORPH_Lpol_A0ABM1BHV8',
         'ORPH_Lpol_A0ABM1BZC6', 'ORPH_Lpol_A0ABM1S0K4']


def readfa(p):
    d, k = {}, None
    for line in open(p):
        if line.startswith('>'):
            k = line[1:].strip().split()[0]; d[k] = ''
        elif k:
            d[k] += line.strip()
    return d


for tag in ('gt0.1', 'gappyout'):
    aln = f'{TREE}/panel.{tag}.faa'
    ml = f'{TREE}/{tag}.treefile'
    if not (os.path.exists(aln) and os.path.exists(ml)):
        print(f'{tag}: inputs missing, skipped'); continue
    taxa = set(readfa(aln))
    print(f'\n=== {tag}  ({len(taxa)} taxa) ===')

    trees = [('ml', ml)]
    for label, group in (('claim', CLAIM), ('calib_mel', MEL), ('calib_proto', PROTO)):
        g = [t for t in group if t in taxa]
        if len(g) != len(group):
            missing = [t for t in group if t not in taxa]
            raise SystemExit(f'ABORT: {label} names taxa not in {tag}: {missing}')
        # EXCLUSIVE constraint: every taxon named, the group bracketed together
        rest = sorted(taxa - set(g))
        con = f'{OUT}/{tag}_{label}.constraint.tree'
        with open(con, 'w') as fh:
            fh.write('((' + ','.join(g) + '),' + ','.join(rest) + ');\n')
        pref = f'{OUT}/{tag}_{label}'
        done = f'{pref}.treefile'
        # Resumable: each constrained search is a full ML search on 258 taxa and takes
        # ~20 min, and this script has already been killed mid-run once. A finished
        # treefile is reused rather than recomputed.
        if os.path.exists(done) and os.path.getsize(done) > 0:
            print(f'  constrained search: {label} - already done, reusing', flush=True)
        else:
            print(f'  constrained search: {label} ({len(g)} taxa forced monophyletic)...',
                  flush=True)
            subprocess.run([IQ, '-s', aln, '-m', MODEL, '-g', con, '-T', str(THREADS),
                            '--prefix', pref, '-redo', '--quiet'], check=True)
        trees.append((label, done))

    # IQ-TREE writes <prefix>.trees itself, so the INPUT must not be named that -
    # an earlier round of this project truncated its own input to 0 bytes that way.
    cat = f'{OUT}/{tag}_candidates_input.trees'
    with open(cat, 'w') as fh:
        for _, p in trees:
            fh.write(open(p).read().strip() + '\n')
    au_done = f'{OUT}/{tag}_au.iqtree'
    if os.path.exists(au_done) and os.path.getsize(au_done) > 0:
        print(f'  AU test already done, reusing {au_done}', flush=True)
        continue
    print(f'  AU test over {len(trees)} trees: {", ".join(l for l, _ in trees)}',
          flush=True)
    subprocess.run([IQ, '-s', aln, '-m', MODEL, '-z', cat, '-zb', '10000', '-au',
                    '-n', '0', '-T', str(THREADS),
                    '--prefix', f'{OUT}/{tag}_au', '-redo', '--quiet'], check=True)
    print(f'  -> {OUT}/{tag}_au.iqtree')

print('\ndone')
