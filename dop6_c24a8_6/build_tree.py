#!/usr/bin/env python3
"""build_tree.py — adrenergic/aminergic GPCR tree including dop-6, and a test of whether the tree
can predict Galpha coupling at all.

Why the second part exists. Inferring a receptor's G protein from its nearest characterised
relative is how SER-5 was assigned Gs (Harris et al. 2009, J Neurosci 29:1446, PMID 19193891;
cited for that purpose in Koelle 2018 WormBook), and it is the obvious thing to do with DOP-6.
But it is homology-only inference, of exactly the kind that has to be validated rather than
assumed. Koelle's Table 1 gives experimentally determined coupling for 15 C. elegans receptors,
so the method can be tested on them: if coupling carries phylogenetic signal in this family, a
prediction for DOP-6 is worth stating as a hypothesis; if it does not, tree position says nothing
and the honest answer is that the question needs an assay.

Panel additions: the complete set of "poorly-studied putative small-molecule GPCRs" from Koelle
2018 Table 1 — T02E9.3 (= dop-5), C24A8.1 (= dop-6), C24A8.6, F35H10.10, F59D12.1, T21B4.4 —
since these are the uncharacterised candidates the project's receptor hunt is actually choosing
between.

Writes output/tree/aln/*, output/tree/trees/*, output/tree/coupling.tsv, output/tree/tree_report.txt
"""
import os, random, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__)); T = os.path.join(HERE, 'output', 'tree')
ALN, TRE = f'{T}/aln', f'{T}/trees'
os.makedirs(ALN, exist_ok=True); os.makedirs(TRE, exist_ok=True)
IQ = '/Users/mikeodonnell/git/receptor_phylogeny/scripts/bin/iqtree3-mt'
PEP = ('/private/tmp/claude-501/-Users-mikeodonnell-git-Malaiwong-Oglu-2026/'
       '93015e40-e331-495a-8362-20ac7ea07c75/scratchpad/cel_pep.fa')
THREADS = '4'
CONTACTS = [57, 120, 121, 124, 125, 128, 129, 131, 132, 133, 206, 210, 213,
            467, 470, 471, 474, 475, 534, 535]
OUTGROUP = ['hCHRM1', 'hCHRM3', 'cel_gar-3']
# Cja-dop-6: 326 aa of a 692 aa protein, 67 % gaps in the alignment — a fragment (see section 1
# of the phylogenetics checklist: full-length only; a fragment can attract to the wrong place
# because it is missing the divergent region that would have placed it correctly).
# dmTAR3: the original query returned the SAME accession as dmTAR2 (IQ-TREE flagged them as
# identical). Drosophila actually has three: Oct-TyrR/CG7485 P22270, TyrR/CG7431 Q9VEG1 and
# TyrRII/CG16766 Q9VEG2, so the third is added below under an accession-tagged label.
DROP = {'Cja-dop-6_CJA41806b.1', 'dmTAR3'}
RELABEL = {'dmTAR1': 'dmOctTyrR_P22270', 'dmTAR2': 'dmTyrRII_Q9VEG2'}
ADD_UP = {'dmTyrR_Q9VEG1': 'Q9VEG1'}
ADD = {'cel_F35H10.10': 'F35H10.10', 'cel_F59D12.1': 'F59D12.1', 'cel_T21B4.4': 'T21B4.4'}

# Galpha coupling, verbatim from Koelle 2018 WormBook Table 1 (NBK535514).
# '' = the table's coupling cell is empty, i.e. no coupling data exists.
COUPLING = {
 'cel_ser-1': 'Gq', 'cel_ser-4': 'Gi/o', 'cel_ser-5': 'Gs*', 'cel_ser-7': 'Gs',
 'cel_dop-1': 'Gq', 'cel_dop-2': 'Gi/o', 'cel_dop-3': 'Go', 'cel_dop-4': 'Gs',
 'cel_octr-1': 'Go', 'cel_ser-3': 'Gq', 'cel_ser-6': 'Gs',
 'cel_ser-2': 'Gi/o', 'cel_tyra-2': 'Gi/o', 'cel_tyra-3': 'Gq',
 'cel_gar-3': 'Gq+Gi/o',
 'cel_dop-5': '', 'cel_dop-6': '', 'cel_F35H10.10': '', 'cel_F59D12.1': '', 'cel_T21B4.4': '',
 # vertebrate/insect references, textbook coupling
 'hADRB1': 'Gs', 'hADRB2': 'Gs', 'hADRA2A': 'Gi/o', 'hADRA2C': 'Gi/o',
 'hADRA1A': 'Gq', 'hADRA1B': 'Gq', 'hADRA1D': 'Gq',
 'hDRD1': 'Gs', 'hDRD2': 'Gi/o', 'hHTR1A': 'Gi/o', 'hHTR2A': 'Gq',
 'hCHRM1': 'Gq', 'hCHRM3': 'Gq',
 'dmOAMB': 'Gq', 'dmOctA2R': 'Gi/o', 'dmOctB1R': 'Gs', 'dmOctB2R': 'Gs', 'dmOctB3R': 'Gs',
 'dmOctTyrR_P22270': 'Gi/o', 'dmTyrRII_Q9VEG2': 'Gi/o', 'dmTyrR_Q9VEG1': 'Gi/o',
 'dmDop1R1': 'Gs', 'dmDop2R': 'Gi/o',
}
rep = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); rep.append(s)

# ── panel ────────────────────────────────────────────────────────────────────
seqs = {}
for r in open(f'{T}/panel_curated.faa').read().split('>')[1:]:
    h, s = r.split('\n', 1); lab = h.strip().split('|')[0]
    if lab in DROP:
        say(f'excluded {lab}: fragment of a 692 aa protein'); continue
    seqs[RELABEL.get(lab, lab)] = s.replace('\n', '')
# the third Drosophila tyramine receptor, fetched by accession
import urllib.request as _u
for lab, acc in ADD_UP.items():
    with _u.urlopen(f'https://rest.uniprot.org/uniprotkb/{acc}.fasta', timeout=60) as _r:
        seqs[lab] = ''.join(_r.read().decode().split('\n')[1:])
    say(f'added {lab:<16} {acc:<16} {len(seqs[lab]):>4} aa  (third Drosophila tyramine receptor)')
# add the remaining Koelle orphans, longest isoform per gene, from the Ensembl proteome
pep = {}
for r in open(PEP).read().split('>')[1:]:
    h, s = r.split('\n', 1); pep[h.split()[0]] = (h, s.replace('\n', ''))
for lab, sn in ADD.items():
    cands = [(k, v) for k, v in pep.items() if re.match(rf'^{re.escape(sn)}([a-z])?\.\d+$', k)]
    if not cands: say(f'{lab}: no protein found for {sn}'); continue
    k, (h, s) = max(cands, key=lambda kv: len(kv[1][1]))
    seqs[lab] = s
    say(f'added {lab:<16} {k:<16} {len(s):>4} aa  (Koelle 2018 "poorly-studied putative small-molecule GPCR")')
fa = f'{ALN}/panel.faa'
with open(fa, 'w') as f:
    for lab in sorted(seqs):
        f.write(f'>{lab}\n')
        for i in range(0, len(seqs[lab]), 60): f.write(seqs[lab][i:i+60] + '\n')
say(f'\n{len(seqs)} sequences into the alignment')
with open(f'{T}/coupling.tsv', 'w') as f:
    f.write('tip\tgalpha_coupling\tsource\n')
    for lab in sorted(seqs):
        c = COUPLING.get(lab, '?')
        src = ('Koelle 2018 WormBook Table 1' if lab.startswith('cel_') else
               'textbook reference receptor' if c else '')
        f.write(f'{lab}\t{c}\t{src}\n')

# ── alignment + trimming ─────────────────────────────────────────────────────
aln = f'{ALN}/panel.linsi.faa'
say('running MAFFT L-INS-i ...')
with open(aln, 'w') as o:
    r = subprocess.run(['mafft', '--localpair', '--maxiterate', '1000', '--thread', THREADS,
                        '--quiet', fa], stdout=o, stderr=subprocess.PIPE)
if r.returncode: say('MAFFT failed'); sys.exit(1)
A = {}
for r_ in open(aln).read().split('>')[1:]:
    h, s = r_.split('\n', 1); A[h.strip()] = s.replace('\n', '')
W = len(next(iter(A.values())))
d6 = A['cel_dop-6']; col = {}; i = 0
for c, ch in enumerate(d6):
    if ch != '-': i += 1; col[i] = c
say(f'alignment {len(A)} x {W} columns\n')
say(f'trim level      columns   Galpha-contact columns surviving (of {len(CONTACTS)})')
trims = []
for name, args in [('gappyout', ['-gappyout']), ('gt0.1', ['-gt', '0.1']),
                   ('gt0.2', ['-gt', '0.2']), ('gt0.3', ['-gt', '0.3'])]:
    out = f'{ALN}/panel.{name}.faa'
    r = subprocess.run(['trimal', '-in', aln, '-out', out, '-colnumbering'] + args,
                       capture_output=True, text=True)
    if r.returncode: say(f'  {name}: trimal failed'); continue
    kept = set()
    for ln in r.stdout.split('\n'):
        if ln.strip().startswith('#ColumnsMap'):
            kept = {int(x) for x in ln.split('\t')[1].replace(' ', '').split(',') if x}
    B = {}
    for rr in open(out).read().split('>')[1:]:
        h, s = rr.split('\n', 1); B[h.strip()] = s.replace('\n', '')
    w = len(next(iter(B.values())))
    surv = sum(1 for p in CONTACTS if p in col and col[p] in kept)
    say(f'  {name:<12} {w:>7}   {surv}')
    trims.append((name, out, w))

# ── trees ────────────────────────────────────────────────────────────────────
say('\nIQ-TREE (ModelFinder, 1000 UFBoot, 1000 SH-aLRT), rooted on the muscarinic outgroup ...')
for name, path, w in [t for t in trims if t[0] == 'gt0.1'] + [t for t in trims if t[0] != 'gt0.1']:
    pre = f'{TRE}/{name}'
    cmd = [IQ, '-s', path, '-m', 'MFP', '-B', '1000', '-alrt', '1000', '-T', THREADS,
           '--prefix', pre, '-o', ','.join(OUTGROUP), '--quiet', '-redo']
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        say(f'  {name}: FAILED -> {(r.stderr or r.stdout)[-300:]}')
    else:
        model = ''
        for ln in open(pre + '.iqtree'):
            if ln.startswith('Best-fit model'): model = ln.split(':')[1].strip(); break
        say(f'  {name}: done, best model {model}')
open(f'{T}/tree_report.txt', 'w').write('\n'.join(rep) + '\n')
say('\n-> output/tree/trees/, output/tree/coupling.tsv')
