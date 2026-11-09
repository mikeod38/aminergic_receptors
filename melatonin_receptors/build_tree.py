#!/usr/bin/env python3
"""build_tree.py — does PCDR-1/F59D12.1 fall inside the melatonin receptor clade?

Panel, 118 taxa:
  - the curated aminergic backbone from ../dop6_c24a8_6 b120 (69 human/worm/fly
    reference receptors, hMTNR1A/1B/hGPR50 among them)
  - the 32 invertebrate melatonin reciprocal-best-hits found in this sub-project
  - PCDR-1 plus the 10 closest relatives of PCDR-1 across the proteome set
  - 6 human neuropeptide/orphan receptors, because PCDR-1's relatives return those
    as their own best human hits; without them the tree cannot place PCDR-1 where
    the BLAST evidence says it belongs

Model is fixed to the one b120 converged on (LG+F+I+R7) rather than re-selected -
this is an exploratory round, and model selection dominates runtime.

Writes output/tree/.
"""
import os, re, glob, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DOP6 = os.path.join(HERE, '..', 'dop6_c24a8_6')
PROT = '/Users/mikeodonnell/git/receptor_phylogeny/analysis/intron_mapping_2026-09-02/proteomes'
IQ   = '/Users/mikeodonnell/git/receptor_phylogeny/scripts/bin/iqtree3-mt'
OUT  = os.path.join(HERE, 'output', 'tree2')
os.makedirs(OUT, exist_ok=True)

MODEL = 'LG+F+I+R7'
NJOBS = 2
CORES = os.sysconf('SC_NPROCESSORS_ONLN')
THREADS = max(1, CORES // NJOBS)

# human neuropeptide/orphan receptors: the families PCDR-1's own relatives hit
NEUROPEPTIDE = {'P49146': 'hNPY2R', 'Q9Y5X5': 'hNPFFR2', 'P29371': 'hTACR3',
                'Q96P65': 'hQRFPR', 'Q9Y2T5': 'hGPR52', 'Q99679': 'hGPR21'}
PCDR_FAMILY = {'Dmel': 'A0A0B4K6Y4', 'Pcau': 'A0ABM1EQA9', 'Cvir': 'A0A8B8D717',
               'Lana': 'A0A1S3KFU7', 'Bflo': 'A0A9J7MKC7', 'Nvec': 'A7S5H8',
               'Drer': 'A0A8M1Q794', 'Spur': 'A0A7M7GHW5', 'Cint': 'H2XQ43',
               'Myes': 'A0A210QRP2'}
PCDR_FAMILY['Dmel'] = 'Q9VZD1'


def readfa(path):
    d, k = {}, None
    for line in open(path):
        if line.startswith('>'):
            k = line[1:].strip(); d[k] = ''
        elif k:
            d[k] += line.strip()
    return d


def from_proteome(species, accs):
    """Pull sequences by accession from one proteome, keyed by accession."""
    out, keep = {}, None
    for line in open(glob.glob(f'{PROT}/{species}_*.fasta')[0]):
        if line.startswith('>'):
            acc = line.split('|')[1] if '|' in line else None
            keep = acc if acc in accs else None
            if keep:
                out[keep] = ''
        elif keep:
            out[keep] += line.strip()
    return out


panel = {}

# 1. aminergic backbone, already curated and already named
b120 = readfa(f'{DOP6}/output/tree/b120/panel.faa')
for name, seq in b120.items():
    if re.match(r'^(h|cel_|dm)', name):
        panel[name] = seq
n_backbone = len(panel)

# 2. the melatonin RBH set from this sub-project (the 3 human refs are in b120)
for name, seq in readfa(f'{HERE}/data/mtnr_profile_seqs.faa').items():
    if not name.startswith('h'):
        panel['MEL_' + name] = seq
n_mel = len(panel) - n_backbone

# 3. PCDR-1 and its family
panel['cel_pcdr-1_F59D12.1'] = readfa(f'{HERE}/data/pcdr1.faa').popitem()[1]
for sp, acc in PCDR_FAMILY.items():
    got = from_proteome(sp, {acc})
    if acc in got:
        panel[f'PCDRfam_{sp}_{acc}'] = got[acc]
    else:
        print(f'  WARNING: {sp} {acc} not found in proteome', file=sys.stderr)
n_pcdr = len(panel) - n_backbone - n_mel

# 3b. every C. elegans protein that outscores PCDR-1 against any melatonin query.
#     Without these the tree cannot answer "is some other worm gene closer to the
#     melatonin receptors than PCDR-1 is" - it can only place PCDR-1 itself.
have = {v for v in panel.values()}
n_cel = 0
for name, seq in readfa(f'{HERE}/data/cele_candidates.faa').items():
    if seq not in have:
        panel[name] = seq
        n_cel += 1

# 4. neuropeptide / orphan outgroup
got = from_proteome('Hsap', set(NEUROPEPTIDE))
for acc, label in NEUROPEPTIDE.items():
    if acc in got:
        panel[label] = got[acc]
    else:
        print(f'  WARNING: {label} {acc} not found', file=sys.stderr)
n_np = len(panel) - n_backbone - n_mel - n_pcdr - n_cel

print(f'panel: {len(panel)} taxa  '
      f'({n_backbone} aminergic backbone, {n_mel} melatonin BLAST hits, '
      f'{n_pcdr} PCDR-1 nearest neighbours, {n_cel} C. elegans candidates, '
      f'{n_np} neuropeptide outgroup)')

raw = f'{OUT}/panel.faa'
with open(raw, 'w') as fh:
    for k, v in panel.items():
        fh.write(f'>{k}\n{v}\n')

# ---- align, trim, build -----------------------------------------------------
aln = f'{OUT}/panel.linsi.faa'
if not os.path.exists(aln):
    print('aligning (mafft L-INS-i)...')
    with open(aln, 'w') as fh:
        subprocess.run(['mafft', '--quiet', '--localpair', '--maxiterate', '1000',
                        '--thread', str(CORES), raw], stdout=fh, check=True)

TRIMS = {'gt0.1': ['-gt', '0.1'], 'gappyout': ['-gappyout']}
procs = []
for tag, flags in TRIMS.items():
    trimmed = f'{OUT}/panel.{tag}.faa'
    subprocess.run(['trimal', '-in', aln, '-out', trimmed, *flags], check=True)
    ncol = len(next(iter(readfa(trimmed).values())))
    print(f'  {tag}: {ncol} columns')
    procs.append(subprocess.Popen(
        [IQ, '-s', trimmed, '-m', MODEL, '-B', '1000', '--alrt', '1000',
         '-T', str(THREADS), '--prefix', f'{OUT}/{tag}', '-redo', '--quiet']))
for p in procs:
    p.wait()
print('trees written to', OUT)
