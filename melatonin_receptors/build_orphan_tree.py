#!/usr/bin/env python3
"""build_orphan_tree.py — is the GPR52/GPR21 family a real orthology group, and how old?

orphan_family.py found 93 sequences above the non-member floor across 21 metazoan
proteomes, PCDR-1 the single C. elegans one. But that floor is set by human ADRB2,
and other species' aminergic receptors can clear it, so the 93 are an upper bound
rather than a membership list. Only a tree separates a real member from aminergic
bleed-through, and it can only do so if the bleed-through has somewhere else to go -
hence the full curated aminergic backbone, which carries species-matched aminergic
receptors for the same lineages the hits came from.

Panel:
  - the b120 aminergic backbone, MINUS its nematode neuropeptide expansion
    (TEST_tkr-*, TEST_npr-* etc). Those are a nematode-specific radiation of
    neuropeptide receptors; they were in the melatonin panel only because that panel
    took every worm gene outscoring PCDR-1 against a melatonin query, and none of
    them comes near this family's floor.
  - all 93 above-floor hits, unfiltered. Membership is decided here, by the tree,
    not upstream by what a record happens to be named: "Probable GPR52" is not
    evidence for a sequence any more than "Oamb-like" is evidence against it.
  - PCDR-1.
  - the melatonin family: the 8 verified clade members plus the 2 lamprey orthologs.
  - 3 human neuropeptide receptors to root on.

The melatonin clade earns its place twice over. It is the sister-group context the
question needs, and its monophyly is established independently of this panel (100/100
in output/tree2, and both lamprey members are reciprocal best hits to MTNR1A), which
makes it the section 7.2 POSITIVE CONTROL this project had been missing: a group that
independent evidence says must be monophyletic. If it fragments here, this alignment
cannot adjudicate at this depth and nothing else the tree says should be believed.

Model fixed to LG+F+I+R7, as in build_tree.py, for comparability with that round.

Writes output/orphan/tree/.
"""
import glob, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DOP6 = os.path.join(HERE, '..', 'dop6_c24a8_6')
PROT = '/Users/mikeodonnell/git/receptor_phylogeny/analysis/intron_mapping_2026-09-02/proteomes'
IQ = '/Users/mikeodonnell/git/receptor_phylogeny/scripts/bin/iqtree3-mt'
OUT = os.path.join(HERE, 'output', 'orphan', 'tree')
os.makedirs(OUT, exist_ok=True)

MODEL = 'LG+F+I+R7'
NJOBS = 2
CORES = os.sysconf('SC_NPROCESSORS_ONLN')
THREADS = max(1, CORES // NJOBS)

NEUROPEPTIDE = {'P49146': 'hNPY2R', 'Q9Y5X5': 'hNPFFR2', 'P29371': 'hTACR3'}


def readfa(path):
    d, k = {}, None
    for line in open(path):
        if line.startswith('>'):
            k = line[1:].strip().split()[0]
            d[k] = ''
        elif k:
            d[k] += line.strip()
    return d


def from_proteome(species, accs):
    out, keep = {}, None
    hits = [p for p in sorted(glob.glob(f'{PROT}/{species}_UP*.fasta'))
            if os.path.getsize(p) > 0]
    if not hits:
        return out
    for line in open(hits[0]):
        if line.startswith('>'):
            acc = line.split('|')[1] if line.count('|') >= 2 else None
            keep = acc if acc in accs else None
            if keep:
                out[keep] = ''
        elif keep:
            out[keep] += line.strip()
    return out


panel, seen = {}, {}


def add(name, seq):
    """Add unless this exact sequence is already in the panel under another label."""
    if seq in seen:
        return False
    seen[seq] = name
    panel[name] = seq
    return True


# 1. aminergic backbone, minus the nematode neuropeptide expansion
n_back = 0
for name, seq in readfa(f'{DOP6}/output/tree/b120/panel.faa').items():
    if name.startswith('TEST_'):
        continue
    n_back += add(f'AM_{name}' if not name.startswith(('h', 'cel_', 'dm')) else name, seq)

# 2. every above-floor hit, unfiltered
n_orph = 0
for name, seq in readfa(f'{HERE}/output/orphan/members.faa').items():
    n_orph += add(f'ORPH_{name}', seq)

# 3. PCDR-1 (may already be in as the Cele orphan hit - add() dedupes by sequence)
pc = readfa(f'{HERE}/data/pcdr1.faa')
n_pcdr = add('cel_pcdr-1_F59D12.1', next(iter(pc.values())))

# 4. the melatonin family, including lamprey
n_mel = 0
for name, seq in readfa(f'{HERE}/data/mel_plus_lamprey.faa').items():
    # mel_plus_lamprey.faa already prefixes its zebrafish entries with MEL_; adding
    # another produced MEL_MEL_Drer_* and broke the positive-control label match.
    lab = name if name.startswith(('h', 'MEL_')) else f'MEL_{name}'
    n_mel += add(lab, seq)

# 5. root
n_np = 0
got = from_proteome('Hsap', set(NEUROPEPTIDE))
for acc, label in NEUROPEPTIDE.items():
    if acc in got:
        n_np += add(label, got[acc])
    else:
        print(f'  WARNING: {label} {acc} not found', file=sys.stderr)

print(f'panel: {len(panel)} taxa  ({n_back} aminergic backbone, {n_orph} above-floor '
      f'orphan-family hits, {n_pcdr} PCDR-1, {n_mel} melatonin family, {n_np} root)')
dropped = [n for n in ('cel_pcdr-1_F59D12.1',) if n not in panel]
if dropped:
    print(f'  (PCDR-1 already present as its orphan hit under '
          f'{seen[next(iter(pc.values()))]})')

raw = f'{OUT}/panel.faa'
with open(raw, 'w') as fh:
    for k, v in panel.items():
        fh.write(f'>{k}\n{v}\n')

aln = f'{OUT}/panel.linsi.faa'
if not os.path.exists(aln):
    print('aligning (mafft L-INS-i)...', flush=True)
    with open(aln, 'w') as fh:
        subprocess.run(['mafft', '--quiet', '--localpair', '--maxiterate', '1000',
                        '--thread', str(CORES), raw], stdout=fh, check=True)

procs = []
for tag, flags in {'gt0.1': ['-gt', '0.1'], 'gappyout': ['-gappyout']}.items():
    trimmed = f'{OUT}/panel.{tag}.faa'
    subprocess.run(['trimal', '-in', aln, '-out', trimmed, *flags], check=True)
    print(f'  {tag}: {len(next(iter(readfa(trimmed).values())))} columns', flush=True)
    procs.append(subprocess.Popen(
        [IQ, '-s', trimmed, '-m', MODEL, '-B', '1000', '--alrt', '1000',
         '-T', str(THREADS), '--prefix', f'{OUT}/{tag}', '-redo', '--quiet']))
for p in procs:
    p.wait()
print('trees written to', OUT)
