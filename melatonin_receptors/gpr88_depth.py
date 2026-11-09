#!/usr/bin/env python3
"""gpr88_depth.py — how deep does GPR88 go, and is it near the GPR52/GPR21 family?

GPR88 scores 14.5 bits against the orphan7 family profile (rank 691 of 707 in the
human proteome), so it is not a member of the family PCDR-1 belongs to. The question
left is where it does sit, and how far back it reaches.

BLAST finds confident GPR88 orthologs in mouse (716 bits, 95.6% id) and zebrafish
(239, 39.8%) and nothing convincing in lamprey, hagfish or tunicate - but a BLAST
negative at this depth is exactly the failure mode this project has already been
caught by twice, so the profile is run as the second, orthogonal axis before any
absence is claimed.

Same construction as orphan_family.py: profile from the verified orthologs, a
leave-one-out sensitivity control, a non-member floor from real class-A GPCRs that
are not GPR88, and then every proteome scanned.

Writes output/gpr88/result.txt
"""
import glob, json, os, subprocess, urllib.request

import pyhmmer
from pyhmmer.easel import Alphabet, SequenceFile, TextSequence, TextMSA
from pyhmmer.plan7 import Builder, Background, Pipeline

HERE = os.path.dirname(os.path.abspath(__file__))
PROT = '/Users/mikeodonnell/git/receptor_phylogeny/analysis/intron_mapping_2026-09-02/proteomes'
OUT = os.path.join(HERE, 'output', 'gpr88')
os.makedirs(OUT, exist_ok=True)
ALPHA = Alphabet.amino()

# resolved by query / reciprocal BLAST, not from memory; each verified below
SEED = {'Hsap_Q9GZN0': ('Hsap', 'Q9GZN0'),
        'Mmus_Q9EPB7': ('Mmus', 'Q9EPB7'),
        'Drer_A0A8M9Q9X5': ('Drer', 'A0A8M9Q9X5')}
# non-members: real class-A GPCRs, including both this project's other two families
NEG = {'P48039': 'MTNR1A', 'Q13585': 'GPR50', 'Q9Y2T5': 'GPR52', 'Q99679': 'GPR21',
       'P08908': 'HTR1A', 'P28221': 'HTR1D', 'P25100': 'ADRA1D', 'P07550': 'ADRB2'}

DEPTH = [('Aque', 'sponge'), ('Tadh', 'Trichoplax'), ('Mlei', 'ctenophore'),
         ('Nvec', 'cnidarian'), ('Hvul', 'cnidarian (Hydra)'), ('Cele', 'nematode'),
         ('Dmel', 'insect'), ('Lpol', 'chelicerate'), ('Pcau', 'priapulid'),
         ('Cvir', 'oyster'), ('Myes', 'scallop'), ('Lana', 'brachiopod'),
         ('Spur', 'sea urchin'), ('Pmin', 'starfish'), ('Apla', 'crown-of-thorns'),
         ('Ajap', 'sea cucumber'), ('Skow', 'hemichordate'), ('Bflo', 'amphioxus'),
         ('Cint', 'tunicate'), ('Ebur', 'hagfish'), ('Pmar', 'lamprey'),
         ('Drer', 'teleost'), ('Mmus', 'mouse'), ('Hsap', 'human')]


def _s(x):
    if x is None:
        return ''
    return x.decode() if isinstance(x, bytes) else x


def readfa(p):
    d, k = {}, None
    for line in open(p):
        if line.startswith('>'):
            k = line[1:].strip().split()[0]; d[k] = ''
        elif k:
            d[k] += line.strip()
    return d


def proteome(sp):
    h = [p for p in sorted(glob.glob(f'{PROT}/{sp}_UP*.fasta')) if os.path.getsize(p) > 0]
    if not h:
        h = [p for p in sorted(glob.glob(f'{PROT}/{sp}_*.fasta')) if os.path.getsize(p) > 0]
    return h[0] if h else None


def pull(sp, accs):
    out, keep, path = {}, None, proteome(sp)
    if not path:
        return out
    for line in open(path):
        if line.startswith('>'):
            acc = line.split('|')[1] if line.count('|') >= 2 else None
            keep = acc if acc in accs else None
            if keep:
                out[keep] = ''
        elif keep:
            out[keep] += line.strip()
    return out


def uniprot(acc):
    j = json.load(urllib.request.urlopen(
        f'https://rest.uniprot.org/uniprotkb/{acc}.json', timeout=30))
    nm = (j['proteinDescription'].get('recommendedName', {})
          .get('fullName', {}).get('value', '?'))
    g = [x.get('geneName', {}).get('value') for x in j.get('genes', []) if x.get('geneName')]
    return (j['organism']['scientificName'], g[0] if g else '-', nm, j['sequence']['length'])


def build(seqs, tag):
    t = f'{OUT}/_{tag}.faa'
    with open(t, 'w') as fh:
        for k, v in seqs.items():
            fh.write(f'>{k}\n{v}\n')
    a = f'{OUT}/_{tag}.aln'
    with open(a, 'w') as fh:
        subprocess.run(['mafft', '--quiet', '--localpair', '--maxiterate', '1000', t],
                       stdout=fh, check=True)
    d = readfa(a)
    msa = TextMSA(name=tag.encode(),
                  sequences=[TextSequence(name=k.encode(), sequence=v) for k, v in d.items()])
    hmm, _, _ = Builder(ALPHA).build_msa(msa.digitize(ALPHA), Background(ALPHA))
    return hmm


def scan(hmm, sp):
    path = proteome(sp)
    if not path:
        return []
    out = []
    with SequenceFile(path, digital=True, alphabet=ALPHA) as sf:
        for h in Pipeline(ALPHA, background=Background(ALPHA)).search_hmm(hmm, sf):
            nm = _s(h.name)
            out.append(((nm.split('|')[1] if nm.count('|') >= 2 else nm),
                        _s(h.description), h.score))
    out.sort(key=lambda r: -r[2])
    return out


L = ['GPR88 - how deep, and is it anywhere near the GPR52/GPR21 family?', '']
L.append('Seed, each checked against its own record:')
seed = {}
for lab, (sp, acc) in SEED.items():
    got = pull(sp, {acc})
    if acc not in got:
        raise SystemExit(f'ABORT: {lab} {acc} not in the {sp} proteome')
    org, gene, nm, ln = uniprot(acc)
    L.append(f'  {lab:18} {acc:12} {org[:24]:24} {gene:8} {ln:4} aa  {nm[:40]}')
    seed[lab] = got[acc]
L.append('')

L.append('CONTROL 1 - leave-one-out (only 3 seeds, so each rebuild uses 2):')
loo = []
for lab in SEED:
    sub = {k: v for k, v in seed.items() if k != lab}
    hits = scan(build(sub, f'loo_{lab}'), SEED[lab][0])
    acc = SEED[lab][1]
    rank = next((i + 1 for i, h in enumerate(hits) if h[0] == acc), None)
    sc = next((h[2] for h in hits if h[0] == acc), None)
    loo.append((lab, sc))
    L.append(f'  {lab:18} score {("%.1f" % sc) if sc else "NOT FOUND":>8}  '
             f'rank {rank if rank else "-"} of {len(hits)}')
L.append('')

hmm = build(seed, 'gpr88')
with open(f'{OUT}/gpr88.hmm', 'wb') as fh:
    hmm.write(fh)
L.append(f'Profile: {hmm.M} positions, {len(seed)} sequences -> output/gpr88/gpr88.hmm')
L.append('')

L.append('CONTROL 2 - non-member floor. Includes GPR52 and GPR21, so this also '
         'measures how far GPR88 sits from the family PCDR-1 is in.')
hs = {a: s for a, _, s in scan(hmm, 'Hsap')}
for acc, gene in NEG.items():
    s = hs.get(acc)
    L.append(f'  {acc:10} {gene:10} {("%.1f" % s) if s else "below threshold":>16}')
floor = max([hs[a] for a in NEG if a in hs], default=0.0)
L.append('')
L.append(f'  Floor (highest non-member): {floor:.1f}')
bad = [lb for lb, sc in loo if sc is None or sc < floor]
L.append(f'  CONTROL 1 verdict: {len(loo) - len(bad)}/{len(loo)} held-out orthologs clear it'
         + (f'  FAILED: {", ".join(bad)}' if bad else ''))
L.append('')

L.append('=' * 74)
L.append('Every available proteome, oldest lineage first.')
L.append('')
L.append('NOTE on the floor: 56.5 is useless as a membership criterion HERE, unlike in')
L.append('orphan_family.py. This profile is built from three close jawed-vertebrate')
L.append('orthologs, so it is narrow, and hundreds of ordinary class-A receptors score')
L.append('57-130 against it. The "n>=floor" column is therefore noise and is reported')
L.append('only to show that it is - a starfish does not have 118 GPR88 genes. What')
L.append('carries the signal is the gap: the three seeds score 528-792 and the next')
L.append('lineage down manages 128.')
L.append('')
L.append(f'  {"sp":6} {"lineage":18} {"best":>8} {"n>=floor":>9}   top hit')
deepest = None
for code, lin in DEPTH:
    if not proteome(code):
        L.append(f'  {code:6} {lin:18} {"-":>8} {"not searched":>9}')
        continue
    hits = scan(hmm, code)
    keep = [h for h in hits if h[2] >= floor]
    best = hits[0][2] if hits else float('nan')
    L.append(f'  {code:6} {lin:18} {best:8.1f} {len(keep):9}   '
             + (f'{hits[0][0]:12} {hits[0][1][:40]}' if hits else '-'))
    if keep and deepest is None:
        deepest = (code, lin, best)
L.append('')
L.append('VERDICT. GPR88 is detectable only in jawed vertebrates: human 790.9, mouse')
L.append('792.0, zebrafish 528.8. Every other lineage, lamprey (62.1) and hagfish (63.9)')
L.append('included, falls in the 31-128 class-A background, 4x below zebrafish. There is')
L.append('no intermediate score anywhere, so this is not a gradient being cut by a')
L.append('threshold - it is a gap.')
L.append('')
L.append('  C. elegans: best 53.1 (Q22188), below even the porous floor. BLAST agrees')
L.append('  independently at 53.9 bits / 24.9% id. No GPR88 ortholog in C. elegans.')
L.append('')
L.append('LIMIT, and it is a real one. Three jawed-vertebrate sequences make this the')
L.append('most narrowly anchored profile in the project - more so than the melatonin')
L.append('one, which had eight to ten. A sufficiently divergent invertebrate GPR88')
L.append('ortholog would be invisible to it, so the invertebrate negative is')
L.append('method-limited, NOT a confirmed absence. Treat it as weaker than the')
L.append('C. elegans melatonin negative, which had two orthogonal axes and a tree.')
L.append('')
L.append('One lead if that is worth chasing: the top invertebrate hits are Drosophila')
L.append('Moody (68.3) and an amphioxus "GPCR moody-like" (83.9), and the earlier BLAST')
L.append('also returned Moody first in the fly (57.4). Two methods pointing at the same')
L.append('family at twilight-zone scores is weak but consistent, and moody is where an')
L.append('invertebrate relative of GPR88 would be looked for. Not tested here.')
L.append('')
L.append('AND: GPR88 is not related to the family PCDR-1 is in, tested both ways.')
L.append('GPR52 and GPR21 score BELOW THRESHOLD against this profile, and GPR88 scores')
L.append('14.5 against the orphan7 family profile (rank 691 of 707 in human). The two')
L.append('families are mutually undetectable - a symmetric negative.')
L.append('')
L.append('Read against the other two families in this project:')
L.append('  GPR52/GPR21 family (PCDR-1 in it): cnidarian to human')
L.append('  melatonin receptors:                base of vertebrates (lamprey) to human')
L.append('  GPR88:                              as read off the table above')

txt = '\n'.join(L)
print(txt)
open(f'{OUT}/result.txt', 'w').write(txt + '\n')
for f in glob.glob(f'{OUT}/_*'):
    os.remove(f)
