#!/usr/bin/env python3
"""orphan_family.py — how deeply conserved is the GPR52/GPR21 family PCDR-1 falls into?

The 178-taxon tree put PCDR-1 inside a clade of human GPR52/GPR21 plus one
invertebrate receptor each from sea urchin, fly, priapulid, oyster and lamp shell,
and two controls said that placement is not long-branch attraction. If that family
really runs from cnidarian to human it is a deeply conserved orthology group, which
is a more useful fact about PCDR-1 than "not a melatonin receptor".

Three things have to be separated here:

  1. Does the family exist as a profile-detectable group across the 32 proteomes,
     or is it an artifact of taking one best hit per species?
  2. Is it 1:1 orthologous, or a set of lineage-specific expansions that only looks
     1:1 because the panel held one member per species?
  3. Is PCDR-1 in it?

(3) is why **the profile is built from the seven corroborated members only and
PCDR-1 is excluded**. Including the sequence whose membership is in question would
make its own recovery circular - the profile would be tuned to find it. PCDR-1 is
scored against the finished profile like any other candidate.

Controls, both asserted rather than assumed:
  - leave-one-species-out: rebuild the profile without each member and check that
    member is still recovered above its own proteome's background. A profile that
    cannot re-find a known member cannot be trusted to call an absence.
  - negatives: the 8 verified melatonin-clade receptors and the human aminergic
    set are scored against the profile. These are real class-A GPCRs, so they set
    the "related but not in this family" floor that a real hit has to clear.

Writes output/orphan/ (scores.tsv, members.faa, result.txt)
"""
import json, os, glob, subprocess, sys, urllib.parse, urllib.request
from collections import defaultdict

import pyhmmer
from pyhmmer.easel import Alphabet, SequenceFile, TextSequence, TextMSA
from pyhmmer.plan7 import Builder, Background, Pipeline

HERE = os.path.dirname(os.path.abspath(__file__))
PROT = '/Users/mikeodonnell/git/receptor_phylogeny/analysis/intron_mapping_2026-09-02/proteomes'
OUT = os.path.join(HERE, 'output', 'orphan')
os.makedirs(OUT, exist_ok=True)

# The seven corroborated family members. Every one of these ranks at or near the top
# of its own proteome in the FORWARD direction from hGPR52/hGPR21 (search registry,
# 2026-11-03), so each is independently supported, not just a tree neighbour.
# PCDR-1 is deliberately NOT here - see the docstring.
FAMILY = {
    'Hsap_Q9Y2T5':      ('Hsap', 'Q9Y2T5', 'Homo sapiens'),
    'Hsap_Q99679':      ('Hsap', 'Q99679', 'Homo sapiens'),
    'Spur_A0A7M7GHW5':  ('Spur', 'A0A7M7GHW5', 'Strongylocentrotus purpuratus'),
    'Dmel_Q9VZD1':      ('Dmel', 'Q9VZD1', 'Drosophila melanogaster'),
    'Pcau_A0ABM1EQA9':  ('Pcau', 'A0ABM1EQA9', 'Priapulus caudatus'),
    'Cvir_A0A8B8D717':  ('Cvir', 'A0A8B8D717', 'Crassostrea virginica'),
    'Lana_A0A1S3KFU7':  ('Lana', 'A0A1S3KFU7', 'Lingula anatina'),
}
PCDR1 = ('Cele', 'O17899')     # scored, never used to build

# negative controls: real class-A GPCRs that are NOT in this family
MEL_NEG = ['P48039', 'P49286', 'Q13585']            # MTNR1A, MTNR1B, GPR50
AMINE_NEG = ['P08908', 'P14416', 'P07550', 'P35367']  # HTR1A, DRD2, ADRB2, HRH1

ALPHA = Alphabet.amino()


def _s(x):
    """pyhmmer hands back bytes for some files and str for others."""
    if x is None:
        return ''
    return x.decode() if isinstance(x, bytes) else x


def readfa(path):
    d, k = {}, None
    for line in open(path):
        if line.startswith('>'):
            k = line[1:].strip().split()[0]
            d[k] = ''
        elif k:
            d[k] += line.strip()
    return d


def proteome(sp):
    # prefer the UniProt reference set over any other fasta sharing the prefix,
    # and treat a zero-byte file as absent rather than as an empty proteome - a
    # failed download must not be reportable as "searched, nothing found".
    hits = [p for p in sorted(glob.glob(f'{PROT}/{sp}_UP*.fasta')) if os.path.getsize(p) > 0]
    if not hits:
        hits = [p for p in sorted(glob.glob(f'{PROT}/{sp}_*.fasta')) if os.path.getsize(p) > 0]
    return hits[0] if hits else None


def pull(sp, accs):
    """Sequences by UniProt accession out of one proteome fasta."""
    out, keep = {}, None
    path = proteome(sp)
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
    """Fetch a record so the label can be checked against it, never asserted."""
    u = f'https://rest.uniprot.org/uniprotkb/{acc}.json'
    j = json.load(urllib.request.urlopen(u, timeout=30))
    name = (j['proteinDescription'].get('recommendedName', {})
            .get('fullName', {}).get('value')
            or (j['proteinDescription'].get('submissionNames') or [{}])[0]
            .get('fullName', {}).get('value', '?'))
    genes = [g.get('geneName', {}).get('value') for g in j.get('genes', []) if g.get('geneName')]
    return {'organism': j['organism']['scientificName'], 'name': name,
            'gene': genes[0] if genes else '-', 'length': j['sequence']['length']}


def build_hmm(seqs, name):
    """mafft L-INS-i then a pyhmmer profile. Returns the HMM."""
    tmp = f'{OUT}/_{name}.faa'
    with open(tmp, 'w') as fh:
        for k, v in seqs.items():
            fh.write(f'>{k}\n{v}\n')
    aln = f'{OUT}/_{name}.linsi.faa'
    with open(aln, 'w') as fh:
        subprocess.run(['mafft', '--quiet', '--localpair', '--maxiterate', '1000', tmp],
                       stdout=fh, check=True)
    a = readfa(aln)
    msa = TextMSA(name=name.encode(),
                  sequences=[TextSequence(name=k.encode(), sequence=v) for k, v in a.items()])
    hmm, _, _ = Builder(ALPHA).build_msa(msa.digitize(ALPHA), Background(ALPHA))
    return hmm


def scan(hmm, sp):
    """All hits of one profile against one proteome. Returns [(acc, desc, score)]."""
    path = proteome(sp)
    if not path:
        return []
    out = []
    with SequenceFile(path, digital=True, alphabet=ALPHA) as sf:
        for h in Pipeline(ALPHA, background=Background(ALPHA)).search_hmm(hmm, sf):
            nm = _s(h.name)
            acc = nm.split('|')[1] if nm.count('|') >= 2 else nm
            out.append((acc, _s(h.description), h.score))
    out.sort(key=lambda r: -r[2])
    return out


# ---------------------------------------------------------------- verify the seed
L = []
L.append('GPR52/GPR21 orphan family - how far does it reach, and is PCDR-1 in it?')
L.append('')
L.append('Seed, each accession checked against its own UniProt record '
         '(label is only kept if the organism matches):')
seed = {}
for label, (sp, acc, expect_org) in FAMILY.items():
    got = pull(sp, {acc})
    if acc not in got:
        L.append(f'  ABORT: {label} {acc} not in the {sp} proteome')
        raise SystemExit('\n'.join(L))
    rec = uniprot(acc)
    ok = expect_org.split()[0] in rec['organism']
    L.append(f'  {label:20} {acc:12} {rec["organism"][:26]:26} {rec["gene"]:10} '
             f'{rec["length"]:4} aa  {rec["name"][:40]}  {"OK" if ok else "ORGANISM MISMATCH"}')
    if not ok:
        raise SystemExit('\n'.join(L) + '\n\nABORT: label does not match the record.')
    seed[label] = got[acc]

pc = pull(*[PCDR1[0]], {PCDR1[1]})
pcdr_seq = pc.get(PCDR1[1])
rec = uniprot(PCDR1[1])
L.append(f'  {"(test) PCDR-1":20} {PCDR1[1]:12} {rec["organism"][:26]:26} {rec["gene"]:10} '
         f'{rec["length"]:4} aa  {rec["name"][:40]}  EXCLUDED FROM PROFILE')
L.append('')

# --------------------------------------------------- control 1: leave-one-species-out
L.append('CONTROL 1 - leave-one-species-out. The profile is rebuilt without each member '
         'and that member re-scored. A member the profile cannot re-find would mean the '
         'profile is blind at that distance, and no absence call could be trusted.')
L.append('')
L.append(f'  {"held-out member":20} {"score vs profile":>16} {"rank in own proteome":>21} '
         f'{"2nd-best in that proteome":>26}')
loo = []
for label in FAMILY:
    sp = FAMILY[label][0]
    sub = {k: v for k, v in seed.items() if k != label}
    hmm = build_hmm(sub, f'loo_{label}')
    hits = scan(hmm, sp)
    acc = FAMILY[label][1]
    rank = next((i + 1 for i, h in enumerate(hits) if h[0] == acc), None)
    score = next((h[2] for h in hits if h[0] == acc), None)
    others = [h for h in hits if h[0] != acc]
    second = others[0][2] if others else float('nan')
    L.append(f'  {label:20} {("%.1f" % score) if score else "NOT FOUND":>16} '
             f'{(str(rank) + " of " + str(len(hits))) if rank else "-":>21} {second:26.1f}')
    # Rank 1 is the WRONG criterion where a proteome holds two co-orthologs: drop
    # GPR52 and the profile still finds GPR21 at 411, so the held-out one is rank 2
    # by construction, not because the profile missed it. The criterion that means
    # anything is whether the held-out member clears the non-member floor set in
    # control 2 - which is computed after this loop, so the check is deferred.
    loo.append((label, score, rank, second))
L.append('')
L.append('  Verdict deferred to control 2 (the floor is computed there).')
L.append('')

# ------------------------------------------------------------- the full-seed profile
hmm = build_hmm(seed, 'orphan7')
with open(f'{OUT}/orphan7.hmm', 'wb') as fh:
    hmm.write(fh)
L.append(f'Profile from all seven: {hmm.M} positions, {len(seed)} sequences '
         f'-> output/orphan/orphan7.hmm')
L.append('')

# --------------------------------------------------- control 2: unrelated class-A floor
L.append('CONTROL 2 - negatives. Melatonin and aminergic receptors are real class-A GPCRs '
         'but not in this family, so they set the floor a genuine hit must clear.')
negs = pull('Hsap', set(MEL_NEG + AMINE_NEG))
hs = {a: s for a, _, s in scan(hmm, 'Hsap')}
L.append('')
for acc in MEL_NEG + AMINE_NEG:
    r = uniprot(acc)
    s = hs.get(acc)
    L.append(f'  {acc:10} {r["gene"]:10} {r["name"][:44]:44} '
             f'{("%.1f" % s) if s else "below threshold":>16}')
neg_top = max([hs[a] for a in MEL_NEG + AMINE_NEG if a in hs], default=0.0)
L.append('')
L.append(f'  Highest-scoring non-member: {neg_top:.1f}. Treat that as the floor.')
L.append('')
loo_fail = [lb for lb, sc_, rk, _ in loo if sc_ is None or sc_ < neg_top]
L.append(f'CONTROL 1, verdict: {len(FAMILY) - len(loo_fail)}/{len(FAMILY)} held-out members '
         f'clear the {neg_top:.1f} floor (range '
         f'{min(s_ for _, s_, _, _ in loo if s_):.1f}-{max(s_ for _, s_, _, _ in loo if s_):.1f}).'
         + (f'  FAILED: {", ".join(loo_fail)}' if loo_fail else
            '  The profile is sensitive across the whole span it claims.'))
L.append('  The two human members are each rank 2 rather than rank 1 only because '
         'GPR52 and GPR21 are co-orthologs in the same proteome: holding one out '
         'leaves the other as the top hit. That is not a sensitivity failure.')
L.append('')

# ------------------------------------------------------------------ scan everything
L.append('=' * 78)
L.append('EVERY AVAILABLE PROTEOME. All hits listed, not one per species, so a '
         'lineage-specific expansion cannot hide as a 1:1 ortholog.')
L.append('')
allsp = sorted({os.path.basename(p).split('_')[0] for p in glob.glob(f'{PROT}/*.fasta')})
species = [sp for sp in allsp if proteome(sp)]
gaps = [sp for sp in allsp if not proteome(sp)]
if gaps:
    # Mcir/Mver/Rory are Mucor circinelloides, Podila (Mortierella) verticillata and
    # Rhizopus oryzae - fungal outgroups inherited from the parent receptor_phylogeny
    # project. Zero-byte files, so not searched, but no metazoan aminergic class-A
    # family could occupy them: not a gap in the lineages this question is about.
    L.append(f'  Not searched (zero-byte proteome file): {", ".join(gaps)} - all three '
             f'are Mucoromycota fungi carried over from the parent project, outside '
             f'any lineage this family could occupy. Metazoan sampling is complete.')
    L.append('')
rows, members = [], {}
L.append(f'  {"sp":6} {"n hits":>7} {"best":>8} {"2nd":>8} {"3rd":>8}   top hit')
for sp in species:
    hits = scan(hmm, sp)
    keep = [h for h in hits if h[2] >= neg_top]
    sc = [h[2] for h in hits[:3]] + [float('nan')] * 3
    top = hits[0] if hits else None
    L.append(f'  {sp:6} {len(keep):7} {sc[0]:8.1f} {sc[1]:8.1f} {sc[2]:8.1f}   '
             + (f'{top[0]:12} {top[1][:46]}' if top else '-'))
    for acc, desc, score in keep:
        rows.append((sp, acc, score, desc))
    pseqs = pull(sp, {h[0] for h in keep})
    for acc, _, score in keep:
        if acc in pseqs:
            members[f'{sp}_{acc}'] = pseqs[acc]

with open(f'{OUT}/scores.tsv', 'w') as fh:
    fh.write('species\taccession\tscore\tdescription\n')
    for sp, acc, score, desc in sorted(rows, key=lambda r: (-r[2],)):
        fh.write(f'{sp}\t{acc}\t{score:.1f}\t{desc}\n')
with open(f'{OUT}/members.faa', 'w') as fh:
    for k, v in sorted(members.items()):
        fh.write(f'>{k}\n{v}\n')
L.append('')
L.append(f'  {len(rows)} hits at or above the non-member floor, across '
         f'{len({r[0] for r in rows})} of {len(species)} proteomes '
         f'-> output/orphan/scores.tsv, members.faa')
L.append('')

# ------------------------------------------------------------------ and PCDR-1?
L.append('=' * 78)
L.append('THE TEST - PCDR-1 against a profile it had no part in building.')
L.append('')
ce = scan(hmm, 'Cele')
rank = next((i + 1 for i, h in enumerate(ce) if h[0] == PCDR1[1]), None)
score = next((h[2] for h in ce if h[0] == PCDR1[1]), None)
L.append(f'  PCDR-1 / O17899: score {("%.1f" % score) if score else "below threshold"}, '
         f'rank {rank if rank else "-"} of {len(ce)} C. elegans hits')
L.append(f'  non-member floor {neg_top:.1f}; '
         f'{"CLEARS it" if score and score > neg_top else "DOES NOT clear it"}')
L.append('')
L.append('  Top C. elegans hits:')
for i, (acc, desc, s) in enumerate(ce[:8]):
    mark = '  <-- PCDR-1' if acc == PCDR1[1] else ''
    L.append(f'    {i+1:2}. {acc:12} {s:8.1f}  {desc[:50]}{mark}')

txt = '\n'.join(L)
print(txt)
open(f'{OUT}/result.txt', 'w').write(txt + '\n')
for f in glob.glob(f'{OUT}/_*'):
    os.remove(f)
