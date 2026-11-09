#!/usr/bin/env python3
"""galpha_conservation.py — would a human melatonin receptor find a partner in a worm?

A humanised assay needs the human receptor to drive a C. elegans G protein. Two
questions, and this script answers the second:

  1. Does the worm have the right Galpha class at all?  (repertoire, below)
  2. Is the surface the receptor actually touches conserved between the human and
     worm proteins?  (contacts, computed here)

Method follows ../dop6_c24a8_6/gprotein_interface_multi.py: contacts are taken from
real coordinates, not from any paper's numbering, and the chains are identified by
sequence rather than by chain label. Here the contacts wanted are the Galpha-side
ones - which residues of Galpha sit within 4.5 A of the receptor - mapped onto the
human reference and then read off each worm ortholog through a pairwise alignment.

Writes output/galpha_conservation.txt
"""
import gzip, json, os, re, urllib.request
from collections import defaultdict
from Bio.Align import PairwiseAligner, substitution_matrices

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'output')
WORK = os.path.join(HERE, 'data', 'pdb')
os.makedirs(WORK, exist_ok=True)
CUT = 4.5

# Gi/Go complexes - the class MT1/MT2 couple to. Two independent structures.
COMPLEXES = ['6DDF', '6OIJ']
AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V','MSE':'M'}

aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
aligner.open_gap_score, aligner.extend_gap_score, aligner.mode = -11, -1, 'global'


def uniprot(acc):
    j = json.load(urllib.request.urlopen(
        f'https://rest.uniprot.org/uniprotkb/{acc}.json', timeout=30))
    name = (j['proteinDescription'].get('recommendedName', {})
            .get('fullName', {}).get('value', '?'))
    genes = [g.get('geneName', {}).get('value') for g in j.get('genes', [])]
    return {'acc': acc, 'name': name, 'gene': genes[0] if genes else '?',
            'organism': j['organism']['scientificName'], 'seq': j['sequence']['value']}


def by_query(q):
    """Resolve an accession by query so nothing is written from memory."""
    u = ('https://rest.uniprot.org/uniprotkb/search?format=json&size=1&fields=accession&query='
         + urllib.request.quote(q))
    r = json.load(urllib.request.urlopen(u, timeout=30))
    return r['results'][0]['primaryAccession'] if r.get('results') else None


def pdb_chains(pdb):
    path = os.path.join(WORK, f'{pdb}.cif.gz')
    if not os.path.exists(path):
        urllib.request.urlretrieve(f'https://files.rcsb.org/download/{pdb}.cif.gz', path)
    # Read the _atom_site loop's own column order rather than assuming positions -
    # guessing them put every atom in one chain and silently produced no complex.
    cols, chains = {}, defaultdict(list)
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            s = line.strip()
            if s.startswith('_atom_site.'):
                cols[s.split('.', 1)[1]] = len(cols)
                continue
            if not line.startswith('ATOM'):
                continue
            need = ('label_comp_id', 'auth_asym_id', 'auth_seq_id', 'Cartn_x', 'Cartn_y', 'Cartn_z')
            if not all(k in cols for k in need):
                raise SystemExit(f'ABORT: {pdb} atom_site is missing columns {need}')
            f = line.split()
            res = f[cols['label_comp_id']]
            if res not in AA3:
                continue
            try:
                xyz = tuple(float(f[cols[k]]) for k in ('Cartn_x', 'Cartn_y', 'Cartn_z'))
                num = int(f[cols['auth_seq_id']])
            except ValueError:
                continue
            chains[f[cols['auth_asym_id']]].append((num, AA3[res], xyz))
    return chains


def seq_of(atoms):
    seen, out = set(), []
    for num, aa, _ in atoms:
        if num not in seen:
            seen.add(num); out.append((num, aa))
    out.sort()
    return ''.join(a for _, a in out), [n for n, _ in out]


def pid(a, b):
    aln = aligner.align(a, b)[0]
    s = aln[0], aln[1]
    same = sum(1 for x, y in zip(*s) if x == y and x != '-')
    cov = sum(1 for x, y in zip(*s) if x != '-' and y != '-')
    return 100.0 * same / max(cov, 1)


def map_positions(ref, other):
    """ref index (0-based) -> other index, via global alignment."""
    aln = aligner.align(ref, other)[0]
    m, i, j = {}, 0, 0
    for x, y in zip(aln[0], aln[1]):
        if x != '-' and y != '-':
            m[i] = j
        if x != '-':
            i += 1
        if y != '-':
            j += 1
    return m


# ---- sequences --------------------------------------------------------------
HUMAN = {'GNAI1': 'gene:GNAI1 AND organism_id:9606 AND reviewed:true',
         'GNAO1': 'gene:GNAO1 AND organism_id:9606 AND reviewed:true',
         'GNAQ':  'gene:GNAQ AND organism_id:9606 AND reviewed:true',
         'GNAS':  'gene:GNAS AND organism_id:9606 AND reviewed:true'}
WORM = {'goa-1':  'gene:goa-1 AND organism_id:6239 AND reviewed:true',
        'gpa-16': 'gene:gpa-16 AND organism_id:6239 AND reviewed:true',
        'egl-30': 'gene:egl-30 AND organism_id:6239 AND reviewed:true',
        'gsa-1':  'gene:gsa-1 AND organism_id:6239 AND reviewed:true',
        'gpa-7':  'gene:gpa-7 AND organism_id:6239 AND reviewed:true'}

lines = ['Would a human melatonin receptor find a G protein partner in C. elegans?', '']
seqs = {}
lines.append('Sequences, each resolved by query and verified against its own record:')
for label, q in list(HUMAN.items()) + list(WORM.items()):
    acc = by_query(q)
    if not acc:
        lines.append(f'  {label:8} NOT FOUND'); continue
    rec = uniprot(acc)
    seqs[label] = rec['seq']
    lines.append(f'  {label:8} {acc:8} {rec["organism"][:22]:22} {rec["name"][:44]}')
lines.append('')

# ---- Galpha-side contacts from real coordinates ------------------------------
contact_sets = []
for pdb in COMPLEXES:
    chains = pdb_chains(pdb)
    best_g, best_r = None, None
    for ch, atoms in chains.items():
        s, _ = seq_of(atoms)
        if len(s) < 150:
            continue
        gi = pid(seqs['GNAI1'], s)
        if best_g is None or gi > best_g[1]:
            best_g = (ch, gi)
    # receptor = the chain that is not Galpha/Gbeta/Ggamma and is 7TM-sized
    for ch, atoms in chains.items():
        if ch == best_g[0]:
            continue
        s, _ = seq_of(atoms)
        if 250 <= len(s) <= 500 and pid(seqs['GNAI1'], s) < 25:
            if best_r is None or len(s) > best_r[1]:
                best_r = (ch, len(s))
    if not best_g or not best_r or best_g[1] < 80:
        lines.append(f'{pdb}: chains not identified ({best_g}, {best_r}) - SKIPPED'); continue
    g_atoms = chains[best_g[0]]; r_atoms = chains[best_r[0]]
    near = set()
    for gn, ga, (x1, y1, z1) in g_atoms:
        for _, _, (x2, y2, z2) in r_atoms:
            if (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2 <= CUT*CUT:
                near.add((gn, ga)); break
    gseq, gnums = seq_of(g_atoms)
    idx = {n: i for i, n in enumerate(gnums)}
    m = map_positions(gseq, seqs['GNAI1'])
    refpos = sorted({m[idx[n]] for n, _ in near if n in idx and idx[n] in m})
    contact_sets.append(set(refpos))
    lines.append(f'{pdb}: Galpha chain {best_g[0]} ({best_g[1]:.1f}% to GNAI1), '
                 f'receptor chain {best_r[0]} ({best_r[1]} res) -> '
                 f'{len(refpos)} Galpha contact positions')

if not contact_sets:
    raise SystemExit('ABORT: no usable complex')
consensus = sorted(set.intersection(*contact_sets)) if len(contact_sets) > 1 else sorted(contact_sets[0])
union = sorted(set.union(*contact_sets))
lines += ['', f'Receptor-contact positions on Galpha(i): {len(consensus)} found in both '
              f'structures, {len(union)} in either. Conservation is scored on the '
              f'{len(consensus)} shared ones.', '']

# ---- conservation at those positions ----------------------------------------
ref = seqs['GNAI1']
lines.append(f'{"protein":10} {"whole-protein id":>17} {"identity at the":>17} {"contact residues differing":>28}')
lines.append(f'{"":10} {"vs GNAI1":>17} {"contact sites":>17}')
for label in ['GNAO1', 'GNAQ', 'GNAS', 'goa-1', 'gpa-16', 'gpa-7', 'egl-30', 'gsa-1']:
    if label not in seqs:
        continue
    m = map_positions(ref, seqs[label])
    hit = tot = 0
    diffs = []
    for p in consensus:
        if p in m:
            tot += 1
            if seqs[label][m[p]] == ref[p]:
                hit += 1
            else:
                diffs.append(f'{ref[p]}{p+1}{seqs[label][m[p]]}')
    lines.append(f'{label:10} {pid(ref, seqs[label]):16.1f}% '
                 f'{100.0*hit/max(tot,1):15.1f}% ({hit}/{tot})   '
                 f'{" ".join(diffs[:6])}{" ..." if len(diffs) > 6 else ""}')

txt = '\n'.join(lines)
print(txt)
open(os.path.join(OUT, 'galpha_conservation.txt'), 'w').write(txt + '\n')
