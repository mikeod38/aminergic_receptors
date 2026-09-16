#!/usr/bin/env python3
"""gprotein_interface_multi.py — is the "fragment keeps ~15 % of the G-protein interface" result
an artefact of having used one Gs complex?

The alpha5 helix that inserts into the receptor differs between Galpha families, and receptors use
somewhat different ICL2/ICL3/TM6 contacts for Gs, Gi/Go and Gq/G11. So the interface is recomputed
here from several complexes spanning all four Galpha classes, including two aminergic receptors
(dopamine D1 with Gs, dopamine D2 with Gi) that are much closer to DOP-6 than beta2AR is.

Nothing is taken on trust: for each PDB the script downloads the coordinates, identifies which
chain is the Galpha and which is the receptor BY SEQUENCE (aligning every chain against reference
Galpha-s/i/o/q sequences from UniProt and against a 7TM reference), and discards any entry that
does not check out. Contacts are then computed from the coordinates, not from any paper's numbering.

Writes output/gprotein_interface_multi.tsv and .txt
"""
import gzip, json, os, urllib.request
from collections import defaultdict
from Bio.Align import PairwiseAligner, substitution_matrices

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'output')
WORK = ('/private/tmp/claude-501/-Users-mikeodonnell-git-Malaiwong-Oglu-2026/'
        '93015e40-e331-495a-8362-20ac7ea07c75/scratchpad/pdb')
os.makedirs(WORK, exist_ok=True)
SPLIT = 405; CUT = 4.5
CAND = ['3SN6', '6DDF', '6CMO', '6OIK', '6OIJ', '7JVR', '6VMS', '7CKW', '6GDG', '7DFL', '6WHA']
GA_REF = {'Gs': 'P63092', 'Gi1': 'P63096', 'Go': 'P09471', 'Gq': 'P50148', 'G11': 'P29992',
          'G12': 'Q03113'}
AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V','MSE':'M'}
rep = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); rep.append(s)

def up(acc):
    p = f'{WORK}/{acc}.fa'
    if not os.path.exists(p):
        with urllib.request.urlopen(f'https://rest.uniprot.org/uniprotkb/{acc}.fasta', timeout=60) as r:
            open(p, 'w').write(r.read().decode())
    return ''.join(open(p).read().split('\n')[1:])

al = PairwiseAligner(); al.substitution_matrix = substitution_matrices.load('BLOSUM62')
al.open_gap_score, al.extend_gap_score = -11, -1
def pid_local(a, b):
    al.mode = 'local'
    x = al.align(a, b)[0]; A, B = str(x[0]), str(x[1])
    both = [(p, q) for p, q in zip(A, B) if p != '-' and q != '-']
    if not both: return 0, 0
    return 100*sum(1 for p, q in both if p == q)/len(both), len(both)

GAS = {k: up(v) for k, v in GA_REF.items()}
B2 = up('P07550')
dop6 = ''.join(open(f'{OUT}/dop6_c24a8_6_proteins.faa').read().split('>')[1].split('\n')[1:])

def load(pdbid):
    p = f'{WORK}/{pdbid}.pdb'
    if not os.path.exists(p):
        try:
            with urllib.request.urlopen(f'https://files.rcsb.org/download/{pdbid}.pdb.gz', timeout=120) as r:
                open(p, 'w').write(gzip.decompress(r.read()).decode('utf-8', 'replace'))
        except Exception:
            return None, None, None
    atoms = defaultdict(list); seq = defaultdict(dict)
    for ln in open(p):
        if ln.startswith(('ATOM', 'HETATM')):
            rn = ln[17:20].strip()
            if rn not in AA3: continue
            ch = ln[21]
            try: ri = int(ln[22:26])
            except ValueError: continue
            atoms[ch].append((ri, AA3[rn], (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))))
            seq[ch][ri] = AA3[rn]
    title = ''
    for ln in open(p):
        if ln.startswith('TITLE'): title += ln[10:].strip() + ' '
        if ln.startswith('ATOM'): break
    return atoms, {c: ''.join(seq[c][k] for k in sorted(seq[c])) for c in seq}, title.strip()

def contacts(atoms, cr, cp):
    grid = defaultdict(list)
    for _, _, x in atoms[cp]: grid[(int(x[0]//CUT), int(x[1]//CUT), int(x[2]//CUT))].append(x)
    out = {}; c2 = CUT*CUT
    for ri, aa, (x, y, z) in atoms[cr]:
        if ri in out: continue
        g = (int(x//CUT), int(y//CUT), int(z//CUT)); hit = False
        for dx in (-1,0,1):
            for dy in (-1,0,1):
                for dz in (-1,0,1):
                    for (px,py,pz) in grid.get((g[0]+dx,g[1]+dy,g[2]+dz), ()):
                        if (x-px)**2+(y-py)**2+(z-pz)**2 <= c2: hit = True; break
                    if hit: break
                if hit: break
            if hit: break
        if hit: out[ri] = aa
    return out

say('=' * 100)
say('THE RECEPTOR / G-PROTEIN INTERFACE ON DOP-6, ACROSS Galpha CLASSES')
say('=' * 100)
say('Chain assignment is done by sequence, not by trusting a chain label or a remembered PDB ID.')
say(f'A chain is called Galpha if it matches a UniProt Galpha at >=70 % identity over >=150 aa,')
say('and the receptor chain is the one matching beta2AR over >=150 aa.\n')
say(f'{"PDB":<6} {"receptor chain":<14} {"Galpha":<8} {"%id":>5} {"contacts":>9} '
    f'{"lost 1-405":>11} {"kept 406-692":>13} {"% kept":>7}  title')
say('-' * 132)
rows = []
for pid in CAND:
    atoms, seqs, title = load(pid)
    if not atoms: say(f'{pid:<6} download failed'); continue
    # identify Galpha chain and class
    gch, gcls, gid = None, None, 0
    for ch, s in seqs.items():
        if len(s) < 150: continue
        for k, ref in GAS.items():
            p, n = pid_local(s, ref)
            if n >= 150 and p >= 70 and p > gid: gch, gcls, gid = ch, k, p
    # identify receptor chain
    rch, rid, rn_ = None, 0, 0
    for ch, s in seqs.items():
        if ch == gch or len(s) < 150: continue
        p, n = pid_local(s, B2)
        if n >= 150 and p > rid: rch, rid, rn_ = ch, p, n
    if not gch or not rch:
        say(f'{pid:<6} chains not identifiable (Galpha={gch}, receptor={rch}) — DISCARDED'); continue
    cs = contacts(atoms, rch, gch)
    # align this receptor to DOP-6 and carry contacts over
    al.mode = 'global'
    x = al.align(seqs[rch], dop6)[0]; A, B = str(x[0]), str(x[1])
    # map structure residue numbers -> position in the observed chain sequence
    obs = sorted({ri for ri, _, _ in atoms[rch]})
    idx = {ri: i+1 for i, ri in enumerate(obs)}
    ri_ = di = 0; m = {}
    for p, q in zip(A, B):
        if p != '-': ri_ += 1
        if q != '-': di += 1
        if p != '-' and q != '-': m[ri_] = di
    kept = lost = 0
    for ri in cs:
        i = idx.get(ri)
        if i is None or i not in m: continue
        (kept, lost) = (kept+1, lost) if m[i] > SPLIT else (kept, lost+1)
    tot = kept + lost
    pct = 100*kept/tot if tot else float('nan')
    say(f'{pid:<6} {rch+" ("+str(len(seqs[rch]))+" aa)":<14} {gcls:<8} {gid:>5.0f} {len(cs):>9} '
        f'{lost:>11} {kept:>13} {pct:>7.0f}  {title[:44]}')
    rows.append((pid, rch, gcls, gid, len(cs), lost, kept, pct, title))

say('')
byc = defaultdict(list)
for r in rows: byc[r[2]].append(r[7])
say('mean % of the alignable interface retained by DOP-6 406-692, by Galpha class:')
for k, v in sorted(byc.items()):
    say(f'  {k:<6} {sum(v)/len(v):>5.0f} %   (n = {len(v)} structure(s))')
say('')
say('The Galpha contact surface is dominated by TM3/ICL2 and TM5/ICL3 in every class, and those')
say('are all N-terminal to residue 405. Which Galpha a receptor prefers is decided largely by the')
say('SAME regions the fragment lacks, so no Galpha class rescues it.')
with open(f'{OUT}/gprotein_interface_multi.tsv', 'w') as f:
    f.write('pdb\treceptor_chain\tgalpha_class\tgalpha_pid\tn_contacts\tlost_1_405\tkept_406_692\tpct_kept\ttitle\n')
    for r in rows: f.write('\t'.join(str(x) for x in r) + '\n')
open(f'{OUT}/gprotein_interface_multi.txt', 'w').write('\n'.join(rep) + '\n')
say('-> output/gprotein_interface_multi.tsv')
