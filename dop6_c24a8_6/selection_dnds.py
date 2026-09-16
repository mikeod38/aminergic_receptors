#!/usr/bin/env python3
"""selection_dnds.py — how constrained is dop-6, relative to the other C. elegans aminergic GPCRs?

For each receptor, take the C. briggsae ortholog (Ensembl Compara), build a codon alignment from
the protein alignment, and compute NG86 dN, dS and dN/dS with Jukes-Cantor correction. dop-6 is
additionally split at residue 405 into the unique half and the half duplicated as C24A8.6.

Writes output/dnds_vs_briggsae.tsv and output/dnds_vs_briggsae.txt
"""
import itertools, json, math, os, sys, time, urllib.request
from Bio.Align import PairwiseAligner, substitution_matrices

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'output')
os.makedirs(OUT, exist_ok=True)
ENS = 'https://rest.ensembl.org'
TARGET = 'caenorhabditis_briggsae'
RECEPTORS = ['dop-6', 'dop-1', 'dop-2', 'dop-3', 'dop-4', 'dop-5', 'ser-1', 'ser-2', 'ser-3',
             'ser-4', 'ser-5', 'ser-6', 'ser-7', 'octr-1', 'tyra-2', 'tyra-3', 'gar-3']
DOP6_SPLIT = 405

def get(u, c='application/json', tries=6):
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(u, headers={'Content-Type': c}), timeout=90) as r:
                t = r.read().decode()
            return json.loads(t) if c == 'application/json' else t.strip()
        except Exception:
            if k == tries - 1: raise
            time.sleep(2 * (k + 1))

CO = {}; i = 0
for b1 in 'TCAG':
    for b2 in 'TCAG':
        for b3 in 'TCAG':
            CO[b1+b2+b3] = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'[i]; i += 1

def syn_sites(c):
    if c not in CO or CO[c] == '*': return 0.0
    s = 0.0
    for k in range(3):
        for b in 'TCAG':
            if b == c[k]: continue
            if CO.get(c[:k]+b+c[k+1:]) == CO[c]: s += 1/3
    return s

def ng86_pair(c1, c2):
    diff = [k for k in range(3) if c1[k] != c2[k]]
    if not diff: return 0.0, 0.0
    sd = nd = 0.0; np_ = 0
    for order in itertools.permutations(diff):
        cur = c1; ok = True; s = n = 0
        for k in order:
            nxt = cur[:k] + c2[k] + cur[k+1:]
            if CO.get(cur) == '*' or CO.get(nxt) == '*': ok = False; break
            if CO.get(cur) == CO.get(nxt): s += 1
            else: n += 1
            cur = nxt
        if ok: sd += s; nd += n; np_ += 1
    if not np_: return 0.0, float(len(diff))
    return sd/np_, nd/np_

def jc(p):
    if p is None or p >= 0.75: return None
    return -0.75 * math.log(1 - 4*p/3)

pa = PairwiseAligner(); pa.substitution_matrix = substitution_matrices.load('BLOSUM62')
pa.open_gap_score, pa.extend_gap_score, pa.mode = -11, -1, 'global'

def canon(gid):
    g = get(f'{ENS}/lookup/id/{gid}?expand=1;content-type=application/json')
    t = ([x for x in g['Transcript'] if x.get('is_canonical')] or g['Transcript'])[0]
    p = get(f'{ENS}/sequence/id/{t["id"]}?type=protein;content-type=text/plain', 'text/plain')
    c = get(f'{ENS}/sequence/id/{t["id"]}?type=cds;content-type=text/plain', 'text/plain')
    return p, c, t['id'], g

def codon_stats(pA, cA, pB, cB, lo=None, hi=None):
    """NG86 over the codon alignment, optionally restricted to Cel residues lo..hi (1-based)."""
    aln = pa.align(pA, pB)[0]; A, B = str(aln[0]), str(aln[1])
    ia = ib = 0; S = N = Sd = Nd = 0.0; ident = naln = 0
    for x, y in zip(A, B):
        if x != '-': ia += 1
        if y != '-': ib += 1
        if x == '-' or y == '-': continue
        if lo is not None and not (lo <= ia <= hi): continue
        c1, c2 = cA[3*(ia-1):3*ia], cB[3*(ib-1):3*ib]
        if len(c1) < 3 or len(c2) < 3 or 'N' in c1+c2: continue
        if c1 not in CO or c2 not in CO: continue
        naln += 1; ident += (x == y)
        s1, s2 = syn_sites(c1), syn_sites(c2)
        S += (s1+s2)/2; N += (3-s1 + 3-s2)/2
        sd, nd = ng86_pair(c1, c2); Sd += sd; Nd += nd
    pS = Sd/S if S else None; pN = Nd/N if N else None
    dS, dN = jc(pS), jc(pN)
    return dict(codons=naln, ident=ident, pid=100*ident/naln if naln else 0,
                S=S, N=N, Sd=Sd, Nd=Nd, dS=dS, dN=dN,
                omega=(dN/dS if (dN is not None and dS not in (None, 0)) else None))

rows = []
log = open(f'{OUT}/dnds_vs_briggsae.txt', 'w')
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); log.write(s + '\n'); log.flush()

say(f'NG86 dN/dS, C. elegans vs C. briggsae, canonical transcripts (Ensembl/WormBase)')
say(f'{"gene":<10} {"Cbr ortholog":<18} {"codons":>6} {"%id_aa":>7} {"dN":>7} {"dS":>7} {"dN/dS":>7}')
say('-' * 70)
for name in RECEPTORS:
    try:
        g = get(f'{ENS}/lookup/symbol/caenorhabditis_elegans/{name}?content-type=application/json')
        gid = g['id']
        h = get(f'{ENS}/homology/id/caenorhabditis_elegans/{gid}?content-type=application/json'
                f';target_species={TARGET};type=orthologues')
        hs = h['data'][0]['homologies']
        if not hs: say(f'{name:<10} (no C. briggsae ortholog)'); continue
        best = max(hs, key=lambda x: len(x['target']['align_seq'].replace('-', '')))
        pA, cA, tA, _ = canon(gid)
        pB, cB, tB, _ = canon(best['target']['id'])
        st = codon_stats(pA, cA, pB, cB)
        say(f'{name:<10} {best["target"]["id"]:<18} {st["codons"]:>6} {st["pid"]:>7.1f} '
            f'{st["dN"] if st["dN"] is not None else float("nan"):>7.4f} '
            f'{st["dS"] if st["dS"] is not None else float("nan"):>7.4f} '
            f'{st["omega"] if st["omega"] is not None else float("nan"):>7.4f}')
        rows.append((name, best['target']['id'], 'whole', st))
        if name == 'dop-6':
            for lab, lo, hi in [('1-405 unique', 1, DOP6_SPLIT), ('406-692 duplicated', DOP6_SPLIT+1, len(pA))]:
                s2 = codon_stats(pA, cA, pB, cB, lo, hi)
                say(f'{"  "+lab:<29} {s2["codons"]:>6} {s2["pid"]:>7.1f} '
                    f'{s2["dN"] if s2["dN"] is not None else float("nan"):>7.4f} '
                    f'{s2["dS"] if s2["dS"] is not None else float("nan"):>7.4f} '
                    f'{s2["omega"] if s2["omega"] is not None else float("nan"):>7.4f}')
                rows.append((name, best['target']['id'], lab, s2))
    except Exception as e:
        say(f'{name:<10} FAILED: {type(e).__name__} {e}')

with open(f'{OUT}/dnds_vs_briggsae.tsv', 'w') as f:
    f.write('gene\tcbr_ortholog\tsegment\tcodons\tident\tpid_aa\tS\tN\tSd\tNd\tdN\tdS\tomega\n')
    for n, o, seg, s in rows:
        f.write(f'{n}\t{o}\t{seg}\t{s["codons"]}\t{s["ident"]}\t{s["pid"]:.2f}\t{s["S"]:.1f}\t{s["N"]:.1f}\t'
                f'{s["Sd"]:.1f}\t{s["Nd"]:.1f}\t'
                f'{"" if s["dN"] is None else f"{s['dN']:.4f}"}\t{"" if s["dS"] is None else f"{s['dS']:.4f}"}\t'
                f'{"" if s["omega"] is None else f"{s['omega']:.4f}"}\n')
say(f'\n-> output/dnds_vs_briggsae.tsv ({len(rows)} rows)')
log.close()
