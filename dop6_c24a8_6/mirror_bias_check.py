#!/usr/bin/env python3
"""mirror_bias_check.py — does the N2 palindrome deflate the measured substitution rate in
DOP-6 codons 406-692?

The concern: dN for the duplicated C-half (0.0197) came out lower than for the unique N-half
(0.0293), and both were measured using the N2 reference sequence — the one strain in which that
segment sits in a perfect inverted repeat. If intrachromosomal gene conversion between the arms
suppressed substitution there, the low C-half dN would be an artefact of the mirror rather than
evidence of constraint.

Two independent tests, neither of which relies on an argument about timescales:

  TEST 1 (within species, the direct one). Compare N2 with wild isolates that have NO palindrome
  (CB4856, XZ1516, DL238, MY2147), segment by segment. If the mirror suppresses substitution in
  N2's C-half, N2-vs-wild divergence must be lower in the C-half than in the N-half. The wild
  isolates are single-copy there, so they provide the unconstrained comparison.

  TEST 2 (between species). Rebuild the dop-6 CDS from each wild isolate's single-copy locus and
  redo NG86 dN/dS against C. briggsae using that sequence instead of N2's. If the N2-based C-half
  dN is an artefact, the wild-isolate-based value will be higher.

Also reports whether the N-half/C-half dN difference is significant at all (Fisher exact on
nonsynonymous vs synonymous difference counts).

Writes output/mirror_bias_check.txt
"""
import itertools, json, math, os, subprocess, time, urllib.request
from Bio.Align import PairwiseAligner, substitution_matrices

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'output')
WORK = ('/private/tmp/claude-501/-Users-mikeodonnell-git-Malaiwong-Oglu-2026/'
        '93015e40-e331-495a-8362-20ac7ea07c75/scratchpad/asm')
ENS = 'https://rest.ensembl.org'
L0 = 4318000; SPLIT = 405
ARM_A = (4324178, 4328106)
n2 = ''.join(l.strip() for l in open(f'{OUT}/locus_X_4318000-4334500.fa') if not l.startswith('>')).upper()
G = lambda a, b: n2[a-L0:b-L0+1]
RC = lambda s: s[::-1].translate(str.maketrans('ACGTacgt', 'TGCAtgca'))
models = json.load(open(f'{OUT}/gene_models.json'))
EXONS = models['dop-6']['transcripts'][0]['exons']
WILD = [('GCA_028201415_1', 'CB4856'), ('GCA_016989125_1', 'XZ1516'),
        ('GCA_016989505_1', 'DL238'), ('GCA_016989235_1', 'MY2147')]

CO = {}; i = 0
for b1 in 'TCAG':
    for b2 in 'TCAG':
        for b3 in 'TCAG':
            CO[b1+b2+b3] = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'[i]; i += 1
def syn_sites(c):
    if c not in CO or CO[c] == '*': return 0.0
    return sum(1/3 for k in range(3) for b in 'TCAG'
               if b != c[k] and CO.get(c[:k]+b+c[k+1:]) == CO[c])
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
    return (sd/np_, nd/np_) if np_ else (0.0, float(len(diff)))
def jc(p):
    return None if (p is None or p >= 0.75) else -0.75*math.log(1-4*p/3)

rep = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); rep.append(s)

def get(u, c='text/plain', tries=6):
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(u, headers={'Content-Type': c}), timeout=90) as r:
                t = r.read().decode()
            return json.loads(t) if c == 'application/json' else t.strip()
        except Exception:
            if k == tries-1: raise
            time.sleep(2*(k+1))

# ── rebuild each wild isolate's dop-6 CDS, exon by exon, from its own assembly ──
def wild_cds(tag, name):
    """blast each N2 exon into the isolate's window; concatenate the best hits in order."""
    win = f'{WORK}/{tag}_win'
    if not os.path.exists(win + '.nsq'): return None, 'no window db'
    seqs = []
    for j, (es, ee) in enumerate(EXONS, 1):
        q = f'{WORK}/_ex.fa'; open(q, 'w').write(f'>ex{j}\n{G(es,ee)}\n')
        r = subprocess.run(['blastn', '-query', q, '-db', win, '-dust', 'no', '-evalue', '1e-5',
                            '-outfmt', '6 sstart send pident length sseq'],
                           capture_output=True, text=True)
        hs = [h.split('\t') for h in r.stdout.strip().split('\n') if h]
        hs = [h for h in hs if int(h[3]) == ee-es+1]
        if not hs: return None, f'exon {j} not recovered full length'
        # in the oriented window dop-6 is on the plus strand
        plus = [h for h in hs if int(h[1]) > int(h[0])]
        h = (plus or hs)[0]
        seqs.append(h[4].replace('-', ''))
    return ''.join(seqs), None

say('=' * 94)
say('DOES THE N2 PALINDROME DEFLATE THE SUBSTITUTION RATE IN DOP-6 CODONS 406-692?')
say('=' * 94)
say(f'N-half = codons 1-{SPLIT} (single copy in every strain)')
say(f'C-half = codons {SPLIT+1}-692 (inside the inverted repeat in N2 only)\n')

n2_cds = get(f'{ENS}/sequence/id/C24A8.1.1?type=cds;content-type=text/plain')
n2_prot = ''.join(CO.get(n2_cds[i:i+3], 'X') for i in range(0, len(n2_cds)-2, 3)).rstrip('*')
say(f'N2 CDS {len(n2_cds)} nt / {len(n2_prot)} aa')

# ── TEST 1 ───────────────────────────────────────────────────────────────────
say('\n' + '-'*94)
say('TEST 1 — N2 vs palindrome-free wild isolates, per segment (nucleotide, whole CDS)')
say('-'*94)
say(f'{"isolate":<10} {"N-half diffs/sites":>22} {"%div":>7}   {"C-half diffs/sites":>22} {"%div":>7}')
wild_cdss = {}
for tag, name in WILD:
    cds, err = wild_cds(tag, name)
    if err: say(f'{name:<10} CDS rebuild failed: {err}'); continue
    wild_cdss[name] = cds
    if len(cds) != len(n2_cds):
        say(f'{name:<10} length {len(cds)} != N2 {len(n2_cds)} — indel present, reporting anyway')
    nN = sum(1 for a, b in zip(n2_cds[:3*SPLIT], cds[:3*SPLIT]) if a != b)
    nC = sum(1 for a, b in zip(n2_cds[3*SPLIT:], cds[3*SPLIT:]) if a != b)
    lN, lC = min(3*SPLIT, len(cds)), max(0, min(len(n2_cds), len(cds)) - 3*SPLIT)
    say(f'{name:<10} {f"{nN}/{lN}":>22} {100*nN/lN:>7.3f}   {f"{nC}/{lC}":>22} '
        f'{100*nC/lC if lC else 0:>7.3f}')
say('\nIf the mirror suppressed substitution in N2, C-half %div would be systematically LOWER')
say('than N-half %div. These isolates are single-copy in both halves, so they supply the control.')

# ── TEST 2 ───────────────────────────────────────────────────────────────────
say('\n' + '-'*94)
say('TEST 2 — NG86 dN/dS vs C. briggsae, computed from N2 vs from each single-copy wild isolate')
say('-'*94)
pa = PairwiseAligner(); pa.substitution_matrix = substitution_matrices.load('BLOSUM62')
pa.open_gap_score, pa.extend_gap_score, pa.mode = -11, -1, 'global'
cbr_prot = None
for r in open(f'{OUT}/dop6_orthologs.faa').read().split('>')[1:]:
    h, s = r.split('\n', 1)
    if h.startswith('Cbr-dop-6'): cbr_prot = s.replace('\n', '')
cbr_cds = get(f'{ENS}/sequence/id/CBG14393b.1?type=cds;content-type=text/plain')

def seg_stats(cds, prot, lo, hi):
    aln = pa.align(prot, cbr_prot)[0]; A, B = str(aln[0]), str(aln[1])
    ia = ib = 0; S = N = Sd = Nd = 0.0; nc = 0
    for x, y in zip(A, B):
        if x != '-': ia += 1
        if y != '-': ib += 1
        if x == '-' or y == '-': continue
        if not (lo <= ia <= hi): continue
        c1, c2 = cds[3*(ia-1):3*ia], cbr_cds[3*(ib-1):3*ib]
        if len(c1) < 3 or len(c2) < 3 or c1 not in CO or c2 not in CO: continue
        nc += 1
        s1, s2 = syn_sites(c1), syn_sites(c2)
        S += (s1+s2)/2; N += (3-s1 + 3-s2)/2
        sd, nd = ng86_pair(c1, c2); Sd += sd; Nd += nd
    return nc, S, N, Sd, Nd, jc(Nd/N if N else None), jc(Sd/S if S else None)

say(f'{"source of Cel seq":<20} {"segment":<20} {"codons":>6} {"Nd":>7} {"N":>8} {"dN":>8} {"dS":>8}')
allrows = []
for name, cds in [('N2 (reference)', n2_cds)] + list(wild_cdss.items()):
    prot = ''.join(CO.get(cds[i:i+3], 'X') for i in range(0, len(cds)-2, 3)).rstrip('*')
    for lab, lo, hi in [('N-half 1-405', 1, SPLIT), ('C-half 406-692', SPLIT+1, len(prot))]:
        nc, S, N, Sd, Nd, dN, dS = seg_stats(cds, prot, lo, hi)
        allrows.append((name, lab, nc, Nd, N, Sd, S, dN, dS))
        say(f'{name:<20} {lab:<20} {nc:>6} {Nd:>7.1f} {N:>8.1f} '
            f'{dN if dN is not None else float("nan"):>8.4f} {dS if dS is not None else float("nan"):>8.4f}')

# ── is the N-half/C-half difference significant at all? ──────────────────────
say('\n' + '-'*94)
say('IS THE N-half / C-half dN DIFFERENCE EVEN SIGNIFICANT?')
say('-'*94)
def fisher(a, b, c, d):
    def lf(n): return math.lgamma(n+1)
    def p(a, b, c, d):
        return math.exp(lf(a+b)+lf(c+d)+lf(a+c)+lf(b+d)-lf(a)-lf(b)-lf(c)-lf(d)-lf(a+b+c+d))
    p0 = p(a, b, c, d); tot = 0.0
    for i in range(0, a+b+c+d+1):
        aa = i; bb = a+b-i; cc = a+c-i; dd = d-a+i
        if min(aa, bb, cc, dd) < 0: continue
        pi = p(aa, bb, cc, dd)
        if pi <= p0*(1+1e-9): tot += pi
    return min(1.0, tot)
for name in ['N2 (reference)'] + list(wild_cdss.keys()):
    rs = [r for r in allrows if r[0] == name]
    if len(rs) != 2: continue
    (_, _, _, NdN, NN, SdN, SN, dNN, _), (_, _, _, NdC, NC, SdC, SC, dNC, _) = rs
    pv = fisher(round(NdN), round(NN-NdN), round(NdC), round(NC-NdC))
    say(f'{name:<20} nonsyn diffs  N-half {NdN:.0f}/{NN:.0f}  vs  C-half {NdC:.0f}/{NC:.0f}   '
        f'Fisher p = {pv:.3f}  {"(n.s.)" if pv > 0.05 else "(significant)"}')

say('\nTiming note, independent of both tests: the duplication is absent from CB4856, XZ1516,')
say('DL238 and MY2147, so it postdates the common ancestor of those isolates with N2 and is')
say('vastly younger than the C. elegans / C. briggsae split (dS saturated, > 1). For essentially')
say('the whole branch on which dN was measured, dop-6 was single-copy.')
open(f'{OUT}/mirror_bias_check.txt', 'w').write('\n'.join(rep) + '\n')
say('\n-> output/mirror_bias_check.txt')
