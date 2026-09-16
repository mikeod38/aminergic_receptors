#!/usr/bin/env python3
"""polymorphism_from_assemblies.py — real dop-6 diversity, without mapping to a duplicated reference.

Short-read variant calls inside the palindrome are missing because reads cannot be placed against
a reference that carries two identical copies — the sequence exists in the data, it just has
nowhere unique to go. Long-read assemblies do not have that problem: each isolate's own contig
carries its own copy. So this computes diversity for dop-6 straight from the per-isolate CDSs that
survey_assemblies.py rebuilt, and contrasts it with what the CaeNDR VCF reports over the same
codons, which quantifies the size of the mapping artefact.

Per segment (N-half codons 1-405, single copy everywhere; C-half codons 406-692, duplicated in
carriers) it reports segregating sites, pi, piN, piS, and the same counts split by whether the
isolate carries the duplication.

Writes output/dop6_polymorphism.tsv and output/dop6_polymorphism.txt
"""
import itertools, json, math, os, subprocess
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'output')
SPLIT = 405
VCF = ('https://caendr-open-access-data-bucket.s3.us-east-2.amazonaws.com/dataset_release/'
       'c_elegans/20250625/variation/WI.20250625.hard-filter.isotype.vcf.gz')
models = json.load(open(f'{OUT}/gene_models.json'))
EXONS = models['dop-6']['transcripts'][0]['exons']          # ascending, + strand gene

CO = {}; i = 0
for b1 in 'TCAG':
    for b2 in 'TCAG':
        for b3 in 'TCAG':
            CO[b1+b2+b3] = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'[i]; i += 1
def syn_sites(c):
    if c not in CO or CO[c] == '*': return 0.0
    return sum(1/3 for k in range(3) for b in 'TCAG'
               if b != c[k] and CO.get(c[:k]+b+c[k+1:]) == CO[c])

rep = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); rep.append(s)

# ── load the per-isolate CDSs and the carrier status ─────────────────────────
seqs = {}
for r in open(f'{OUT}/dop6_cds_isolates.fna').read().split('>')[1:]:
    h, s = r.split('\n', 1); seqs[h.strip()] = s.replace('\n', '').upper()
carrier = {}
if os.path.exists(f'{OUT}/assembly_survey.tsv'):
    for ln in open(f'{OUT}/assembly_survey.tsv').read().strip().split('\n')[1:]:
        p = ln.split('\t')
        if len(p) >= 4 and p[3].isdigit(): carrier[p[0].replace(' ', '_')] = int(p[3])

L = Counter(len(s) for s in seqs.values())
say('=' * 92)
say('dop-6 DIVERSITY COMPUTED FROM ASSEMBLIES, NOT FROM READS MAPPED TO WBcel235')
say('=' * 92)
say(f'{len(seqs)} CDSs; lengths {dict(L)}')
ref = seqs.get('N2_WBcel235_reference')
usable = {k: v for k, v in seqs.items() if len(v) == len(ref)}
say(f'{len(usable)} of them are the reference length ({len(ref)} nt) and enter the analysis')
drop = [k for k in seqs if k not in usable]
if drop: say(f'excluded (indel relative to N2, needs codon-aware alignment): {", ".join(drop)}')

# ── segregating sites and diversity, per segment ─────────────────────────────
def per_seg(keys):
    """segregating sites, pi, and nonsyn/syn split per segment"""
    n = len(keys)
    # cols[j] is alignment column j across the isolates (zip transposes the sequences)
    cols = list(zip(*[usable[k] for k in keys]))
    out = {'N': dict(S=0, pi_num=0.0, sites=0, nsyn=0, syn=0, Nsit=0.0, Ssit=0.0),
           'C': dict(S=0, pi_num=0.0, sites=0, nsyn=0, syn=0, Nsit=0.0, Ssit=0.0)}
    for ci in range(0, len(ref)-2, 3):
        cr = ref[ci:ci+3]
        if cr not in CO: continue
        seg = 'N' if ci//3 < SPLIT else 'C'
        out[seg]['Ssit'] += syn_sites(cr); out[seg]['Nsit'] += 3 - syn_sites(cr)
        for k in range(3):
            if ci + k >= len(cols): continue
            al = Counter(x for x in cols[ci+k] if x in 'ACGT')
            m = sum(al.values()); out[seg]['sites'] += 1
            if m < 2 or len(al) < 2: continue
            out[seg]['S'] += 1
            # pi contribution (per site, unbiased)
            h = 1 - sum((v/m)**2 for v in al.values())
            out[seg]['pi_num'] += h * m/(m-1)
            # classify by translating the codon with each allele
            majors = sorted(al, key=lambda x: -al[x])
            a1, a2 = majors[0], majors[1]
            c1 = cr[:k] + a1 + cr[k+1:]; c2 = cr[:k] + a2 + cr[k+1:]
            if CO.get(c1) == CO.get(c2): out[seg]['syn'] += 1
            else: out[seg]['nsyn'] += 1
    return out

groups = [('all isolates', list(usable)),
          ('carriers (2 copies)', [k for k in usable if carrier.get(k) == 2]),
          ('non-carriers (1 copy)', [k for k in usable if carrier.get(k) == 1])]
say('')
say(f'{"group":<24} {"n":>4} {"segment":<16} {"S":>4} {"pi":>9} {"nonsyn":>7} {"syn":>5} '
    f'{"piN":>8} {"piS":>8} {"piN/piS":>8}')
say('-' * 104)
tsv = [('group', 'n', 'segment', 'S', 'pi', 'nonsyn', 'syn', 'piN', 'piS', 'piN_piS')]
for gname, keys in groups:
    if len(keys) < 2:
        say(f'{gname:<24} {len(keys):>4}  (too few to analyse)'); continue
    o = per_seg(keys)
    for seg, lab in (('N', 'N-half 1-405'), ('C', 'C-half 406-692')):
        d = o[seg]
        pi = d['pi_num']/d['sites'] if d['sites'] else 0
        piN = d['nsyn']/d['Nsit'] if d['Nsit'] else 0
        piS = d['syn']/d['Ssit'] if d['Ssit'] else 0
        ratio = piN/piS if piS else float('nan')
        say(f'{gname if seg=="N" else "":<24} {len(keys) if seg=="N" else "":>4} {lab:<16} '
            f'{d["S"]:>4} {pi:>9.5f} {d["nsyn"]:>7} {d["syn"]:>5} {piN:>8.5f} {piS:>8.5f} {ratio:>8.3f}')
        tsv.append((gname, len(keys), lab, d['S'], f'{pi:.6f}', d['nsyn'], d['syn'],
                    f'{piN:.6f}', f'{piS:.6f}', f'{ratio:.4f}'))
    say('')

# ── what the mapping-based call set says over the same codons ────────────────
say('-' * 104)
say('THE SAME CODONS, AS SEEN BY SHORT READS MAPPED TO THE DUPLICATED WBcel235 REFERENCE')
say('-' * 104)
regs = ','.join(f'X:{a}-{b}' for a, b in EXONS)
r = subprocess.run(['bcftools', 'query', '-r', regs, '-f', '%POS\n', VCF],
                   capture_output=True, text=True, timeout=900)
if r.returncode == 0:
    pos = [int(x) for x in r.stdout.strip().split('\n') if x]
    nN = sum(1 for p in pos if p < 4324178); nC = len(pos) - nN
    exN = sum(b-a+1 for a, b in EXONS if b < 4324178); exC = sum(b-a+1 for a, b in EXONS if a > 4324178)
    say(f'CaeNDR 20250625 hard-filter isotype VCF, dop-6 exons only, 684 isotypes:')
    say(f'  N-half exons ({exN} bp, single copy in every strain): {nN} variant sites')
    say(f'  C-half exons ({exC} bp, duplicated in carriers):      {nC} variant sites')
    tsv.append(('CaeNDR VCF (684 isotypes)', 684, 'N-half exons', nN, '', '', '', '', '', ''))
    tsv.append(('CaeNDR VCF (684 isotypes)', 684, 'C-half exons', nC, '', '', '', '', '', ''))
else:
    say(f'bcftools query failed: {r.stderr[-300:]}')
say('')
say('The gap between these two blocks is the artefact: the assemblies recover real variation in')
say('the C-half from the same natural isolates, because each isolate is compared in its own')
say('coordinate system instead of being forced onto a reference that carries two identical copies.')

with open(f'{OUT}/dop6_polymorphism.tsv', 'w') as f:
    for t in tsv: f.write('\t'.join(str(x) for x in t) + '\n')
open(f'{OUT}/dop6_polymorphism.txt', 'w').write('\n'.join(rep) + '\n')
say('-> output/dop6_polymorphism.tsv, output/dop6_polymorphism.txt')
