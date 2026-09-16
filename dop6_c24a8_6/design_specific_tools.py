#!/usr/bin/env python3
"""design_specific_tools.py — reagents that hit dop-6 without touching C24A8.6, and vice versa.

The palindrome (X:4,324,178-4,328,106 == revcomp X:4,329,303-4,333,231, 0 mismatches) means every
site in dop-6 exons 10-17 has a second perfect copy in C24A8.6. Only three parts of the locus are
single-copy: dop-6 exons 1-9 + introns (X:4,320,863-4,324,177), the spacer (X:4,328,107-4,329,302),
and the flanks. This script enumerates reagents in those regions and verifies each against the whole
WBcel235 genome by BLAST.

  1. SpCas9 guides (N20-NGG) in dop-6 unique exons, genome-unique, scored for position in the ORF
  2. RT-PCR / qPCR primer pairs whose amplicon lies wholly in dop-6 unique exons and spans an intron
  3. diagnostic PCRs across the four palindrome junctions, to confirm the structure in a real strain
  4. an assessment of whether C24A8.6 can be deleted on its own

Needs a BLAST db of WBcel235 (path in GENOME below). Writes output/specific_tools.txt and
output/dop6_guides.tsv
"""
import json, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'output')
GENOME = os.environ.get('CE_BLASTDB',
    '/private/tmp/claude-501/-Users-mikeodonnell-git-Malaiwong-Oglu-2026/'
    '93015e40-e331-495a-8362-20ac7ea07c75/scratchpad/genome/ce')
L0 = 4318000
ARM_A = (4324178, 4328106); ARM_B = (4329303, 4333231)
SPACER = (4328107, 4329302)
seq = ''.join(l.strip() for l in open(f'{OUT}/locus_X_4318000-4334500.fa') if not l.startswith('>'))
G = lambda a, b: seq[a - L0:b - L0 + 1].upper()
RC = lambda s: s[::-1].translate(str.maketrans('ACGTacgt', 'TGCAtgca'))
models = json.load(open(f'{OUT}/gene_models.json'))
dop6 = models['dop-6']['transcripts'][0]
EXONS = dop6['exons']                                    # sorted ascending, + strand gene
UNIQUE_EXONS = [(s, e) for s, e in EXONS if e < ARM_A[0]]   # exons 1-9

rep = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); rep.append(s)

def blast_count(q, name='q'):
    """perfect (100 % id, full-length) hits of q in the genome, both strands"""
    fa = f'/tmp/{name}.fa'
    open(fa, 'w').write(f'>{name}\n{q}\n')
    r = subprocess.run(['blastn', '-query', fa, '-db', GENOME, '-task', 'blastn-short' if len(q) < 30 else 'blastn',
                        '-dust', 'no', '-evalue', '1000', '-perc_identity', '100', '-outfmt',
                        '6 sseqid sstart send pident length mismatch'], capture_output=True, text=True)
    hits = [l.split('\t') for l in r.stdout.strip().split('\n') if l]
    full = [h for h in hits if int(h[4]) == len(q) and int(h[5]) == 0]
    return full

say('=' * 92)
say('REAGENTS THAT DISTINGUISH dop-6 FROM C24A8.6')
say('=' * 92)
say(f'single-copy regions of the locus:')
say(f'  dop-6 5\' half   X:{models["dop-6"]["start"]}-{ARM_A[0]-1}   '
    f'{ARM_A[0]-models["dop-6"]["start"]} bp   (exons 1-9, codons 1-405)')
say(f'  spacer          X:{SPACER[0]}-{SPACER[1]}   {SPACER[1]-SPACER[0]+1} bp   (C24A8.5)')
say(f'everything in X:{ARM_A[0]}-{ARM_A[1]} and X:{ARM_B[0]}-{ARM_B[1]} exists in two perfect copies.\n')

# ── 1. guides ────────────────────────────────────────────────────────────────
say('-' * 92); say('1. SpCas9 GUIDES IN dop-6 UNIQUE EXONS (protospacer + PAM, genome-verified)')
say('-' * 92)
cum = 0; guides = []
for i, (es, ee) in enumerate(UNIQUE_EXONS, 1):
    ex = G(es, ee)
    for strand, s in (('+', ex), ('-', RC(ex))):
        for m in re.finditer(r'(?=([ACGT]{21})([AG]G))', s):
            proto = m.group(1)[1:]                      # 20 nt
            pam = m.group(1)[-1] + m.group(2)           # NGG
            if proto.count('G') + proto.count('C') < 6: continue
            off = m.start() + 1
            gpos = es + off if strand == '+' else ee - off - 21
            # codon of the cut site (3 bp 5' of PAM)
            cut = gpos + 17 if strand == '+' else gpos + 5
            ncod = sum(min(e, cut) - s2 + 1 for s2, e in EXONS if s2 <= cut) // 3
            guides.append(dict(exon=i, strand=strand, proto=proto, pam=pam, start=gpos,
                               cut=cut, codon=ncod, gc=100*(proto.count('G')+proto.count('C'))/20))
say(f'{len(guides)} candidate protospacers in exons 1-9; BLASTing each against WBcel235 ...')
keep = []
for g in sorted(guides, key=lambda x: (x['exon'], x['start'])):
    hits = blast_count(g['proto'], 'g')
    g['nhits'] = len(hits)
    if len(hits) == 1: keep.append(g)
say(f'{len(keep)} are unique in the genome (1 perfect full-length hit).\n')
say(f'{"exon":>4} {"codon":>6} {"str":>4} {"protospacer (20 nt)":<22} {"PAM":<4} {"GC%":>5} {"cut at X:":>12}')
best = [g for g in keep if 2 <= g['exon'] <= 6 and 35 <= g['gc'] <= 75]
for g in (best or keep)[:14]:
    say(f'{g["exon"]:>4} {g["codon"]:>6} {g["strand"]:>4} {g["proto"]:<22} {g["pam"]:<4} '
        f'{g["gc"]:>5.0f} {g["cut"]:>12}')
with open(f'{OUT}/dop6_guides.tsv', 'w') as f:
    f.write('exon\tcodon\tstrand\tprotospacer\tpam\tgc_pct\tstart_X\tcut_X\tgenome_perfect_hits\n')
    for g in keep:
        f.write(f'{g["exon"]}\t{g["codon"]}\t{g["strand"]}\t{g["proto"]}\t{g["pam"]}\t{g["gc"]:.0f}\t'
                f'{g["start"]}\t{g["cut"]}\t{g["nhits"]}\n')
say(f'\n-> all {len(keep)} unique guides in output/dop6_guides.tsv')
say('NOTE: a guide anywhere in dop-6 exons 10-17 has TWO perfect targets (the second in C24A8.6);')
say('      two cuts 9 kb apart in inverted orientation delete or invert the intervening DNA.')

# how many guides would a naive exon 10-17 design have?
n_amb = 0
for (es, ee) in [(s, e) for s, e in EXONS if s >= ARM_A[0]]:
    for strand, s in (('+', G(es, ee)), ('-', RC(G(es, ee)))):
        n_amb += len(re.findall(r'(?=[ACGT]{21}[AG]G)', s))
say(f'      (for reference, exons 10-17 contain {n_amb} protospacers, every one of them doubled)')

# ── 2. RT-PCR / qPCR primers ─────────────────────────────────────────────────
say('\n' + '-' * 92); say('2. dop-6-SPECIFIC RT-PCR / qPCR PRIMERS (amplicon inside exons 1-9, spans an intron)')
say('-' * 92)
def tm(p):
    """nearest-neighbour Tm, 50 mM Na+, 250 nM oligo (SantaLucia 1998 unified)"""
    dH = {'AA':-7.9,'TT':-7.9,'AT':-7.2,'TA':-7.2,'CA':-8.5,'TG':-8.5,'GT':-8.4,'AC':-8.4,
          'CT':-7.8,'AG':-7.8,'GA':-8.2,'TC':-8.2,'CG':-10.6,'GC':-9.8,'GG':-8.0,'CC':-8.0}
    dS = {'AA':-22.2,'TT':-22.2,'AT':-20.4,'TA':-21.3,'CA':-22.7,'TG':-22.7,'GT':-22.4,'AC':-22.4,
          'CT':-21.0,'AG':-21.0,'GA':-22.2,'TC':-22.2,'CG':-27.2,'GC':-24.4,'GG':-19.9,'CC':-19.9}
    h = 0.1 if p[0] in 'GC' else 2.3; s = -2.8 if p[0] in 'GC' else 4.1
    h += 0.1 if p[-1] in 'GC' else 2.3; s += -2.8 if p[-1] in 'GC' else 4.1
    for i in range(len(p)-1): h += dH[p[i:i+2]]; s += dS[p[i:i+2]]
    import math
    return (1000*h)/(s + 1.987*math.log(250e-9/4)) - 273.15 + 16.6*math.log10(0.05)
# build the spliced mRNA of exons 1-9 with a map back to genome
mseq, mpos = '', []
for es, ee in UNIQUE_EXONS:
    mseq += G(es, ee); mpos += list(range(es, ee+1))
bounds = []; c = 0
for es, ee in UNIQUE_EXONS: c += ee-es+1; bounds.append(c)       # cumulative exon ends in mRNA
say(f'exons 1-9 spliced = {len(mseq)} nt; exon junctions in mRNA at {bounds[:-1]}')
pairs = []
for fs in range(0, len(mseq)-20):
    F = mseq[fs:fs+22]
    if not (40 <= 100*(F.count('G')+F.count('C'))/len(F) <= 60): continue
    tF = tm(F)
    if not (58 <= tF <= 62): continue
    for rs in range(fs+80, min(fs+260, len(mseq)-20)):
        R = RC(mseq[rs:rs+22]); tR = tm(R)
        if not (58 <= tR <= 62) or abs(tF-tR) > 1.5: continue
        amp = rs+22-fs
        if not any(fs < b < rs for b in bounds[:-1]): continue     # must span a junction
        pairs.append((abs(tF-tR), fs, rs, F, R, tF, tR, amp))
pairs.sort()
say(f'{len(pairs)} candidate pairs; checking the best few for genome uniqueness\n')
shown = 0
for d, fs, rs, F, R, tF, tR, amp in pairs:
    if shown >= 3: break
    if len(blast_count(F, 'F')) != 1 or len(blast_count(R, 'R')) != 1: continue
    j = [k+1 for k, b in enumerate(bounds[:-1]) if fs < b < rs]
    say(f'  pair {shown+1}: amplicon {amp} bp, spans exon junction(s) {j}')
    say(f'     F  5\'-{F}-3\'  Tm {tF:.1f} C   X:{mpos[fs]}')
    say(f'     R  5\'-{R}-3\'  Tm {tR:.1f} C   X:{mpos[rs+21]}')
    say(f'     genomic product would be larger by the introns spanned -> cDNA-specific')
    shown += 1

# ── 3. palindrome junction PCRs ──────────────────────────────────────────────
say('\n' + '-' * 92); say('3. DIAGNOSTIC PCRs ACROSS THE PALINDROME JUNCTIONS')
say('-' * 92)
say('The palindrome is a 3.9 kb inverted repeat with a 1.2 kb loop: it will not amplify or clone')
say('easily as one piece. Each junction is best tested with one primer in unique flanking DNA and')
say('one in the spacer, giving a product that exists only if the arm is present and in that orientation.\n')
for lab, a, b in [('proximal (dop-6 side, arm A start)', ARM_A[0]-400, ARM_A[0]+400),
                  ('armA -> spacer',                     ARM_A[1]-400, ARM_A[1]+400),
                  ('spacer -> armB',                     ARM_B[0]-400, ARM_B[0]+400),
                  ('distal (cst-2 side, arm B end)',     ARM_B[1]-400, ARM_B[1]+400)]:
    say(f'  {lab:<38} junction at X:{(a+b)//2}')
say('\n  Unique anchors available: X:<4,324,178 (dop-6 exons 1-9 side) and X:4,328,107-4,329,302')
say('  (spacer) and X:>4,333,231 (cst-2 side). Arm-internal primers are ambiguous by construction.')

# ── 4. can C24A8.6 be removed on its own? ────────────────────────────────────
say('\n' + '-' * 92); say('4. CAN C24A8.6 BE DELETED WITHOUT TOUCHING dop-6?')
say('-' * 92)
cst2 = models['cst-2']
ov = min(cst2['end'], ARM_B[1]) - max(cst2['start'], ARM_B[0]) + 1
say(f'C24A8.6 X:{models["C24A8.6"]["start"]}-{models["C24A8.6"]["end"]} lies wholly inside arm B, so no')
say('guide inside it is specific. The only specific cut sites flanking arm B are in the spacer')
say(f'(X:{SPACER[0]}-{SPACER[1]}, unique) and distal to X:{ARM_B[1]}.')
say(f'BUT arm B overlaps cst-2 by {ov} bp: cst-2 exon 9 runs X:{cst2["transcripts"][0]["exons"][0][0]}-'
    f'{cst2["transcripts"][0]["exons"][0][1]} and its 5\' {ARM_B[1]-cst2["start"]+1} bp are inside the arm.')
say('So a spacer-to-distal deletion removes all of C24A8.6 and leaves dop-6 intact, but it also')
say('truncates cst-2 unless the distal cut is placed inside the overlap and cst-2 is restored by HDR.')
say('')
say('Practical consequence: the decisive experiment is not a C24A8.6 null. It is a dop-6 null made')
say('in the unique 5\' half (leaving both palindrome arms untouched) plus cDNA rescue. That assigns')
say('the phenotype to dop-6 without ever having to separate the two copies.')

open(f'{OUT}/specific_tools.txt', 'w').write('\n'.join(rep) + '\n')
say('\n-> output/specific_tools.txt, output/dop6_guides.tsv')
