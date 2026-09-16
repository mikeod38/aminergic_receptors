#!/usr/bin/env python3
"""crossspecies_palindrome.py — is the dop-6 inverted duplication C. elegans-specific?

For each Caenorhabditis / outgroup nematode with a dop-6 ortholog (Ensembl Compara, see
output/homology_dop6.json), pull the ortholog's protein and its surrounding genomic region,
self-align the region to look for the same inverted repeat, and verify the reported identity
to C. elegans DOP-6 with an independent local alignment (Compara %id is asymmetric and
coverage-dependent, so it is not trusted here).

Writes output/crossspecies.txt, output/dop6_orthologs.faa.
"""
import json, os, subprocess, time, urllib.request
from Bio.Align import PairwiseAligner, substitution_matrices

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'output')
ENS = 'https://rest.ensembl.org'
FLANK = 12000
ORTHOLOGS = [  # species, gene id (from Compara), label
    ('caenorhabditis_briggsae',  'WBGene00034912', 'Cbr-dop-6'),
    ('caenorhabditis_brenneri',  'WBGene00158432', 'Cbn-dop-6'),
    ('caenorhabditis_remanei_gca010183535v1rs', 'GCK72_024233', 'Cre-dop-6'),
    ('caenorhabditis_japonica',  'WBGene00217654', 'Cja-dop-6'),
    ('pristionchus_pacificus',   'WBGene00089991', 'Ppa-dop-6'),
]
def get(url, ctype='application/json', tries=6):
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={'Content-Type': ctype}), timeout=90) as r:
                t = r.read().decode()
            return json.loads(t) if ctype == 'application/json' else t.strip()
        except Exception:
            if k == tries - 1: raise
            time.sleep(3 * (k + 1))
RC = lambda s: s[::-1].translate(str.maketrans('ACGTacgt', 'TGCAtgca'))

rep = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); rep.append(s)

dop6 = ''.join(open(f'{OUT}/dop6_c24a8_6_proteins.faa').read().split('>')[1].split('\n')[1:])
al = PairwiseAligner(); al.substitution_matrix = substitution_matrices.load('BLOSUM62')
al.open_gap_score, al.extend_gap_score, al.mode = -11, -1, 'global'
faa = [('Cel-dop-6|C24A8.1', dop6)]

say('=' * 90)
say('CROSS-SPECIES CHECK: dop-6 ortholog identity (independently realigned) and locus palindrome')
say('=' * 90)
for sp, gid, label in ORTHOLOGS:
    try:
        g = get(f'{ENS}/lookup/id/{gid}?expand=1;content-type=application/json')
    except Exception as e:
        say(f'\n{label:<10} {sp}: LOOKUP FAILED ({e})'); continue
    canon = [t for t in g.get('Transcript', []) if t.get('is_canonical')] or g.get('Transcript', [])
    t = canon[0]
    prot = get(f'{ENS}/sequence/id/{t["id"]}?type=protein;content-type=text/plain', 'text/plain')
    faa.append((f'{label}|{t["id"]}', prot))
    aln = al.align(dop6, prot)[0]
    A, B = str(aln[0]), str(aln[1])
    both = [(x, y) for x, y in zip(A, B) if x != '-' and y != '-']
    idn = sum(1 for x, y in both if x == y)
    say(f'\n{label:<10} {sp}')
    say(f'   {g["seq_region_name"]}:{g["start"]}-{g["end"]} strand {g["strand"]:+d}  '
        f'{t["id"]} {len(prot)} aa, {len(t.get("Exon",[]))} exons')
    say(f'   global realignment vs Cel DOP-6 (692 aa): {idn}/{len(both)} aligned positions identical '
        f'= {100*idn/len(both):.1f} %  (Compara reported a coverage-weighted number)')

    # locus self-alignment: same inverted repeat?
    a0, b0 = max(1, g['start'] - FLANK), g['end'] + FLANK
    try:
        seq = get(f'{ENS}/sequence/region/{sp}/{g["seq_region_name"]}:{a0}..{b0}:1?content-type=text/plain', 'text/plain')
    except Exception as e:
        say(f'   locus fetch FAILED ({e})'); continue
    fp = f'/tmp/{label}_locus.fa'
    with open(fp, 'w') as f:
        f.write(f'>{label} {g["seq_region_name"]}:{a0}-{b0}\n')
        for i in range(0, len(seq), 60): f.write(seq[i:i+60] + '\n')
    subprocess.run(['makeblastdb', '-in', fp, '-dbtype', 'nucl', '-out', f'/tmp/{label}db'],
                   check=True, capture_output=True)
    r = subprocess.run(['blastn', '-query', fp, '-db', f'/tmp/{label}db', '-dust', 'no', '-evalue', '1e-5',
                        '-outfmt', '6 qstart qend sstart send pident length sstrand'],
                       check=True, capture_output=True, text=True)
    hits = [l.split('\t') for l in r.stdout.strip().split('\n') if l]
    rep_hits = [h for h in hits if int(h[5]) > 200 and not (h[0] == '1' and int(h[1]) == len(seq))]
    say(f'   locus {g["seq_region_name"]}:{a0}-{b0} ({len(seq)} bp) self-blast: '
        f'{len(rep_hits)} repeat hit(s) >200 bp')
    for h in sorted(rep_hits, key=lambda x: -int(x[5]))[:6]:
        say(f'      q{h[0]}-{h[1]}  s{h[2]}-{h[3]}  {h[4]} % over {h[5]} bp  ({h[6]})')
    if not rep_hits: say('      -> NO internal repeat: the dop-6 locus is single-copy in this species')

with open(f'{OUT}/dop6_orthologs.faa', 'w') as f:
    for h, s in faa:
        f.write(f'>{h}\n')
        for i in range(0, len(s), 60): f.write(s[i:i+60] + '\n')
open(f'{OUT}/crossspecies.txt', 'w').write('\n'.join(rep) + '\n')
say(f'\n-> output/crossspecies.txt, output/dop6_orthologs.faa')
