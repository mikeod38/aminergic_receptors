#!/usr/bin/env python3
"""fetch_locus.py — pull the dop-6 / C24A8.6 locus (X:4,318,000-4,334,000) and every
transcript model in it from Ensembl (WBcel235 / WormBase), and cache to ./output.

Outputs: output/locus_X_4.318-4.334Mb.fa, output/gene_models.json,
         output/dop6_c24a8_6_proteins.faa, output/dop6_c24a8_6_cds.fna
"""
import json, os, sys, time, urllib.request, urllib.error

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'output')
os.makedirs(OUT, exist_ok=True)
ENS = 'https://rest.ensembl.org'
GENES = {'dop-6': 'WBGene00016037', 'C24A8.5': 'WBGene00016039', 'C24A8.6': 'WBGene00044004',
         'cst-2': 'WBGene00016038'}
REG = ('X', 4318000, 4334500)

def get(url, ctype='application/json', tries=5):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={'Content-Type': ctype})
            with urllib.request.urlopen(req, timeout=60) as r:
                b = r.read().decode()
            return json.loads(b) if ctype == 'application/json' else b.strip()
        except Exception as e:
            if k == tries - 1: raise
            time.sleep(2 * (k + 1))

# ── locus genomic sequence ───────────────────────────────────────────────────
chrom, a, b = REG
seq = get(f'{ENS}/sequence/region/caenorhabditis_elegans/{chrom}:{a}..{b}:1?content-type=text/plain', 'text/plain')
with open(f'{OUT}/locus_X_{a}-{b}.fa', 'w') as f:
    f.write(f'>X:{a}-{b} WBcel235 plus strand\n')
    for i in range(0, len(seq), 60): f.write(seq[i:i+60] + '\n')
print(f'locus {chrom}:{a}-{b}  {len(seq)} nt -> output/locus_X_{a}-{b}.fa')

models = {}
prots, cdss = [], []
for name, gid in GENES.items():
    g = get(f'{ENS}/lookup/id/{gid}?expand=1;content-type=application/json')
    rec = {'name': name, 'id': gid, 'chr': g['seq_region_name'], 'start': g['start'],
           'end': g['end'], 'strand': g['strand'], 'biotype': g['biotype'],
           'description': g.get('description'), 'transcripts': []}
    for t in g.get('Transcript', []):
        tr = {'id': t['id'], 'biotype': t['biotype'], 'start': t['start'], 'end': t['end'],
              'canonical': bool(t.get('is_canonical')),
              'exons': sorted([(e['start'], e['end']) for e in t.get('Exon', [])]),
              'protein_length': (t.get('Translation') or {}).get('length')}
        if t['biotype'] == 'protein_coding':
            p = get(f'{ENS}/sequence/id/{t["id"]}?type=protein;content-type=text/plain', 'text/plain')
            c = get(f'{ENS}/sequence/id/{t["id"]}?type=cds;content-type=text/plain', 'text/plain')
            tr['protein'] = p; tr['cds_len'] = len(c)
            prots.append((f'{name}|{t["id"]}', p)); cdss.append((f'{name}|{t["id"]}', c))
        rec['transcripts'].append(tr)
    models[name] = rec
    print(f'{name:<9} {gid} {rec["chr"]}:{rec["start"]}-{rec["end"]} strand {rec["strand"]:+d} '
          f'{len(rec["transcripts"])} transcript(s)')
    for t in rec['transcripts']:
        print(f'   {t["id"]:<16} {t["biotype"]:<16} {len(t["exons"])} exons  '
              f'aa={t["protein_length"]}  {"CANONICAL" if t["canonical"] else ""}')

json.dump(models, open(f'{OUT}/gene_models.json', 'w'), indent=1)
def wfa(path, recs):
    with open(path, 'w') as f:
        for h, s in recs:
            f.write(f'>{h}\n')
            for i in range(0, len(s), 60): f.write(s[i:i+60] + '\n')
wfa(f'{OUT}/dop6_c24a8_6_proteins.faa', prots)
wfa(f'{OUT}/dop6_c24a8_6_cds.fna', cdss)
print(f'\n{len(prots)} proteins -> output/dop6_c24a8_6_proteins.faa')
