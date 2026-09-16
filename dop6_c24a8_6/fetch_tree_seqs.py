#!/usr/bin/env python3
"""fetch_tree_seqs.py — assemble the adrenergic-family GPCR panel for the dop-6 tree.

Every sequence is retrieved by a logged UniProt query (reviewed/Swiss-Prot preferred), never by
a remembered accession. Writes:
  output/tree/panel.faa            the sequences
  output/tree/provenance.tsv       id, organism, group, accession, query used, length, note
  output/tree/search_registry.tsv  every query submitted, with its result count
Groups follow the ligand classes the tree has to separate: adrenergic (vertebrate), octopamine,
tyramine, dopamine, serotonin, muscarinic (outgroup), plus the C. elegans orphans under test.
"""
import json, os, time, urllib.parse, urllib.request
from datetime import date

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'output', 'tree')
os.makedirs(OUT, exist_ok=True)
UP = 'https://rest.uniprot.org/uniprotkb/search'

# (label, group, uniprot query) — query is the provenance, not an accession
PANEL = [
    # ── vertebrate adrenergic ────────────────────────────────────────────────
    ('hADRA1A', 'adrenergic_a1', 'gene_exact:ADRA1A AND organism_id:9606 AND reviewed:true'),
    ('hADRA1B', 'adrenergic_a1', 'gene_exact:ADRA1B AND organism_id:9606 AND reviewed:true'),
    ('hADRA2A', 'adrenergic_a2', 'gene_exact:ADRA2A AND organism_id:9606 AND reviewed:true'),
    ('hADRA2C', 'adrenergic_a2', 'gene_exact:ADRA2C AND organism_id:9606 AND reviewed:true'),
    ('hADRB1',  'adrenergic_b',  'gene_exact:ADRB1 AND organism_id:9606 AND reviewed:true'),
    ('hADRB2',  'adrenergic_b',  'gene_exact:ADRB2 AND organism_id:9606 AND reviewed:true'),
    # ── vertebrate dopamine ──────────────────────────────────────────────────
    ('hDRD1', 'dopamine_v', 'gene_exact:DRD1 AND organism_id:9606 AND reviewed:true'),
    ('hDRD2', 'dopamine_v', 'gene_exact:DRD2 AND organism_id:9606 AND reviewed:true'),
    # ── vertebrate serotonin + muscarinic outgroup ───────────────────────────
    ('hHTR1A', 'serotonin_v',  'gene_exact:HTR1A AND organism_id:9606 AND reviewed:true'),
    ('hHTR2A', 'serotonin_v',  'gene_exact:HTR2A AND organism_id:9606 AND reviewed:true'),
    ('hCHRM1', 'muscarinic',   'gene_exact:CHRM1 AND organism_id:9606 AND reviewed:true'),
    ('hCHRM3', 'muscarinic',   'gene_exact:CHRM3 AND organism_id:9606 AND reviewed:true'),
    # ── Drosophila octopamine / tyramine (the invertebrate adrenergic set) ───
    ('dmOAMB',     'octopamine_a', 'gene_exact:Oamb AND organism_id:7227 AND reviewed:true'),
    ('dmOctA2R',   'octopamine_a', 'gene:Octalpha2R AND organism_id:7227'),
    ('dmOctB1R',   'octopamine_b', 'gene:Octbeta1R AND organism_id:7227'),
    ('dmOctB2R',   'octopamine_b', 'gene:Octbeta2R AND organism_id:7227'),
    ('dmOctB3R',   'octopamine_b', 'gene:Octbeta3R AND organism_id:7227'),
    ('dmTAR1',     'tyramine',     'gene_exact:TyrR AND organism_id:7227 AND reviewed:true'),
    ('dmTAR2',     'tyramine',     'gene:TyrRII AND organism_id:7227'),
    ('dmTAR3',     'tyramine',     'gene:TyrRIII AND organism_id:7227'),
    # ── Drosophila dopamine ──────────────────────────────────────────────────
    ('dmDop1R1', 'dopamine_i', 'gene_exact:Dop1R1 AND organism_id:7227 AND reviewed:true'),
    ('dmDop2R',  'dopamine_i', 'gene_exact:Dop2R AND organism_id:7227 AND reviewed:true'),
    # ── other protostome adrenergic (Bauknecht & Jekely-style breadth) ───────
    ('plADRA',  'adrenergic_inv', 'gene:adra AND organism_name:"Platynereis dumerilii"'),
    ('plADRB',  'adrenergic_inv', 'gene:adrb AND organism_name:"Platynereis dumerilii"'),
]
# C. elegans receptors come from Ensembl (WormBase) instead, by gene symbol
CEL = ['dop-1', 'dop-2', 'dop-3', 'dop-4', 'dop-5', 'dop-6', 'ser-1', 'ser-2', 'ser-3', 'ser-4',
       'ser-5', 'ser-6', 'ser-7', 'octr-1', 'tyra-2', 'tyra-3', 'gar-3']
ENS = 'https://rest.ensembl.org'

def get(u, c='application/json', tries=5):
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(u, headers={'Content-Type': c}), timeout=90) as r:
                t = r.read().decode()
            return json.loads(t) if c == 'application/json' else t.strip()
        except Exception:
            if k == tries - 1: raise
            time.sleep(2 * (k + 1))

today = date.today().isoformat()
prov, reg, seqs = [], [], []
def log(*a): print(' '.join(str(x) for x in a), flush=True)

for label, group, q in PANEL:
    url = f'{UP}?query={urllib.parse.quote(q)}&format=json&size=5&fields=accession,id,protein_name,organism_name,length,sequence,reviewed'
    try:
        d = get(url)
    except Exception as e:
        reg.append((today, 'uniprot', q, 'ERROR', str(e)[:60])); log(f'{label:<10} QUERY FAILED'); continue
    res = d.get('results', [])
    reg.append((today, 'uniprot', q, str(len(res)), ''))
    if not res:
        log(f'{label:<10} NO HIT   query="{q}"'); continue
    rev = [r for r in res if r.get('entryType', '').startswith('UniProtKB reviewed')]
    r = (rev or res)[0]
    acc = r['primaryAccession']; sq = r['sequence']['value']
    org = r['organism']['scientificName']
    pname = (r.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value')
             or r.get('proteinDescription', {}).get('submissionNames', [{}])[0].get('fullName', {}).get('value', '?'))
    seqs.append((f'{label}|{acc}', sq))
    prov.append((label, org, group, acc, 'uniprot:' + q, len(sq),
                 ('reviewed' if rev else 'unreviewed') + f'; {pname}'))
    log(f'{label:<10} {acc:<10} {len(sq):>4} aa  {org[:28]:<28} {"SwissProt" if rev else "TrEMBL"}')

for name in CEL:
    try:
        g = get(f'{ENS}/lookup/symbol/caenorhabditis_elegans/{name}?expand=1;content-type=application/json')
        t = ([x for x in g['Transcript'] if x.get('is_canonical')] or g['Transcript'])[0]
        p = get(f'{ENS}/sequence/id/{t["id"]}?type=protein;content-type=text/plain', 'text/plain')
    except Exception as e:
        reg.append((today, 'ensembl', f'symbol:{name}', 'ERROR', str(e)[:60])); log(f'cel_{name:<8} FAILED'); continue
    reg.append((today, 'ensembl', f'symbol:{name}', '1', ''))
    grp = {'dop-6': 'orphan_under_test'}.get(name,
          'dopamine_cel' if name.startswith('dop') else
          'muscarinic' if name.startswith('gar') else
          'tyramine' if name.startswith('tyra') else
          'octopamine_cel' if name == 'octr-1' else 'serotonin_cel')
    seqs.append((f'cel_{name}|{t["id"]}', p))
    prov.append((f'cel_{name}', 'Caenorhabditis elegans', grp, g['id'], f'ensembl symbol:{name}',
                 len(p), f'canonical transcript {t["id"]}'))
    log(f'cel_{name:<8} {g["id"]:<18} {len(p):>4} aa')

# C24A8.6 and the nematode dop-6 orthologs already fetched by earlier scripts
extra = os.path.join(HERE, 'output', 'dop6_orthologs.faa')
if os.path.exists(extra):
    for r in open(extra).read().split('>')[1:]:
        h, s = r.split('\n', 1); s = s.replace('\n', '')
        lab = h.strip().split('|')[0]
        if lab == 'Cel-dop-6': continue
        seqs.append((h.strip().replace('|', '_'), s))
        prov.append((lab, lab.split('-')[0], 'dop6_ortholog', h.strip().split('|')[1],
                     'ensembl compara ortholog of WBGene00016037', len(s), 'see crossspecies.txt'))
p6 = os.path.join(HERE, 'output', 'dop6_c24a8_6_proteins.faa')
for r in open(p6).read().split('>')[1:]:
    h, s = r.split('\n', 1)
    if h.startswith('C24A8.6|C24A8.6b'):
        seqs.append(('cel_C24A8.6b|C24A8.6b.1', s.replace('\n', '')))
        prov.append(('cel_C24A8.6b', 'Caenorhabditis elegans', 'duplicate_fragment', 'WBGene00044004',
                     'ensembl WBGene00044004', len(s.replace('\n', '')),
                     'DOP-6 406-692 + 17 aa cryptic exon; TM6-TM7 fragment, no ligand pocket'))

with open(f'{OUT}/panel.faa', 'w') as f:
    for h, s in seqs:
        f.write(f'>{h}\n')
        for i in range(0, len(s), 60): f.write(s[i:i+60] + '\n')
with open(f'{OUT}/provenance.tsv', 'w') as f:
    f.write('id\torganism\tgroup\taccession\tsource_query\tlength\tnote\n')
    for r in prov: f.write('\t'.join(str(x) for x in r) + '\n')
with open(f'{OUT}/search_registry.tsv', 'w') as f:
    f.write('date\tdatabase\tquery\tn_results\terror\n')
    for r in reg: f.write('\t'.join(str(x) for x in r) + '\n')
log(f'\n{len(seqs)} sequences -> output/tree/panel.faa')
log(f'{len(prov)} provenance rows, {len(reg)} logged queries')
