#!/usr/bin/env python3
"""fetch_invert_breadth.py — invertebrate and nematode depth for the amine receptor panel.

v2 was 36 human, 14 Drosophila and 3 others, which is the wrong shape for placing a nematode
orphan: the reference clusters were vertebrate, and DOP-6's only close relatives were its own
orthologs. Two additions fix that.

  1. NEMATODE DEPTH — the most important. For each characterised C. elegans aminergic receptor,
     pull its orthologs in other nematodes from Ensembl Compara, the same route already used for
     dop-6. This lets nematode-specific groups resolve on their own rather than each worm gene
     hanging alone off a vertebrate scaffold.
  2. PROTOSTOME BREADTH — insects beyond Drosophila, plus molluscs and annelids, from UniProt.

Orthologs are filtered on length (within the family's own range) and on alignment coverage, so
fragments and mis-annotated models do not enter, per the curation rules.

Writes output/tree3/panel.faa, provenance.tsv, search_registry.tsv
"""
import json, os, time, urllib.parse, urllib.request
from datetime import date

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'output', 'tree3'); os.makedirs(OUT, exist_ok=True)
V2 = os.path.join(HERE, 'output', 'tree2')
ENS = 'https://rest.ensembl.org'; UP = 'https://rest.uniprot.org/uniprotkb'

# nematode species to accept orthologs from (well-assembled, spanning the phylum)
NEMATODES = {'caenorhabditis_briggsae':'Cbr', 'caenorhabditis_brenneri':'Cbn',
             'caenorhabditis_remanei_gca010183535v1rs':'Cre', 'pristionchus_pacificus':'Ppa',
             'ascaris_suum':'Asu', 'brugia_malayi':'Bma', 'strongyloides_ratti':'Sra',
             'haemonchus_contortus':'Hco', 'trichinella_spiralis':'Tsp',
             'necator_americanus_gca031761385v1rs':'Nam'}
# C. elegans receptors to seed the ortholog search with (characterised ones + the orphans)
SEEDS = ['dop-1','dop-2','dop-3','dop-4','dop-5','ser-1','ser-2','ser-3','ser-4','ser-6','ser-7',
         'octr-1','tyra-2','tyra-3']
MAXPER = 3          # at most this many orthologs per C. elegans gene, best coverage first
LEN_RANGE = (330, 900)

PROTOSTOME = [  # (label, uniprot query, what it should be)
 ('amOA1',   'organism_name:"Apis mellifera" AND protein_name:"Octopamine receptor" AND reviewed:true', 'Octopamine'),
 ('amTAR1',  'organism_name:"Apis mellifera" AND protein_name:"Tyramine receptor"', 'Tyramine'),
 ('am5HT1A', 'organism_name:"Apis mellifera" AND protein_name:"5-hydroxytryptamine receptor"', 'hydroxytryptamine'),
 ('amDop1',  'organism_name:"Apis mellifera" AND protein_name:"Dopamine receptor"', 'Dopamine'),
 ('bmTAR',   'organism_name:"Bombyx mori" AND protein_name:"Tyramine receptor"', 'Tyramine'),
 ('tcOA1',   'organism_name:"Tribolium castaneum" AND protein_name:"Octopamine receptor"', 'ctopamine'),
 ('apOA',    'organism_name:"Aplysia californica" AND protein_name:"Octopamine receptor"', 'ctopamine'),
 ('ap5HT',   'organism_name:"Aplysia californica" AND protein_name:"5-hydroxytryptamine receptor"', 'hydroxytryptamine'),
 ('lsOA',    'organism_name:"Lymnaea stagnalis" AND protein_name:"receptor"', 'receptor'),
 ('plADRA2', 'organism_name:"Platynereis dumerilii" AND protein_name:"adrenergic receptor"', 'adrenergic'),
 ('smTAR',   'organism_name:"Schistosoma mansoni" AND protein_name:"receptor" AND reviewed:true', 'receptor'),
]
def get(u, js=True, tries=5):
    for k in range(tries):
        try:
            hdr = {'Accept':'application/json'} if js else {'Content-Type':'text/plain'}
            with urllib.request.urlopen(urllib.request.Request(u, headers=hdr), timeout=90) as r:
                t = r.read().decode()
            return json.loads(t) if js else t.strip()
        except Exception:
            if k == tries-1: raise
            time.sleep(2*(k+1))
def pname(r):
    d = r.get('proteinDescription', {})
    return (d.get('recommendedName',{}).get('fullName',{}).get('value')
            or (d.get('submissionNames') or [{}])[0].get('fullName',{}).get('value','?'))

today = date.today().isoformat()
seqs, prov, reg = {}, [], []
def log(*a): print(' '.join(str(x) for x in a), flush=True)

# carry v2 over
for r in open(f'{V2}/panel.faa').read().split('>')[1:]:
    h, s = r.split('\n',1); seqs[h.strip()] = s.replace('\n','')
for ln in open(f'{V2}/provenance.tsv').read().strip().split('\n')[1:]:
    prov.append(ln.split('\t'))
log(f'carried {len(seqs)} sequences over from v2\n')

log('--- nematode orthologs via Ensembl Compara ---')
for g in SEEDS:
    try:
        gg = get(f'{ENS}/lookup/symbol/caenorhabditis_elegans/{g}?content-type=application/json')
        h = get(f'{ENS}/homology/id/caenorhabditis_elegans/{gg["id"]}'
                f'?content-type=application/json;type=orthologues')
    except Exception as e:
        reg.append((today,'ensembl',f'homology:{g}','ERROR','')); log(f'  {g:<8} FAILED'); continue
    hs = h['data'][0]['homologies']
    reg.append((today,'ensembl',f'homology:{g}',str(len(hs)),''))
    cands = []
    for x in hs:
        sp = x['target']['species']
        if sp not in NEMATODES: continue
        tl = len(x['target']['align_seq'].replace('-',''))
        if not (LEN_RANGE[0] <= tl <= LEN_RANGE[1]): continue
        a, b = x['source']['align_seq'], x['target']['align_seq']
        both = sum(1 for p,q in zip(a,b) if p!='-' and q!='-')
        cov = both / tl
        if cov < 0.70: continue
        cands.append((cov, sp, x['target']['id'], tl))
    cands.sort(reverse=True)
    taken = set(); n = 0
    for cov, sp, tid, tl in cands:
        if sp in taken or n >= MAXPER: continue
        try:
            p = get(f'{ENS}/sequence/id/{tid}?type=protein;content-type=text/plain', js=False)
        except Exception: continue
        if not (LEN_RANGE[0] <= len(p) <= LEN_RANGE[1]): continue
        lab = f'{NEMATODES[sp]}_{g}'
        if lab in seqs: continue
        seqs[lab] = p; taken.add(sp); n += 1
        prov.append([lab, sp, '', tid, f'ensembl compara ortholog of {g}', str(len(p)),
                     f'coverage {cov:.2f}'])
        log(f'  {lab:<16} {tid:<20} {len(p):>4} aa  cov={cov:.2f}')
    if n == 0: log(f'  {g:<8} no ortholog passed the length/coverage filter')

log('\n--- protostome breadth via UniProt ---')
for lab, q, kw in PROTOSTOME:
    try:
        d = get(f'{UP}/search?query={urllib.parse.quote(q)}&format=json&size=6'
                f'&fields=accession,protein_name,organism_name,length,sequence')
    except Exception:
        reg.append((today,'uniprot',q,'ERROR','')); log(f'  {lab:<10} FAILED'); continue
    res = d.get('results',[]); reg.append((today,'uniprot',q,str(len(res)),''))
    hit = next((r for r in res if kw.lower() in pname(r).lower()
                and LEN_RANGE[0] <= r['sequence']['length'] <= LEN_RANGE[1]), None)
    if not hit: log(f'  {lab:<10} no hit passing name+length filter'); continue
    if lab in seqs: continue
    seqs[lab] = hit['sequence']['value']
    prov.append([lab, hit['organism']['scientificName'], '', hit['primaryAccession'],
                 'uniprot:'+q, str(hit['sequence']['length']), pname(hit)])
    log(f'  {lab:<10} {hit["primaryAccession"]:<11} {hit["sequence"]["length"]:>4} aa  '
        f'{hit["organism"]["scientificName"][:22]:<22} {pname(hit)[:36]}')

with open(f'{OUT}/panel.faa','w') as f:
    for k in sorted(seqs):
        f.write(f'>{k}\n')
        for i in range(0,len(seqs[k]),60): f.write(seqs[k][i:i+60]+'\n')
with open(f'{OUT}/provenance.tsv','w') as f:
    f.write('id\torganism\tligand\taccession\tsource_query\tlength\tnote\n')
    for r in prov: f.write('\t'.join(str(x) for x in r)+'\n')
with open(f'{OUT}/search_registry.tsv','w') as f:
    f.write('date\tdatabase\tquery\tn_results\terror\n')
    for r in reg: f.write('\t'.join(str(x) for x in r)+'\n')
from collections import Counter
c = Counter('nematode' if any(k.startswith(v+'_') for v in NEMATODES.values()) else
            'human' if k.startswith('h') else
            'C. elegans' if k.startswith('cel_') else 'other invertebrate' for k in seqs)
log(f'\n{len(seqs)} sequences -> output/tree3/panel.faa')
for k,v in sorted(c.items()): log(f'   {k:<20} {v}')
