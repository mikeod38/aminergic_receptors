#!/usr/bin/env python3
"""fetch_tree_seqs_v2.py — a panel dense enough to resolve the aminergic clusters.

v1 was too thin on vertebrates to answer the question it was being asked. The canonical amine
receptor topology splits largely by G protein — ADRA1 + 5-HT2 (Gq), ADRA2 + DRD2-like + 5-HT1
(Gi/o), ADRB + DRD1-like + 5-HT4/6/7 (Gs) — and v1 carried 2 of 13 human 5-HT receptors, 2 of 5
dopamine receptors, no histamine receptors and no TAAR. With one 5-HT2 and three alpha-1 there is
no way to resolve the Gq cluster, which is exactly the node DOP-6's placement was being read from.

v2 takes the complete human amine receptor complement plus the rest of the C. elegans set, so the
reference clusters can form on their own before any orphan is placed against them. HTR6 matters
particularly: Harris 2009 identified SER-5 as 5-HT6-like and v1 had no HTR6 at all.

Writes output/tree2/panel.faa, provenance.tsv, search_registry.tsv
"""
import json, os, time, urllib.parse, urllib.request
from datetime import date

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'output', 'tree2')
os.makedirs(OUT, exist_ok=True)
UP = 'https://rest.uniprot.org/uniprotkb'
ENS = 'https://rest.ensembl.org'
PEP = ('/private/tmp/claude-501/-Users-mikeodonnell-git-Malaiwong-Oglu-2026/'
       '93015e40-e331-495a-8362-20ac7ea07c75/scratchpad/cel_pep.fa')

HUMAN = {  # the complete human aminergic complement, by gene symbol -> ligand class
 'HTR1A':'5-HT','HTR1B':'5-HT','HTR1D':'5-HT','HTR1E':'5-HT','HTR1F':'5-HT',
 'HTR2A':'5-HT','HTR2B':'5-HT','HTR2C':'5-HT','HTR4':'5-HT','HTR5A':'5-HT',
 'HTR6':'5-HT','HTR7':'5-HT',
 'DRD1':'DA','DRD2':'DA','DRD3':'DA','DRD4':'DA','DRD5':'DA',
 'ADRA1A':'NE','ADRA1B':'NE','ADRA1D':'NE','ADRA2A':'NE','ADRA2B':'NE','ADRA2C':'NE',
 'ADRB1':'NE','ADRB2':'NE','ADRB3':'NE',
 'HRH1':'HA','HRH2':'HA','HRH3':'HA','HRH4':'HA',
 'TAAR1':'TA',
 # melatonin receptors: an N-acetylated indoleamine branch of the same cluster. Included because
 # DOP-5 — DOP-6's sister in the v1 tree — carries G rather than the amine-binding D at position
 # 3.32, which is the melatonin-receptor state, so this is a hypothesis-driven addition rather
 # than generic breadth. GPR50 is the melatonin-related orphan.
 'MTNR1A':'MEL','MTNR1B':'MEL','GPR50':'MEL',
 'CHRM1':'ACh','CHRM2':'ACh','CHRM3':'ACh','CHRM5':'ACh',
}
# Accessions below were each verified against their own UniProt record (organism + protein name)
# by verify_accessions.py. Three earlier entries were recalled rather than retrieved and were
# badly wrong - P31356 is a squid rhodopsin, P20906 a Pseudomonas decarboxylase and Q9VNJ0 a
# Drosophila ubiquitin ligase, all three labelled as fly amine receptors. They were caught by the
# position-3.32 positive control, not by any naming or length check. Nothing here is hardcoded
# from memory any more: the remaining fly serotonin receptors are fetched by query and filtered on
# organism AND on their own protein name.
FLY = {
 'dmOAMB':'Q7JQF1','dmOctA2R':'A0A0B4KHD0','dmOctB1R':'Q9VCZ3','dmOctB2R':'Q4LBB9',
 'dmOctB3R':'Q4LBB6','dmOctTyrR':'P22270','dmTyrR':'Q9VEG1','dmTyrRII':'Q9VEG2',
 'dmDop1R1':'P41596','dmDop2R':'Q8IS44',
 # NB no fly 5-HT receptors hardcoded here - they come from FLY_5HT_QUERY below, labelled from
 # their own gene names, because the stale 'dm5HT1A':'P20905' entry duplicated the 5-HT7 record.
}
# Drosophila 5-HT receptor nomenclature does not map cleanly onto UniProt entry names (a
# gene_exact:5-HT1B query returns an entry whose own name is "receptor 2B"; 5-HT7 and 5-HT1A
# queries return the same accession). So rather than asserting an identity per label and then
# discovering the record disagrees, all Drosophila serotonin receptors are pulled in one query and
# each is LABELLED FROM ITS OWN RECORD - gene name plus accession. The panel then says what the
# database says, with no claim of mine layered on top.
FLY_5HT_QUERY = 'organism_id:7227 AND protein_name:"5-hydroxytryptamine"'
OTHER = {'bmOAR1':'Q17232', 'plADRA1':'A0A1J0F5X6'}
CEL = ['dop-1','dop-2','dop-3','dop-4','dop-5','dop-6','ser-1','ser-2','ser-3','ser-4','ser-5',
       'ser-6','ser-7','octr-1','tyra-2','tyra-3','gar-1','gar-2','gar-3']
CEL_SEQNAME = {'cel_F35H10.10':'F35H10.10','cel_F59D12.1':'F59D12.1','cel_T21B4.4':'T21B4.4'}

def get(u, tries=5, js=True):
    for k in range(tries):
        try:
            hdr = {'Accept': 'application/json'} if js else {'Content-Type': 'text/plain'}
            with urllib.request.urlopen(urllib.request.Request(u, headers=hdr), timeout=90) as r:
                t = r.read().decode()
            return json.loads(t) if js else t.strip()
        except Exception:
            if k == tries-1: raise
            time.sleep(2*(k+1))
def pname(r):
    d = r.get('proteinDescription', {})
    return (d.get('recommendedName', {}).get('fullName', {}).get('value')
            or (d.get('submissionNames') or [{}])[0].get('fullName', {}).get('value', '?'))

today = date.today().isoformat(); seqs = {}; prov = []; reg = []
def log(*a): print(' '.join(str(x) for x in a), flush=True)

for sym, lig in HUMAN.items():
    q = f'gene_exact:{sym} AND organism_id:9606 AND reviewed:true'
    try:
        d = get(f'{UP}/search?query={urllib.parse.quote(q)}&format=json&size=5'
                f'&fields=accession,protein_name,gene_names,length,sequence')
    except Exception:
        reg.append((today,'uniprot',q,'ERROR','')); log(f'h{sym:<8} QUERY FAILED'); continue
    res = d.get('results', []); reg.append((today,'uniprot',q,str(len(res)),''))
    hit = next((r for r in res if any(g.get('geneName',{}).get('value','').upper()==sym
                                      for g in r.get('genes',[]))), res[0] if res else None)
    if not hit: log(f'h{sym:<8} NO HIT'); continue
    lab = 'h'+sym; seqs[lab] = hit['sequence']['value']
    prov.append((lab,'Homo sapiens',lig,hit['primaryAccession'],'uniprot:'+q,
                 str(hit['sequence']['length']),pname(hit)))
    log(f'{lab:<10} {hit["primaryAccession"]:<11} {hit["sequence"]["length"]:>4} aa  {pname(hit)[:44]}')

try:
    d = get(f'{UP}/search?query={urllib.parse.quote(FLY_5HT_QUERY)}&format=json&size=40'
            f'&fields=accession,protein_name,organism_name,length,sequence,gene_names')
    res = d.get('results', []); reg.append((today,'uniprot',FLY_5HT_QUERY,str(len(res)),''))
    seen=set()
    for r_ in res:
        acc=r_['primaryAccession']; L=r_['sequence']['length']
        if acc in seen or r_['organism'].get('taxonId')!=7227: continue
        if not (330 <= L <= 900): 
            log(f'  skip {acc} ({L} aa, outside the family length range): {pname(r_)[:44]}'); continue
        gn=''
        for g in r_.get('genes',[]):
            gn=g.get('geneName',{}).get('value',''); break
        lab='dm'+(gn or acc).replace('-','')
        if lab in seqs: continue
        seen.add(acc); seqs[lab]=r_['sequence']['value']
        prov.append((lab, r_['organism']['scientificName'], '5-HT', acc,
                     'uniprot:'+FLY_5HT_QUERY, str(L), pname(r_)))
        log(f'{lab:<12} {acc:<11} {L:>4} aa  gene={gn:<8} {pname(r_)[:40]}')
except Exception as e:
    log(f'fly 5-HT query failed: {e}')

for lab, acc in {**FLY, **OTHER}.items():
    try: r = get(f'{UP}/{acc}.json')
    except Exception: log(f'{lab:<10} FAILED'); continue
    seqs[lab] = r['sequence']['value']
    prov.append((lab, r['organism']['scientificName'], '', acc, f'uniprot accession {acc}',
                 str(r['sequence']['length']), pname(r)))
    reg.append((today,'uniprot',f'accession:{acc}','1',''))
    log(f'{lab:<10} {acc:<11} {r["sequence"]["length"]:>4} aa  {pname(r)[:44]}')

for g in CEL:
    try:
        gg = get(f'{ENS}/lookup/symbol/caenorhabditis_elegans/{g}?expand=1;content-type=application/json')
        t = ([x for x in gg['Transcript'] if x.get('is_canonical')] or gg['Transcript'])[0]
        p = get(f'{ENS}/sequence/id/{t["id"]}?type=protein;content-type=text/plain', js=False)
    except Exception: log(f'cel_{g:<7} FAILED'); continue
    seqs[f'cel_{g}'] = p
    prov.append((f'cel_{g}','Caenorhabditis elegans','',gg['id'],f'ensembl symbol:{g}',str(len(p)),
                 f'canonical {t["id"]}'))
    reg.append((today,'ensembl',f'symbol:{g}','1',''))
    log(f'cel_{g:<7} {gg["id"]:<16} {len(p):>4} aa')

import re
pep = {}
for r in open(PEP).read().split('>')[1:]:
    h, s = r.split('\n', 1); pep[h.split()[0]] = s.replace('\n','')
for lab, sn in CEL_SEQNAME.items():
    c = [(k,v) for k,v in pep.items() if re.match(rf'^{re.escape(sn)}([a-z])?\.\d+$', k)]
    if not c: log(f'{lab:<10} not found'); continue
    k, v = max(c, key=lambda kv: len(kv[1])); seqs[lab] = v
    prov.append((lab,'Caenorhabditis elegans','',sn,'Ensembl WBcel235 proteome',str(len(v)),
                 'Koelle 2018 orphan'))
    log(f'{lab:<10} {k:<16} {len(v):>4} aa')
# dop-6 orthologs carried over from v1 (minus the C. japonica fragment)
v1 = os.path.join(HERE,'output','dop6_orthologs.faa')
for r in open(v1).read().split('>')[1:]:
    h, s = r.split('\n',1); lab = h.strip().split('|')[0]
    if lab in ('Cel-dop-6','Cja-dop-6'): continue
    seqs[lab.replace('-','_')] = s.replace('\n','')
    prov.append((lab.replace('-','_'), lab.split('-')[0], '', h.strip().split('|')[1],
                 'Ensembl Compara ortholog of WBGene00016037', str(len(s.replace('\n',''))), ''))

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
log(f'\n{len(seqs)} sequences -> output/tree2/panel.faa  (v1 had 49)')
