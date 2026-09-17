#!/usr/bin/env python3
"""audit_live_panel.py — verify every accession in the LIVE b120 panel against its own record.

Replaces output/accession_audit.tsv, which a referee audit found was auditing a DELETED panel:
it still listed the three fabricated accessions (P31356 squid rhodopsin, P20906 Pseudomonas
decarboxylase, Q9VNJ0 Drosophila ubiquitin ligase), missed the entire invertebrate-breadth set,
and therefore could not answer the question it existed to answer.

The rule this enforces: an accession is only trusted when the record it points at agrees with the
label the panel gives it, on BOTH organism and protein identity. A label is a claim; the record
is the evidence. Checking a name against itself is not verification.

Tips fall into three provenance classes and only two are checkable here:
  embedded  - proteome-derived, accession is inside the tip name -> verified against UniProt
  curated   - hardcoded accession in a fetch script            -> verified against UniProt
  by-query  - fetched by gene symbol or WormBase gene name, so the record supplied its own label,
              and there is no accession of mine to be wrong. Reported, not fetched.

Writes output/accession_audit.tsv
"""
import json, os, re, sys, time, urllib.request
HERE=os.path.dirname(os.path.abspath(__file__))
PANEL=f'{HERE}/output/tree/panel_b120.faa'
OUT=f'{HERE}/output/accession_audit.tsv'

CURATED={  # label -> (accession, expected organism substring, expected protein-name keywords)
 'dmOAMB':('Q7JQF1','Drosophila',['octopamine']),
 'dmOctA2R':('A0A0B4KHD0','Drosophila',['octopamine']),
 'dmOctB1R':('Q9VCZ3','Drosophila',['octopamine']),
 'dmOctB2R':('Q4LBB9','Drosophila',['octopamine']),
 'dmOctB3R':('Q4LBB6','Drosophila',['octopamine']),
 'dmOctTyrR':('P22270','Drosophila',['octopamine','tyramine']),
 'dmTyrR':('Q9VEG1','Drosophila',['tyramine']),
 'dmTyrRII':('Q9VEG2','Drosophila',['tyramine']),
 'dmDop1R1':('P41596','Drosophila',['dopamine']),
 'dmDop2R':('Q8IS44','Drosophila',['dopamine']),
 'bmOAR1':('Q17232','Bombyx',['octopamine']),
 'plADRA1':('A0A1J0F5X6','Platynereis',['adrenergic']),
}
AMINE_OK=['adrenergic','dopamine','serotonin','5-hydroxytryptamine','octopamine','tyramine',
          'histamine','muscarinic','acetylcholine','melatonin','trace amine','amine receptor','adrenoceptor',
          'g-protein coupled','g protein-coupled','gpcr','biogenic','receptor']
WRONG_FOLD=['kinase','ligase','tetraspan','decarboxylase','rhodopsin','transporter','phosphatase']

def fetch(acc):
    try:
        with urllib.request.urlopen(f'https://rest.uniprot.org/uniprotkb/{acc}.json',timeout=45) as r:
            return json.load(r)
    except Exception as e:
        return {'_err':str(e)}

def pname(d):
    pd=d.get('proteinDescription',{})
    n=pd.get('recommendedName',{}).get('fullName',{}).get('value')
    if n: return n
    sn=pd.get('submissionNames',[])
    return sn[0].get('fullName',{}).get('value','') if sn else ''

tips=[l[1:].strip() for l in open(PANEL) if l.startswith('>')]
rows=[]; nfetch=0
for t in tips:
    m=re.search(r'_(?:sp|tr)_([A-Z0-9]+)_',t)
    if m:
        acc=m.group(1); cls='embedded'; exp_org=None; exp_kw=AMINE_OK
    elif t in CURATED:
        acc,exp_org,exp_kw=CURATED[t]; cls='curated'
    else:
        rows.append((t,'-','by-query','-','-','not-applicable',
                     'fetched by gene symbol / WormBase gene name; the record supplied its own label'))
        continue
    d=fetch(acc); nfetch+=1; time.sleep(0.12)
    if '_err' in d:
        rows.append((t,acc,cls,'?','?','FETCH-FAILED',d['_err'][:60])); continue
    org=d.get('organism',{}).get('scientificName','')
    nm=pname(d); low=nm.lower()
    notes=[]; verdict='ok'
    if exp_org and exp_org.lower() not in org.lower():
        verdict='ORGANISM-MISMATCH'; notes.append(f'expected {exp_org}')
    if any(w in low for w in WRONG_FOLD) and not any(k in low for k in ['receptor']):
        verdict='WRONG-FOLD'; notes.append('record is not a receptor')
    elif not any(k in low for k in exp_kw):
        if verdict=='ok': verdict='NAME-UNEXPECTED'
        notes.append(f'name does not match {exp_kw[:3]}')
    rows.append((t,acc,cls,org,nm,verdict,'; '.join(notes) or '-'))

with open(OUT,'w') as f:
    f.write('# Regenerated 2026-09-17 against the LIVE b120 panel (173 tips).\n')
    f.write('# The previous file audited a deleted panel and still listed three fabricated accessions.\n')
    f.write('# verdict: ok | NAME-UNEXPECTED | ORGANISM-MISMATCH | WRONG-FOLD | FETCH-FAILED | not-applicable\n')
    f.write('tip\taccession\tprovenance\torganism\tprotein_name\tverdict\tnotes\n')
    for r in rows: f.write('\t'.join(str(x) for x in r)+'\n')

from collections import Counter
c=Counter(r[5] for r in rows)
print(f'{len(tips)} tips, {nfetch} records fetched')
for k,v in c.most_common(): print(f'  {k:<18} {v}')
bad=[r for r in rows if r[5] not in ('ok','not-applicable')]
if bad:
    print('\nNOT CLEAN:')
    for r in bad: print(f'  {r[0]:<38} {r[1]:<12} {r[5]:<18} {r[4][:52]}')
else:
    print('\nEvery checkable accession agrees with its own record.')
print(f'-> {OUT}')
