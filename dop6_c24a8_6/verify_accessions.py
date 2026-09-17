#!/usr/bin/env python3
"""verify_accessions.py — check every UniProt accession in the panel against its own record.

Two accessions in the v2 panel were hardcoded from memory rather than retrieved by a logged query,
and both were wrong in a way that no length or naming check would have caught:
  P31356, labelled "Drosophila 5-HT2A", is Todarodes pacificus RHODOPSIN
  P20906, labelled "Drosophila 5-HT7",  is Pseudomonas putida benzoylformate decarboxylase
A squid opsin and a bacterial enzyme were sitting in an amine receptor panel. The position-3.32
analysis flagged both as lacking the amine-binding aspartate, which read as a method failure until
the accessions were checked - the method was right and the data were wrong.

So: every accession is re-fetched and its own organism and protein name are checked against what
the label claims. Nothing is trusted because it was in a dict.

Writes output/accession_audit.tsv
"""
import json, os, time, urllib.request

HERE=os.path.dirname(os.path.abspath(__file__)); OUT=os.path.join(HERE,'output')
UP='https://rest.uniprot.org/uniprotkb'
PROV=f'{OUT}/tree2/provenance.tsv'
EXPECT_ORG={'h':'Homo sapiens','dm':'Drosophila melanogaster','bm':'Bombyx mori',
            'pl':'Platynereis dumerilii','am':'Apis mellifera','tc':'Tribolium castaneum',
            'ap':'Aplysia californica'}
# what the label claims the protein should be, as a keyword its own name must contain
EXPECT_KW={'OAMB':'ctopamine','Oct':'ctopamine','Tyr':'yramine','Dop':'opamine','5HT':'hydroxytrypt',
           'HTR':'hydroxytrypt','DRD':'opamine','ADR':'drenergic','CHRM':'uscarinic',
           'HRH':'istamine','TAAR':'race amine','MTNR':'elatonin','GPR50':'elatonin',
           'OAR':'ctopamine','HisCl':'istamine','ADRA':'drenergic'}
def get(u,tries=4):
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(u,headers={'Accept':'application/json'}),timeout=60) as r:
                return json.loads(r.read().decode())
        except Exception:
            if k==tries-1: return None
            time.sleep(2*(k+1))
def pname(r):
    d=r.get('proteinDescription',{})
    return (d.get('recommendedName',{}).get('fullName',{}).get('value')
            or (d.get('submissionNames') or [{}])[0].get('fullName',{}).get('value','?'))
rows=[]; bad=[]
print(f'{"label":<14} {"acc":<12} {"organism":<26} {"protein name":<44} verdict')
print('-'*118)
for ln in open(PROV).read().strip().split('\n')[1:]:
    p=ln.split('\t'); lab,acc=p[0],p[3]
    if not (len(acc) in (6,10) and acc[0].isalpha() and acc[1].isdigit()): continue
    r=get(f'{UP}/{acc}.json')
    if r is None: print(f'{lab:<14} {acc:<12} LOOKUP FAILED'); continue
    org=r['organism']['scientificName']; nm=pname(r)
    exp_org=''
    for pre,o in sorted(EXPECT_ORG.items(),key=lambda x:-len(x[0])):
        if lab.startswith(pre): exp_org=o; break
    kw=''
    for k,v in EXPECT_KW.items():
        if k.lower() in lab.lower(): kw=v; break
    org_ok = (not exp_org) or org.startswith(exp_org.split()[0])
    nm_ok  = (not kw) or kw.lower() in nm.lower()
    verdict='ok' if (org_ok and nm_ok) else ('WRONG ORGANISM' if not org_ok else 'NAME MISMATCH')
    if verdict!='ok': bad.append((lab,acc,org,nm,verdict))
    print(f'{lab:<14} {acc:<12} {org[:26]:<26} {nm[:44]:<44} {verdict}')
    rows.append((lab,acc,org,nm,verdict))
print()
print(f'{len(bad)} problem(s):')
for b in bad: print('   ',b[0],b[1],'->',b[2],'|',b[3],'|',b[4])
with open(f'{OUT}/accession_audit.tsv','w') as f:
    f.write('label\taccession\torganism\tprotein_name\tverdict\n')
    for r_ in rows: f.write('\t'.join(r_)+'\n')
print(f'\n-> output/accession_audit.tsv')
