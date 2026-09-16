#!/usr/bin/env python3
"""curate_tree_panel.py — verify and repair the adrenergic-family panel before any tree is built.

fetch_tree_seqs.py collected sequences by UniProt query and labelled their ligand class from the
gene-name prefix. Both steps have known failure modes and both produced errors here:

  * gene_exact:ADRA1A returned P25100, whose own protein name is "Alpha-1D adrenergic receptor"
  * TyrRII and TyrRIII queries both returned the same accession Q9VEG2
  * Octbeta3R returned a 1256 aa entry, far outside the family's length range
  * the C. elegans "ser-*" genes were labelled serotonin by prefix, but SER-2 is a tyramine
    receptor and SER-3 and SER-6 are octopamine receptors by direct assay

So this script (a) re-fetches each entry and checks its OWN protein name against the ligand class
claimed for it, (b) replaces the broken queries, (c) applies ligand classes from the curated
pharmacology table rather than from gene names, with the evidence tier attached, and (d) applies
the length and fragment rules: full-length only, within the family's own length range.

Writes output/tree/panel_curated.faa, output/tree/provenance_curated.tsv, output/tree/curation_log.txt
"""
import json, os, time, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'output', 'tree')
UP = 'https://rest.uniprot.org/uniprotkb'

# ligand class from the lab's curated evidence file (WormPooledCloning cel_monoamine_evidence.md),
# not from the gene name. tier: IDA = direct assay on this gene; IBA/none = inference only.
CEL_LIGAND = {
 'dop-1': ('dopamine', 'IDA', 'GTPgS + [125I]iodo-LSD binding, PMID 14739932 / 11814642'),
 'dop-2': ('dopamine', 'IDA', '[125I]iodo-LSD binding + FSK cAMP, PMID 12887685'),
 'dop-3': ('dopamine', 'IDA', 'inhibition of FSK-induced cAMP, PMID 16001968'),
 'dop-4': ('dopamine', 'IDA', 'cAMP activation in 293 cells, PMID 16001968'),
 'dop-5': ('ORPHAN', 'none', 'no pharmacology; named by homology only, PMID 15378064'),
 'dop-6': ('ORPHAN', 'none', 'no pharmacology; named by homology only, PMID 15378064'),
 'ser-1': ('serotonin', 'IDA', 'PMID 10098838'),
 'ser-2': ('tyramine', 'IDA', 'PMID 12354282 — despite the ser- name'),
 'ser-3': ('octopamine', 'IDA', 'PMID 18033297'),
 'ser-4': ('serotonin', 'IDA', 'PMID 18033297'),
 'ser-5': ('uncertain', 'IBA', 'ligand and coupling IBA only'),
 'ser-6': ('octopamine', 'IDA', 'PMID 22124329'),
 'ser-7': ('serotonin', 'IBA', 'no direct evidence in WormBase'),
 'octr-1': ('octopamine', 'IDA', 'PMID 22124329; strongest OA evidence in the genome'),
 'tyra-2': ('tyramine', 'IBA', 'WormBase "amine" IBA only'),
 'tyra-3': ('tyramine/OA', 'IBA', 'WormBase annotates OA IBA; ligand conflicting'),
 'gar-3': ('acetylcholine', 'IDA', 'muscarinic — outgroup, not a monoamine receptor'),
}
# replacements / additions, as (label, expected ligand class, keyword that MUST appear in the
# entry's own protein name, uniprot query)
FIX = [
 ('hADRA1A', 'adrenergic_a1', 'Alpha-1A', 'gene_exact:ADRA1A AND organism_id:9606 AND reviewed:true AND protein_name:"Alpha-1A adrenergic receptor"'),
 ('hADRA1D', 'adrenergic_a1', 'Alpha-1D', 'accession:P25100'),
 ('dmTAR2',  'tyramine',      'Tyramine', 'gene_exact:Tyr-II AND organism_id:7227'),
 ('dmTAR3',  'tyramine',      'Tyramine', 'gene_exact:TyrR-III AND organism_id:7227'),
 ('dmOctB3R','octopamine_b',  'Octopamine', 'accession:Q4LBB6'),
 # extra breadth for the octopamine / tyramine clade
 ('amOA1',   'octopamine_a',  'Octopamine', 'gene:Oa1 AND organism_name:"Apis mellifera" AND reviewed:true'),
 ('bmOAR1',  'octopamine_a',  'Octopamine', 'organism_name:"Bombyx mori" AND protein_name:"Octopamine receptor"'),
 ('plADRA1', 'adrenergic_inv','adrenergic','organism_name:"Platynereis dumerilii" AND protein_name:"adrenergic receptor"'),
]
def get(u, tries=5):
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(u, headers={'Accept': 'application/json'}), timeout=60) as r:
                return json.loads(r.read().decode())
        except Exception:
            if k == tries-1: raise
            time.sleep(2*(k+1))
def pname(r):
    d = r.get('proteinDescription', {})
    return (d.get('recommendedName', {}).get('fullName', {}).get('value')
            or (d.get('submissionNames') or [{}])[0].get('fullName', {}).get('value', '?'))

log = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); log.append(s)

seqs, prov = {}, {}
for r in open(f'{OUT}/panel.faa').read().split('>')[1:]:
    h, s = r.split('\n', 1); seqs[h.strip().split('|')[0]] = (h.strip(), s.replace('\n', ''))
for ln in open(f'{OUT}/provenance.tsv').read().strip().split('\n')[1:]:
    p = ln.split('\t'); prov[p[0]] = p

say('=' * 100); say('PANEL CURATION'); say('=' * 100)

# ── 1. verify every UniProt accession against its own protein name ───────────
say('\n1. VERIFYING EACH UniProt ENTRY AGAINST ITS OWN PROTEIN NAME')
say('-' * 100)
EXPECT = {'adrenergic': 'adrenergic', 'octopamine': 'ctopamine', 'tyramine': 'yramine',
          'dopamine': 'opamine', 'serotonin': 'hydroxytryptamine', 'muscarinic': 'uscarinic'}
bad = []
for lab, p in sorted(prov.items()):
    acc = p[3]
    if not (len(acc) in (6, 10) and acc[0].isalpha() and acc[1].isdigit()): continue
    try: r = get(f'{UP}/{acc}.json')
    except Exception as e: say(f'  {lab:<10} {acc:<12} LOOKUP FAILED'); continue
    nm = pname(r); org = r['organism']['scientificName']
    grp = p[2]; key = grp.split('_')[0]
    ok = EXPECT.get(key, '') .lower() in nm.lower() if key in EXPECT else True
    flag = '' if ok else '   <-- NAME DOES NOT MATCH THE CLAIMED CLASS'
    if not ok: bad.append(lab)
    say(f'  {lab:<10} {acc:<12} {org[:24]:<24} {nm[:46]:<46}{flag}')
say(f'\n  {len(bad)} mismatch(es): {", ".join(bad) if bad else "none"}')

# ── 2. apply the fixes ───────────────────────────────────────────────────────
say('\n2. REPAIRING BROKEN / DUPLICATED QUERIES AND ADDING BREADTH')
say('-' * 100)
for lab, grp, kw, q in FIX:
    url = (f'{UP}/search?query={urllib.parse.quote(q)}&format=json&size=5'
           f'&fields=accession,protein_name,organism_name,length,sequence')
    try: d = get(url)
    except Exception as e: say(f'  {lab:<10} QUERY FAILED'); continue
    res = d.get('results', [])
    hit = next((r for r in res if kw.lower() in pname(r).lower()), None)
    if not hit:
        say(f'  {lab:<10} no entry whose own name contains "{kw}"  (query returned {len(res)})'); continue
    acc = hit['primaryAccession']; sq = hit['sequence']['value']
    seqs[lab] = (f'{lab}|{acc}', sq)
    prov[lab] = [lab, hit['organism']['scientificName'], grp, acc, 'uniprot:' + q, str(len(sq)),
                 f'verified: {pname(hit)}']
    say(f'  {lab:<10} {acc:<12} {len(sq):>4} aa  {pname(hit)[:52]}')

# ── 3. ligand classes from pharmacology, not gene names ──────────────────────
say('\n3. RE-LABELLING C. elegans RECEPTORS BY CURATED PHARMACOLOGY')
say('-' * 100)
for g, (lig, tier, note) in CEL_LIGAND.items():
    lab = f'cel_{g}'
    if lab not in prov: continue
    old = prov[lab][2]
    prov[lab][2] = f'{lig}({tier})'
    prov[lab][6] = note
    mark = '  CHANGED' if old.split('_')[0] not in lig else ''
    say(f'  {lab:<12} {old:<16} -> {lig} [{tier}]{mark}')

# ── 4. length / fragment rules ───────────────────────────────────────────────
say('\n4. LENGTH AND FRAGMENT RULES')
say('-' * 100)
L = sorted(len(s) for _, s in seqs.values())
import statistics
med = statistics.median(L)
say(f'  panel length range {L[0]}-{L[-1]}, median {med:.0f}')
DROP = {}
for lab, (h, s) in list(seqs.items()):
    if lab == 'cel_C24A8.6b':
        DROP[lab] = ('fragment by construction — it IS DOP-6 406-692; kept out of the tree so it '
                     'cannot destabilise it, and reported separately as an exact substring')
    elif len(s) < 0.55 * med:
        DROP[lab] = f'fragment: {len(s)} aa, below 55 % of the panel median ({med:.0f})'
for lab, why in DROP.items():
    say(f'  DROP {lab:<14} {why}')
    seqs.pop(lab, None); prov.pop(lab, None)
say(f'  retained {len(seqs)} sequences')

with open(f'{OUT}/panel_curated.faa', 'w') as f:
    for lab in sorted(seqs):
        h, s = seqs[lab]; f.write(f'>{h}\n')
        for i in range(0, len(s), 60): f.write(s[i:i+60] + '\n')
with open(f'{OUT}/provenance_curated.tsv', 'w') as f:
    f.write('id\torganism\tligand_class\taccession\tsource_query\tlength\tnote\n')
    for lab in sorted(prov): f.write('\t'.join(prov[lab]) + '\n')
open(f'{OUT}/curation_log.txt', 'w').write('\n'.join(log) + '\n')
say(f'\n-> output/tree/panel_curated.faa ({len(seqs)} seqs), provenance_curated.tsv, curation_log.txt')
