#!/usr/bin/env python3
"""nematode_placement.py — date the dop-3/dop-5/dop-6 radiation by PLACEMENT, not best hits.

Replaces nematode_dating2.py, whose reciprocal-best-hit matrix was retracted: RBH cannot resolve
members of a recent paralogue radiation, because the reverse hit lands on whichever paralogue
scores highest and a gene that is genuinely present is recorded as absent.

The method here makes no orthology call at all. It scans every available nematode proteome with
the monoamine HMM, takes EVERY candidate at or above the panel's own boundary, and lets the tree
assign paralogy. A clade is then dated by which nematode lineages have a member inside it.

Sampling spans the phylum, which the previous 5-species panel did not:
  Dorylaimia   (clade I,   earliest-diverging)  Trichinella spiralis, Trichuris muris
  Spirurina    (clade III)                      Brugia malayi
  Tylenchina   (clade IV)                       Strongyloides ratti
  Rhabditina   (clade V)                        Haemonchus contortus, Necator americanus,
                                                Pristionchus pacificus, Caenorhabditis briggsae

Writes output/tree/placement/ (candidates, panel, report).
"""
import os, re, sys, subprocess, csv
from collections import defaultdict
import pyhmmer

HERE=os.path.dirname(os.path.abspath(__file__))
STORE=os.path.expanduser('~/data/proteomes')
HMM='/Users/mikeodonnell/git/WormPooledCloning/data/raw/MonoamineHMM_Model.hmm'
OUT=f'{HERE}/output/tree/placement'; os.makedirs(OUT,exist_ok=True)
CUT=120.0

LINEAGE={'Tsp':('Trichinella spiralis','I   Dorylaimia (earliest-diverging)'),
         'Tmu':('Trichuris muris','I   Dorylaimia'),
         'Bma':('Brugia malayi','III Spirurina'),
         'Sra':('Strongyloides ratti','IV  Tylenchina'),
         'Hco':('Haemonchus contortus','V   Rhabditina (Strongylida)'),
         'Nam':('Necator americanus','V   Rhabditina (Strongylida)'),
         'Ppa':('Pristionchus pacificus','V   Rhabditina (Diplogasteromorpha)'),
         'Cbr':('Caenorhabditis briggsae','V   Rhabditina (Caenorhabditis)')}

rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)

man={r['tag']:r for r in csv.DictReader(open(f'{STORE}/MANIFEST.tsv'),delimiter='\t')}
say('='*100); say('DATING THE dop-3/dop-5/dop-6 RADIATION BY PHYLOGENETIC PLACEMENT'); say('='*100)
say(f'monoamine HMM, inclusion at score >= {CUT} (the b120 panel boundary), no orthology call made.\n')

with pyhmmer.plan7.HMMFile(HMM) as fh: hmm=next(fh)
alpha=hmm.alphabet
cands=[]
for tag in LINEAGE:
    row=man.get(tag)
    if not row: say(f'  {tag}: NOT IN STORE - skipped'); continue
    if row['complete']!='yes':
        say(f'  {tag}: manifest says INCOMPLETE ({row["n_proteins"]}/{row["expected_n"]}) - REFUSED'); continue
    path=row['file']
    with pyhmmer.easel.SequenceFile(path,digital=True,alphabet=alpha) as sf:
        seqs=list(sf)
    hits=list(pyhmmer.hmmsearch([hmm],seqs,cpus=6))[0]
    keep=[h for h in hits if h.score>=CUT]
    sp,cl=LINEAGE[tag]
    say(f'  {tag}  {sp:<26} {len(seqs):>6} proteins -> {len(keep):>3} candidates >= {CUT}')
    def _n(x): return x.decode() if isinstance(x,bytes) else x
    byname={_n(s.name):s for s in seqs}
    for h in keep:
        nm=_n(h.name)
        acc=nm.split('|')[1] if '|' in nm else nm
        cands.append({'tag':tag,'species':sp,'clade':cl,'acc':acc,'score':h.score,
                      'seq':byname[nm].textize().sequence})

say(f'\n{len(cands)} candidates total across {len(set(c["tag"] for c in cands))} proteomes')
say('\nper lineage:')
byc=defaultdict(int)
for c in cands: byc[c['clade']]+=1
for k in sorted(byc): say(f'   {k:<42} {byc[k]}')

with open(f'{OUT}/candidates.faa','w') as f:
    for c in cands:
        f.write(f'>{c["tag"]}_{c["acc"]}\n{c["seq"]}\n')
with open(f'{OUT}/candidates.tsv','w') as f:
    f.write('tip\ttag\tspecies\tclade\taccession\thmm_score\tlength\n')
    for c in cands:
        f.write(f'{c["tag"]}_{c["acc"]}\t{c["tag"]}\t{c["species"]}\t{c["clade"]}\t{c["acc"]}\t{c["score"]:.1f}\t{len(c["seq"])}\n')
say(f'\n-> {OUT}/candidates.faa  and  candidates.tsv')
open(f'{OUT}/scan.txt','w').write('\n'.join(rep)+'\n')
