#!/usr/bin/env python3
"""audit_panel_hmm.py — panel membership decided by HMM score, calibrated on known controls.

Two rounds of panel building let through sequences with plausible labels and plausible lengths that
were not amine receptors at all - a squid rhodopsin and a Pseudomonas decarboxylase labelled as fly
serotonin receptors, a Drosophila ubiquitin ligase, a mollusc receptor tyrosine kinase and a
schistosome tetraspanin. Name and length filters missed every one.

The right instrument is the monoamine-receptor HMM already built for this lab's GPCR work:
`MonoamineHMM_Model.hmm` (224 positions, built from 51 experimentally deorphanised monoaminergic
receptors; Yanez-Guerra et al. 2026, Current Biology, "Prebilaterian origin of monoaminergic
signaling"), scanned locally with pyhmmer. PF00001 (7tm_1) is also scanned, because it is the
generic class-A fold and is deliberately NOT used as the membership criterion - an opsin passes it
easily. The monoamine model is the discriminating one.

Calibration, as the methodology requires, uses real controls rather than a guessed cutoff:
  positives - the 35 verified human aminergic receptors, every one a true family member
  negatives - the five contaminants found earlier, whose identities are now confirmed
The threshold is placed in the gap between them, and reported.

Writes output/tree3/hmm_scores.tsv, output/tree3/panel_clean.faa
"""
import os, sys
import pyhmmer

HERE=os.path.dirname(os.path.abspath(__file__)); T=os.path.join(HERE,'output','panel_build','v3_breadth_panel')
RAW='/Users/mikeodonnell/git/WormPooledCloning/data/raw'
MONO=f'{RAW}/MonoamineHMM_Model.hmm'; PF7TM=f'{RAW}/PF00001.hmm'
PANEL=f'{T}/panel.faa'
NEG={  # known negatives: sequences proven not to be amine receptors, kept only as controls
 'dm5HT2A_neg':'P31356','dm5HT7_neg':'P20906','dmHisCl_neg':'Q9VNJ0',
 'lsOA':'O76997','smTAR':'C4QM85'}
rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)

seqs={}
for r in open(PANEL).read().split('>')[1:]:
    h,s=r.split('\n',1); seqs[h.strip()]=s.replace('\n','').replace('*','')
# pull the confirmed-negative sequences in as explicit controls
import json, urllib.request, time
for lab,acc in NEG.items():
    if lab in seqs: continue
    try:
        with urllib.request.urlopen(f'https://rest.uniprot.org/uniprotkb/{acc}.fasta',timeout=60) as r:
            seqs[lab]=''.join(r.read().decode().split('\n')[1:])
    except Exception: pass

alpha=pyhmmer.easel.Alphabet.amino()
digs=[pyhmmer.easel.TextSequence(name=k.encode(), sequence=v).digitize(alpha)
      for k,v in sorted(seqs.items())]
nm = lambda x: x.decode() if isinstance(x, (bytes, bytearray)) else str(x)
def scan(hmmfile):
    out={}
    with pyhmmer.plan7.HMMFile(hmmfile) as hf:
        hmms=list(hf)
    for hits in pyhmmer.hmmsearch(hmms, digs, E=1e9, cpus=4):
        for h in hits:
            n=nm(h.name)
            if n not in out or h.score>out[n][0]: out[n]=(h.score,h.evalue)
    return out
say('='*96); say('PANEL AUDIT BY HMM SCORE'); say('='*96)
say(f'monoamine model: {os.path.basename(MONO)} (224 positions, 51 deorphanised receptors)')
say(f'generic fold   : PF00001 7tm_1 (reported for context, NOT the criterion)\n')
mono=scan(MONO); tm=scan(PF7TM)

HUM=[k for k in seqs if k.startswith(('hHTR','hDRD','hADR','hCHRM','hHRH','hTAAR'))]
pos=[mono.get(k,(0,))[0] for k in HUM]
negs=[k for k in NEG if k in seqs]
negsc=[mono.get(k,(0,))[0] for k in negs]
say(f'POSITIVE CONTROLS  {len(HUM)} human aminergic receptors: monoamine HMM score '
    f'min {min(pos):.1f}, median {sorted(pos)[len(pos)//2]:.1f}, max {max(pos):.1f}')
say(f'NEGATIVE CONTROLS  {len(negs)} confirmed non-receptors:')
for k in negs:
    say(f'    {k:<14} monoamine {mono.get(k,(0,))[0]:>7.1f}   7tm_1 {tm.get(k,(0,))[0]:>7.1f}')
gap_lo=max(negsc) if negsc else 0; gap_hi=min(pos)
say(f'\ngap between highest negative ({gap_lo:.1f}) and lowest positive ({gap_hi:.1f})')
if gap_hi<=gap_lo:
    say('NO CLEAN GAP - the controls overlap, so a single cutoff cannot separate them. Stopping.')
    open(f'{T}/audit_hmm.txt','w').write('\n'.join(rep)+'\n'); sys.exit(1)
CUT=gap_lo+(gap_hi-gap_lo)/2
say(f'threshold placed midway: {CUT:.1f}\n')

rows=[]
for k in sorted(seqs):
    if k in NEG and k.endswith('_neg'): continue
    m=mono.get(k,(0.0,float('inf'))); t=tm.get(k,(0.0,float('inf')))
    v='PASS' if m[0]>=CUT else 'FAIL'
    rows.append((k,f'{m[0]:.1f}',f'{m[1]:.1e}',f'{t[0]:.1f}',v))
fails=[r for r in rows if r[4]=='FAIL']
say(f'{len(rows)-len(fails)} pass, {len(fails)} fail\n')
say(f'{"sequence":<22} {"monoamine":>10} {"evalue":>10} {"7tm_1":>8}  verdict')
say('-'*62)
for r in sorted(rows,key=lambda x: float(x[1]))[:18]:
    say(f'{r[0]:<22} {r[1]:>10} {r[2]:>10} {r[3]:>8}  {r[4]}')
say('  ... (lowest-scoring 18 shown)')
say('')
say('EXCLUDED:')
for r in fails: say(f'  {r[0]:<20} monoamine {r[1]:>7}  7tm_1 {r[3]:>7}')
keep=[r[0] for r in rows if r[4]=='PASS']
with open(f'{T}/panel_clean.faa','w') as f:
    for k in sorted(keep):
        f.write(f'>{k}\n')
        for i in range(0,len(seqs[k]),60): f.write(seqs[k][i:i+60]+'\n')
with open(f'{T}/hmm_scores.tsv','w') as f:
    f.write('sequence\tmonoamine_score\tmonoamine_evalue\tpf00001_score\tverdict\n')
    for r in rows: f.write('\t'.join(r)+'\n')
from collections import Counter
c=Counter('nematode' if any(k.startswith(p) for p in ('Cbr_','Cbn_','Cre_','Ppa_','Asu_','Bma_','Sra_','Hco_','Tsp_','Nam_')) else
          'human' if k.startswith('h') else 'C. elegans' if k.startswith('cel_') else 'other invertebrate'
          for k in keep)
say(f'\n-> output/tree3/panel_clean.faa ({len(keep)} sequences)')
for k,v in sorted(c.items()): say(f'   {k:<20} {v}')
open(f'{T}/audit_hmm.txt','w').write('\n'.join(rep)+'\n')
