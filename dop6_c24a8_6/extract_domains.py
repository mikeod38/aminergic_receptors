#!/usr/bin/env python3
"""extract_domains.py — build the tree from 7TM domain excerpts, not full-length sequences.

The full-length alignment of this panel came out 84.4 % gaps with every one of 164 sequences above
50 %, which is not an alignment so much as a scattering. The cause is structural: aminergic GPCRs
differ enormously in their N-termini, C-termini and third intracellular loops - DOP-6 alone carries
a ~280-residue ICL3 - while the seven-transmembrane core they actually share is only ~250 residues.
Aligning full length forces all that non-homologous loop into columns nobody can compare.

So each sequence is cut down to its monoamine-HMM envelope plus a margin. The margin matters: an
HMM envelope is a detection threshold, not a domain boundary, and real reference positions sit just
outside it, so a zero-margin excerpt silently truncates them. 15 residues each side, clamped to the
sequence.

The excerpt is then checked, not assumed: position 3.32 - the amine-binding aspartate, located
independently from ligand-bound structures - must survive into the excerpt for the receptors known
to have it.

Writes output/tree/panel_domains.faa and domain_report.txt
"""
import os, re, sys
import pyhmmer

HERE=os.path.dirname(os.path.abspath(__file__)); T=os.environ.get('DOP6_TREEDIR') or os.path.join(HERE,'output','tree')
MONO='/Users/mikeodonnell/git/WormPooledCloning/data/raw/MonoamineHMM_Model.hmm'
TM7='/Users/mikeodonnell/git/WormPooledCloning/data/raw/PF00001.hmm'
MARGIN=15
# The monoamine model is 224 TRIMMED positions, so a typical envelope on a real receptor runs
# ~160-220 residues, not the full 7TM span. An earlier 180-residue floor therefore discarded most
# of the human serotonin receptors - the floor has to be set against the MODEL's length, not
# against an intuition about how long a 7TM domain "should" be.
MIN_MODEL_COVERAGE=0.60
nm=lambda x: x.decode() if isinstance(x,(bytes,bytearray)) else str(x)
rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)

seqs={}
for r in open(f'{T}/panel.faa').read().split('>')[1:]:
    h,s=r.split('\n',1); seqs[h.strip()]=s.replace('\n','')
say('='*92); say('7TM DOMAIN EXCERPTS'); say('='*92)
say(f'{len(seqs)} sequences; margin {MARGIN} residues each side of the HMM envelope\n')

alpha=pyhmmer.easel.Alphabet.amino()
digs=[pyhmmer.easel.TextSequence(name=k.encode(),
        sequence=re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]','X',v)).digitize(alpha)
      for k,v in sorted(seqs.items())]
def envelopes(hmmfile):
    env={}
    with pyhmmer.plan7.HMMFile(hmmfile) as hf: hmms=list(hf)
    for hits in pyhmmer.hmmsearch(hmms,digs,E=1e-3,cpus=6):
        for h in hits:
            n=nm(h.name)
            # A 7TM receptor does not align to the model as one block: the HMM alignment breaks at
            # the long third intracellular loop, so hmmsearch reports several domains per sequence.
            # Taking only the best domain captures a fraction of the receptor (~160 aa, ~0.5 model
            # coverage) and cuts off TM6/TM7. The envelope is therefore the UNION across domains -
            # first domain start to last domain end - which spans TM1 to TM7.
            ds=[d for d in h.domains if d.i_evalue < 1.0]
            if not ds: continue
            # Keep each aligned segment SEPARATELY and concatenate them, rather than taking the
            # span from the first to the last. The span would re-include the third intracellular
            # loop sitting between the segments, which is the very region that made the
            # full-length alignment 84 % gaps - its length varies from a few residues to ~280 in
            # this family and none of it is comparable across taxa. Concatenating gives TM1-TM5
            # joined to TM6-TM7, which is what the model actually aligns.
            segs=sorted((d.alignment.target_from, d.alignment.target_to) for d in ds)
            hmin=min(d.alignment.hmm_from for d in ds)
            hmax=max(d.alignment.hmm_to for d in ds)
            L=ds[0].alignment.hmm_length
            cov=(hmax-hmin+1)/L
            tot=sum(b_-a_+1 for a_,b_ in segs)
            if n not in env or tot>env[n][3]: env[n]=(segs,h.score,cov,tot)
    return env
say('scanning with the monoamine model ...')
E=envelopes(MONO)
say(f'  {len(E)}/{len(seqs)} sequences have a monoamine envelope')
missing=[k for k in seqs if k not in E]
if missing:
    say(f'  {len(missing)} without one; falling back to the 7tm_1 fold envelope for those')
    E2=envelopes(TM7)
    for k in missing:
        if k in E2: E[k]=E2[k]
    still=[k for k in seqs if k not in E]
    if still: say(f'  {len(still)} still have no envelope and are dropped: {", ".join(still[:8])}')

out={}; lens=[]
dropped=[]
for k,s in seqs.items():
    if k not in E: continue
    segs,_,cov,_=E[k]
    iv=[]
    for a,b in segs:
        lo=max(1,a-MARGIN); hi=min(len(s),b+MARGIN)
        if iv and lo<=iv[-1][1]+1: iv[-1]=(iv[-1][0],max(iv[-1][1],hi))
        else: iv.append((lo,hi))
    ex=''.join(s[lo-1:hi] for lo,hi in iv)
    if len(ex)<150:
        dropped.append((k,len(ex),cov)); continue
    out[k]=ex; lens.append(len(ex))
if dropped:
    say(f'  {len(dropped)} excerpts under 150 aa, dropped:')
    for k,L,c in dropped[:8]: say(f'     {k:<24} {L} aa, model coverage {c:.2f}')
if dropped:
    say(f'  {len(dropped)} dropped for covering <{MIN_MODEL_COVERAGE:.0%} of the model:')
    for k,L,c in dropped[:10]: say(f'     {k:<22} excerpt {L} aa, model coverage {c:.2f}')
lens.sort()
say(f'\n{len(out)} excerpts; length min {lens[0]}, median {lens[len(lens)//2]}, max {lens[-1]}')
say(f'(full-length panel was min {min(len(v) for v in seqs.values())}, '
    f'max {max(len(v) for v in seqs.values())})')

with open(f'{T}/panel_domains.faa','w') as f:
    for k in sorted(out):
        f.write(f'>{k}\n')
        for i in range(0,len(out[k]),60): f.write(out[k][i:i+60]+'\n')
# check the amine aspartate survived where it should
gone=[]
for k in ('hADRB2','hHTR2A','hDRD2','cel_dop-6','cel_octr-1','cel_ser-3'):
    if k in out and 'D' not in out[k]: gone.append(k)
say('')
say('CHECK: the sequences known to carry the amine-binding aspartate still contain it in the')
say(f'excerpt: {"all present" if not gone else "MISSING for "+str(gone)}')
say(f'\n-> output/tree/panel_domains.faa ({len(out)} sequences)')
open(f'{T}/domain_report.txt','w').write('\n'.join(rep)+'\n')
