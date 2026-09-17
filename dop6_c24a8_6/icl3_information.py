#!/usr/bin/env python3
"""icl3_information.py — is the discarded ICL3 carrying phylogenetic signal we are throwing away?

The 7TM excerpting removes the third intracellular loop because its length varies from ~15 to ~820
residues across this panel and it cannot be aligned across the family. But "cannot be aligned across
the family" is not the same as "carries no information": an insertion can be a perfectly good
lineage marker at shallow depth even when it is unalignable at depth. Worth testing rather than
assuming, and worth testing in two separate ways, because the two possible kinds of signal fail
differently.

  TEST 1 - is ICL3 LENGTH phylogenetically clustered? Length is a trait like any other. Compared
           against the ML tree: the mean absolute length difference between sister tips, versus a
           null from shuffling lengths across tips. If close relatives have similar ICL3 lengths,
           length is conserved and carries signal; if not, it turns over too fast to be a marker.

  TEST 2 - is ICL3 SEQUENCE alignable at the depths that matter? Percent identity of ICL3 between
           orthologues (same gene, different nematode species), between paralogues within
           C. elegans, and across the family - compared with the same three numbers for the 7TM
           core from the same sequence pairs. If ICL3 identity collapses to background while the
           core stays alignable, the loop is unusable outside very shallow comparisons.

Writes output/icl3_information.txt
"""
import os, re, random
from collections import defaultdict
import pyhmmer
from Bio.Align import PairwiseAligner, substitution_matrices

HERE=os.path.dirname(os.path.abspath(__file__)); T=os.path.join(HERE,'output','tree')
MONO='/Users/mikeodonnell/git/WormPooledCloning/data/raw/MonoamineHMM_Model.hmm'
random.seed(0)
nm=lambda x: x.decode() if isinstance(x,(bytes,bytearray)) else str(x)
rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)

seqs={}
for r in open(f'{T}/panel.faa').read().split('>')[1:]:
    h,s=r.split('\n',1); seqs[h.strip()]=s.replace('\n','')
alpha=pyhmmer.easel.Alphabet.amino()
digs=[pyhmmer.easel.TextSequence(name=k.encode(),
        sequence=re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]','X',v)).digitize(alpha)
      for k,v in sorted(seqs.items())]
with pyhmmer.plan7.HMMFile(MONO) as hf: hmms=list(hf)
blocks={}
for hits in pyhmmer.hmmsearch(hmms,digs,E=1e-3,cpus=6):
    for h in hits:
        ds=sorted([d for d in h.domains if d.i_evalue<1.0],
                  key=lambda d:d.alignment.target_from)
        if len(ds)>=2: blocks[nm(h.name)]=[(d.alignment.target_from,d.alignment.target_to) for d in ds]
# ICL3 = the largest skipped stretch between consecutive blocks
icl3={}; icl3seq={}
for k,bs in blocks.items():
    gaps=[(bs[i+1][0]-bs[i][1]-1, bs[i][1], bs[i+1][0]) for i in range(len(bs)-1)]
    g,a,b=max(gaps)
    if g<=0: continue
    icl3[k]=g; icl3seq[k]=seqs[k][a:b-1]
say('='*92); say('IS THE DISCARDED ICL3 CARRYING SIGNAL?'); say('='*92)
L=sorted(icl3.values())
say(f'{len(icl3)} sequences have a measurable ICL3')
say(f'length: min {L[0]}, median {L[len(L)//2]}, max {L[-1]}, '
    f'{sum(1 for x in L if x>100)} over 100 aa\n')

# ── TEST 1: phylogenetic signal in ICL3 length ───────────────────────────────
tf=f'{T}/trees/gt0.1.treefile'
if os.path.exists(tf):
    def parse(nw):
        nw=nw.strip().rstrip(';'); pos=[0]
        def node():
            if nw[pos[0]]=='(':
                pos[0]+=1; kids=[]
                while True:
                    kids.append(node())
                    if nw[pos[0]]==',': pos[0]+=1; continue
                    if nw[pos[0]]==')': pos[0]+=1; break
                m=re.match(r'[^,\)\(:]*',nw[pos[0]:]); pos[0]+=len(m.group(0))
                if pos[0]<len(nw) and nw[pos[0]]==':':
                    m2=re.match(r':[-\d.eE+]*',nw[pos[0]:]); pos[0]+=len(m2.group(0))
                return ('int',kids)
            m=re.match(r'[^,\)\(:]*',nw[pos[0]:]); lab=m.group(0); pos[0]+=len(lab)
            if pos[0]<len(nw) and nw[pos[0]]==':':
                m2=re.match(r':[-\d.eE+]*',nw[pos[0]:]); pos[0]+=len(m2.group(0))
            return ('tip',lab)
        return node()
    t=parse(open(tf).read())
    pairs=[]
    def walk(n):
        if n[0]=='tip': return [n[1]]
        tl=[]
        kids=[walk(k) for k in n[1]]
        if len(kids)==2 and len(kids[0])==1 and len(kids[1])==1:
            pairs.append((kids[0][0],kids[1][0]))
        for k in kids: tl+=k
        return tl
    walk(t)
    use=[(a,b) for a,b in pairs if a in icl3 and b in icl3]
    if use:
        obs=sum(abs(icl3[a]-icl3[b]) for a,b in use)/len(use)
        keys=[k for k in icl3]
        null=[]
        for _ in range(10000):
            vals=[icl3[k] for k in keys]; random.shuffle(vals)
            sh=dict(zip(keys,vals))
            null.append(sum(abs(sh[a]-sh[b]) for a,b in use)/len(use))
        p=(sum(1 for x in null if x<=obs)+1)/10001
        say('TEST 1 - phylogenetic signal in ICL3 LENGTH')
        say(f'  {len(use)} sister pairs on the ML tree')
        say(f'  mean |ICL3 length difference| between sisters: {obs:.1f} aa')
        say(f'  same statistic with lengths shuffled across tips: {sum(null)/len(null):.1f} aa')
        say(f'  p = {p:.4f}  -> {"length IS conserved between close relatives" if p<0.05 else "no detectable signal"}')
else:
    say('TEST 1 skipped: no tree yet')

# ── TEST 2: ICL3 identity vs core identity, at three depths ──────────────────
al=PairwiseAligner(); al.substitution_matrix=substitution_matrices.load('BLOSUM62')
al.open_gap_score,al.extend_gap_score,al.mode=-11,-1,'global'
core={}
for r in open(f'{T}/panel_domains.faa').read().split('>')[1:]:
    h,s=r.split('\n',1); core[h.strip()]=s.replace('\n','')
def pid(a,b):
    if not a or not b: return None
    x=al.align(a,b)[0]; A,B=str(x[0]),str(x[1])
    both=[(p,q) for p,q in zip(A,B) if p!='-' and q!='-']
    return 100*sum(1 for p,q in both if p==q)/len(both) if both else None
NEM=('Cbr_','Cbn_','Cre_','Ppa_','Asu_','Bma_','Sra_','Hco_','Tsp_','Nam_')
orth=[]; para=[]; far=[]
cel=[k for k in icl3 if k.startswith('cel_')]
for k in icl3:
    if not k.startswith(NEM): continue
    g=k.split('_',1)[1]; partner=f'cel_{g}'
    if partner in icl3: orth.append((k,partner))
for i in range(len(cel)):
    for j in range(i+1,len(cel)): para.append((cel[i],cel[j]))
hum=[k for k in icl3 if k.startswith('h')]
for c in cel[:12]:
    for hh in hum[:12]: far.append((c,hh))
say('')
say('TEST 2 - ICL3 sequence identity vs 7TM core identity, at three depths')
say(f'  {"comparison":<34} {"n":>5} {"ICL3 %id":>10} {"core %id":>10}')
for lab,pl in (('orthologues (same gene, other nematode)',orth[:60]),
               ('paralogues within C. elegans',random.sample(para,min(60,len(para)))),
               ('C. elegans vs human',random.sample(far,min(60,len(far))))):
    iv=[]; cv=[]
    for a,b in pl:
        x=pid(icl3seq.get(a),icl3seq.get(b));  y=pid(core.get(a),core.get(b))
        if x is not None: iv.append(x)
        if y is not None: cv.append(y)
    if iv: say(f'  {lab:<34} {len(iv):>5} {sum(iv)/len(iv):>9.1f}% {sum(cv)/len(cv) if cv else 0:>9.1f}%')
say('')
say('Reading: if ICL3 identity is near the core at shallow depth but collapses toward background')
say('(~10-15 % for unrelated sequences of this composition) between paralogues and across phyla,')
say('then the loop is a usable marker only within orthologue groups - far shallower than the nodes')
say('this tree has to resolve - and excising it costs the deep topology nothing.')
open(f'{HERE}/output/icl3_information.txt','w').write('\n'.join(rep)+'\n')
say('\n-> output/icl3_information.txt')
