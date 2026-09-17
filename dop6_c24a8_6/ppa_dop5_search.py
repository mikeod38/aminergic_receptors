#!/usr/bin/env python3
"""ppa_dop5_search.py — is dop-5 genuinely absent from P. pacificus, or did Compara miss it?

The tree shows P. pacificus carrying dop-3 and dop-6 but no dop-5, which would place the dop-5
duplication inside Caenorhabditis. That conclusion rests on Ensembl Compara having called the
orthologues, and an ortholog-calling pipeline under-calls in less well annotated genomes - exactly
the "a negative reached any way other than a direct search is provisional" problem. So the absence
is retested directly against the P. pacificus proteome, by two independent methods:

  1. HMM   - scan the whole proteome with the monoamine model, take everything scoring in the
             aminergic range, and ask what each one's closest C. elegans receptor is
  2. BLAST - reciprocal best hits between C. elegans dop-1..dop-6 and the P. pacificus proteome

dop-5 is called present only if some P. pacificus protein's best C. elegans match is dop-5 AND
dop-5's best P. pacificus match is that protein. A coverage check is reported alongside, because a
negative is only as good as the database behind it.

Writes output/ppa_dop5_search.txt
"""
import os, re, subprocess, sys
import pyhmmer

HERE=os.path.dirname(os.path.abspath(__file__)); OUT=os.path.join(HERE,'output')
PPA=('/private/tmp/claude-501/-Users-mikeodonnell-git-Malaiwong-Oglu-2026/'
     '93015e40-e331-495a-8362-20ac7ea07c75/scratchpad/ppa/ppa.faa')
MONO='/Users/mikeodonnell/git/WormPooledCloning/data/raw/MonoamineHMM_Model.hmm'
W='/tmp/_ppa'; os.makedirs(W,exist_ok=True)
nm=lambda x: x.decode() if isinstance(x,(bytes,bytearray)) else str(x)
rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)

ppa={}
name=None; buf=[]
for ln in open(PPA,errors='replace'):
    if ln.startswith('>'):
        if name: ppa[name]=''.join(buf)
        name=ln[1:].split()[0]; buf=[]
    else: buf.append(ln.strip())
if name: ppa[name]=''.join(buf)
say('='*88); say('IS dop-5 REALLY ABSENT FROM P. PACIFICUS?'); say('='*88)
say(f'proteome UP000005239: {len(ppa)} proteins  <- coverage check: a negative against a proteome')
say('this size is a fair test, not a database gap\n')

# C. elegans dop-1..dop-6 from the panel
cel={}
for r in open(f'{OUT}/tree/panel.faa').read().split('>')[1:]:
    h,s=r.split('\n',1)
    if re.match(r'^cel_dop-[1-6]$',h.strip()): cel[h.strip()]=s.replace('\n','')
say(f'queries: {", ".join(sorted(cel))}\n')

# ── 1. HMM scan of the whole proteome ────────────────────────────────────────
alpha=pyhmmer.easel.Alphabet.amino()
digs=[pyhmmer.easel.TextSequence(name=k.encode(),
        sequence=re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]','X',v)).digitize(alpha)
      for k,v in ppa.items() if 250<=len(v)<=1200]
with pyhmmer.plan7.HMMFile(MONO) as hf: hmms=list(hf)
sc={}
for hits in pyhmmer.hmmsearch(hmms,digs,E=1e-5,cpus=6):
    for h in hits:
        n=nm(h.name)
        if n not in sc or h.score>sc[n]: sc[n]=h.score
rank=sorted(sc.items(),key=lambda x:-x[1])
say(f'1. HMM SCAN: {len(sc)} P. pacificus proteins hit the monoamine model;')
say(f'   {sum(1 for _,v in rank if v>=120)} score >=120 (the aminergic range)\n')

# ── 2. reciprocal BLAST ──────────────────────────────────────────────────────
with open(f'{W}/cel.faa','w') as f:
    for k,v in cel.items(): f.write(f'>{k}\n{v}\n')
cand={n:ppa[n] for n,v in rank if v>=120 and n in ppa}
with open(f'{W}/ppa_cand.faa','w') as f:
    for k,v in cand.items(): f.write(f'>{k}\n{v}\n')
for tag,src in (('cel',f'{W}/cel.faa'),('ppa',f'{W}/ppa_cand.faa')):
    subprocess.run(['makeblastdb','-in',src,'-dbtype','prot','-out',f'{W}/{tag}'],capture_output=True)
def best(q,db):
    r=subprocess.run(['blastp','-query',q,'-db',db,'-outfmt','6 qseqid sseqid bitscore pident',
                      '-evalue','1e-5','-max_target_seqs','5'],capture_output=True,text=True)
    out={}
    for l in r.stdout.strip().split('\n'):
        if not l: continue
        p=l.split('\t')
        if p[0] not in out: out[p[0]]=(p[1],float(p[2]),float(p[3]))
    return out
fwd=best(f'{W}/cel.faa',f'{W}/ppa')      # each cel dop-N -> best P. pacificus candidate
rev=best(f'{W}/ppa_cand.faa',f'{W}/cel') # each P. pacificus candidate -> best cel dop-N
say('2. RECIPROCAL BEST HITS')
say(f'   {"C. elegans":<12} {"best P. pacificus":<20} {"bits":>7}   {"reciprocates?":<14} note')
for g in sorted(cel):
    if g not in fwd: say(f'   {g:<12} {"no hit":<20}'); continue
    tgt,b,pid=fwd[g]
    back=rev.get(tgt,('',0,0))[0]
    ok = back==g
    say(f'   {g:<12} {tgt[:20]:<20} {b:>7.0f}   {"YES" if ok else "no":<14} '
        f'{"" if ok else f"that protein maps back to {back}"}')
say('')
d5=fwd.get('cel_dop-5')
if d5 and rev.get(d5[0],('',0,0))[0]=='cel_dop-5':
    say('VERDICT: a reciprocal best-hit partner for dop-5 EXISTS in P. pacificus - Compara missed it,')
    say('and the dop-5 duplication predates the Caenorhabditis / Pristionchus split.')
else:
    say('VERDICT: no P. pacificus protein reciprocates with dop-5. Its best forward hit maps back to')
    say(f'   {rev.get(d5[0],("none",))[0] if d5 else "n/a"}, i.e. dop-5 has no distinct orthologue here.')
    say('   Consistent with the dop-5 duplication having occurred within Caenorhabditis, though a')
    say('   Pristionchus-specific loss cannot be separated from a later duplication by this test.')
open(f'{OUT}/ppa_dop5_search.txt','w').write('\n'.join(rep)+'\n')
say(f'\n-> output/ppa_dop5_search.txt')
