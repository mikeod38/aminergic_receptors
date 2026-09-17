#!/usr/bin/env python3
"""nematode_dating.py — bracket the dop-3 / dop-5 / dop-6 duplications across nematode phylogeny.

The three genes form an exclusively nematode clade (UFboot 91 with the rogue removed), sister to the
dop-2 lineage which is shared with Drosophila. Which nematode clades carry each of the three genes
brackets when each duplication happened. Trichinella is the decisive taxon: it sits in Dorylaimia,
the earliest-diverging nematode clade, so a dop-3/5/6 member there would place the radiation before
the nematode crown, and its absence would confine it to Chromadorea.

Compara cannot answer this - it missed dop-5 in P. pacificus entirely - so each proteome is searched
directly, by HMM scan followed by reciprocal best hits against C. elegans dop-1..dop-6. A gene is
called present only when the reciprocity holds in both directions, and the proteome size is reported
with every call so a negative can be weighed against the database behind it.

Writes output/nematode_dating.tsv and .txt
"""
import gzip, json, os, re, subprocess, sys, time, urllib.parse, urllib.request
import pyhmmer

HERE=os.path.dirname(os.path.abspath(__file__)); OUT=os.path.join(HERE,'output')
W=('/private/tmp/claude-501/-Users-mikeodonnell-git-Malaiwong-Oglu-2026/'
   '93015e40-e331-495a-8362-20ac7ea07c75/scratchpad/nem')
os.makedirs(W,exist_ok=True)
MONO='/Users/mikeodonnell/git/WormPooledCloning/data/raw/MonoamineHMM_Model.hmm'
# clade assignment follows the standard nematode phylogeny (Blaxter clades I-V)
SPECIES=[('Tspi','UP000054776','Trichinella spiralis','I  Dorylaimia (earliest-diverging)'),
         ('Tmur','UP000046395','Trichuris muris',      'I  Dorylaimia'),
         ('Bmal','UP000006672','Brugia malayi',        'III Spirurina'),
         ('Srat','UP000035682','Strongyloides ratti',  'IV Tylenchina'),
         ('Hcon','UP000025227','Haemonchus contortus', 'V  Rhabditina (Strongylida)'),
         ('Name','UP000053676','Necator americanus',   'V  Rhabditina (Strongylida)')]
nm=lambda x: x.decode() if isinstance(x,(bytes,bytearray)) else str(x)
rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)

cel={}
for r in open(f'{OUT}/tree/panel.faa').read().split('>')[1:]:
    h,s=r.split('\n',1)
    if re.match(r'^cel_dop-[1-6]$',h.strip()): cel[h.strip()]=s.replace('\n','')
with open(f'{W}/cel.faa','w') as f:
    for k,v in cel.items(): f.write(f'>{k}\n{v}\n')
subprocess.run(['makeblastdb','-in',f'{W}/cel.faa','-dbtype','prot','-out',f'{W}/cel'],capture_output=True)
alpha=pyhmmer.easel.Alphabet.amino()
with pyhmmer.plan7.HMMFile(MONO) as hf: HMMS=list(hf)

say('='*100); say('WHICH NEMATODE CLADES CARRY dop-3 / dop-5 / dop-6?'); say('='*100)
say('presence called only on reciprocal best hits in both directions; proteome size reported so a')
say('negative can be weighed against its database\n')
rows=[]
for tag,upid,sci,clade in SPECIES:
    fa=f'{W}/{tag}.faa'
    if not os.path.exists(fa) or os.path.getsize(fa)<500_000:
        # UniProt's stream endpoint truncates on large queries - urlretrieve reports success on a
        # partial file and IncompleteRead on others. curl with --retry and an explicit size check is
        # the only reliable form; a truncated proteome would turn a real gene into a false absence.
        url=(f'https://rest.uniprot.org/uniprotkb/stream?query={urllib.parse.quote("proteome:"+upid)}'
             f'&format=fasta')
        ok=False
        for attempt in range(4):
            say(f'{tag}: downloading {upid} (attempt {attempt+1}) ...')
            r=subprocess.run(['curl','-s','--retry','3','--retry-delay','5','--max-time','900',
                              '-o',fa,url],capture_output=True)
            n=sum(1 for ln in open(fa,errors='replace') if ln.startswith('>')) if os.path.exists(fa) else 0
            if r.returncode==0 and n>3000: ok=True; say(f'   got {n} proteins'); break
            say(f'   incomplete ({n} proteins, curl rc={r.returncode}), retrying')
            time.sleep(8)
        if not ok:
            say(f'{tag}: DOWNLOAD FAILED after 4 attempts - no call made for this species'); continue
    seqs={}; name=None; buf=[]
    for ln in open(fa,errors='replace'):
        if ln.startswith('>'):
            if name: seqs[name]=''.join(buf)
            name=ln[1:].split()[0]; buf=[]
        else: buf.append(ln.strip())
    if name: seqs[name]=''.join(buf)
    say(f'{tag} ({sci}): {len(seqs)} proteins')
    if len(seqs)<3000:
        say(f'   PROTEOME LOOKS TRUNCATED - not calling absence from this'); rows.append((tag,sci,clade,len(seqs),{},'TRUNCATED')); continue
    digs=[pyhmmer.easel.TextSequence(name=k.encode(),
            sequence=re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]','X',v)).digitize(alpha)
          for k,v in seqs.items() if 250<=len(v)<=1400]
    sc={}
    for hits in pyhmmer.hmmsearch(HMMS,digs,E=1e-5,cpus=6):
        for h in hits:
            n=nm(h.name)
            if n not in sc or h.score>sc[n]: sc[n]=h.score
    cand={k:seqs[k] for k,v in sc.items() if v>=100}
    say(f'   {len(sc)} hit the monoamine model, {len(cand)} score >=100')
    with open(f'{W}/{tag}_cand.faa','w') as f:
        for k,v in cand.items(): f.write(f'>{k}\n{v}\n')
    subprocess.run(['makeblastdb','-in',f'{W}/{tag}_cand.faa','-dbtype','prot','-out',f'{W}/{tag}db'],
                   capture_output=True)
    def best(q,db):
        r=subprocess.run(['blastp','-query',q,'-db',db,'-outfmt','6 qseqid sseqid bitscore',
                          '-evalue','1e-5','-max_target_seqs','5'],capture_output=True,text=True)
        o={}
        for l in r.stdout.strip().split('\n'):
            if not l: continue
            p=l.split('\t')
            if p[0] not in o: o[p[0]]=(p[1],float(p[2]))
        return o
    fwd=best(f'{W}/cel.faa',f'{W}/{tag}db'); rev=best(f'{W}/{tag}_cand.faa',f'{W}/cel')
    calls={}
    for g in sorted(cel):
        if g not in fwd: calls[g]=None; continue
        t,b=fwd[g]
        calls[g]=(t,b) if rev.get(t,('',0))[0]==g else None
    rows.append((tag,sci,clade,len(seqs),calls,''))
    got=[g[4:] for g in sorted(calls) if calls[g]]
    say(f'   reciprocal: {", ".join(got) if got else "none"}')

say('')
say('='*100)
say(f'{"clade":<34} {"species":<24} {"proteins":>9}  dop-1 dop-2 dop-3 dop-4 dop-5 dop-6')
say('-'*100)
for tag,sci,clade,n,calls,note in sorted(rows,key=lambda r:r[2]):
    marks=' '.join(('  +  ' if calls.get(f'cel_dop-{i}') else '  -  ') for i in range(1,7))
    say(f'{clade:<34} {sci:<24} {n:>9}  {marks} {note}')
say(f'{"":<34} {"C. elegans (reference)":<24} {"":>9}    +     +     +     +     +     +')
with open(f'{OUT}/nematode_dating.tsv','w') as f:
    f.write('tag\tspecies\tclade\tproteins\t'+'\t'.join(f'dop-{i}' for i in range(1,7))+'\n')
    for tag,sci,clade,n,calls,note in rows:
        f.write(f'{tag}\t{sci}\t{clade}\t{n}\t'+
                '\t'.join(('present' if calls.get(f'cel_dop-{i}') else 'absent') for i in range(1,7))+'\n')
open(f'{OUT}/nematode_dating.txt','w').write('\n'.join(rep)+'\n')
say('\n-> output/nematode_dating.tsv')
