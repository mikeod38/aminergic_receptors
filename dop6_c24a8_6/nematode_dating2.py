#!/usr/bin/env python3
"""nematode_dating2.py — RETRACTED 2026-09-17. Its output is a method artefact; do not reuse.

*** THIS SCRIPT CANNOT DO WHAT ITS NAME SAYS. Its results were withdrawn on referee audit. ***

The intent was to bracket the dop-3/dop-5/dop-6 duplications by asking which nematode clades carry
each gene, with Trichinella (Dorylaimia, earliest-diverging) as the decisive taxon. The method is
reciprocal best hits in both directions after an HMM scan of each whole proteome.

WHY IT FAILS: RBH cannot resolve members of a RECENT PARALOGUE RADIATION, which is exactly what
dop-3/dop-5/dop-6 are. When several paralogues are mutually more similar than any is to the query,
the reverse hit lands on whichever paralogue happens to score highest, so a gene that is genuinely
present gets recorded as absent. Diagnostic that settled it: in Haemonchus one protein reciprocates
with dop-5 at 976 bits while reverse hits pile onto dop-4 and dop-1 indiscriminately. Every
"absent" this script emits means "reciprocated to a different paralogue".

Confirmed counterexamples, all already in the shared proteome store, all called absent here:
Necator dop-6 (RBH 1110/1121 bits, 81.7%), Strongyloides dop-6, Haemonchus dop-5.

WHAT SURVIVES: the dop-2 column only - present in all six species and shared with Drosophila.

THE CORRECT METHOD is phylogenetic placement: add every candidate scoring above threshold from each
proteome to the panel and rebuild the tree, letting the topology assign paralogy. Best-hit
assignment cannot. See reviews/2026-09-16_dop6_conclusions/DISPOSITION.md.

The one design decision worth keeping: **absence is only called from a proteome the manifest marks
complete**, because a truncated download silently converts a present gene into a confident absence,
and two files in this project already did exactly that.

Every call is a reciprocal best hit in both directions, after an HMM scan of the whole proteome -
Compara cannot be used here, having missed dop-5 in P. pacificus outright. Proteomes come from the
shared store and **absence is only called from a proteome the manifest marks complete**, because a
truncated download silently converts a present gene into a confident absence; two files in this
project already did exactly that.

Writes output/nematode_dating.tsv and .txt
"""
import os, re, subprocess, sys
import pyhmmer

HERE=os.path.dirname(os.path.abspath(__file__)); OUT=os.path.join(HERE,'output')
STORE=os.path.expanduser('~/data/proteomes')
MONO='/Users/mikeodonnell/git/WormPooledCloning/data/raw/MonoamineHMM_Model.hmm'
W='/tmp/_nemdate'; os.makedirs(W,exist_ok=True)
CLADE={'Tspi':'I   Dorylaimia (earliest-diverging)','Tmur':'I   Dorylaimia',
       'Bmal':'III Spirurina','Srat':'IV  Tylenchina','Hcon':'V   Rhabditina (Strongylida)',
       'Name':'V   Rhabditina (Strongylida)','Ppac':'V   Rhabditina (Diplogasteromorpha)'}
nm=lambda x: x.decode() if isinstance(x,(bytes,bytearray)) else str(x)
rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)

man={}
for ln in open(os.path.join(STORE,'MANIFEST.tsv')).read().strip().split('\n')[1:]:
    p=ln.split('\t')
    if len(p)>=8: man[p[0]]=dict(species=p[1],n=p[4],exp=p[5],complete=p[6],path=p[7],src=p[8])
cel={}
for r in open(f'{OUT}/tree/panel.faa').read().split('>')[1:]:
    h,s=r.split('\n',1)
    if re.match(r'^cel_dop-[1-6]$',h.strip()): cel[h.strip()]=s.replace('\n','')
with open(f'{W}/cel.faa','w') as f:
    for k,v in cel.items(): f.write(f'>{k}\n{v}\n')
subprocess.run(['makeblastdb','-in',f'{W}/cel.faa','-dbtype','prot','-out',f'{W}/cel'],capture_output=True)
alpha=pyhmmer.easel.Alphabet.amino()
with pyhmmer.plan7.HMMFile(MONO) as hf: HMMS=list(hf)

say('='*104); say('WHEN DID dop-3 / dop-5 / dop-6 ARISE? PRESENCE ACROSS NEMATODE CLADES'); say('='*104)
say('reciprocal best hits after a whole-proteome HMM scan; ABSENCE CALLED ONLY FROM A PROTEOME')
say('THE MANIFEST MARKS COMPLETE\n')
rows=[]
for tag in ('Tspi','Tmur','Bmal','Srat','Hcon','Name','Ppac'):
    if tag not in man: say(f'{tag}: not in the store yet'); continue
    m=man[tag]
    if m['complete']!='yes':
        say(f'{tag} ({m["species"]}): INCOMPLETE ({m["n"]}/{m["exp"]}) - skipped, no calls made'); continue
    seqs={}; name=None; buf=[]
    for ln in open(m['path'],errors='replace'):
        if ln.startswith('>'):
            if name: seqs[name]=''.join(buf)
            name=ln[1:].split()[0]; buf=[]
        else: buf.append(ln.strip())
    if name: seqs[name]=''.join(buf)
    digs=[pyhmmer.easel.TextSequence(name=k.encode(),
            sequence=re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]','X',v)).digitize(alpha)
          for k,v in seqs.items() if 250<=len(v)<=1400]
    sc={}
    for hits in pyhmmer.hmmsearch(HMMS,digs,E=1e-5,cpus=6):
        for h in hits:
            n=nm(h.name)
            if n not in sc or h.score>sc[n]: sc[n]=h.score
    cand={k:seqs[k] for k,v in sc.items() if v>=100}
    with open(f'{W}/{tag}.faa','w') as f:
        for k,v in cand.items(): f.write(f'>{k}\n{v}\n')
    subprocess.run(['makeblastdb','-in',f'{W}/{tag}.faa','-dbtype','prot','-out',f'{W}/{tag}db'],
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
    fwd=best(f'{W}/cel.faa',f'{W}/{tag}db'); rev=best(f'{W}/{tag}.faa',f'{W}/cel')
    calls={}
    for g in sorted(cel):
        t=fwd.get(g)
        calls[g]= (t[0],t[1]) if (t and rev.get(t[0],('',0))[0]==g) else None
    rows.append((tag,m['species'],CLADE.get(tag,'?'),int(m['n']),calls))
    got=[g[4:] for g in sorted(calls) if calls[g]]
    say(f'{tag} {m["species"]:<24} {m["n"]:>6} proteins  {len(sc):>3} HMM hits, {len(cand):>2} >=100'
        f'   reciprocal: {", ".join(got) if got else "none"}')
say('')
say('='*104)
say(f'{"clade":<36} {"species":<24} {"n":>6}  ' + '  '.join(f'dop-{i}' for i in range(1,7)))
say('-'*104)
for tag,sp,cl,n,calls in sorted(rows,key=lambda r:r[2]):
    marks='  '.join((' +   ' if calls.get(f'cel_dop-{i}') else ' -   ') for i in range(1,7))
    say(f'{cl:<36} {sp:<24} {n:>6}  {marks}')
say(f'{"V   Rhabditina (Caenorhabditis)":<36} {"C. elegans (reference)":<24} {"":>6}   +     +     +     +     +     +')
say('')
say('READING: dop-2 present in Dorylaimia and shared with Drosophila = the ancestral lineage.')
say('The clade at which dop-3/dop-5/dop-6 first appear brackets the radiation.')
with open(f'{OUT}/nematode_dating.tsv','w') as f:
    f.write('tag\tspecies\tclade\tproteins\t'+'\t'.join(f'dop-{i}' for i in range(1,7))+'\n')
    for tag,sp,cl,n,calls in rows:
        f.write(f'{tag}\t{sp}\t{cl}\t{n}\t'+'\t'.join(
            ('present' if calls.get(f'cel_dop-{i}') else 'absent') for i in range(1,7))+'\n')
open(f'{OUT}/nematode_dating.txt','w').write('\n'.join(rep)+'\n')
say(f'-> output/nematode_dating.tsv')
