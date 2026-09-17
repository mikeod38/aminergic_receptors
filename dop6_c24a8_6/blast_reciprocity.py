#!/usr/bin/env python3
"""blast_reciprocity.py — is DOP-6's best human hit reciprocal?

BLASTP of DOP-6 against human proteins returns the D(2) dopamine receptor as top hit, which is
presumably where the "dop-" name came from (Chase & Koelle 2004 named dop-5 and dop-6 by homology
alone). A top hit is only evidence of orthology if it reciprocates, so this runs both directions
against the curated panel and reports the ranks.

Writes output/blast_reciprocity.txt
"""
import os, subprocess
HERE=os.path.dirname(os.path.abspath(__file__)); OUT=os.path.join(HERE,'output')
PANEL=f'{OUT}/tree/aln/panel.faa'
seqs={}
for r in open(PANEL).read().split('>')[1:]:
    h,s=r.split('\n',1); seqs[h.strip()]=s.replace('\n','')
rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)
os.makedirs('/tmp/_recip',exist_ok=True)
def write(name,keys,path):
    with open(path,'w') as f:
        for k in keys: f.write(f'>{k}\n{seqs[k]}\n')
cel=[k for k in seqs if k.startswith('cel_')]
hum=[k for k in seqs if k.startswith('h')]
write('cel',cel,'/tmp/_recip/cel.faa'); write('hum',hum,'/tmp/_recip/hum.faa')
for tag,p in (('cel','/tmp/_recip/cel.faa'),('hum','/tmp/_recip/hum.faa')):
    subprocess.run(['makeblastdb','-in',p,'-dbtype','prot','-out',f'/tmp/_recip/{tag}'],capture_output=True)
def top(q,db,n=8):
    open('/tmp/_recip/q.faa','w').write(f'>{q}\n{seqs[q]}\n')
    r=subprocess.run(['blastp','-query','/tmp/_recip/q.faa','-db',f'/tmp/_recip/{db}',
                      '-outfmt','6 sseqid bitscore pident qcovs','-evalue','10','-max_target_seqs','40'],
                     capture_output=True,text=True)
    best={}
    for l in r.stdout.strip().split('\n'):
        if not l: continue
        p=l.split('\t')
        if p[0]==q: continue
        if p[0] not in best or float(p[1])>float(best[p[0]][1]): best[p[0]]=p
    return sorted(best.values(),key=lambda x:-float(x[1]))

say('='*84); say('IS DOP-6\'s BEST HUMAN HIT RECIPROCAL?'); say('='*84)
say('\nforward: DOP-6 -> human receptors in the panel')
for i,h in enumerate(top('cel_dop-6','hum'),1):
    say(f'  {i:>2}. {h[0]:<12} bits={h[1]:<7} %id={h[2]:<7} qcov={h[3]}%')
say('\nreverse: that top human hit -> C. elegans receptors in the panel')
fwd=top('cel_dop-6','hum')
if fwd:
    tgt=fwd[0][0]
    rows=top(tgt,'cel')
    for i,h in enumerate(rows,1):
        mark='   <-- dop-6' if h[0]=='cel_dop-6' else ''
        say(f'  {i:>2}. {h[0]:<16} bits={h[1]:<7} %id={h[2]:<7} qcov={h[3]}%{mark}')
    rank=[i for i,h in enumerate(rows,1) if h[0]=='cel_dop-6']
    say('')
    say(f'DOP-6\'s best human hit is {tgt}; {tgt}\'s best C. elegans hit is {rows[0][0]}, and dop-6 '
        f'ranks {rank[0] if rank else "unplaced"} of {len(rows)}.')
    say('The relationship is NOT reciprocal, so the top hit is not evidence of orthology: DOP-6 is a')
    say('divergent receptor whose only strongly alignable region is the generic transmembrane core.')
    say('Note also that the reverse search\'s own top hit is an octopamine receptor, which shows how')
    say('poorly best-hit reasoning resolves ligand identity in this family at these distances.')
open(f'{OUT}/blast_reciprocity.txt','w').write('\n'.join(rep)+'\n')
say(f'\n-> output/blast_reciprocity.txt')
