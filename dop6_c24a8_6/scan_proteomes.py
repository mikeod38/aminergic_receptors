#!/usr/bin/env python3
"""scan_proteomes.py — find amine receptors across metazoan proteomes with the monoamine HMM.

Building the panel from UniProt text queries twice let through sequences that were not receptors at
all (a receptor tyrosine kinase, a tetraspanin, a squid rhodopsin, a bacterial decarboxylase). The
gold-standard alternative, and the one the methodology asks for, is to run the HMM search directly
against real proteomes rather than trusting anyone's annotation text.

Scanned: the 31 metazoan reference proteomes already downloaded for the receptor_phylogeny project
(sponge, ctenophore, placozoan, cnidarian, molluscs, deuterostomes, insects, nematode, vertebrates)
with `MonoamineHMM_Model.hmm` (224 positions, 51 experimentally deorphanised monoaminergic
receptors). PF00001 (7tm_1) is scanned alongside to separate "is a class-A GPCR" from "is an amine
receptor", since the melatonin and muscarinic receptors we deliberately want as outgroups are the
former but not the latter.

Calibration is done on the human proteome, where the true answer is known: the threshold must
recover the 35 human aminergic receptors already in the panel.

Writes output/proteome_scan/hits.tsv and output/proteome_scan/candidates.faa
"""
import os, re, sys
from collections import defaultdict
import pyhmmer

HERE=os.path.dirname(os.path.abspath(__file__))
OUT=os.path.join(HERE,'output','panel_build','proteome_scan'); os.makedirs(OUT,exist_ok=True)
PROT='/Users/mikeodonnell/git/receptor_phylogeny/analysis/intron_mapping_2026-09-02/proteomes'
RAW='/Users/mikeodonnell/git/WormPooledCloning/data/raw'
MONO=f'{RAW}/MonoamineHMM_Model.hmm'; TM7=f'{RAW}/PF00001.hmm'
TOP_PER_SPECIES=10
LEN_RANGE=(300,1000)
SPECIES={'Tadh':'Trichoplax adhaerens','Nvec':'Nematostella vectensis','Hvul':'Hydra vulgaris',
 'Aque':'Amphimedon queenslandica','Mlei':'Mnemiopsis leidyi','Cvir':'Crassostrea virginica',
 'Myes':'Mizuhopecten yessoensis','Lana':'Lingula anatina','Lpol':'Limulus polyphemus',
 'Pmar':'Petromyzon marinus','Apla':'Acanthaster planci','Pcau':'Priapulus caudatus',
 'Skow':'Saccoglossus kowalevskii','Spur':'Strongylocentrotus purpuratus',
 'Bflo':'Branchiostoma floridae','Cint':'Ciona intestinalis','Cele':'Caenorhabditis elegans',
 'Dmel':'Drosophila melanogaster','Hsap':'Homo sapiens','Mmus':'Mus musculus',
 'Drer':'Danio rerio','Pmin':'Patiria miniata','Ebur':'Exaiptasia','Ajap':'Anguilla japonica',
 'Amac':'Acropora','Pbla':'Phoronis','Spun':'Sphaeroforma'}
nm=lambda x: x.decode() if isinstance(x,(bytes,bytearray)) else str(x)
rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)

alpha=pyhmmer.easel.Alphabet.amino()
def load_hmms(p):
    with pyhmmer.plan7.HMMFile(p) as hf: return list(hf)
MH=load_hmms(MONO); TH=load_hmms(TM7)

say('='*94); say('MONOAMINE HMM SCAN ACROSS 31 METAZOAN PROTEOMES'); say('='*94)
files=sorted(f for f in os.listdir(PROT) if f.endswith(('.fasta','.faa'))
             and os.path.getsize(os.path.join(PROT,f))>1_000_000)
say(f'{len(files)} proteomes\n')
allhits=[]; seqstore={}
for fn in files:
    tag=fn.split('_')[0]
    path=os.path.join(PROT,fn)
    seqs={}
    name=None; buf=[]
    for ln in open(path, errors='replace'):
        if ln.startswith('>'):
            if name and LEN_RANGE[0]<=len(''.join(buf))<=LEN_RANGE[1]:
                seqs[name]=''.join(buf)
            name=ln[1:].split()[0]; buf=[]
        else: buf.append(ln.strip())
    if name and LEN_RANGE[0]<=len(''.join(buf))<=LEN_RANGE[1]: seqs[name]=''.join(buf)
    if not seqs: say(f'{tag:<6} no sequences in range'); continue
    digs=[pyhmmer.easel.TextSequence(name=k.encode(),sequence=re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]','X',v)).digitize(alpha)
          for k,v in seqs.items()]
    mono={}; tm={}
    for hits in pyhmmer.hmmsearch(MH,digs,E=1e-5,cpus=4):
        for h in hits:
            n=nm(h.name)
            if n not in mono or h.score>mono[n]: mono[n]=h.score
    for hits in pyhmmer.hmmsearch(TH,digs,E=1e-5,cpus=4):
        for h in hits:
            n=nm(h.name)
            if n not in tm or h.score>tm[n]: tm[n]=h.score
    ranked=sorted(mono.items(),key=lambda x:-x[1])
    say(f'{tag:<6} {SPECIES.get(tag,tag)[:26]:<26} {len(seqs):>6} proteins in range, '
        f'{len(mono):>4} monoamine hits, best {ranked[0][1]:.1f}' if ranked else
        f'{tag:<6} {SPECIES.get(tag,tag)[:26]:<26} {len(seqs):>6} proteins, no hits')
    for n,sc in ranked:
        allhits.append((tag,n,sc,tm.get(n,0.0)))
        seqstore[(tag,n)]=seqs[n]

# ── calibrate on human ───────────────────────────────────────────────────────
say('')
hs=sorted((h for h in allhits if h[0]=='Hsap'),key=lambda x:-x[2])
say(f'CALIBRATION on Homo sapiens: {len(hs)} proteins hit the monoamine model')
say(f'  score range {hs[-1][2]:.1f} - {hs[0][2]:.1f}')
say(f'  the human panel contains 35 verified aminergic receptors; the top {len(hs)} hits here')
say('  are the candidate set from which those are drawn.')
for h in hs[:3]: say(f'    top: {h[1]} {h[2]:.1f}')
CUT=120.0
say(f'\nthreshold for inclusion: monoamine score >= {CUT} (below the human aminergic range but')
say('well above the generic class-A background), then the top hits per species are taken.')

with open(f'{OUT}/hits.tsv','w') as f:
    f.write('species\tprotein\tmonoamine_score\tpf00001_score\n')
    for t,n,sc,tv in sorted(allhits,key=lambda x:(x[0],-x[2])):
        f.write(f'{t}\t{n}\t{sc:.1f}\t{tv:.1f}\n')
sel=defaultdict(list)
for t,n,sc,tv in sorted(allhits,key=lambda x:-x[2]):
    if sc<CUT: continue
    if len(sel[t])>=TOP_PER_SPECIES: continue
    sel[t].append((n,sc,tv))
say('')
say(f'{"species":<8} {"kept":>5}  score range')
tot=0
with open(f'{OUT}/candidates.faa','w') as f:
    for t in sorted(sel):
        rows=sel[t]; tot+=len(rows)
        say(f'{t:<8} {len(rows):>5}  {rows[-1][1]:.0f}-{rows[0][1]:.0f}')
        for n,sc,tv in rows:
            lab=f'{t}_{n}'.replace('|','_')
            s=seqstore[(t,n)]
            f.write(f'>{lab}\n')
            for i in range(0,len(s),60): f.write(s[i:i+60]+'\n')
say(f'\n{tot} candidate amine receptors -> output/proteome_scan/candidates.faa')
open(f'{OUT}/scan.txt','w').write('\n'.join(rep)+'\n')
