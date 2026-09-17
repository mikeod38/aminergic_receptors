#!/usr/bin/env python3
"""nematode_expansion.py — how far does the clade holding dop-3/dop-5/dop-6 extend?

CAUTION (2026-09-17, referee audit). This script reports the deepest clade around dop-6 with no
non-nematode member, and that output was previously over-read in two ways. Both are now printed in
the report itself, but read them here before reusing this script:

  - Every species the clade actually returns is RHABDITINA. "No non-nematode member" is a statement
    about the panel, not about the phylum. Whether the radiation is Rhabditina-restricted or
    nematode-wide is unresolved.
  - The taxa that appear absent are absent FROM THE PANEL. Real orthologues exist in the proteome
    store and were never added (Cbr dop-3 = Q61H86, Ppa dop-5 = A0A2A6CA79, plus Necator and
    Strongyloides dop-6 and Haemonchus dop-5). The non-Rhabditina nematodes returned nothing because
    reciprocal best hits cannot resolve this radiation - the retracted dating artefact.

So this script CANNOT date the duplications and must not be used to. Dating needs phylogenetic
placement of the unplaced sequences above.

Reports, for each nested clade around dop-6: taxonomic composition, the deepest clade with no
outside member, and which lineages in the panel carry each gene.

Writes output/tree/nematode_expansion.txt
"""
import os, re
from collections import Counter, defaultdict

HERE=os.path.dirname(os.path.abspath(__file__)); T=os.environ.get('DOP6_TREEDIR') or os.path.join(HERE,'output','tree')
NEM={'Cbr':'C. briggsae','Cbn':'C. brenneri','Cre':'C. remanei','Ppa':'P. pacificus',
     'Asu':'A. suum','Bma':'B. malayi','Sra':'S. ratti','Hco':'H. contortus',
     'Tsp':'T. spiralis','Nam':'N. americanus','cel':'C. elegans'}
# rough phylogenetic placement, for reading duplication timing
CLADE={'C. elegans':'Rhabditina (Caenorhabditis)','C. briggsae':'Rhabditina (Caenorhabditis)',
 'C. brenneri':'Rhabditina (Caenorhabditis)','C. remanei':'Rhabditina (Caenorhabditis)',
 'P. pacificus':'Rhabditina (Diplogastero.)','N. americanus':'Rhabditina (Strongylida)',
 'H. contortus':'Rhabditina (Strongylida)','S. ratti':'Tylenchina','A. suum':'Spirurina',
 'B. malayi':'Spirurina','T. spiralis':'Dorylaimia (basal)'}
rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)
def parse(nw):
    nw=nw.strip().rstrip(';'); pos=[0]
    def node():
        if nw[pos[0]]=='(':
            pos[0]+=1; k=[]
            while True:
                k.append(node())
                if nw[pos[0]]==',': pos[0]+=1; continue
                if nw[pos[0]]==')': pos[0]+=1; break
            m=re.match(r'[^,\)\(:]*',nw[pos[0]:]); lab=m.group(0); pos[0]+=len(lab)
            if pos[0]<len(nw) and nw[pos[0]]==':':
                m2=re.match(r':[-\d.eE+]*',nw[pos[0]:]); pos[0]+=len(m2.group(0))
            return {'tip':False,'lab':lab,'kids':k}
        m=re.match(r'[^,\)\(:]*',nw[pos[0]:]); lab=m.group(0); pos[0]+=len(lab)
        if pos[0]<len(nw) and nw[pos[0]]==':':
            m2=re.match(r':[-\d.eE+]*',nw[pos[0]:]); pos[0]+=len(m2.group(0))
        return {'tip':True,'lab':lab,'kids':[]}
    return node()
def tips(n): return [n['lab']] if n['tip'] else [x for k in n['kids'] for x in tips(k)]
def sp_of(t):
    p=t.split('_',1)[0]
    return NEM.get(p)
def is_nem(t): return sp_of(t) is not None

say('='*96); say('HOW FAR DOES THE CLADE HOLDING dop-3 / dop-5 / dop-6 EXTEND?'); say('='*96)
for trim in ('gt0.1','gappyout'):
    tf=f'{T}/trees/{trim}.treefile'
    if not os.path.exists(tf): continue
    t=parse(open(tf).read()); path=[]
    def find(n,anc):
        if n['tip']:
            if n['lab']=='cel_dop-6': path.extend(anc)
            return
        for k in n['kids']: find(k,anc+[n])
    find(t,[])
    say(f'\n--- {trim}')
    say(f'{"tips":>5} {"UFboot":>7} {"nematode":>9} {"other":>6}   first non-nematode members')
    deepest=None
    for n in reversed(path):
        tl=tips(n)
        if len(tl)>40: break
        m=re.match(r'([\d.]+)/([\d.]+)',n['lab'] or '')
        sup=f'{float(m.group(2)):.0f}' if m else '-'
        nm=[x for x in tl if is_nem(x)]; ot=[x for x in tl if not is_nem(x)]
        if not ot: deepest=(len(tl),sup,tl)
        say(f'{len(tl):>5} {sup:>7} {len(nm):>9} {len(ot):>6}   {", ".join(sorted(ot)[:4]) or "-"}')
    if deepest:
        sz,sup,tl=deepest
        genes=sorted({('cel_'+x.split("_",1)[1]) if not x.startswith('cel_') else x for x in tl})
        say(f'\n  deepest clade containing dop-6 with no outside member (SEE CAUTION BELOW): {sz} tips, UFboot {sup}')
        say(f'  genes in it: {", ".join(g[4:] for g in genes)}')
        bysp=defaultdict(list)
        for x in tl:
            g=x.split("_",1)[1] if not x.startswith('cel_') else x[4:]
            bysp[sp_of(x)].append(g)
        say('')
        say(f'  {"species":<24} {"clade":<28} genes present')
        for s in sorted(bysp,key=lambda z: CLADE.get(z,'zz')):
            say(f'  {s:<24} {CLADE.get(s,"?"):<28} {", ".join(sorted(bysp[s]))}')
say('')
say('-'*96)
say('Reading (rewritten 2026-09-17 after referee audit - the previous text overreached twice):')
say('')
say('  1. EVERY SPECIES IN THIS CLADE IS RHABDITINA. The clade column above says so. A clade with')
say('     no non-nematode member is NOT the same as a nematode-wide radiation, and this one is')
say('     restricted to Caenorhabditis + Pristionchus. Do not write "nematode-specific".')
say('  2. THE ABSENCES ARE SAMPLING, NOT BIOLOGY. Real orthologues exist and are not in the panel:')
say('     Cbr dop-3 = Q61H86 (1092 bits, and already used in dnds_vs_briggsae.tsv), Ppa dop-5 =')
say('     A0A2A6CA79, plus clean RBH dop-6 in Necator (1110 bits) and Strongyloides and dop-5 in')
say('     Haemonchus, all sitting unused in ~/data/proteomes/. The species table below is Compara')
say('     coverage, not a survey. The non-Rhabditina nematodes in the panel returned no dop-3/5/6')
say('     orthologue because reciprocal best hits cannot resolve this radiation - the retracted')
say('     dating artefact - not because none exists.')
say('  3. NO ANCESTRAL STATE WAS RECONSTRUCTED. "From a D2-like ancestor" was asserted from the')
say('     sister relationship alone. Placement inside the D2-like radiation is supported; the')
say('     ancestral ligand is not.')
say('')
say('What survives: the characterised evidence near DOP-6 is thin - its first characterised')
say('neighbour is dop-3, four genes out. Whether that is structural or a sampling failure is OPEN')
say('until the sequences named above are placed phylogenetically.')
open(f'{T}/nematode_expansion.txt','w').write('\n'.join(rep)+'\n')
say('\n-> output/tree/nematode_expansion.txt')
