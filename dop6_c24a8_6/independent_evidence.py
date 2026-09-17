#!/usr/bin/env python3
"""independent_evidence.py — how much INDEPENDENT evidence sits near DOP-6?

Clade sizes in tips overstate the evidence badly here, because ~40 of the 163 tips are Compara
orthologues of C. elegans genes. Five orthologues of dop-3 are one gene observed five times, not
five independent data points about what dop-3's clade does: they inherit ligand and coupling by
descent, exactly the pseudo-replication problem that inflates significance in shared-character
analyses. The orthologues earn their place by resolving the topology; they must not be counted as
evidence.

So this recounts each nested clade around DOP-6 in three currencies:
  tips                - what the tree draws
  independent genes   - orthologue groups collapsed to one
  characterised genes - of those, how many have an experimentally determined ligand or coupling

and reports the evidence tier of each characterised one, because "DOP-3 is Go" is an inference from
a mutant phenotype, not a binding assay, and that distinction should not vanish into a clade count.

Writes output/tree/independent_evidence.txt
"""
import os, re
from collections import Counter

HERE=os.path.dirname(os.path.abspath(__file__)); T=os.environ.get('DOP6_TREEDIR') or os.path.join(HERE,'output','tree')
# experimentally determined, with evidence tier, from Koelle 2018 WormBook Table 1 and the lab's
# curated pharmacology table. IDA = direct assay; IMP = inferred from mutant phenotype.
CHAR={
 'cel_dop-1':('DA','Gq','ligand IDA (GTPgS, iodo-LSD binding); coupling IDA'),
 'cel_dop-2':('DA','Gi/o','ligand IDA (iodo-LSD, FSK cAMP); coupling IDA'),
 'cel_dop-3':('DA','Go','ligand IDA (FSK cAMP); coupling IMP, from mutant phenotype'),
 'cel_dop-4':('DA','Gs','ligand IDA (cAMP in 293 cells); coupling inferred'),
 'hDRD2':('DA','Gi/o','vertebrate reference'),'hDRD3':('DA','Gi/o','vertebrate reference'),
 'hDRD4':('DA','Gi/o','vertebrate reference'),'hDRD1':('DA','Gs','vertebrate reference'),
 'hDRD5':('DA','Gs','vertebrate reference'),
 'dmDop2R':('DA','Gi/o','insect reference'),'dmDop1R1':('DA','Gs','insect reference'),
 'cel_ser-1':('5HT','Gq','ligand IDA'),'cel_ser-4':('5HT','Gi/o','ligand IDA'),
 'cel_ser-7':('5HT','Gs','coupling IDA'),'cel_ser-2':('TA','Gi/o','ligand IDA'),
 'cel_ser-3':('OA','Gq','ligand IDA'),'cel_ser-6':('OA','Gs','ligand IDA'),
 'cel_octr-1':('OA','Gi/o','ligand IDA; strongest OA evidence in the genome'),
 'cel_tyra-2':('TA','Gi/o','IBA only'),'cel_tyra-3':('TA/OA','Gq','IBA only'),
}
for k in ['hHTR1A','hHTR1B','hHTR1D','hHTR1E','hHTR1F']: CHAR[k]=('5HT','Gi/o','vertebrate reference')
for k in ['hHTR2A','hHTR2B','hHTR2C']: CHAR[k]=('5HT','Gq','vertebrate reference')
for k in ['hADRA2A','hADRA2B','hADRA2C']: CHAR[k]=('NE','Gi/o','vertebrate reference')
NEMPRE=('Cbr_','Cbn_','Cre_','Ppa_','Asu_','Bma_','Sra_','Hco_','Tsp_','Nam_')
def gene_of(t):
    """collapse a tip to the independent gene it represents"""
    if t.startswith(NEMPRE):
        # normalise '_' to '-' in the gene part: some orthologue tips are named
        # Cbn_dop_6 while the C. elegans tip is cel_dop-6, and without this the two
        # never match, so dop-6 gets counted as two independent genes at every level.
        return 'cel_'+t.split('_',1)[1].replace('_','-')
    return t
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

say('='*98); say('HOW MUCH INDEPENDENT EVIDENCE IS NEAR DOP-6?'); say('='*98)
say('Orthologue groups are collapsed to one gene: five orthologues of dop-3 are ONE observation')
say('about what dop-6\'s neighbourhood does, not five.\n')
t=parse(open(f'{T}/trees/gt0.1.treefile').read()); path=[]
def find(n,anc):
    if n['tip']:
        if n['lab']=='cel_dop-6': path.extend(anc)
        return
    for k in n['kids']: find(k,anc+[n])
find(t,[])
say(f'{"tips":>5} {"UFboot":>7} {"indep genes":>12} {"characterised":>14}   which')
say('-'*98)
for n in reversed(path):
    tl=tips(n)
    if len(tl)>45: break
    m=re.match(r'([\d.]+)/([\d.]+)',n['lab'] or '')
    sup=f'{float(m.group(2)):.0f}' if m else '-'
    genes=sorted(set(gene_of(x) for x in tl))
    ch=[g for g in genes if g in CHAR]
    lst=', '.join(f'{g[4:] if g.startswith("cel_") else g}({CHAR[g][1]})' for g in ch)
    say(f'{len(tl):>5} {sup:>7} {len(genes):>12} {len(ch):>14}   {lst or "none"}')
say('')
say('EVIDENCE TIER of each characterised gene that appears in the nested clades above:')
seen=set()
for n in reversed(path):
    tl=tips(n)
    if len(tl)>45: break
    for g in sorted(set(gene_of(x) for x in tl)):
        if g in CHAR and g not in seen:
            seen.add(g)
            lg,cp,note=CHAR[g]
            say(f'   {(g[4:] if g.startswith("cel_") else g):<12} ligand {lg:<6} coupling {cp:<6} {note}')
open(f'{T}/independent_evidence.txt','w').write('\n'.join(rep)+'\n')
say('\n-> output/tree/independent_evidence.txt')
