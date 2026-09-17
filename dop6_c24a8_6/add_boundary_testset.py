#!/usr/bin/env python3
"""add_boundary_testset.py — convert an annotation-based panel boundary into a tested one.

The panel excluded every gene labelled a neuropeptide receptor and included cel_F59D12.1, even
though 36 of the excluded genes score HIGHER on the monoamine HMM than F59D12.1 does (52.8). The
boundary was therefore drawn by nomenclature, not by evidence - the practice the methodology
forbids in both directions.

So the highest-scoring excluded genes are added to the panel as an explicitly labelled TEST SET. If
they form a coherent group outside the aminergic radiation, the exclusion was right and the score at
which they separate becomes a defensible, tested boundary. If any falls inside the radiation, it was
mis-annotated - which matters directly, since the project is looking for an unidentified octopamine
receptor and these genes have never been excluded by anything except their names.

Writes output/tree/panel_boundary.faa and boundary_testset.tsv
"""
import os, re
HERE=os.path.dirname(os.path.abspath(__file__)); T=os.path.join(HERE,'output','tree')
SC='/Users/mikeodonnell/git/WormPooledCloning/data/processed/gpcr/cel_gpcr_monoamine_scores.tsv'
PEP=('/private/tmp/claude-501/-Users-mikeodonnell-git-Malaiwong-Oglu-2026/'
     '93015e40-e331-495a-8362-20ac7ea07c75/scratchpad/cel_pep.fa')
N_TEST=12
rep=[]
def say(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True); rep.append(s)

panel={}
for r in open(f'{T}/panel.faa').read().split('>')[1:]:
    h,s=r.split('\n',1); panel[h.strip()]=s.replace('\n','')
have={k[4:] for k in panel if k.startswith('cel_')}
rows=[]
for l in open(SC).read().strip().split('\n')[1:]:
    p=l.split('\t')
    try: sc=float(p[4])
    except (ValueError,IndexError): continue          # 'no_hit' rows
    rows.append((p[1],p[2],p[3],sc))
rows.sort(key=lambda x:-x[3])
excluded=[r for r in rows if r[1] not in have]
say('='*92); say('TESTING THE PANEL BOUNDARY'); say('='*92)
say(f'F59D12.1 (in the panel) scores 52.8; {sum(1 for r in excluded if r[3]>52.8)} excluded genes score higher.\n')
test=excluded[:N_TEST]
say(f'{"gene":<14} {"WBGene":<18} {"their label":<24} {"score":>7}')
for wb,g,fam,sc in test: say(f'{g:<14} {wb:<18} {fam:<24} {sc:>7.1f}')

# proteins, longest isoform per gene, from the Ensembl proteome
byg={}
for r in open(PEP).read().split('>')[1:]:
    h,s=r.split('\n',1); seq=s.replace('\n','')
    m=re.search(r'gene:(\S+)',h)
    if not m: continue
    g=m.group(1)
    if g not in byg or len(seq)>len(byg[g][1]): byg[g]=(h.split()[0],seq)
add={}
for wb,g,fam,sc in test:
    if wb not in byg: say(f'  {g}: no protein found for {wb}'); continue
    tid,seq=byg[wb]
    add[f'test_{g}']=seq
say(f'\n{len(add)} test sequences recovered')
with open(f'{T}/panel_boundary.faa','w') as f:
    for k in sorted(panel):
        f.write(f'>{k}\n')
        for i in range(0,len(panel[k]),60): f.write(panel[k][i:i+60]+'\n')
    for k in sorted(add):
        f.write(f'>{k}\n')
        for i in range(0,len(add[k]),60): f.write(add[k][i:i+60]+'\n')
with open(f'{T}/boundary_testset.tsv','w') as f:
    f.write('label\tgene\twbgene\ttheir_label\tmonoamine_score\tlength\n')
    for wb,g,fam,sc in test:
        if f'test_{g}' in add:
            f.write(f'test_{g}\t{g}\t{wb}\t{fam}\t{sc}\t{len(add[f"test_{g}"])}\n')
say(f'-> output/tree/panel_boundary.faa ({len(panel)+len(add)} sequences)')
open(f'{T}/boundary_testset.txt','w').write('\n'.join(rep)+'\n')
