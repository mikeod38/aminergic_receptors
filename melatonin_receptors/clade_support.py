import re,sys

ORPH = {'hGPR52','hGPR21','PCDRfam_Spur_A0A7M7GHW5','PCDRfam_Dmel_Q9VZD1',
        'PCDRfam_Pcau_A0ABM1EQA9','PCDRfam_Cvir_A0A8B8D717','PCDRfam_Lana_A0A1S3KFU7'}
MEL = {'hMTNR1A','hMTNR1B','hGPR50','MEL_Drer_A0A8M1NRI5','MEL_Drer_A0A8M9PN98',
       'MEL_Drer_A2BGT9','MEL_Drer_P51046','MEL_Drer_P51049'}
PC = 'cel_pcdr-1_F59D12.1'

class N:
    def __init__(s,par=None): s.par,s.ch,s.name,s.sup=par,[],None,None

def parse(t):
    t=t.strip().rstrip(';'); root=N(); cur=root; i=0; buf=''
    def flush():
        nonlocal buf
        if buf:
            # name[:len] or sup:len
            parts=buf.split(':')
            lab=parts[0]
            if cur.ch: cur.sup=lab or None
            else: cur.name=lab or None
        buf=''
    while i<len(t):
        c=t[i]
        if c=='(':
            n=N(cur); cur.ch.append(n); cur=n
        elif c==',':
            flush(); n=N(cur.par); cur.par.ch.append(n); cur=n
        elif c==')':
            flush(); cur=cur.par
        else: buf+=c
        i+=1
    flush()
    return root

def tips(n):
    if not n.ch: return {n.name}
    s=set()
    for c in n.ch: s|=tips(c)
    return s

def walk(n,out):
    out.append(n)
    for c in n.ch: walk(c,out)
    return out

def sup(n):
    if n.sup: return n.sup
    if n.par is None: return 'root'
    sib=[c for c in n.par.ch if c is not n]
    if n.par.par is None and len(sib)==1: return sib[0].sup or 'unlabelled'
    return 'unlabelled'

def mrca(root,target):
    best=None
    for n in walk(root,[]):
        tp=tips(n)
        if target<=tp and (best is None or len(tp)<len(tips(best))): best=n
    return best

for tag,path in [('ML, PCDR-1 present','output/tree2/gt0.1.treefile'),
                 ('PCDR-1 REMOVED','output/lba/no_pcdr1.treefile'),
                 ('LG+C20+F+G (PMSF)','output/lba/pmsf.treefile')]:
    root=parse(open(path).read())
    all_t=tips(root)
    print(f'=== {tag} ===')
    for lbl,grp in [('orphan family',ORPH),('melatonin clade',MEL)]:
        g=grp & all_t
        n=mrca(root,g)
        tp=tips(n)
        extra=sorted(tp-g)
        print(f'  {lbl}: {len(tp)} tips, support {sup(n)}'
              + (f'  | also inside: {", ".join(extra)}' if extra else '  | exactly monophyletic'))
    if PC in all_t:
        m=mrca(root,MEL & all_t)
        print(f'  pcdr-1 inside melatonin clade? {PC in tips(m)}')
    print()
