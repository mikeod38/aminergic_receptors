#!/usr/bin/env python3
"""export_annotations.py — hand the tree to FigTree and iTOL with this project's annotations attached.

FigTree v1.4.4 is already installed and is better than a hand-written viewer at everything generic:
zooming, rerooting, ordering, exporting. What it cannot know is the project-specific layer - the
residue at position 3.32, the Galpha coupling from Koelle 2018 Table 1, the ligand evidence tier, and
which tips are boundary-test genes. So rather than maintain a renderer, emit those as annotations the
standard tools consume:

  <out>.figtree.nexus   NEXUS with a per-tip [&...] annotation block; open directly in FigTree and
                        colour/filter by any field from the Tip Labels and Appearance panels
  <out>_itol_*.txt      iTOL datasets: ligand colour strip, coupling text, 3.32 residue, and a
                        binary marker for the genes under test

Usage: export_annotations.py <treefile> <outprefix>
"""
import os, re, sys
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import importlib.util
spec=importlib.util.spec_from_file_location('mtv',
    os.path.join(os.path.dirname(os.path.abspath(__file__)),'make_tree_viewer.py'))
tf, pre = sys.argv[1], sys.argv[2]
HERE=os.path.dirname(os.path.abspath(__file__))

# reuse the annotation tables by importing them out of the viewer generator
src=open(os.path.join(HERE,'make_tree_viewer.py')).read()
ns={'__name__':'notmain','__file__':os.path.join(HERE,'make_tree_viewer.py'),
    'sys':sys,'os':os,'re':re,'json':__import__('json')}
body=src[:src.find('data=json.dumps')]
body=body.replace("tf   = sys.argv[1]",f"tf   = {tf!r}").replace("out  = sys.argv[2]","out  = '/dev/null'")
body=body.replace("title= sys.argv[3] if len(sys.argv)>3 else os.path.basename(tf)","title=''")
exec(body,ns)
LIG,COUP,POCK,pretty = ns['LIG'],ns['COUP'],ns['POCK'],ns['pretty']
nw_probe=open(tf).read()
BASIS={}
best=(-1,None)
for cand in ('output/tree/basis_b120.tsv','output/tree/basis.tsv'):
    fp=os.path.join(HERE,cand)
    if not os.path.exists(fp): continue
    d={}
    for ln in open(fp).read().strip().split('\n')[1:]:
        q=ln.split('\t'); d[q[0]]=q[1].split(' (')[0]
    # choose the basis table whose taxa actually appear in this tree - the wrong one silently
    # produces an empty dataset, which is how the test-set file came out with no DATA rows
    hit=sum(1 for k in d if k in nw_probe)
    if hit>best[0]: best=(hit,d)
if best[1]: BASIS=best[1]
print(f'   basis table matched {best[0]} taxa in this tree')
nw=open(tf).read().strip()
tipnames=re.findall(r'[\(,]([^\(\),:]+):',nw)
tipnames=[t for t in tipnames if t and not re.match(r'^[\d.]+/[\d.]+$',t)]
LC={'OA':'#C15C3E','TA':'#E0A458','TA/OA':'#D4813F','NE':'#8E3E29','DA':'#1A2C7A',
    '5-HT':'#3C8DAD','5-HT?':'#9EC5D4','ACh':'#7A7A7A','MEL':'#6B4E8E','HA':'#4F8A5B'}

# ── FigTree NEXUS ────────────────────────────────────────────────────────────
ann=nw
for t in sorted(set(tipnames),key=len,reverse=True):
    sp,g,kind=pretty(t)
    a=[f'ligand="{LIG.get(t,"unknown")}"',f'coupling="{COUP.get(t,"unknown")}"',
       f'res332="{POCK.get(t,"?")}"',f'basis="{BASIS.get(t,kind)}"',
       f'species="{sp or "?"}"',f'gene="{g}"',
       f'ligcolour="{LC.get(LIG.get(t,""),"#FFFFFF")}"']
    ann=re.sub(r'(?<![\w.])'+re.escape(t)+r'(?=:)', t+'[&'+','.join(a)+']', ann, count=1)
nexus=('#NEXUS\nbegin taxa;\n\tdimensions ntax=%d;\n\ttaxlabels\n%s\n;\nend;\n\n'
       'begin trees;\n\ttree tree_1 = [&R] %s\nend;\n'
       % (len(set(tipnames)), '\n'.join('\t'+t for t in sorted(set(tipnames))), ann))
open(pre+'.figtree.nexus','w').write(nexus)

# ── iTOL datasets ────────────────────────────────────────────────────────────
def itol(name, header, rows):
    open(f'{pre}_itol_{name}.txt','w').write(header+'\n'+'\n'.join(rows)+'\n')
itol('ligand',
 'DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\tligand\nCOLOR\t#333333\n'
 'LEGEND_TITLE\texperimentally determined ligand\n'
 'LEGEND_SHAPES\t'+'\t'.join('1' for _ in LC)+'\n'
 'LEGEND_COLORS\t'+'\t'.join(LC.values())+'\n'
 'LEGEND_LABELS\t'+'\t'.join(LC.keys())+'\nDATA',
 [f'{t}\t{LC[LIG[t]]}\t{LIG[t]}' for t in sorted(set(tipnames)) if LIG.get(t) in LC])
itol('coupling','DATASET_TEXT\nSEPARATOR TAB\nDATASET_LABEL\tGalpha coupling\nCOLOR\t#555555\nDATA',
 [f'{t}\t{COUP[t]}\t-1\t#555555\tnormal\t1\t0' for t in sorted(set(tipnames)) if COUP.get(t)])
itol('res332','DATASET_TEXT\nSEPARATOR TAB\nDATASET_LABEL\tresidue 3.32\nCOLOR\t#B03030\nDATA',
 [f'{t}\t{POCK[t]}\t-1\t{"#B03030" if POCK[t]!="D" else "#333333"}\tbold\t1\t0'
  for t in sorted(set(tipnames)) if POCK.get(t)])
subj=[t for t in sorted(set(tipnames))
      if BASIS.get(t) in ('subject','boundary test') or t.startswith('TEST_')
      or t in ('cel_dop-6','cel_dop-5')]
itol('testset','DATASET_BINARY\nSEPARATOR TAB\nDATASET_LABEL\tunder test\nCOLOR\t#B03030\n'
     'FIELD_SHAPES\t2\nFIELD_LABELS\tunder test\nFIELD_COLORS\t#B03030\nDATA',
     [f'{t}\t1' for t in subj])
if not subj: print('   WARNING: test-set dataset is empty - no taxa matched')
# ── tidy table, the natural input for ggtree's %<+% join and for any other tool ──────────
with open(pre+'.annotations.tsv','w') as f:
    f.write('label\tgene\tspecies\tligand\tcoupling\tres332\tbasis\tkind\tligcolour\n')
    for t in sorted(set(tipnames)):
        sp,g,kind=pretty(t)
        f.write('\t'.join([t,g,sp or 'NA',LIG.get(t,'NA'),COUP.get(t,'NA'),POCK.get(t,'NA'),
                            BASIS.get(t,kind),kind,LC.get(LIG.get(t,''),'NA')])+'\n')
print(f'-> {pre}.annotations.tsv      (tidy per-tip table; ggtree, R, pandas)')
print(f'-> {pre}.figtree.nexus        (open in FigTree; annotations on the Tip Labels panel)')
for n in ('ligand','coupling','res332','testset'):
    print(f'-> {pre}_itol_{n}.txt')
