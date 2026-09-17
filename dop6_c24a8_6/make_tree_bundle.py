#!/usr/bin/env python3
"""make_tree_bundle.py — ship a tree and its annotations as one versioned, validated unit.

Annotation files and trees drift apart. iTOL accepts a dataset whose IDs match no tip and renders
nothing; FigTree does the same with an unmatched NEXUS block. Both fail silently, so a stale overlay
looks like a real result with a feature missing. That happened here: a test-set dataset built against
one taxon set was loaded onto another and came out completely empty.

A bundle fixes it by construction. Each is a directory holding the tree, every annotation file built
FROM that tree, and a manifest recording the tree's SHA-256, its tip count, and - the part that
matters - the fraction of each dataset's IDs that actually matched a tip. Anything below 100 % is
flagged in the manifest rather than discovered on screen.

Usage: make_tree_bundle.py <treefile> <bundle-name> [--aln FILE] [--note TEXT]
"""
import hashlib, os, re, shutil, subprocess, sys, datetime

HERE=os.path.dirname(os.path.abspath(__file__))
tf=sys.argv[1]; name=sys.argv[2]
aln=None; note=''
for i,a in enumerate(sys.argv):
    if a=='--aln' and i+1<len(sys.argv): aln=sys.argv[i+1]
    if a=='--note' and i+1<len(sys.argv): note=sys.argv[i+1]
BD=os.path.join(HERE,'output','tree','bundles',name)
os.makedirs(BD,exist_ok=True)
nw=open(tf).read()
sha=hashlib.sha256(nw.encode()).hexdigest()[:16]
tips=[t for t in re.findall(r'[\(,]([^\(\),:]+):',nw) if t and not re.match(r'^[\d.]+/[\d.]+$',t)]
tips=sorted(set(tips))
shutil.copy(tf, os.path.join(BD,'tree.newick'))

# build annotations from THIS tree
r=subprocess.run([sys.executable,os.path.join(HERE,'export_annotations.py'),tf,os.path.join(BD,'ann')],
                 capture_output=True,text=True)
if r.returncode: print(r.stdout+r.stderr); sys.exit(1)

# validate: what fraction of each dataset's IDs are real tips of this tree?
rows=[]
for f in sorted(os.listdir(BD)):
    if not f.startswith('ann_itol_'): continue
    path=os.path.join(BD,f)
    lines=open(path).read().split('\n')
    try: di=lines.index('DATA')
    except ValueError: rows.append((f,0,0,'no DATA block')); continue
    ids=[l.split('\t')[0] for l in lines[di+1:] if l.strip()]
    hit=sum(1 for i in ids if i in tips)
    pct=100*hit/len(ids) if ids else 0
    flag='' if ids and hit==len(ids) else ('EMPTY' if not ids else f'{len(ids)-hit} unmatched')
    rows.append((f,len(ids),pct,flag))
nex=os.path.join(BD,'ann.figtree.nexus')
nexhit=sum(1 for t in tips if t+'[&' in open(nex).read()) if os.path.exists(nex) else 0

model=''
iq=tf.replace('.treefile','.iqtree')
if os.path.exists(iq):
    for ln in open(iq):
        if ln.startswith('Best-fit model'): model=ln.split(':')[1].strip(); break
man=[f'# Tree bundle: {name}','',
     f'Built {datetime.date.today().isoformat()}.'+(f'  {note}' if note else ''),'',
     '| | |','|---|---|',
     f'| tree | `tree.newick` |',f'| source | `{os.path.relpath(tf,HERE)}` |',
     f'| sha256 (first 16) | `{sha}` |',f'| tips | {len(tips)} |',
     f'| model | {model or "n/a"} |']
if aln: man.append(f'| alignment | `{os.path.relpath(aln,HERE)}` |')
man += ['','## Annotation datasets, validated against THIS tree','',
        '| file | rows | IDs matching a tip | status |','|---|---|---|---|']
for f,n,pct,flag in rows:
    man.append(f'| `{f}` | {n} | {pct:.0f} % | {flag or "ok"} |')
man.append(f'| `ann.figtree.nexus` | {nexhit} tips annotated | '
           f'{100*nexhit/len(tips):.0f} % | {"ok" if nexhit==len(tips) else "incomplete"} |')
man += ['','## Use','',
 '**iTOL** — upload `tree.newick`, then drag the `ann_itol_*.txt` files onto the rendered tree.',
 '**FigTree** — open `ann.figtree.nexus`; annotations appear on the Tip Labels panel.','',
 'Do not mix these files with another tree. Every dataset here was generated from the tree whose',
 'sha256 is recorded above; loading them onto a different taxon set fails silently, showing an',
 'empty or partial overlay with no error.','']
open(os.path.join(BD,'BUNDLE.md'),'w').write('\n'.join(man)+'\n')
bad=[r for r in rows if r[3]]
print(f'-> output/tree/bundles/{name}/   ({len(tips)} tips, sha {sha})')
for f,n,pct,flag in rows: print(f'   {f:<26} {n:>4} rows  {pct:>3.0f}% matched  {flag}')
if bad: print(f'   WARNING: {len(bad)} dataset(s) do not fully match this tree')
