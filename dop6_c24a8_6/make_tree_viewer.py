#!/usr/bin/env python3
"""make_tree_viewer.py — emit a self-contained interactive HTML viewer for a tree.

A 170-tip PDF is not browsable. This writes a single HTML file (no network, no CDN) with the tree
and its annotations embedded, offering the things a static figure cannot: collapse a clade to a
triangle, search and jump to a taxon, filter to C. elegans or to the genes under test, toggle
support values, and colour by ligand or by Galpha coupling.

Usage: make_tree_viewer.py <treefile> <out.html> [title]
"""
import json, os, re, sys

tf   = sys.argv[1]
out  = sys.argv[2]
title= sys.argv[3] if len(sys.argv)>3 else os.path.basename(tf)
HERE = os.path.dirname(os.path.abspath(__file__))

def parse(nw):
    nw=nw.strip().rstrip(';'); pos=[0]; uid=[0]
    def node():
        if nw[pos[0]]=='(':
            pos[0]+=1; kids=[]
            while True:
                kids.append(node())
                if nw[pos[0]]==',': pos[0]+=1; continue
                if nw[pos[0]]==')': pos[0]+=1; break
            m=re.match(r'[^,\)\(:]*',nw[pos[0]:]); lab=m.group(0); pos[0]+=len(lab); bl=0.0
            if pos[0]<len(nw) and nw[pos[0]]==':':
                m2=re.match(r':([-\d.eE+]*)',nw[pos[0]:]); pos[0]+=len(m2.group(0)); bl=float(m2.group(1) or 0)
            uid[0]+=1
            return {'id':uid[0],'name':'','sup':lab,'bl':bl,'kids':kids}
        m=re.match(r'[^,\)\(:]*',nw[pos[0]:]); lab=m.group(0); pos[0]+=len(lab); bl=0.0
        if pos[0]<len(nw) and nw[pos[0]]==':':
            m2=re.match(r':([-\d.eE+]*)',nw[pos[0]:]); pos[0]+=len(m2.group(0)); bl=float(m2.group(1) or 0)
        uid[0]+=1
        return {'id':uid[0],'name':lab,'sup':'','bl':bl,'kids':[]}
    return node()

tree=parse(open(tf).read())
def tipcount(n): return 1 if not n['kids'] else sum(tipcount(k) for k in n['kids'])

# ── annotations ─────────────────────────────────────────────────────────────
LIG={'cel_dop-1':'DA','cel_dop-2':'DA','cel_dop-3':'DA','cel_dop-4':'DA','hDRD1':'DA','hDRD2':'DA',
 'hDRD3':'DA','hDRD4':'DA','hDRD5':'DA','dmDop1R1':'DA','dmDop2R':'DA','amDop1':'DA',
 'cel_ser-1':'5-HT','cel_ser-4':'5-HT','cel_ser-7':'5-HT','dm5HT1A':'5-HT','dm5HT1B':'5-HT',
 'dm5HT2A':'5-HT','dm5HT7':'5-HT','ap5HT':'5-HT','am5HT1A':'5-HT','cel_ser-5':'5-HT?',
 'cel_ser-3':'OA','cel_ser-6':'OA','cel_octr-1':'OA','dmOAMB':'OA','dmOctA2R':'OA','dmOctB1R':'OA',
 'dmOctB2R':'OA','dmOctB3R':'OA','bmOAR1':'OA','amOA1':'OA','tcOA1':'OA','apOA':'OA',
 'cel_ser-2':'TA','cel_tyra-2':'TA','dmTyrR':'TA','dmTyrRII':'TA','amTAR1':'TA','bmTAR':'TA',
 'cel_tyra-3':'TA/OA','dmOctTyrR':'TA/OA',
 'hCHRM1':'ACh','hCHRM2':'ACh','hCHRM3':'ACh','hCHRM5':'ACh','cel_gar-1':'ACh','cel_gar-2':'ACh',
 'cel_gar-3':'ACh','hMTNR1A':'MEL','hMTNR1B':'MEL','hGPR50':'MEL',
 'hHRH1':'HA','hHRH2':'HA','hHRH3':'HA','hHRH4':'HA','hTAAR1':'TA'}
for k in ['hHTR1A','hHTR1B','hHTR1D','hHTR1E','hHTR1F','hHTR2A','hHTR2B','hHTR2C','hHTR4','hHTR5A',
          'hHTR6','hHTR7']: LIG[k]='5-HT'
for k in ['hADRA1A','hADRA1B','hADRA1D','hADRA2A','hADRA2B','hADRA2C','hADRB1','hADRB2','hADRB3',
          'plADRA1','plADRA2']: LIG[k]='NE'
COUP={'cel_ser-1':'Gq','cel_ser-4':'Gi/o','cel_ser-7':'Gs','cel_dop-1':'Gq','cel_dop-2':'Gi/o',
 'cel_dop-3':'Go','cel_dop-4':'Gs','cel_octr-1':'Go','cel_ser-3':'Gq','cel_ser-6':'Gs',
 'cel_ser-2':'Gi/o','cel_tyra-2':'Gi/o','cel_tyra-3':'Gq','hDRD1':'Gs','hDRD5':'Gs','hDRD2':'Gi/o',
 'hDRD3':'Gi/o','hDRD4':'Gi/o','hADRB1':'Gs','hADRB2':'Gs','hADRB3':'Gs','hADRA1A':'Gq',
 'hADRA1B':'Gq','hADRA1D':'Gq','hADRA2A':'Gi/o','hADRA2B':'Gi/o','hADRA2C':'Gi/o','hHTR1A':'Gi/o',
 'hHTR1B':'Gi/o','hHTR1D':'Gi/o','hHTR1E':'Gi/o','hHTR1F':'Gi/o','hHTR2A':'Gq','hHTR2B':'Gq',
 'hHTR2C':'Gq','hHTR4':'Gs','hHTR6':'Gs','hHTR7':'Gs','hCHRM1':'Gq','hCHRM3':'Gq','hCHRM5':'Gq',
 'hCHRM2':'Gi/o','dmDop1R1':'Gs','dmDop2R':'Gi/o','dmOctB1R':'Gs','dmOctB2R':'Gs','dmOctB3R':'Gs',
 'dmOAMB':'Gq','dmOctA2R':'Gi/o','dmOctTyrR':'Gi/o','dmTyrR':'Gi/o','dmTyrRII':'Gi/o'}
POCK={}
# Read the annotation that belongs to THIS panel, not a hardcoded global one. The viewer used to
# load output/tree/pocket.tsv unconditionally - the v1 table - so after the tips were renamed to
# hyphens it silently failed to match Cbn_dop-6/Cbr_dop-6/Cre_dop-6/Ppa_dop-6 and showed them with
# no residue at 3.32, while cel_dop-6 resolved because it always had a hyphen. Look next to the
# treefile first; fall back to the global table only if the panel has none.
_pdir=os.path.dirname(os.path.dirname(os.path.abspath(tf)))   # <panel>/trees/x.treefile -> <panel>
for cand in (os.path.join(_pdir,'pocket.tsv'),
             os.path.join(HERE,'output','tree','pocket.tsv')):
    if os.path.exists(cand):
        for ln in open(cand).read().strip().split('\n')[1:]:
            p=ln.split('\t'); POCK[p[0]]=p[1]
        break
SPCODE={'Ajap':'A. japonica','Apla':'A. planci','Bflo':'B. floridae','Cint':'C. intestinalis',
 'Cvir':'C. virginica','Drer':'D. rerio','Ebur':'E. diaphana','Hvul':'H. vulgaris','Lana':'L. anatina',
 'Lpol':'L. polyphemus','Mlei':'M. leidyi','Myes':'M. yessoensis','Nvec':'N. vectensis',
 'Pcau':'P. caudatus','Pmar':'P. marinus','Pmin':'P. miniata','Skow':'S. kowalevskii',
 'Spur':'S. purpuratus','Tadh':'T. adhaerens','Aque':'A. queenslandica','Cbr':'C. briggsae',
 'Cbn':'C. brenneri','Cre':'C. remanei','Ppa':'P. pacificus','Asu':'A. suum','Bma':'B. malayi',
 'Sra':'S. ratti','Hco':'H. contortus','Tsp':'T. spiralis','Nam':'N. americanus'}
SCAN=re.compile(r'^([A-Z][a-z]{3})_(?:tr|sp)_([A-Za-z0-9]+)_')
NEMR=re.compile(r'^([A-Z][a-z]{2})_(.+)$')
def pretty(lab):
    if lab.startswith('TEST_'): return 'C. elegans', lab[5:], 'test'
    if lab.startswith('cel_'):  return 'C. elegans', lab[4:], 'cel'
    m=SCAN.match(lab)
    if m:
        a=m.group(2); a=a[3:] if a.startswith('A0A') else a
        return SPCODE.get(m.group(1),m.group(1)), a, 'scan'
    m=NEMR.match(lab)
    if m and m.group(1) in SPCODE: return SPCODE[m.group(1)], m.group(2).replace('_','-'), 'nem'
    for pre,sp in (('h','H. sapiens'),('dm','D. melanogaster'),('am','A. mellifera'),
                   ('bm','B. mori'),('tc','T. castaneum'),('ap','A. californica'),('pl','P. dumerilii')):
        if lab.startswith(pre) and len(lab)>len(pre): return sp, lab[len(pre):], 'ref'
    return '', lab, 'ref'
def deco(n):
    if n['kids']:
        for k in n['kids']: deco(k)
        return
    sp,g,kind=pretty(n['name'])
    n['sp']=sp; n['g']=g; n['kind']=kind
    n['lig']=LIG.get(n['name'],''); n['cp']=COUP.get(n['name'],''); n['p332']=POCK.get(n['name'],'')
deco(tree)
data=json.dumps(tree,separators=(',',':'))

HTML = """<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title><style>
:root{--bg:#fff;--fg:#1a1a1a;--mut:#8a8a8a;--line:#555;--pan:#f6f6f6;--bd:#e0e0e0;--hi:#ffe98a;--sub:#b03030;--nd:#7a7a7a}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#15161a;--fg:#e8e8e8;--mut:#8a8a8a;--line:#c0c0c0;--pan:#1e1f24;--bd:#33343a;--hi:#6a5600;--sub:#ff8a80;--nd:#9a9a9a}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:13px/1.4 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif}
header{position:sticky;top:0;z-index:5;background:var(--pan);border-bottom:1px solid var(--bd);padding:8px 16px;display:flex;flex-wrap:wrap;gap:10px;align-items:center}
h1{font-size:14px;margin:0 10px 0 0;font-weight:600}
input[type=search]{padding:5px 9px;border:1px solid var(--bd);border-radius:6px;background:var(--bg);color:var(--fg);min-width:180px;font-size:13px}
button,select{padding:5px 10px;border:1px solid var(--bd);border-radius:6px;background:var(--bg);color:var(--fg);cursor:pointer;font-size:12px}
button:hover{background:var(--bd)}
button.on{background:var(--fg);color:var(--bg);border-color:var(--fg)}
/* Floating legend: position:fixed so it stays put while a 268-tip tree scrolls past. */
#legend{position:fixed;right:14px;bottom:14px;z-index:20;background:var(--pan);border:1px solid var(--bd);
 border-radius:8px;padding:8px 10px;font-size:11px;box-shadow:0 2px 10px rgba(0,0,0,.18);max-height:70vh;overflow:auto}
#legend.hide{display:none}
#legend h4{margin:0 0 6px;font-size:11px;font-weight:600;color:var(--mut);text-transform:uppercase;letter-spacing:.04em}
#legend .row{display:flex;align-items:center;gap:7px;padding:1px 0;white-space:nowrap}
#legend .sw{width:11px;height:11px;border-radius:2px;border:1px solid rgba(128,128,128,.55);flex:none}
#legend .n{color:var(--mut);margin-left:auto;padding-left:10px;font-variant-numeric:tabular-nums}
#legend .ft{margin-top:6px;padding-top:5px;border-top:1px solid var(--bd);color:var(--mut);font-size:10px;max-width:190px;white-space:normal}
#legToggle{position:fixed;right:14px;bottom:14px;z-index:21}
label{font-size:12px;color:var(--mut);display:flex;align-items:center;gap:4px;cursor:pointer}
#stat{margin-left:auto;color:var(--mut);font-size:11.5px}
#scroll{overflow:auto;height:calc(100vh - 47px)}
.branch{stroke:var(--line);stroke-width:1.1;fill:none}
.node{fill:var(--bg);stroke:var(--nd);stroke-width:1.2;cursor:pointer}
.node:hover{stroke:var(--fg);stroke-width:2}
.node.col{fill:var(--nd)}
.nodehit{fill:transparent;cursor:pointer}
.tri{fill:var(--nd);opacity:.45;cursor:pointer}
.tri:hover{opacity:.75}
text{font-size:11px;fill:var(--fg)}
text.sp{fill:var(--mut);font-size:10px;font-style:italic}
text.gene{font-style:italic}
text.gene.up{font-style:normal}
text.sup{fill:var(--mut);font-size:8.5px}
text.cp{fill:var(--mut);font-size:10px}
text.p3{font-family:ui-monospace,Menlo,monospace;font-size:10px}
text.p3.bad{fill:var(--sub);font-weight:700}
text.sub{fill:var(--sub);font-weight:700}
text.cladelab{fill:var(--mut);font-size:10.5px}
rect.hl{fill:var(--hi);opacity:.55}
</style></head><body>
<header>
<h1>__TITLE__</h1>
<input type="search" id="q" placeholder="search gene / species">
<button id="collapseAll">collapse to depth 3</button><button id="expandAll">expand all</button>
<button id="onlyCel">C. elegans only</button><button id="onlyTest">test set</button>
<label><input type="checkbox" id="showSup" checked> support</label>
<label><input type="checkbox" id="useBL"> branch lengths</label>
<label>colour <select id="colorBy"><option value="lig">ligand</option><option value="cp">G&#945; coupling</option><option value="none">none</option></select></label>
<span id="stat"></span>
</header><div id="scroll"><svg id="svg" xmlns="http://www.w3.org/2000/svg"></svg></div>
<div id="legend"></div>
<button id="legToggle" style="display:none" onclick="document.getElementById('legend').classList.remove('hide');this.style.display='none'">legend</button>
<script>
const TREE=__DATA__;
const LC={'OA':'#C15C3E','TA':'#E0A458','TA/OA':'#D4813F','NE':'#8E3E29','DA':'#1A2C7A','5-HT':'#3C8DAD','5-HT?':'#9EC5D4','ACh':'#7A7A7A','MEL':'#6B4E8E','HA':'#4F8A5B','':'none'};
const CC={'Gs':'#C15C3E','Gq':'#1A2C7A','Gi/o':'#3C8DAD','Go':'#3C8DAD','':'none'};
const NS='http://www.w3.org/2000/svg';
const collapsed=new Set(); let onlyCel=false,onlyTest=false;
const ROW=15, XSTEP=17, LPAD=14;
function tips(n){return n.kids.length?n.kids.flatMap(tips):[n];}
function passes(n){ if(onlyTest) return n.kind==='test'; if(onlyCel) return n.kind==='cel'||n.kind==='test'; return true; }
function el(t,a){const e=document.createElementNS(NS,t);for(const k in a)e.setAttribute(k,a[k]);return e;}
function render(){
  const svg=document.getElementById('svg'); while(svg.firstChild)svg.removeChild(svg.firstChild);
  const showSup=document.getElementById('showSup').checked;
  const useBL=document.getElementById('useBL').checked;
  const cb=document.getElementById('colorBy').value;
  const q=document.getElementById('q').value.trim().toLowerCase();
  let y=0, maxx=0;
  // Branch-length mode needs an adaptive scale: a fixed px-per-substitution factor throws labels
  // off-screen on a deep tree. Measure the largest root-to-tip distance first and fit it to BLW px.
  const BLW=560;
  let maxlen=0;
  (function ml(n,x){ const d=x+(n.bl||0); if(!n.kids.length||collapsed.has(n.id)){ if(d>maxlen)maxlen=d; return;} n.kids.forEach(k=>ml(k,d)); })(TREE,0);
  const BLS = maxlen>0 ? BLW/maxlen : 1;
  // assign coordinates; a node is a leaf of the drawing if collapsed or terminal
  function lay(n,depth,xlen){
    const drawLeaf = !n.kids.length || collapsed.has(n.id);
    n._x = useBL ? LPAD + xlen*BLS : LPAD + depth*XSTEP;
    if(n._x>maxx) maxx=n._x;
    if(drawLeaf){
      if(!n.kids.length && !passes(n)){ n._y=null; return false; }
      if(n.kids.length && !tips(n).some(passes)){ n._y=null; return false; }
      n._y = (y++)*ROW + 12; n._leafdraw=true; return true;
    }
    n._leafdraw=false;
    const ok=n.kids.map(k=>lay(k,depth+1,xlen+(k.bl||0))).filter(Boolean);
    const ys=n.kids.filter(k=>k._y!=null).map(k=>k._y);
    if(!ys.length){ n._y=null; return false; }
    n._y=(Math.min(...ys)+Math.max(...ys))/2; return true;
  }
  lay(TREE,0,0);
  const H=y*ROW+(useBL?46:30), LABX=Math.max(maxx+26, 210);
  const W=LABX+430;
  svg.setAttribute('width',W); svg.setAttribute('height',H);
  const g=el('g',{});
  // branches
  (function draw(n){
    if(n._y==null) return;
    if(!n._leafdraw){
      const kids=n.kids.filter(k=>k._y!=null);
      if(kids.length){
        const y0=Math.min(...kids.map(k=>k._y)), y1=Math.max(...kids.map(k=>k._y));
        g.appendChild(el('path',{d:`M${n._x} ${y0} L${n._x} ${y1}`,class:'branch'}));
        kids.forEach(k=>{ g.appendChild(el('path',{d:`M${n._x} ${k._y} L${k._x} ${k._y}`,class:'branch'})); draw(k); });
      }
    }
  })(TREE);
  // leaves + collapsed clades + node handles
  let shown=0;
  (function paint(n){
    if(n._y==null) return;
    if(n._leafdraw && n.kids.length){           // collapsed clade
      const lv=tips(n);
      const tri=el('polygon',{points:`${n._x},${n._y} ${n._x+20},${n._y-5} ${n._x+20},${n._y+5}`,class:'tri'});
      tri.onclick=()=>{collapsed.delete(n.id);render();};
      g.appendChild(tri);
      const cels=lv.filter(x=>x.kind==='cel'||x.kind==='test').map(x=>x.g);
      const t=el('text',{x:n._x+25,y:n._y+3.5,class:'cladelab'});
      t.textContent=lv.length+' tips'+(cels.length?'  ·  '+cels.slice(0,6).join(', ')+(cels.length>6?' …':''):'');
      g.appendChild(t); shown+=lv.length;
      if(showSup&&n.sup){const s=el('text',{x:n._x-3,y:n._y-4,class:'sup','text-anchor':'end'});s.textContent=n.sup;g.appendChild(s);}
      return;
    }
    if(!n.kids.length){                          // real tip
      const hit=q&&(n.g+' '+n.sp+' '+n.name).toLowerCase().includes(q);
      if(hit) g.appendChild(el('rect',{x:LABX-6,y:n._y-7,width:400,height:14,class:'hl'}));
      g.appendChild(el('path',{d:`M${n._x} ${n._y} L${LABX-10} ${n._y}`,class:'branch','stroke-dasharray':'1.5,2.5',opacity:.45}));
      const col = cb==='lig'?(LC[n.lig]||'none') : cb==='cp'?(CC[n.cp]||'none') : 'none';
      g.appendChild(el('rect',{x:LABX-4,y:n._y-4.5,width:8,height:9,fill:col,stroke:'var(--mut)','stroke-width':.6,rx:1.5}));
      if(n.cp){const c=el('text',{x:LABX+10,y:n._y+3.5,class:'cp'});c.textContent=n.cp;g.appendChild(c);}
      if(n.p332){const p=el('text',{x:LABX+56,y:n._y+3.5,class:'p3'+(n.p332!=='D'?' bad':'')});p.textContent=n.p332;g.appendChild(p);}
      const gt=el('text',{x:LABX+72,y:n._y+3.5,class:'gene'+(n.kind==='scan'?' up':'')+(n.name==='cel_dop-6'?' sub':'')});
      gt.textContent=(n.kind==='test'?'TEST ':'')+n.g; g.appendChild(gt);
      const sp=el('text',{x:LABX+218,y:n._y+3.5,class:'sp'}); sp.textContent=n.sp; g.appendChild(sp);
      shown++; return;
    }
    // internal node: clickable handle sitting ON the node
    const c=el('circle',{cx:n._x,cy:n._y,r:3.6,class:'node'});
    c.onclick=()=>{collapsed.add(n.id);render();};
    const lv=tips(n);
    c.appendChild(el('title',{})).textContent=lv.length+' tips — click to collapse'+(n.sup?'  ('+n.sup+')':'');
    g.appendChild(c);
    if(showSup&&n.sup){const s=el('text',{x:n._x-5,y:n._y-5,class:'sup','text-anchor':'end'});s.textContent=n.sup;g.appendChild(s);}
    n.kids.forEach(paint);
  })(TREE);
  if(useBL && maxlen>0){
    const nice=[0.05,0.1,0.2,0.5,1,2,5].find(v=>v*BLS>60)||1;
    const y0=H-12, x0=LPAD;
    g.appendChild(el('path',{d:`M${x0} ${y0} L${x0+nice*BLS} ${y0}`,class:'branch','stroke-width':1.4}));
    g.appendChild(el('path',{d:`M${x0} ${y0-4} L${x0} ${y0+4}`,class:'branch'}));
    g.appendChild(el('path',{d:`M${x0+nice*BLS} ${y0-4} L${x0+nice*BLS} ${y0+4}`,class:'branch'}));
    const t=el('text',{x:x0+nice*BLS/2,y:y0-7,class:'sup','text-anchor':'middle'});
    t.textContent=nice+' substitutions/site'; g.appendChild(t);
  }
  svg.appendChild(g);
  document.getElementById('stat').textContent=shown+' of '+tips(TREE).length+' tips'+(collapsed.size?'  ·  '+collapsed.size+' collapsed':'');
}
document.getElementById('collapseAll').onclick=()=>{collapsed.clear();
  (function w(n,d){ if(!n.kids.length) return; if(d>=3){collapsed.add(n.id);return;} n.kids.forEach(k=>w(k,d+1)); })(TREE,0); render();};
document.getElementById('expandAll').onclick=()=>{collapsed.clear();render();};
document.getElementById('onlyCel').onclick=e=>{onlyCel=!onlyCel;e.target.classList.toggle('on',onlyCel);render();};
document.getElementById('onlyTest').onclick=e=>{onlyTest=!onlyTest;e.target.classList.toggle('on',onlyTest);render();};
function drawLegend(){
  // Only list values that actually occur in this tree, with counts - a legend showing classes the
  // panel does not contain is worse than none.
  const cb=document.getElementById('colorBy').value;
  const box=document.getElementById('legend'), btn=document.getElementById('legToggle');
  if(cb==='none'){ box.style.display='none'; btn.style.display='none'; return; }
  box.style.display='';
  const key = cb==='lig' ? 'lig' : 'cp';
  const COL = cb==='lig' ? LC : CC;
  const cnt={};
  (function walk(n){ if(!n.kids.length){ const v=n[key]; if(v) cnt[v]=(cnt[v]||0)+1; } n.kids.forEach(walk); })(TREE);
  const order=Object.keys(COL).filter(k=>k&&cnt[k]).sort((a,b)=>cnt[b]-cnt[a]);
  const unl=(function(){let u=0;(function w(n){if(!n.kids.length&&!n[key])u++;n.kids.forEach(w);})(TREE);return u;})();
  let h='<h4>'+(cb==='lig'?'ligand':'G\u03B1 coupling')+'</h4>';
  for(const k of order) h+='<div class="row"><span class="sw" style="background:'+COL[k]+'"></span>'+k+'<span class="n">'+cnt[k]+'</span></div>';
  if(unl) h+='<div class="row"><span class="sw" style="background:transparent"></span>no experimental assignment<span class="n">'+unl+'</span></div>';
  h+='<div class="ft">Colour marks an experimentally determined value only. Uncoloured tips are breadth, not evidence. <a href="#" id="legHide">hide</a></div>';
  box.innerHTML=h;
  // Attach the handler rather than inlining onclick: an inline handler here would need
  // nested quotes inside a quoted JS string, which is how this shipped broken once.
  const hb=document.getElementById('legHide');
  if(hb) hb.addEventListener('click',function(e){e.preventDefault();
    box.classList.add('hide'); btn.style.display='';});
}
const _r=render; render=function(){ _r(); drawLegend(); };
for(const id of ['showSup','useBL','colorBy','q']) document.getElementById(id)[id==='q'?'oninput':'onchange']=render;
render();
</script></body></html>"""

html = HTML.replace('__DATA__',data).replace('__TITLE__',title)
# Syntax-check the embedded script before writing. A JS error anywhere in this block leaves the
# page completely blank with no indication why, which is exactly how the legend patch shipped
# broken: an inline onclick needed nested quotes inside a quoted string and silently terminated it.
def _check(html):
    import re as _re, shutil as _sh, subprocess as _sp, tempfile as _tf
    node=_sh.which('node')
    if not node: return 'skipped (node not found)'
    m=_re.search(r'<script>(.*?)</script>', html, _re.S)
    if not m: return 'no script block found'
    js=m.group(1)
    with _tf.NamedTemporaryFile('w',suffix='.js',delete=False) as f:
        f.write('(function(){'+js+'})();'); path=f.name
    r=_sp.run([node,'--check',path],capture_output=True,text=True)
    return 'OK' if r.returncode==0 else 'SYNTAX ERROR\n'+(r.stderr or r.stdout)[:600]

_chk=_check(html)
if _chk.startswith('SYNTAX ERROR'):
    sys.stderr.write(f'REFUSING TO WRITE {out}: embedded JS is broken\n{_chk}\n'); sys.exit(1)
open(out,'w').write(html)
print(f'-> {out}  ({tipcount(tree)} tips, {os.path.getsize(out)//1024} KB)')
