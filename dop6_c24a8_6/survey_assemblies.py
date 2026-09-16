#!/usr/bin/env python3
"""survey_assemblies.py — the dop-6 palindrome across every available C. elegans assembly,
plus a mapping-free polymorphism dataset for the gene.

Two jobs, both done from assemblies rather than from reads mapped to WBcel235:

  A. COPY NUMBER. For each assembly: locate dop-6 exon 4 (single copy, outside the repeat),
     take a window in the gene's orientation, count copies of the 3,929 bp repeat unit, and — for
     carriers — measure the identity between the two arms. Arm-vs-arm identity is the quantity
     that separates "one recent origin" from "old duplication homogenised by gene conversion".

  B. POLYMORPHISM WITHOUT MAPPING. CaeNDR variant calls vanish inside the repeat because reads
     map to WBcel235 at MAPQ 0. The sequence is not missing, only unplaceable against a duplicated
     reference. So rebuild the dop-6 CDS exon by exon out of each assembly directly and compute
     real diversity per segment from that alignment. Carriers contribute their arm-A copy, which
     is the orthologue of the single copy in non-carriers.

Writes output/assembly_survey.tsv, output/assembly_survey.txt, output/dop6_cds_isolates.fna
"""
import gzip, json, os, re, shutil, subprocess, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'output')
WORK = ('/private/tmp/claude-501/-Users-mikeodonnell-git-Malaiwong-Oglu-2026/'
        '93015e40-e331-495a-8362-20ac7ea07c75/scratchpad/asm')
os.makedirs(WORK, exist_ok=True)
L0 = 4318000; ARM = (4324178, 4328106); SPLIT = 405
n2 = ''.join(l.strip() for l in open(f'{OUT}/locus_X_4318000-4334500.fa') if not l.startswith('>')).upper()
G = lambda a, b: n2[a-L0:b-L0+1]
RC = lambda s: s[::-1].translate(str.maketrans('ACGTacgt', 'TGCAtgca'))
models = json.load(open(f'{OUT}/gene_models.json'))
EXONS = models['dop-6']['transcripts'][0]['exons']
REPEAT = G(*ARM); EX4 = G(4321367, 4321664)
UP, DOWN = 6000, 16000

# every C. elegans assembly at NCBI, annotated with what kind of strain it is.
# lab = N2-derived laboratory strain; wild = natural isolate; ref = the reference itself.
PANEL = [
 ('GCA_000002985.3','WBcel235','N2 reference (Sanger cosmids)','ref'),
 ('GCA_028201515.1','ASM2820151v1','N2 (U. Oregon)','lab'),
 ('GCA_963921025.1','N2','N2 (Wellcome Sanger)','lab'),
 ('GCA_052911875.1','C.elegans_CGC1_1.0','CGC1','lab'),
 ('GCA_039880965.1','ASM3988096v1','UA44','lab'),
 ('GCA_029748435.1','ASM2974843v1','BY250','lab'),
 ('GCA_054508145.1','ASM5450814v1','yA','lab'),
 ('GCA_056506345.1','ASM5650634v1','yA1244','lab'),
 ('GCA_028201415.1','ASM2820141v1','CB4856 (U. Oregon)','wild'),
 ('GCA_020450165.1','ASM2045016v1','CB4856 (Seoul)','wild'),
 ('GCA_963921145.1','CB4856','CB4856 (Sanger)','wild'),
 ('GCA_022984815.1','Ce_DL226','DL226','wild'),
 ('GCA_016989125.1','XZ1516_Canu','XZ1516','wild'),
 ('GCA_016989505.1','DL238_Canu','DL238','wild'),
 ('GCA_963920975.1','DL238','DL238 (Sanger)','wild'),
 ('GCA_016989235.1','MY2147_Canu','MY2147','wild'),
 ('GCA_016989295.1','EG4725_Canu','EG4725','wild'),
 ('GCA_963921005.1','EG4725','EG4725 (Sanger)','wild'),
 ('GCA_016989245.1','JU2600_Canu','JU2600','wild'),
 ('GCA_016989455.1','ECA36_Canu','ECA36','wild'),
 ('GCA_963921125.1','ECA36','ECA36 (Sanger)','wild'),
 ('GCA_016989115.1','NIC526_Canu','NIC526','wild'),
 ('GCA_016989365.1','JU1400_Canu','JU1400','wild'),
 ('GCA_016989285.1','JU2526_Canu','JU2526','wild'),
 ('GCA_016989105.1','MY2693_Canu','MY2693','wild'),
 ('GCA_016989275.1','JU310_Canu','JU310','wild'),
 ('GCA_016989385.1','ECA396_Canu','ECA396','wild'),
 ('GCA_016989095.1','QX1794_Canu','QX1794','wild'),
 ('GCA_016989145.1','NIC2_Canu','NIC2','wild'),
 ('GCA_054661295.1','C.elegans_KR314_1.0','KR314','wild'),
 ('GCA_054661375.1','C.elegans_CB4855_1.0','CB4855','wild'),
 ('GCA_054661335.1','C.elegans_CB4858_1.0','CB4858','wild'),
 ('GCA_054661275.1','C.elegans_RC301_1.0','RC301','wild'),
 ('GCA_045866515.1','C. elegans_AB1_1.0','AB1','wild'),
 ('GCA_054661435.1','C.elegans_AB3_1.0','AB3','wild'),
 ('GCA_054661355.1','C.elegans_CB4857_1.0','CB4857','wild'),
 ('GCA_054661395.1','C.elegans_CB4853_1.0','CB4853','wild'),
 ('GCA_054661315.1','C.elegans_CB4932_1.0','CB4932','wild'),
 ('GCA_054661415.1','C.elegans_CB4852_1.0','CB4852','wild'),
 ('GCA_963921135.1','QX1791','QX1791','wild'),
 ('GCA_963921055.1','GXW1','GXW1','wild'),
 ('GCA_963920995.1','DL200','DL200','wild'),
 ('GCA_963921035.1','PB306','PB306','wild'),
 ('GCA_963921075.1','JU775','JU775','wild'),
 ('GCA_963921105.1','AB1','AB1 (Sanger)','wild'),
 ('GCA_963921115.1','LKC34','LKC34','wild'),
 ('GCA_963921155.1','QX1211','QX1211','wild'),
 ('GCA_963921015.1','PS2025','PS2025','wild'),
 ('GCA_963921045.1','JT11398','JT11398','wild'),
 ('GCA_963921065.1','JU751','JU751','wild'),
 ('GCA_963921085.1','ECA248','ECA248','wild'),
 ('GCA_963920985.1','EG4946','EG4946','wild'),
 ('GCA_963920965.1','JU394','JU394','wild'),
 ('GCA_963921095.1','MY23','MY23','wild'),
 ('GCA_013403715.1','ASM1340371v1','JU1395','wild'),
 ('GCA_027422445.1','ASM2742244v1','ALT1','wild'),
 ('GCA_027422425.1','ASM2742242v1','ALT2','wild'),
]
rep = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); rep.append(s)

def ensure(acc, name):
    tag = acc.replace('.', '_'); fna = f'{WORK}/{tag}.fna'
    if not os.path.exists(fna):
        p = f'{acc[4:7]}/{acc[7:10]}/{acc[10:13]}'
        d = urllib.parse.quote(f'{acc}_{name}')
        url = f'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/{p}/{d}/{d}_genomic.fna.gz'
        try:
            urllib.request.urlretrieve(url, fna + '.gz')
            with gzip.open(fna + '.gz', 'rb') as i, open(fna, 'wb') as o: shutil.copyfileobj(i, o)
            os.remove(fna + '.gz')
        except Exception as e:
            return None, f'download failed: {str(e)[:70]}'
    if not os.path.exists(f'{WORK}/{tag}.nsq'):
        subprocess.run(['makeblastdb', '-in', fna, '-dbtype', 'nucl', '-out', f'{WORK}/{tag}'],
                       capture_output=True)
    return tag, None

def blast(q, db, fmt='6 qseqid sstart send pident length sstrand sseq', ev='1e-10'):
    r = subprocess.run(['blastn', '-query', q, '-db', db, '-dust', 'no', '-evalue', ev,
                        '-outfmt', fmt], capture_output=True, text=True)
    return [h.split('\t') for h in r.stdout.strip().split('\n') if h]

open(f'{WORK}/_ex4.fa', 'w').write(f'>ex4\n{EX4}\n')
open(f'{WORK}/_rep.fa', 'w').write(f'>repeat\n{REPEAT}\n')
rows = []; cdss = []
say('=' * 104)
say('dop-6 PALINDROME ACROSS ALL AVAILABLE C. ELEGANS ASSEMBLIES')
say('=' * 104)
say(f'{"strain":<22} {"kind":<5} {"accession":<18} {"copies":>6} {"arm-arm identity":>18}  note')
say('-' * 104)
for acc, name, strain, kind in PANEL:
    tag, err = ensure(acc, name)
    if err: say(f'{strain:<22} {kind:<5} {acc:<18} {"-":>6} {"-":>18}  {err}'); continue
    hs = [h for h in blast(f'{WORK}/_ex4.fa', f'{WORK}/{tag}', ev='1e-30') if int(h[4]) > 250]
    if len(hs) != 1:
        say(f'{strain:<22} {kind:<5} {acc:<18} {"-":>6} {"-":>18}  exon 4 found {len(hs)}x'); continue
    contig, s0, e0 = hs[0][0], int(hs[0][1]), int(hs[0][2])
    # blast qseqid is the query name; get subject id separately
    r = subprocess.run(['blastn', '-query', f'{WORK}/_ex4.fa', '-db', f'{WORK}/{tag}', '-dust', 'no',
                        '-evalue', '1e-30', '-outfmt', '6 sseqid sstart send length'],
                       capture_output=True, text=True)
    hh = [x.split('\t') for x in r.stdout.strip().split('\n') if x and int(x.split('\t')[3]) > 250]
    contig, s0, e0 = hh[0][0], int(hh[0][1]), int(hh[0][2])
    fwd = e0 > s0
    lo = max(1, (s0 - UP) if fwd else (s0 - DOWN)); hi = (e0 + DOWN) if fwd else (e0 + UP)
    rr = subprocess.run(['seqkit', 'subseq', '--chr', contig, '-r', f'{lo}:{hi}', f'{WORK}/{tag}.fna'],
                        capture_output=True, text=True)
    win = ''.join(l for l in rr.stdout.split('\n')[1:] if not l.startswith('>')).upper()
    if not win: say(f'{strain:<22} {kind:<5} {acc:<18} extraction failed'); continue
    if not fwd: win = RC(win)
    wf = f'{WORK}/{tag}_win.fa'; open(wf, 'w').write(f'>{tag}\n{win}\n')
    subprocess.run(['makeblastdb', '-in', wf, '-dbtype', 'nucl', '-out', f'{WORK}/{tag}_win'],
                   capture_output=True)
    cps = [h for h in blast(f'{WORK}/_rep.fa', f'{WORK}/{tag}_win') if int(h[4]) > 1500]
    ncp = len(cps); armid = ''
    if ncp == 2:
        sf = subprocess.run(['blastn', '-query', wf, '-db', f'{WORK}/{tag}_win', '-dust', 'no',
                             '-evalue', '1e-20', '-outfmt', '6 pident length sstrand'],
                            capture_output=True, text=True)
        inv = [x.split('\t') for x in sf.stdout.strip().split('\n')
               if x and x.split('\t')[2] == 'minus' and int(x.split('\t')[1]) > 1500]
        if inv:
            b = max(inv, key=lambda x: int(x[1])); armid = f'{b[0]} % / {b[1]} bp'
    # rebuild the CDS (arm-A copy for carriers: the first repeat copy in gene orientation)
    cds, bad = [], None
    for j, (es, ee) in enumerate(EXONS, 1):
        open(f'{WORK}/_e.fa', 'w').write(f'>e{j}\n{G(es,ee)}\n')
        eh = [h for h in blast(f'{WORK}/_e.fa', f'{WORK}/{tag}_win', ev='1e-5')
              if int(h[4]) == ee-es+1 and h[5] == 'plus']
        if not eh: bad = f'exon {j} not recovered'; break
        eh.sort(key=lambda h: int(h[1]))
        cds.append(eh[0][6].replace('-', ''))
    if not bad:
        cdss.append((strain.replace(' ', '_'), ''.join(cds)))
    rows.append((strain, kind, acc, ncp, armid, bad or ''))
    say(f'{strain:<22} {kind:<5} {acc:<18} {ncp:>6} {armid:>18}  {bad or ""}')

with open(f'{OUT}/assembly_survey.tsv', 'w') as f:
    f.write('strain\tkind\taccession\trepeat_copies\tarm_arm_identity\tnote\n')
    for r in rows: f.write('\t'.join(str(x) for x in r) + '\n')
with open(f'{OUT}/dop6_cds_isolates.fna', 'w') as f:
    f.write(f'>N2_WBcel235_reference\n{"".join(G(a,b) for a,b in EXONS)}\n')
    for n, s in cdss: f.write(f'>{n}\n{s}\n')
from collections import Counter
say('')
for kind in ('lab', 'wild', 'ref'):
    c = Counter(r[3] for r in rows if r[1] == kind)
    say(f'{kind:<5}: ' + ', '.join(f'{k} copies x{v}' for k, v in sorted(c.items())))
say(f'\n-> output/assembly_survey.tsv ({len(rows)} assemblies), '
    f'output/dop6_cds_isolates.fna ({len(cdss)+1} CDS)')
open(f'{OUT}/assembly_survey.txt', 'w').write('\n'.join(rep) + '\n')
