#!/usr/bin/env python3
"""assembly_copy_number.py — is the dop-6 palindrome in independent long-read assemblies?

The N2 reference (WBcel235) descends from the Sanger cosmid assembly, where a 3.9 kb perfect
inverted repeat is exactly the structure most likely to be mis-assembled, and which is also
unstable in E. coli. This script asks two questions of independent long-read assemblies —
another N2, and wild isolates:

  Q1 (artifact control): does an independent long-read N2 assembly contain the palindrome?
  Q2 (wild isolates):    is the palindrome present outside N2, i.e. how old is it?

Test: BLAST two probes against each assembly.
  probe_arm     500 bp from inside arm A (dop-6 exon 15 region) -> expect 2 hits if palindrome present
  probe_unique  500 bp from dop-6 exon 4 (single copy)          -> expect 1 hit (a control for ploidy/
                                                                   assembly duplication artefacts)
  junctions     4 x 300 bp spanning each arm/flank boundary     -> present only if the arrangement matches
Writes output/assembly_copy_number.txt
"""
import gzip, os, shutil, subprocess, sys, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'output')
WORK = ('/private/tmp/claude-501/-Users-mikeodonnell-git-Malaiwong-Oglu-2026/'
        '93015e40-e331-495a-8362-20ac7ea07c75/scratchpad/asm')
os.makedirs(WORK, exist_ok=True)
L0 = 4318000
ARM_A = (4324178, 4328106); ARM_B = (4329303, 4333231)
seq = ''.join(l.strip() for l in open(f'{OUT}/locus_X_4318000-4334500.fa') if not l.startswith('>'))
G = lambda a, b: seq[a-L0:b-L0+1].upper()

ASSEMBLIES = [
    ('GCA_028201515.1', 'ASM2820151v1', 'N2 (Bristol) long-read, U. Oregon'),
    ('GCA_028201415.1', 'ASM2820141v1', 'CB4856 (Hawaiian) long-read, U. Oregon'),
    ('GCA_052911875.1', 'C.elegans_CGC1_1.0', 'CGC1 complete'),
    ('GCA_016989125.1', 'XZ1516_Canu', 'XZ1516 (divergent isotype)'),
    ('GCA_016989505.1', 'DL238_Canu', 'DL238 (divergent isotype)'),
    ('GCA_016989235.1', 'MY2147_Canu', 'MY2147'),
    ('GCA_016989295.1', 'EG4725_Canu', 'EG4725'),
    ('GCA_016989245.1', 'JU2600_Canu', 'JU2600'),
    ('GCA_016989455.1', 'ECA36_Canu', 'ECA36'),
    ('GCA_016989115.1', 'NIC526_Canu', 'NIC526'),
    ('GCA_016989365.1', 'JU1400_Canu', 'JU1400'),
    ('GCA_016989285.1', 'JU2526_Canu', 'JU2526'),
    ('GCA_016989105.1', 'MY2693_Canu', 'MY2693'),
    ('GCA_016989275.1', 'JU310_Canu', 'JU310'),
    ('GCA_016989385.1', 'ECA396_Canu', 'ECA396'),
    ('GCA_016989095.1', 'QX1794_Canu', 'QX1794'),
    ('GCA_016989145.1', 'NIC2_Canu', 'NIC2'),
    ('GCA_020450165.1', 'ASM2045016v1', 'CB4856 (Seoul, independent assembly)'),
    ('GCA_022984815.1', 'Ce_DL226', 'DL226 (Seoul)'),
    ('GCA_054661295.1', 'C.elegans_KR314_1.0', 'KR314'),
    ('GCA_054661375.1', 'C.elegans_CB4855_1.0', 'CB4855'),
    ('GCA_045866515.1', 'C. elegans_AB1_1.0', 'AB1'),
    ('GCA_054661275.1', 'C.elegans_RC301_1.0', 'RC301'),
]
PROBES = {
    'probe_unique_ex4': G(4321367, 4321664),                 # dop-6 exon 4, outside the palindrome
    'probe_arm_ex15':   G(4326730, 4326834) + G(4326835, 4327134),  # inside arm A
    'junc_armA_start':  G(ARM_A[0]-150, ARM_A[0]+149),
    'junc_armA_end':    G(ARM_A[1]-149, ARM_A[1]+150),
    'junc_armB_start':  G(ARM_B[0]-150, ARM_B[0]+149),
    'junc_armB_end':    G(ARM_B[1]-149, ARM_B[1]+150),
}
pf = f'{WORK}/probes.fa'
with open(pf, 'w') as f:
    for k, v in PROBES.items(): f.write(f'>{k}\n{v}\n')

rep = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); rep.append(s)

say('=' * 100)
say('IS THE dop-6 PALINDROME PRESENT IN INDEPENDENT LONG-READ ASSEMBLIES?')
say('=' * 100)
say('probe_unique_ex4 : dop-6 exon 4, single copy in WBcel235       -> expect 1 hit')
say('probe_arm_ex15   : 405 bp inside arm A, two copies in WBcel235 -> expect 2 hits if palindrome present')
say('junc_*           : 300 bp spanning each repeat boundary        -> expect 1 hit each if arrangement matches')
say('')

for acc, name, desc in ASSEMBLIES:
    tag = acc.replace('.', '_')
    fna = f'{WORK}/{tag}.fna'
    if not os.path.exists(fna):
        p = acc[4:7] + '/' + acc[7:10] + '/' + acc[10:13]
        url = f'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/{p}/{acc}_{name}/{acc}_{name}_genomic.fna.gz'
        try:
            say(f'downloading {acc} ({desc}) ...')
            urllib.request.urlretrieve(url, fna + '.gz')
            with gzip.open(fna + '.gz', 'rb') as i, open(fna, 'wb') as o: shutil.copyfileobj(i, o)
            os.remove(fna + '.gz')
        except Exception as e:
            say(f'  {acc}: download failed ({e})'); continue
    if not os.path.exists(f'{WORK}/{tag}.nsq'):
        subprocess.run(['makeblastdb', '-in', fna, '-dbtype', 'nucl', '-out', f'{WORK}/{tag}'],
                       capture_output=True, check=True)
    r = subprocess.run(['blastn', '-query', pf, '-db', f'{WORK}/{tag}', '-dust', 'no', '-evalue', '1e-20',
                        '-outfmt', '6 qseqid sseqid sstart send pident length'],
                       capture_output=True, text=True)
    hits = [l.split('\t') for l in r.stdout.strip().split('\n') if l]
    say(f'\n{acc}  {desc}')
    size = sum(1 for _ in open(fna)) and os.path.getsize(fna)//10**6
    say(f'  assembly {size} Mb')
    for probe in PROBES:
        hs = [h for h in hits if h[0] == probe and int(h[5]) >= 0.8*len(PROBES[probe])]
        strong = [h for h in hs if float(h[4]) >= 97]
        ids = ', '.join(f'{h[1].split()[0]}:{h[2]}-{h[3]} ({float(h[4]):.1f}%)' for h in strong[:4])
        say(f'  {probe:<18} {len(strong)} hit(s) >=97% over >=80% of probe   {ids}')

say('\nReading: probe_arm_ex15 = 2 means the duplication is present in that assembly; = 1 means that')
say('assembly is single-copy there. The junction probes say whether the arrangement is the same one.')
open(f'{OUT}/assembly_copy_number.txt', 'w').write('\n'.join(rep) + '\n')
say('\n-> output/assembly_copy_number.txt')
