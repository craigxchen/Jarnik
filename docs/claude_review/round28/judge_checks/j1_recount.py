# Judge check J1 (unrefereed): recount boxes.csv (exact integer bookkeeping)
import csv, collections
D='/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28'
rows=list(csv.DictReader(open(D+'/boxes.csv')))
print('boxes',len(rows))
bad=[r['spec'] for r in rows if r['reach_pair_P1']!=r['reach_pair_P2'] or r['H_agree']!='True']
print('P1/P2 disagreements or H disagreement:',bad)
print('zeta not exact:',[r['spec'] for r in rows if r['zeta_all_exact']!='True' or r['zeta_pair_exact']!='True'])
def kind(r): return 'random' if r['family'].startswith('rand') or r['spec'].startswith('rand') else 'structured'
print('families:',collections.Counter(r['family'] for r in rows))
for zname in ('zeta_all','zeta_pair'):
    g=collections.defaultdict(lambda:[0,0,0,0,0])
    for r in rows:
        rp=int(r['reach_pair_P1']); z=int(r[zname]); rho=int(r['rho'])
        key=(int(r['d']),kind(r))
        g[key][0]+=1; g[key][1]+= rp>z; g[key][2]+= rp==z; g[key][3]+= rp<z; g[key][4]+= z<rho
    print('vs',zname,' (d,kind): n, >, =, <, faithful')
    for k in sorted(g): print('  ',k,g[k])
    print('  reach>%s:'%zname,[r['spec'] for r in rows if int(r['reach_pair_P1'])>int(r[zname])])
# best of pair/cube
gt=[r['spec'] for r in rows if max(int(r['reach_pair_P1']),int(r['reach_cube_P1'] or 0))>int(r['zeta_all'])]
print('best(pair,cube) > zeta_all:',gt)
# Sylvester L=3 class
print('Sylvester L>=3 rows:',[(r['spec'],r['d'],r['reach_pair_P1'],r['zeta_all'],r['rho']) for r in rows if 'syl' in r['spec'] and int(r['L'])==3])
# confounding: rho/L>4
hi=[r for r in rows if float(r['rho_over_L'])>4]
print('rho/L>4: n',len(hi),'random',sum(kind(r)=='random' for r in hi),'zonotopal',sum(int(r['zeta_all'])==int(r['rho']) for r in hi),
      'reach>zeta',sum(int(r['reach_pair_P1'])>int(r['zeta_all']) for r in hi))
print('structured max rho/L:',max(float(r['rho_over_L']) for r in rows if kind(r)=='structured'))
# within structured: reach - zeta by L
s=collections.defaultdict(list)
for r in rows:
    if kind(r)=='structured': s[(int(r['d']),int(r['L']))].append(int(r['reach_pair_P1'])-int(r['zeta_all']))
print('structured reach-zeta by (d,L):',dict(sorted(s.items())))
rnd=[r for r in rows if kind(r)=='random']
print('random: n',len(rnd),'zonotopal',sum(int(r['zeta_all'])==int(r['rho']) for r in rnd),'reach<zeta_all',sum(int(r['reach_pair_P1'])<int(r['zeta_all']) for r in rnd))
