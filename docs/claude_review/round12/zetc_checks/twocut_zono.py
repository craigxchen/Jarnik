# certify zeta = rho for two-cut blocks (L>=1): rho_p (mod-p, an upper bound for rho_Q) vs verified zonotope
import re, time
from zono import zeta_lower
from tclib import box_points
tot = ok = 0
for line in open('twocut_rho_out.txt'):
    m = re.match(r'a=(\d+) b=(\d+) c=(\d+) L1=(\d+) L2=(\d+) \|A\|=(\d+) rho=(\d+)', line)
    if not m: continue
    a,b,c,L1,L2,N,rho = map(int, m.groups())
    if L1 == 0 or L2 == 0: continue
    d = a+b+c; S1 = list(range(a+b)); S2 = list(range(a, d))
    pts = box_points(d, [(S1,L1),(S2,L2)])
    z, edges = zeta_lower(pts, [S1,S2], [L1,L2], d)
    tot += 1; ok += (z == rho)
    print('a=%d b=%d c=%d L=(%d,%d) rho_p=%d zeta_verified=%d %s' % (a,b,c,L1,L2,rho,z,'' if z==rho else 'GAP'), flush=True)
print('blocks', tot, 'certified zeta=rho:', ok)
