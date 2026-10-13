import random, itertools
from lamA_direct import *
random.seed(3)
pR, secs = lam_data(good_mob(()))
def pair(p, q):  # Plucker pairing for ordering (01,02,03,12,13,23)
    return p[0]*q[5] - p[1]*q[4] + p[2]*q[3] + p[3]*q[2] - p[4]*q[1] + p[5]*q[0]
print('self-pairing zero?', [pair(s, s).is_zero() for s in secs], pair(pR, pR).is_zero())
for a, b in [(0,1),(0,2),(0,3),(1,2),(4,5)]:
    pp = pair(secs[a], secs[b])
    fac = pp.factor()
    print((a,b), 'meet', meet(secs[a], secs[b]), 'pairing deg', pp.degree(), 'factor degs/mults', [(f.degree(), m) for f, m in fac[1]])
print('Sigma pairing with s0 factor:', [(f.degree(), m) for f, m in pair(pR, secs[0]).factor()[1]])
