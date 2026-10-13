import random, itertools
from lamA_direct import *
random.seed(3)
mob = good_mob(())
pR, secs = lam_data(mob)
def pair(p, q):
    return p[0]*q[5] - p[1]*q[4] + p[2]*q[3] + p[3]*q[2] - p[4]*q[1] + p[5]*q[0]
pp = pair(secs[0], secs[2])
roots = [-int(f[0]) % P for f, m in pp.factor()[1] if f.degree() == 1]
# recompute Q_t at root to test singularity
K = kernel_poly(mob); L1, L2 = K[:4], K[4:]
qbp = [[[C(Qk[i][j]) for j in range(4)] for i in range(4)] for Qk in qb]
rows = [[bil(Qk, u, v) for Qk in qbp] for (u, v) in [(L1, L1), (L1, L2), (L2, L2)]]
cs = []
for j in range(4):
    sub = [[r[k] for k in range(4) if k != j] for r in rows]
    d = det3(sub); cs.append(d if j % 2 == 0 else -d)
cs = prim(cs)
for t0 in roots[:6]:
    Qm = [[sum(int(cs[k](t0)) * qb[k][i][j] for k in range(4)) % P for j in range(4)] for i in range(4)]
    print('t0', t0, 'det Q_t0 =', det_mod(Qm), ' pairing(s0,s2)=0; pairing(s0,s1) at t0:', int(pair(secs[0], secs[1])(t0)), ' s0 == s2?', all(int((secs[0][i]*secs[2][j]-secs[0][j]*secs[2][i])(t0)) == 0 for i in range(6) for j in range(6)))
