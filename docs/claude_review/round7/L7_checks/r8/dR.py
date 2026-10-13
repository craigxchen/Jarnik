# d_R(C) = degree of t -> R_{x(t)} in Plucker P^5 along Mobius test curves; d_R(beta) = sum_T d_R(mu_T) - 32 d_R(mu)
import random, itertools
from lamA_direct import kernel_poly, plucker, good_mob
from star import LAB
random.seed(9)
def dR(T):
    vals = []
    for rep in range(2):
        K = kernel_poly(good_mob(T)); p = plucker(K[:4], K[4:])
        vals.append(max(c.degree() for c in p))
    return vals
d0 = dR(()); print('mu', d0)
tot = 0; per = {}
for T in itertools.combinations(range(7), 3):
    v = dR(T); per[T] = v[0]; tot += v[0]
    print([LAB[a] for a in T], v, flush=True)
print('d_R(beta) =', tot - 32*d0[0])
# also forgetful fibres
from elamp import md, rnd
for j in range(7):
    mob = [[md(rnd()), md(rnd()), 0, 1] if a == j else [0, md(rnd()), 0, 1] for a in range(7)]
    K = kernel_poly(mob); p = plucker(K[:4], K[4:])
    print('f_' + LAB[j], 'd_R =', max(c.degree() for c in p))
