# Exact (over Z) Plucker degree d_R(C) of t -> R_{C(t)} for integer Mobius test curves (cf. dR.py, mod p).
import sys, random, itertools
import lamA_exact as E
from star import LAB
def dR(mob):
    K = E.kernel_poly(mob)
    return max(c.degree() for c in E.content_div(E.wedge(K[:4], K[4:])) if not c.is_zero())
if __name__ == '__main__':
    random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
    d0 = dR(E.int_mob(())); tot = 0; odd = {}
    for T in itertools.combinations(range(7), 3):
        d = dR(E.int_mob(T)); tot += d
        if d != d0: odd[''.join(LAB[a] for a in T)] = d
    print('d_R(mu) =', d0, ' d_R(mu_T) != d_R(mu) for', odd, ' d_R(beta) =', tot - 32 * d0)
