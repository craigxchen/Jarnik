import random
from cox_form import *
rng = random.Random(5)
primes = [5]
shown = 0
for _ in range(4000):
    g = {A: (1, 0) for A in CUTS}
    ch = ([], [])
    pi = gauss_prime(5)
    for _k in range(rng.randint(0, 3)):
        A = rng.choice(CUTS)
        if rng.random() < 0.5:
            g[A] = gmul(g[A], pi); ch[0].append(A)
        else:
            g[A] = gmul(g[A], gconj(pi)); ch[1].append(A)
    zs = zeta_from_gammas((1, 0), g)
    prim = gnorm(ggcd_list(zs)) == 1
    c = chain_ok(ch)
    if prim != c and shown < 6:
        shown += 1
        print('prim', prim, 'chain', c, 'pi-sets', [sorted(A) for A in ch[0]], 'pibar-sets', [sorted(A) for A in ch[1]])
