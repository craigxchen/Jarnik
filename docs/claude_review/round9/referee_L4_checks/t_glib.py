from glib import *
import random
rng = random.Random(1)
# test torus logs at several q^A
for q, A in [(3,3),(5,2),(7,2),(13,1),(23,2),(31,1),(47,1),(79,1),(97,2),(1999,1)]:
    T = Torus(q, A)
    for _ in range(50):
        w = (rng.randrange(T.Q), rng.randrange(T.Q))
        if (w[0]**2+w[1]**2) % q == 0: continue
        h = ratio(w, T.Q)
        x = T.log(h)
        assert mpow(T.g, x, T.Q) == h
    # i has order 4
    assert mpow((0,1), 4, T.Q) == (1,0)
    assert (4*T.iota) % T.m == 0 and (2*T.iota) % T.m != 0
    # parity of log(w/conj w) = [(Norm w / q) = -1]
    for _ in range(50):
        w = (rng.randrange(T.Q), rng.randrange(T.Q))
        n = (w[0]**2+w[1]**2)
        if n % q == 0: continue
        x = T.log(ratio(w, T.Q))
        assert (x % 2) == (1 if legendre(n, q) == -1 else 0), (q, A)
print("glib OK")
