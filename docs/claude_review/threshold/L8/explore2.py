from c5ansatz import *
import random, numpy as np
random.seed(2)

def interp(xs, ys):
    n = len(xs); res = []
    for i in range(n):
        num = [Fr(1)]; den = Fr(1)
        for j in range(n):
            if j != i:
                num = mul(num, linpoly(1, -xs[j])); den *= xs[i] - xs[j]
        res = add(res, scal(ys[i] / den, num))
    return res

def detfun(rho, sigma):
    try:
        return det_frac(matrix(rho, sigma))
    except ZeroDivisionError:
        return None

for trial in range(2):
    rho = {e: Fr(random.randint(-20, 20), random.randint(1, 5)) for e in EDGES}
    sigma = {T: Fr(random.randint(-20, 20), random.randint(1, 5)) for T in TRIPLES}
    for T in TRIPLES:
        # det * prod over rows i in T of A_{iT} * (Vandermonde-ish) is a polynomial in sigma_T
        def N(x):
            s = dict(sigma); s[T] = x
            d = detfun(rho, s)
            D = Fr(1)
            for i in T:
                e1, e2 = nbr_edges(i)
                D *= (x - rho[e1]) * (x - rho[e2])
            return d * D
        xs = [Fr(997 + 13 * k, 7) for k in range(25)]
        p = interp(xs, [N(x) for x in xs])
        rr = rational_roots(p)
        print(T, 'deg', deg(p), 'rational roots', rr, 'other sigmas', sorted(sigma[t] for t in TRIPLES if t != T))

print('---- trial with distinct values, report rho')
random.seed(5)
vals = random.sample(range(-60, 60), 10)
rho = {e: Fr(vals[k], 1) for k, e in enumerate(EDGES)}
sigma = {T: Fr(vals[5 + k], 1) for k, T in enumerate(TRIPLES)}
print('rho', {e: str(v) for e, v in rho.items()})
print('sigma', {T: str(v) for T, v in sigma.items()})
for T in TRIPLES:
    def N(x):
        s = dict(sigma); s[T] = x
        d = detfun(rho, s)
        D = Fr(1)
        for i in T:
            e1, e2 = nbr_edges(i)
            D *= (x - rho[e1]) * (x - rho[e2])
        return d * D
    xs = [Fr(997 + 13 * k, 7) for k in range(25)]
    p = interp(xs, [N(x) for x in xs])
    rr = rational_roots(p)
    q = p
    for r in rr:
        while ev(q, r) == 0:
            q, _ = divmod_(q, linpoly(1, -r))
    print(T, 'deg', deg(p), 'rat roots', [str(r) for r in rr], 'remaining deg', deg(q))
