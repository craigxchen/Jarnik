from c5ansatz import *
import random
random.seed(1)

def interp(xs, ys):
    # Lagrange interpolation exact
    n = len(xs); res = []
    for i in range(n):
        num = [Fr(1)]; den = Fr(1)
        for j in range(n):
            if j != i:
                num = mul(num, linpoly(1, -xs[j])); den *= xs[i] - xs[j]
        res = add(res, scal(ys[i] / den, num))
    return res

for trial in range(3):
    rho = {e: Fr(random.randint(-20, 20), random.randint(1, 5)) for e in EDGES}
    sigma = {T: Fr(random.randint(-20, 20), random.randint(1, 5)) for T in TRIPLES}
    e12 = (1, 2)
    def N(x):
        r = dict(rho); r[e12] = x
        D = prod([linpoly(-1, sigma[T]) for T in triples_of(1)] + [linpoly(-1, sigma[T]) for T in triples_of(2)])
        return det_frac(matrix(r, sigma)) * ev(D, x)
    xs = [Fr(1000 + 7 * k, 3) for k in range(12)]
    p = interp(xs, [N(x) for x in xs])
    print('deg N in rho12:', deg(p))
    print(' sigma(1,2,4)=', sigma[(1,2,4)], ' roots:', rational_roots(p))
    # check multiplicity of root sigma124
    q = p
    m = 0
    while q and ev(q, sigma[(1,2,4)]) == 0:
        q, _ = divmod_(q, linpoly(1, -sigma[(1,2,4)])); m += 1
    print(' multiplicity of sigma124:', m, ' remaining degree', deg(q))
    import numpy as np
    print(' remaining roots (float):', np.roots([float(c) for c in reversed(q)]))
