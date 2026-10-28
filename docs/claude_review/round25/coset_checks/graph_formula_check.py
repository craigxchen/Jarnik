"""Check of the exact step-sum formula of Corollary 4.2 (graph coordinates):
 sum_{v : last(mu_v) = k} y_v = Pi_{J,<k} * sum_x y(x) E_{>k}(x) (Pi_{j_k} - p_{j_k}^{-|f(x)_k|/2}),
 mu_v = v_J + C^T v_{J^c}, f(x) = C^T x.  Brute force over a box |v_i| <= Nbox (tails < 1e-12)."""
import itertools, math
import numpy as np
ps = [13, 17, 29, 37, 41]          # J = first m coords, J^c = rest
m = 3
C = np.array([[1, -2, 0], [0, 1, 1]])   # rows c_i for the 2 dependent coords (i not in J): phi_i = <c_i, phi_J> + g_i
r = len(ps); Nbox = 22
def y_of(v): return math.prod(p ** (-abs(a) / 2) for p, a in zip(ps, v))
brute = [0.0] * m
rng = range(-Nbox, Nbox + 1)
# brute force with pruning: enumerate x in box for dependent coords, v_J in box
for x in itertools.product(rng, repeat=r - m):
    yx = math.prod(p ** (-abs(a) / 2) for p, a in zip(ps[m:], x))
    if yx < 1e-14: continue
    f = np.array(x) @ C
    for vJ in itertools.product(range(-13, 14), repeat=m):
        yv = yx * math.prod(p ** (-abs(a) / 2) for p, a in zip(ps[:m], vJ))
        if yv < 1e-18: continue
        mu = np.array(vJ) + f
        nz = np.nonzero(mu)[0]
        if len(nz) == 0: continue
        brute[nz[-1]] += yv
Pi = [1 + 2 / (math.sqrt(p) - 1) for p in ps]
formula = [0.0] * m
for k in range(m):
    PiJlt = math.prod(Pi[:k])
    s = 0.0
    for x in itertools.product(rng, repeat=r - m):
        yx = math.prod(p ** (-abs(a) / 2) for p, a in zip(ps[m:], x))
        if yx < 1e-14: continue
        f = np.array(x) @ C
        E = math.prod(ps[kk] ** (-abs(f[kk]) / 2) for kk in range(k + 1, m))
        s += yx * E * (Pi[k] - ps[k] ** (-abs(f[k]) / 2))
    formula[k] = PiJlt * s
print('brute  ', brute)
print('formula', formula)
print('max rel diff', max(abs(a - b) / b for a, b in zip(brute, formula)))
