# Independent reproduction of the VII.1 construction with q = p = 3, k = 2 (no constant extension needed),
# in ramified models eta^2 = Delta over F_3[t].
import sys, itertools
import numpy as np
from fp import *

def chi(y):
    y %= 3
    return 0 if y == 0 else (1 if y == 1 else -1)

ALLFORMS = [(1, 0), (0, 1), (1, 1), (1, 2)]   # x0, x1, x0+x1, x0-x1 (mutually nonproportional over F_3)

def build(Delta, C, s, forms, signs, verbose=True, check_sf=True):
    O = Order(Delta)
    k = 2
    Dd = deg(O.D)
    labels = list(itertools.product(range(3), repeat=k))
    # v = eta^3 C^3 (x0 + x1 t^3) = eta * Delta * C^3 * (x0 + x1 t^3)
    C3 = powp(C, 3)
    base = mul(O.D, C3)
    def label(x):
        w = add(const(x[0]), mono(3, x[1]))
        return (np.zeros(0, dtype=np.int64), mul(base, w))
    V = {x: label(x) for x in labels}
    # T
    if s % 3:
        T = add(const(1), mul(C3, mono(s)))
    else:
        T = add(const(1), mul(C3, add(mono(s), mono(1))))
    L = deg(T)
    Tel = (T, np.zeros(0, dtype=np.int64))
    m = len(forms)
    E = {}
    for x in labels:
        E[x] = [1 + chi(sg * (f[0] * x[0] + f[1] * x[1])) for f, sg in zip(forms, signs)]
    gv = {x: min(E[x]) for x in labels}
    g = sum(gv.values())
    n = len(labels)
    # marks
    fac = {x: O.sub(Tel, V[x]) for x in labels}
    zs = []
    for j in range(m):
        z = (const(1), np.zeros(0, dtype=np.int64))
        for x in labels:
            e = E[x][j] - gv[x]
            if e: z = O.mul(z, O.pow(fac[x], e))
        zs.append(z)
    out = {'L': L, 'g': g, 'n': n}
    norms = [O.norm(z) for z in zs]
    out['equal_norms'] = all(np.array_equal(norms[0], N_) for N_ in norms)
    N0 = norms[0]
    out['degN'] = deg(N0)
    out['D2'] = [O.pole2(z) for z in zs]   # twice pole order
    out['distinct'] = len({(tuple(z[0]), tuple(z[1])) for z in zs}) == m
    # pair differences
    w2 = []
    for i in range(m):
        for j in range(i + 1, m):
            w2.append(O.pole2(O.sub(zs[i], zs[j])))
    out['maxdiff2'] = max(w2) if w2 else None
    out['critical'] = (max(w2) <= out['D2'][0] // 2) if w2 else None  # 2*deg(diff) <= D  <=> pole2(diff) <= D2/2
    # R
    R = (const(1), np.zeros(0, dtype=np.int64))
    for x in labels:
        R = O.mul(R, fac[x])
    out['R_in_base'] = len(R[1]) == 0
    Rp = R[0]
    if check_sf:
        out['R_squarefree'] = squarefree(Rp)
        # G(T) conj(G(T)) * N0 == R^2
        Gn = const(1)
        for x in labels:
            e = gv[x]
            if e:
                Gn = mul(Gn, powp(O.norm(fac[x]), e))
        # careful: Norm(G(T)) = prod_x Norm(T - v_x)^{g(x)} -- each norm is (T-v)(T+v)
        out['N0_times_normG_eq_R2'] = np.array_equal(mul(N0, Gn), mul(Rp, Rp))
        me, dec = max_exponent(N0)
        out['N0_max_exponent'] = me
        out['N0_mult_degrees'] = {kk: deg(v) for kk, v in dec.items()}
    out['zs'] = zs
    out['N0'] = N0
    return out

if __name__ == '__main__':
    pass
