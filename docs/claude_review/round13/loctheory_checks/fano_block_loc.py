"""Complete Fano block Y_3(j) (all 7 columns of a 3-dim XOR block, large-n Sylvester image):
which local families generate J(Y_3(j)) in degrees <= D?

Coordinates y_a = k(S_a), a in F_2^3 \\ 0.  Class sums c_z (z in F_2^3), sum c = 0, y_a = sum_{z: a.z=1} c_z, so
c_z = -(1/4) sum_{a != 0} y_a chi_a(z)  (chi_a(z) = (-1)^{a.z}); class-sum directions are integer directions of A.
Local ideal candidates (all generators exact over Q, dims mod p; sandwich as in §1.4):
  W1 = grid + width powers of all directions v = sum_z m_z c_z, m in {-1,0,1}^8 (mod constants)
  T  = pullbacks from the 7 lines (3-circuits) and 7 quadrangles (4-circuits)   [Theorem 3.3 family]
"""
import itertools, sys
from fractions import Fraction
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from lfi import monos, hilbert, leading_forms, ideal_dim_modp, P_DEFAULT
from circuits_loc import section, proj

def chi(a, z):
    return 1 - 2 * (bin(a & z).count('1') % 2)

def poly_power(coeffs, power, nvars):
    """(sum_q coeffs[q] y_q)^power as dict exponent->int (coeffs integers)"""
    poly = {tuple([0] * nvars): 1}
    for _ in range(power):
        new = {}
        for e, c in poly.items():
            for q, a in enumerate(coeffs):
                if a == 0:
                    continue
                f = list(e); f[q] += 1; f = tuple(f)
                new[f] = new.get(f, 0) + c * a
        poly = new
    return {e: c for e, c in poly.items() if c}

def run(j, D):
    F = list(range(1, 8))
    P = section(3, F, j)
    n = 7
    HP = hilbert(P, D, cap=2 * j, modp=True)
    print(f'Y_3({j}): |P|={len(P)} H_P(modp, up to {D}) = {HP}')
    grid = []
    for q in range(n):
        e = [0] * n; e[q] = 2 * j + 1
        grid.append({tuple(e): 1})
    # width powers of class-combination directions, written in y-coordinates (scaled to integers)
    W = []
    seenv = set()
    for m in itertools.product((-1, 0, 1), repeat=8):
        if all(x == 0 for x in m):
            continue
        # v(y) = sum_z m_z c_z = -(1/4) sum_a y_a sum_z m_z chi_a(z)
        coef = [-sum(m[z] * chi(a, z) for z in range(8)) for a in F]   # = 4 v
        if all(c == 0 for c in coef):
            continue
        g = 0
        for c in coef:
            from math import gcd
            g = gcd(g, abs(c))
        prim = tuple(c // g for c in coef)
        if prim[next(i for i, c in enumerate(prim) if c)] < 0:
            prim = tuple(-c for c in prim)
        if prim in seenv:
            continue
        seenv.add(prim)
        vals = set(sum(c * y for c, y in zip(prim, pt)) for pt in P)
        h = len(vals)
        if h <= D:
            W.append(poly_power(list(prim), h, n))
    print(f'  {len(seenv)} class-combination directions, {len(W)} width powers of degree <= {D}')
    T = []
    lines = [s for s in itertools.combinations(range(7), 3) if (F[s[0]] ^ F[s[1]] ^ F[s[2]]) == 0]
    quads = [s for s in itertools.combinations(range(7), 4) if (F[s[0]] ^ F[s[1]] ^ F[s[2]] ^ F[s[3]]) == 0]
    for G in lines + quads:
        PG = proj(P, G)
        for i in range(2 * j + 1, D + 1):
            top, basis = leading_forms(PG, i, cap=2 * j)
            for v in basis:
                g = {}
                for e, c in zip(top, v):
                    if c:
                        ee = [0] * n
                        for gq, eq in zip(G, e):
                            ee[gq] = eq
                        g[tuple(ee)] = c
                T.append(g)
    for name, gens in [('grid', grid), ('grid+W', grid + W), ('grid+T', grid + T), ('grid+T+W', grid + T + W)]:
        defs = []
        for i in range(D + 1):
            Ni = len(monos(n, i, 2 * j))
            dJ = Ni - (HP[i] if i < len(HP) else 0)
            defs.append(dJ - ideal_dim_modp(i, n, gens, cap=2 * j))
        print(f'  {name}: deficits dim J_i - dim L_i, i=0..{D}: {defs}')
        sys.stdout.flush()

if __name__ == '__main__':
    run(2, 8)
