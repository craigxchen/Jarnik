"""Exact certificate: J(Y_3(2))_7 is NOT contained in L_7, L = grid + lines + quadrangles + width powers of all
class-combination directions (m in {-1,0,1}^8).  Everything over Z/Q:
  dim_Q J(P)_e = N_e - H_P(e), H_P from exact fmpz ranks of evaluation matrices (degrees <= 7);
  dim_Q L_e from exact fmpz rank of the generator-multiple matrix in reduced monomials (exponents <= 4).
"""
import itertools, sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import flint
from lfi import monos, evalrow, leading_forms
from circuits_loc import section, proj
from fano_block_loc import chi, poly_power

def exact_h(P, D, cap):
    mons = []
    h = []
    for i in range(D + 1):
        mons += monos(len(P[0]), i, cap)
        M = flint.fmpz_mat([evalrow(pt, mons) for pt in P])
        h.append(M.rank())
        print('   exact h_P(%d) = %d' % (i, h[-1])); sys.stdout.flush()
    return [h[0]] + [h[i] - h[i - 1] for i in range(1, len(h))]

def exact_ideal_dim(i, n, gens, cap):
    top = monos(n, i, cap)
    idx = {e: a for a, e in enumerate(top)}
    rows = []
    for g in gens:
        dg = sum(next(iter(g)))
        if dg > i:
            continue
        for mu in monos(n, i - dg, None):
            row = [0] * len(top)
            for e, c in g.items():
                f = tuple(a + b for a, b in zip(e, mu))
                if max(f) > cap:
                    continue
                row[idx[f]] += c
            if any(row):
                rows.append(row)
    if not rows:
        return 0, len(top)
    return flint.fmpz_mat(rows).rank(), len(top)

if __name__ == '__main__':
    j, D = 2, 7
    F = list(range(1, 8)); n = 7; cap = 2 * j
    P = section(3, F, j)
    HP = exact_h(P, D, cap)
    print('exact H_P up to 7:', HP)
    from math import gcd
    W, seen = [], set()
    for m in itertools.product((-1, 0, 1), repeat=8):
        coef = [-sum(m[z] * chi(a, z) for z in range(8)) for a in F]
        if all(c == 0 for c in coef):
            continue
        g = 0
        for c in coef:
            g = gcd(g, abs(c))
        prim = tuple(c // g for c in coef)
        if prim[next(i for i, c in enumerate(prim) if c)] < 0:
            prim = tuple(-c for c in prim)
        if prim in seen:
            continue
        seen.add(prim)
        h = len(set(sum(c * y for c, y in zip(prim, pt)) for pt in P))
        if h <= D:
            W.append(poly_power(list(prim), h, n))
    T = []
    lines = [s for s in itertools.combinations(range(7), 3) if (F[s[0]] ^ F[s[1]] ^ F[s[2]]) == 0]
    quads = [s for s in itertools.combinations(range(7), 4) if (F[s[0]] ^ F[s[1]] ^ F[s[2]] ^ F[s[3]]) == 0]
    for G in lines + quads:
        PG = proj(P, G)
        for i in range(2 * j + 1, D + 1):
            top, basis = leading_forms(PG, i, cap=cap)
            for v in basis:
                g = {}
                for e, c in zip(top, v):
                    if c:
                        ee = [0] * n
                        for gq, eq in zip(G, e):
                            ee[gq] = eq
                        g[tuple(ee)] = c
                T.append(g)
    for i in range(5, D + 1):
        dL, Ni = exact_ideal_dim(i, n, T + W, cap)
        dJ = Ni - HP[i]
        print(f'degree {i}: reduced monomials {Ni}, dim_Q J(P)_i/grid = {dJ}, dim_Q L_i/grid = {dL}, deficit = {dJ - dL}')
        sys.stdout.flush()
