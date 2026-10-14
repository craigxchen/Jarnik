"""Locality inside Sylvester column-family images P = Lambda_F cap [-j,j]^F.

L_b := J(grid) + sum_{G subset F, |G| <= b} pi_G^* J(pi_G(P))   (coordinate sub-families)
Certification (sandwich): generators of J(pi_G P) computed EXACTLY over Q (fmpz nullspace), then the ideal
dimension mod p (<= rational dim of L_b <= dim_Q J(P)); H_P mod p (h_P^(p) <= h_P).  If the cumulative mod-p counts
agree, L_b = J(P) in those degrees over Q.  If they disagree the deficit is only evidence of non-locality
(unless the gap is certified separately).
"""
import itertools, sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from lfi import monos, hilbert, leading_forms, ideal_dim_modp, P_DEFAULT

def dot(a, z):
    return bin(a & z).count('1') % 2

def section(t, F, j):
    allc = list(range(1, 2 ** t))
    rest = [a for a in allc if a not in F]
    mod = 2 ** (t - 1)
    R = mod
    pts = []
    for yF in itertools.product(range(-j, j + 1), repeat=len(F)):
        for yr in itertools.product(range(0, R), repeat=len(rest)):
            full = dict(zip(F, yF)); full.update(zip(rest, yr))
            if all(sum(full[a] * (1 - 2 * dot(a, z)) for a in allc) % mod == 0 for z in range(2 ** t)):
                pts.append(yF); break
    return pts

def proj(P, G):
    return sorted(set(tuple(y[g] for g in G) for y in P))

def local_gens(P, b, j, Dmax, b_list=None):
    """exact generators of pi_G^* J(pi_G P) for |G| <= b (degrees <= Dmax), embedded in |F| variables."""
    n = len(P[0])
    gens = []
    for q in range(n):  # grid
        e = [0] * n; e[q] = 2 * j + 1
        gens.append({tuple(e): 1})
    for bb in range(2, b + 1):
        for G in itertools.combinations(range(n), bb):
            PG = proj(P, G)
            if len(PG) == (2 * j + 1) ** bb:
                continue  # full grid: nothing new (Lemma: grids invisible)
            for i in range(2 * j + 1, Dmax + 1):
                top, basis = leading_forms(PG, i, cap=2 * j)
                for v in basis:
                    g = {}
                    for e, c in zip(top, v):
                        if c:
                            ee = [0] * n
                            for gq, eq in zip(G, e):
                                ee[gq] = eq
                            g[tuple(ee)] = c
                    gens.append(g)
    return gens

def run(t, F, j, Dmax=None):
    P = section(t, F, j)
    n = len(F)
    HP = hilbert(P, 10 ** 6, cap=2 * j, modp=True)
    rP = len(HP) - 1
    if Dmax is None:
        Dmax = rP + 1
    print(f't={t} F={F} j={j} |P|={len(P)} H_P(modp)={HP} rho<= {rP}')
    for b in range(1, n + 1):
        gens = local_gens(P, b, j, Dmax)
        cum_L, cum_J = 0, 0
        defs = []
        for i in range(Dmax + 1):
            Ni = len(monos(n, i, 2 * j))
            dL = ideal_dim_modp(i, n, gens, cap=2 * j)
            dJ = Ni - (HP[i] if i < len(HP) else 0)
            cum_L += dL; cum_J += dJ
            defs.append(dJ - dL)
        ok = all(x == 0 for x in defs)
        print(f'   b={b}: deficits dim J(P)_i - dim L_b,i (i=0..{Dmax}) = {defs}  {"LOCAL (certified)" if ok else ""}')
        sys.stdout.flush()
        if ok:
            break

if __name__ == '__main__':
    cases = [
        (2, [1, 2, 3], 2), (2, [1, 2, 3], 3),
        (3, [1, 2, 4, 7], 2),            # single 4-circuit
        (3, [1, 2, 3, 4, 5], 2),         # code dim 2 (two 3-circuits sharing 1 + a 4-circuit)
        (3, [1, 2, 3, 4, 5, 6], 2),      # code dim 3
        (3, [1, 2, 3, 4, 5, 6, 7], 2),   # complete Fano block Y_3(2) (non-binary)
        (3, [1, 2, 3, 4, 5, 6, 7], 3),
    ]
    which = sys.argv[1:] and [int(a) for a in sys.argv[1:]] or range(len(cases))
    for k in which:
        run(*cases[k])
