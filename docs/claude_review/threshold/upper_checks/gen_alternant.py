"""Generalised alternants: c = sum_O gamma_O Alt(delta_{v_O}) over S_M-orbits O of points of the shell
S^M(p,n) with pairwise distinct coordinates.  For deg P < C(M,2)+d+1, Alt P = Vandermonde * S with S
symmetric of degree <= d, so c has ord >= C(M,2)+d+1 as soon as
    sum_O gamma_O Vand(v_O) S(v_O) = 0   for all symmetric S of degree <= d   (on sum k = p-n).
Exact rational linear algebra; then an exact moment check of the resulting measure.
Example of interest: M=5, (p,n)=(6,3): products of alternants give only 13, phi_5 = 14."""
import itertools, sys
from fractions import Fraction as F

def distinct_orbit_reps(M, p, n):
    reps = set()
    # choose M distinct integers in [-n, p] with positive sum p and negative sum n
    for comb in itertools.combinations(range(-n, p + 1), M):
        if sum(x for x in comb if x > 0) == p and -sum(x for x in comb if x < 0) == n:
            reps.add(tuple(comb))
    return sorted(reps)

def vand(v):
    r = 1
    for i in range(len(v)):
        for j in range(i + 1, len(v)):
            r *= (v[j] - v[i])
    return r

def esym(v, k):
    return sum(_prod(c) for c in itertools.combinations(v, k))

def _prod(c):
    r = 1
    for x in c:
        r *= x
    return r

def sym_basis_values(v, d):
    """values at v of the monomials e_2^b2 e_3^b3 ... e_M^bM with sum k*b_k <= d (e_1 is constant)."""
    M = len(v)
    e = [None, None] + [esym(v, k) for k in range(2, M + 1)]
    vals = []
    def rec(k, rem, cur):
        if k > M:
            vals.append(cur); return
        b = 0
        while k * b <= rem:
            rec(k + 1, rem - k * b, cur * e[k] ** b)
            b += 1
    rec(2, d, 1)
    return vals

def nullspace(rows, ncols):
    """rational nullspace of a matrix given as list of rows (Fractions)."""
    A = [list(map(F, r)) for r in rows]
    piv = []; r = 0
    for c in range(ncols):
        pr = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if pr is None: continue
        A[r], A[pr] = A[pr], A[r]
        inv = 1 / A[r][c]
        A[r] = [x * inv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
    free = [c for c in range(ncols) if c not in piv]
    basis = []
    for fcol in free:
        vec = [F(0)] * ncols; vec[fcol] = F(1)
        for i, pc in enumerate(piv):
            vec[pc] = -A[i][fcol]
        basis.append(vec)
    return basis

def sign(perm):
    s = 1; seen = [False] * len(perm)
    for i in range(len(perm)):
        if not seen[i]:
            L = 0; j = i
            while not seen[j]:
                seen[j] = True; j = perm[j]; L += 1
            if L % 2 == 0: s = -s
    return s

def build_and_check(M, p, n, d):
    reps = distinct_orbit_reps(M, p, n)
    rows = []
    nb = len(sym_basis_values(reps[0], d)) if reps else 0
    cols = [[vand(v) * x for x in sym_basis_values(v, d)] for v in reps]
    rows = [[cols[o][i] for o in range(len(reps))] for i in range(nb)]
    ns = nullspace(rows, len(reps))
    print(f"M={M} (p,n)=({p},{n}) d={d}: {len(reps)} distinct-coordinate orbits, {nb} symmetric conditions,"
          f" nullspace dim {len(ns)}")
    if not ns:
        return None
    gam = ns[0]
    # clear denominators
    from math import lcm
    den = 1
    for g in gam: den = lcm(den, g.denominator)
    gam = [int(g * den) for g in gam]
    meas = {}
    for g, v in zip(gam, reps):
        if g == 0: continue
        for perm in itertools.permutations(range(M)):
            k = tuple(v[perm[i]] for i in range(M))
            meas[k] = meas.get(k, 0) + g * sign(perm)
    meas = {k: c for k, c in meas.items() if c != 0}
    target = M * (M - 1) // 2 + d + 1
    # exact moment check over all monomials of degree < target in the first M-1 coordinates
    pts = list(meas.items())
    def monos(nv, dd):
        if nv == 1:
            yield (dd,); return
        for a in range(dd, -1, -1):
            for rest in monos(nv - 1, dd - a):
                yield (a,) + rest
    bad = 0
    for dd in range(target):
        for al in monos(M - 1, dd):
            s = 0
            for k, c in pts:
                t = c
                for i, e in enumerate(al):
                    if e: t *= k[i] ** e
                s += t
            if s != 0: bad += 1
    # check that the order is not higher than claimed is not needed (upper bound is Theorem 3.1)
    print(f"   gamma = {gam}; support {len(meas)} points; all moments of degree < {target} vanish: {bad == 0}")
    return bad == 0

if __name__ == '__main__':
    ok = build_and_check(5, 6, 3, 3)
    ok2 = build_and_check(5, 3, 6, 3)
    print("RESULT:", "PASS" if ok and ok2 else "FAIL")
