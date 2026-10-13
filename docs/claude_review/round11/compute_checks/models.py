"""Find regime-faithful model instances.
Conditions on (n, profile, T, L, m), uniform slab |k(S)| <= L for every column S:
 (i)   slabs cut: A != B_T
 (ii)  h_A(beta_p) = 2*width_A(k(S_p)) < m for every column
 (iii) V(A,m) != 0, i.e. m <= rho(A) := max{j : H_A(j-1) < |A|}
 (iv)  regime-faithful: pi(A) < m, i.e. A contains no lattice parallelepiped with sum L_i >= m
       (otherwise ZE holds trivially at every nondegenerate point, compute.md Prop 3.1)
 (v)   full-rank profile (functionals k(S_p) span (R^d)^*), as for nondegenerate clusters.
rho via exact fmpz rank (monotone search).  pi test exhaustive (see has_parallelepiped).
"""
import itertools, math, sys, random
import numpy as np
import flint
from core import *

def rho(A):
    U = len(A)
    j = 0
    while True:
        r = moment_matrix(A, j + 1).rank()   # dim Poly_{<=j}|_A
        if r == U:
            return j   # (S/J)_j != 0 for j <= this value; V(A,m)!=0 iff m <= j
        j += 1

def in_Q(pt, T, slabs):
    if any(abs(x) > T for x in pt):
        return False
    kf = list(pt) + [-sum(pt)]
    return all(abs(sum(kf[x] for x in S)) <= L for S, L in slabs)

def primitive_vectors(d, R):
    out = []
    for v in itertools.product(range(-R, R + 1), repeat=d):
        if all(x == 0 for x in v):
            continue
        if math.gcd(*[abs(x) for x in v]) != 1:
            continue
        # sign normalisation: first nonzero positive
        fnz = next(x for x in v if x != 0)
        if fnz < 0:
            continue
        out.append(v)
    return out

def has_parallelepiped(A, T, slabs, m, d):
    """Exhaustive: is there P = k0 + sum_{i<=j} t_i v_i (0<=t_i<=L_i), v_i lin. indep. integer vectors,
    sum L_i = m, all lattice points in A?  Non-primitive v_i are covered by primitive ones with larger L.
    Since P is inside the box, L_i * |v_i|_inf <= 2T.  P subset A iff its vertices lie in the convex polytope Q."""
    Aset = set(A)
    R = 2 * T
    vecs = primitive_vectors(d, R)
    cap = {v: (2 * T) // max(abs(x) for x in v) for v in vecs}
    vecs = [v for v in vecs if cap[v] >= 1]
    vecs.sort(key=lambda v: -cap[v])
    best = None
    for j in range(1, d + 1):
        for combo in itertools.combinations(vecs, j):
            caps = [cap[v] for v in combo]
            if sum(caps) < m:
                continue
            if np.linalg.matrix_rank(np.array(combo, dtype=float)) < j:
                continue
            # enumerate L with sum = m, L_i <= cap_i, L_i >= 1
            def Ls(i, rem):
                if i == j - 1:
                    if 1 <= rem <= caps[i]:
                        yield (rem,)
                    return
                for a in range(1, min(caps[i], rem - (j - 1 - i)) + 1):
                    for rest in Ls(i + 1, rem - a):
                        yield (a,) + rest
            for Lv in Ls(0, m):
                gens = [tuple(Li * x for x in v) for Li, v in zip(Lv, combo)]
                for k0 in A:
                    ok = True
                    for eps in itertools.product((0, 1), repeat=j):
                        pt = tuple(k0[t] + sum(e * g[t] for e, g in zip(eps, gens)) for t in range(d))
                        if pt not in Aset:
                            ok = False
                            break
                    if ok:
                        return (k0, combo, Lv)
    return None

def fullrank(n, prof):
    d = n - 1
    rows = [[(1 if x in S else 0) - (1 if (n - 1) in S else 0) for x in range(d)] for S in prof]
    return np.linalg.matrix_rank(np.array(rows, dtype=float)) == d

def canon(n, S):
    Sc = tuple(x for x in range(n) if x not in S)
    return min(tuple(S), Sc, key=lambda s: (len(s), s))

def random_profile(n, rng):
    r = rng.randint(n, 3 * n)
    prof = set()
    for _ in range(r):
        c = [rng.choice([1, -1]) for _ in range(n)]
        S = tuple(x for x in range(n) if c[x] == 1)
        if 0 < len(S) < n:
            prof.add(canon(n, S))
    return tuple(sorted(prof))

def analyse(n, prof, T, L, verbose=True):
    slabs = [(S, L) for S in prof]
    A = build_A(n, T, slabs)
    d = n - 1
    if len(A) == (2 * T + 1) ** d:
        return []
    hmax = max(2 * width(A, S) for S in prof)
    rh = rho(A)
    res = []
    for m in range(hmax + 1, rh + 1):
        P = has_parallelepiped(A, T, slabs, m, d)
        res.append((m, P is None, P))
    return dict(A=len(A), hmax=hmax, rho=rh, ms=res)

if __name__ == "__main__":
    n = int(sys.argv[1])
    rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    profs = []
    if n == 4:
        profs.append(((0, 1), (0, 2), (0, 3)))          # Hadamard = skew double core q=3
    while len(profs) < (int(sys.argv[3]) if len(sys.argv) > 3 else 6):
        p = random_profile(n, rng)
        if fullrank(n, p) and p not in profs:
            profs.append(p)
    Ts = {4: [2, 3, 4], 5: [2, 3]}[n]
    for prof in profs:
        for L in [1, 2, 3]:
            for T in Ts:
                if T < L:
                    continue
                if (2 * T + 1) ** (n - 1) > 2500:
                    continue
                r = analyse(n, prof, T, L)
                if not r:
                    continue
                good = [mm for mm, ok, P in r['ms'] if ok]
                print("n=%d prof=%s T=%d L=%d |A|=%d hmax=%d rho=%d faithful m: %s  (pi>=m witnesses: %s)" % (
                    n, prof, T, L, r['A'], r['hmax'], r['rho'], good,
                    [(mm, P[1], P[2]) for mm, ok, P in r['ms'] if not ok][:2]), flush=True)
