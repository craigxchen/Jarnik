"""Search for product certificates prod_i <z_i,u> in span{Phi_c : c in V(A,m)} (= J(A)_m^perp) in small exact models.

usage: python3 prodsearch.py model L m [maxsupp] [maxentry]
Pool: balanced integer z in Z^n_0 (n = d+1 points incl. base 0), primitive, sign-normalised, support 2..maxsupp,
entries in [-maxentry, maxentry].  The linear form is <z', u> with z' = z restricted to the non-base points.
Pruning (all NECESSARY conditions, so no certificate inside the pool is lost):
  * width counts: #{i : <z'_i, v> != 0} <= tau_A(v) for every test direction v in {-2..2}^d with tau_A(v) < m
    (tau_A(v) = #values of <k',v> on A minus 1; refJ Prop 2.1: J_tau subset J(A));
    this includes v = e_x (incidence <= 2) and v = -1 (base incidence <= 2).
Leaves: membership of the product in T_m mod p (search); survivors verified exactly over Q (fmpq solve).
"""
import sys, time, itertools
from math import gcd
import numpy as np
from fplib import *
from flint import fmpq_mat, fmpz_mat


def pool_vectors(n, maxsupp, maxentry):
    out = []
    seen = set()
    for s in range(2, maxsupp + 1):
        for S in itertools.combinations(range(n), s):
            for vals in itertools.product([e for e in range(-maxentry, maxentry + 1) if e], repeat=s):
                if sum(vals) != 0:
                    continue
                g = 0
                for v in vals:
                    g = gcd(g, abs(v))
                if g != 1:
                    continue
                if vals[0] < 0:
                    continue  # sign normalisation: first nonzero entry positive
                z = [0] * n
                for x, v in zip(S, vals):
                    z[x] = v
                t = tuple(z)
                if t not in seen:
                    seen.add(t)
                    out.append(t)
    return out


def main():
    name, L, m = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    maxsupp = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    maxentry = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    M, pats = MODELS[name]()
    A = box(M, pats, L)
    d = M - 1
    n = M
    Aarr = np.array(A, dtype=np.int64)
    t0 = time.time()
    # test directions
    dirs = []
    for v in itertools.product(range(-2, 3) if d <= 8 else range(-1, 2), repeat=d):
        if any(v):
            fnz = next(x for x in v if x)
            if fnz > 0:
                dirs.append(v)
    D = np.array(dirs, dtype=np.int64)
    tau = np.zeros(D.shape[0], dtype=np.int64)
    for s0 in range(0, D.shape[0], 2000):
        vals = Aarr @ D[s0:s0 + 2000].T  # |A| x chunk
        vs = np.sort(vals, axis=0)
        tau[s0:s0 + 2000] = (np.diff(vs, axis=0) != 0).sum(axis=0)
    keep = tau < m
    D = D[keep]
    tau = tau[keep]
    print(f'{name} L={L} m={m} |A|={len(A)} test dirs with tau<m: {len(D)}  ({time.time()-t0:.1f}s)', flush=True)
    pool = pool_vectors(n, maxsupp, maxentry)
    P = np.array(pool, dtype=np.int64)[:, 1:]  # non-base coordinates
    NZ = (P @ D.T) != 0  # pool x dirs
    # quick filter: a generator alone must satisfy counts (always true since tau >= 2 >= 1)
    print(f'pool size {len(pool)}', flush=True)
    countonly = len(sys.argv) > 6 and sys.argv[6] == 'countonly'
    if countonly:
        monos_m, Tb = reduced_monomials(d, 0), np.zeros((0, 1), dtype=np.int64)
    else:
        monos_m, Tb = tangent_space_mod(A, m)
    print(f'dim T_m (mod p) = {Tb.shape[0]}  N_m = {len(monos_m)}', flush=True)
    idx = {a: i for i, a in enumerate(monos_m)}
    Tn = to_nmod(Tb) if Tb.shape[0] else None
    rT = Tb.shape[0]
    cnt = np.zeros(len(D), dtype=np.int64)
    found = []
    leaves = [0]
    chosen = []

    def in_T(Z):
        poly = product_coeffs(Z, m)
        vec = np.zeros(len(monos_m), dtype=object)
        for e, c in poly.items():
            if max(e) > 2:
                return False
            vec[idx[e]] = fmpq(c, multinom(m, e))
        if rT == 0:
            return False
        rows = Tb.tolist() + [[int((x.p * pow(int(x.q), -1, PR)) % PR) if x != 0 else 0 for x in vec]]
        Nn = nmod_mat(len(rows), len(monos_m), [r for row in rows for r in row], PR)
        return Nn.rank() == rT

    supp = [frozenset(x for x in range(n) if pool[i][x]) for i in range(len(pool))]
    inc = [0] * n
    capacity = 2 * n

    def rec(start, used):
        k = len(chosen)
        if k == m:
            leaves[0] += 1
            Z = [list(pool[i][1:]) for i in chosen]
            if countonly:
                if len(found) < 50:
                    found.append(tuple(chosen))
            elif in_T(Z):
                found.append(tuple(chosen))
            return
        if used + 2 * (m - k) > capacity:
            return
        for i in range(start, len(pool)):
            S = supp[i]
            if used + len(S) + 2 * (m - k - 1) > capacity:
                continue
            if any(inc[x] >= 2 for x in S):
                continue
            row = NZ[i]
            cnt[row] += 1
            if (cnt <= tau).all():
                for x in S:
                    inc[x] += 1
                chosen.append(i)
                rec(i, used + len(S))
                chosen.pop()
                for x in S:
                    inc[x] -= 1
            cnt[row] -= 1

    rec(0, 0)
    print(f'leaves (multisets passing all width counts): {leaves[0]};  in T_m mod p: {len(found)}  ({time.time()-t0:.1f}s)',
          flush=True)
    for f in found[:3]:
        print('  ', [pool[i] for i in f])
    # zonotope test for every found certificate: does k0 + sum_i [0,1] z_i' fit in A for some k0?
    if found and not countonly:
        Aset = set(A)
        nz = 0
        for f in found:
            G = np.array([pool[i][1:] for i in f])
            verts = np.array([G[list(I)].sum(axis=0) if I else np.zeros(d, dtype=np.int64)
                              for r in range(len(G) + 1) for I in itertools.combinations(range(len(G)), r)])
            lo = verts.min(axis=0); hi = verts.max(axis=0)
            ok = any(all(tuple(v + np.array(k0)) in Aset for v in verts)
                     for k0 in itertools.product(*[range(-1 - int(l), 2 - int(h)) for l, h in zip(lo, hi)]))
            nz += ok
        print(f'certificates found: {len(found)}; of which genuine zonotopes (translate inside A): {nz}', flush=True)
    # exact verification of the first few
    if found and not countonly:
        low = []
        for e in range(m):
            low += reduced_monomials(d, e)
        Elow = eval_matrix(A, low)
        Em = eval_matrix(A, monos_m)
        Mfull = np.vstack([Elow, Em])
        for f in found[:5]:
            Z = [list(pool[i][1:]) for i in f]
            poly = product_coeffs(Z, m)
            rhs = [fmpq(0)] * len(low) + [fmpq(0)] * len(monos_m)
            for e, c in poly.items():
                rhs[len(low) + idx[e]] = fmpq(c, multinom(m, e))
            Mq = fmpq_mat(Mfull.shape[0], Mfull.shape[1], [int(x) for x in Mfull.flatten()])
            Aug = fmpq_mat(Mfull.shape[0], Mfull.shape[1] + 1,
                           [x for r in range(Mfull.shape[0]) for x in
                            [int(y) for y in Mfull[r]] + [rhs[r]]])
            r1 = fmpz_mat(Mfull.shape[0], Mfull.shape[1], [int(x) for x in Mfull.flatten()]).rank()
            r2 = Aug.rref()[0].rank() if hasattr(Aug.rref()[0], 'rank') else None
            print('   exact: rank M =', r1, ' rank [M|rhs] =', r2, ' => in span over Q:', r1 == r2, flush=True)


if __name__ == '__main__':
    main()
