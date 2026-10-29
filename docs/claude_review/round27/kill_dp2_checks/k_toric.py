# Killer: maximal order at e of sections of t*D on the Losev-Manin toric variety L_k (permutohedral fan in Z^k/Z1),
# D = sum_A d_A E_A. Exact for "no section of order r" via rank mod a prime (rank_p <= rank_Q);
# existence of order r_max via exact integer rank (flint fmpz_mat).
import itertools, sys, time
import flint

def subsets(k):
    for r in range(1, k):
        for A in itertools.combinations(range(k), r):
            yield A

def lattice_points(k, d, t):
    # m in Z^k, sum m = 0, sum_{i in A} m_i >= -t d_A
    lo = [-t*d[(i,)] for i in range(k)]
    hi = [t*d[tuple(j for j in range(k) if j != i)] for i in range(k)]
    pts = []
    for m in itertools.product(*[range(lo[i], hi[i] + 1) for i in range(k - 1)]):
        mk = -sum(m)
        if not (lo[k-1] <= mk <= hi[k-1]):
            continue
        mm = m + (mk,)
        if all(sum(mm[i] for i in A) >= -t*d[A] for A in subsets(k)):
            pts.append(mm)
    return pts

def monomials(nv, deg_lt):
    out = []
    for dd in range(deg_lt):
        for c in itertools.combinations_with_replacement(range(nv), dd):
            e = [0]*nv
            for i in c: e[i] += 1
            out.append(tuple(e))
    return out

def matrix(pts, r, k):
    mons = monomials(k - 1, r)
    rows = []
    for e in mons:
        rows.append([eval_mon(m, e) for m in pts])
    return rows

def eval_mon(m, e):
    v = 1
    for i, ei in enumerate(e):
        if ei: v *= m[i]**ei
    return v

P = 2305843009213693951  # 2^61 - 1, prime
def rank_mod(rows, ncols):
    if not rows: return 0
    M = flint.nmod_mat(len(rows), ncols, [x % P for row in rows for x in row], P)
    return M.rank()

def max_order(k, d, t, rmax=40):
    pts = lattice_points(k, d, t)
    n = len(pts)
    r = 1
    while True:
        rows = matrix(pts, r + 1, k)   # conditions for order >= r+1
        rk = rank_mod(rows, n)
        if rk == n:   # no section of order r+1 (certified)
            break
        r += 1
        if r > rmax: break
    # certify existence of order r exactly
    rows = matrix(pts, r, k)
    if len(rows) < n:
        exist = True
    else:
        Mz = flint.fmpz_mat(rows)
        exist = Mz.rank() < n
    return n, r, exist

if __name__ == '__main__':
    tmax = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    for k in (5, 4):
        subs = list(subsets(k))
        mK = {A: 1 for A in subs}
        cases = [("-K", mK)]
        if k == 5:
            cases.append(("L'=2(-K)-B_{1,4}", {A: (1 if len(A) in (1, 4) else 2) for A in subs}))
        # K + sum Z_ab with orientation rep: Z_ab ~ sum_{A: b in A, a not in A} E_A  (a<b)
        nA = {A: sum(1 for a in range(k) for b in range(a + 1, k) if (b in A and a not in A)) for A in subs}
        cases.append(("K+sum Z_ab", {A: nA[A] - 1 for A in subs}))
        for name, d in cases:
            for t in range(1, tmax + 1):
                t0 = time.time()
                n, r, ex = max_order(k, d, t)
                print(f"k={k} D={name} t={t}: #pts={n}, max order at e = {r} (exists exactly: {ex}), per unit {r/t:.3f}  ({time.time()-t0:.1f}s)", flush=True)
