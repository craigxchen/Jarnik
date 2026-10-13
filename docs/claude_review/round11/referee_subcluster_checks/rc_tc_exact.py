# Referee: EXACT decision of the single tangent cones in the section-6.5 models (author used BFGS on 12 of 280).
# TC(A,m) = {u : k -> <k,u>^m restricted to A lies in Poly_{<m}|_A}.  If no nonzero polynomial of degree <= m
# vanishes on A (evaluation matrix of all monomials of degree <= m has full column rank), then TC = {0}: a form of
# degree m agreeing on A with a polynomial of degree < m would give a nonzero degree-m polynomial vanishing on A.
# Full column rank mod p implies full column rank over Q, so "TC = {0}" below is certified exactly.
import itertools, math, numpy as np
src = open('/home/user/Jarnik/docs/claude_review/round10/referee_determinant_checks/ref_fakes.py').read().split('\nout=[]')[0]
exec(src)   # profile constructors only (read-only reuse)
P = 1000003
def rank_modp(rows):
    A = np.array(rows, dtype=np.int64) % P
    nr, nc = A.shape; r = 0
    for c in range(nc):
        if r == nr: break
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0: continue
        pv = r + nz[0]
        if pv != r: A[[r, pv]] = A[[pv, r]]
        inv = pow(int(A[r, c]), P-2, P); A[r] = (A[r]*inv) % P
        col = A[:, c].copy(); col[r] = 0; idx = np.nonzero(col)[0]
        if len(idx): A[idx] = (A[idx] - (col[idx, None]*A[r][None, :]) % P) % P
        r += 1
    return r
def mons(d, deg):
    return [a for t in range(deg+1) for a in itertools.product(range(t+1), repeat=d) if sum(a) == t]
for name, Aexp in [("skew double core q=7", skew_double_core(7)), ("round-5 Paley-I b=8 q=7", round5_profile(7))]:
    M, r = Aexp.shape; n, T = 4, 3
    beta = (2*Aexp - 1).T
    nsys = 0; trivial = 0; ms = set(); sizes = set()
    for J in itertools.combinations(range(M), n):
        for b in J:
            bi = J.index(b); others = [i for i in range(n) if i != bi]
            box = np.array(list(itertools.product(range(-T, T+1), repeat=n-1)))
            K = np.zeros((len(box), n), dtype=int); K[:, others] = box; K[:, bi] = -box.sum(1)
            X = K @ beta[:, list(J)].T
            mx = np.abs(X).max(axis=1); w = np.sort(mx)[len(mx)//2]
            A = box[mx <= w]
            m = 1
            while 4*math.comb(m+1+n-2, n-1) <= len(A): m += 1
            ms.add(m); sizes.add(len(A))
            mm = mons(n-1, m)
            rows = [[int(np.prod([int(k[i])**a[i] for i in range(n-1)])) for a in mm] for k in A]
            rk = rank_modp(rows)
            nsys += 1
            if rk == len(mm): trivial += 1
    print(f"{name}: systems {nsys}, m values {sorted(ms)}, |A| {min(sizes)}..{max(sizes)}; "
          f"TC = {{0}} certified (no degree-<=m polynomial vanishes on A) in {trivial}/{nsys}", flush=True)
