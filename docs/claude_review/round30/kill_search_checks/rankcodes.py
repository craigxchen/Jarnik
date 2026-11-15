# Independent max size of rank<=r base-form codes with distance >= L+1 and pairwise balance >= h (kill_search).
# Rows: base 0, basis e_1..e_r, extra rows lambda in Lambda_r. Columns (u, s), u in {0,1}^r, s = sign Im q.
import sys, itertools, numpy as np
from fractions import Fraction
from scipy.optimize import milp, LinearConstraint, Bounds
import sympy as sp

def Lam(r):
    U = [u for u in itertools.product([0,1], repeat=r) if any(u)]
    out = set()
    for B in itertools.combinations(U, r):
        Mb = sp.Matrix(B)
        if Mb.det() == 0: continue
        inv = Mb.inv()
        for vals in itertools.product([0,1], repeat=r):
            lam = tuple(inv * sp.Matrix(vals))
            out.add(lam)
    basis = set(tuple(sp.Integer(1) if i == j else sp.Integer(0) for i in range(r)) for j in range(r))
    out -= basis; out.discard(tuple([sp.Integer(0)]*r))
    return sorted(out, key=lambda t: tuple(float(x) for x in t))

def rowvals(rows, r):
    # rows: list of functionals; returns allowed u list and value table
    allu = list(itertools.product([0,1], repeat=r))
    A = []
    for u in allu:
        vals = [sum(l[i]*u[i] for i in range(r)) for l in rows]
        if all(v in (0, 1) for v in vals): A.append((u, [int(v) for v in vals]))
    return A

def feasible(rows, r, D, L, h, balance=True):
    A = rowvals(rows, r)
    nv = 2*len(A)
    if nv == 0: return False
    cons_A, lo = [], []
    R = len(rows)
    for i in range(R):
        for j in range(i+1, R):
            a1 = np.zeros(nv); a2 = np.zeros(nv)
            for idx, (u, v) in enumerate(A):
                if v[i] < v[j]: a1[2*idx] += 1; a2[2*idx+1] += 1
                elif v[i] > v[j]: a2[2*idx] += 1; a1[2*idx+1] += 1
            if balance:
                cons_A += [a1, a2]; lo += [h, h]
            cons_A.append(a1 + a2); lo.append(L+1)
    cons_A.append(np.ones(nv)); lo.append(D)
    Am = np.array(cons_A); lo = np.array(lo, float); hi = np.full(len(lo), np.inf); hi[-1] = D
    res = milp(c=np.zeros(nv), constraints=LinearConstraint(Am, lo, hi), integrality=np.ones(nv),
               bounds=Bounds(0, D))
    return res.status == 0

def maxcode(D, c, r, balance=True, target=None):
    K = (D + 2*c)//2; L = D - 1 - K; h = -(-L//2)
    zero = tuple([sp.Integer(0)]*r)
    basis = [tuple(sp.Integer(1) if i == j else sp.Integer(0) for i in range(r)) for j in range(r)]
    base_rows = [zero] + basis
    if not feasible(base_rows, r, D, L, h, balance): return 0, None
    lam = [l for l in Lam(r) if feasible(base_rows + [l], r, D, L, h, balance)]
    n = len(lam)
    comp = np.zeros((n, n), bool)
    for i in range(n):
        for j in range(i+1, n):
            comp[i, j] = comp[j, i] = feasible(base_rows + [lam[i], lam[j]], r, D, L, h, balance)
    best = [r+1 if True else 0, base_rows]
    if lam: best = [r+2, base_rows + [lam[0]]]
    # orbit representatives under S_r for first extra row
    def canon(l):
        return min(tuple(l[p] for p in perm) for perm in itertools.permutations(range(r)))
    def dfs(chosen, cand):
        size = r + 1 + len(chosen)
        if size > best[0]:
            best[0] = size; best[1] = base_rows + [lam[i] for i in chosen]
        if target and best[0] >= target: return
        for t, i in enumerate(cand):
            if size + 1 + (len(cand) - t - 1) <= best[0]: return
            pass
            rows = base_rows + [lam[k] for k in chosen] + [lam[i]]
            if len(chosen) >= 1 and not feasible(rows, r, D, L, h, balance): continue
            dfs(chosen + [i], [k for k in cand[t+1:] if comp[i, k]])
    dfs([], list(range(n)))
    return best[0], best[1], n

if __name__ == "__main__":
    D, c, r = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    bal = (sys.argv[4] == "1") if len(sys.argv) > 4 else True
    out = maxcode(D, c, r, bal)
    print("D", D, "c", c, "rank<=", r, "balance", bal, "max size", out[0], "|Lambda feasible|", out[2] if len(out) > 2 else None)
    print("rows:", [tuple(str(x) for x in row) for row in out[1]])
