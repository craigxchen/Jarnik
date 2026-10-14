# Exact: Paley-12 (L=2) 8-point sections A_Y (Y = 7 finite points + base infinity), m = 8.
import itertools, random, flint, sympy as sp, time, sys
q = 11; D0 = [0,1,3,4,5,9]; L = 2
cols = [sorted((i + x) % q for x in D0) for i in range(q)]
def section(Y):
    pts = []
    for k in itertools.product((-1,0,1), repeat=7):
        kk = dict(zip(Y, k))
        if all(abs(sum(kk.get(x, 0) for x in c)) <= L for c in cols): pts.append(k)
    return pts
mons = sorted([e for e in itertools.product((0,1,2), repeat=7)], key=sum)
def ev(k, e):
    v = 1
    for a, b in zip(k, e):
        if b: v *= a**b
    return v
random.seed(3)
for Y in [(0,1,2,3,4,5,6), (0,1,2,3,4,5,7), (0,1,2,3,4,5,9)]:
    t0 = time.time()
    A = section(Y)
    M8 = flint.fmpz_mat([[ev(k, e) for e in mons if sum(e) < 8] for k in A])
    M9 = flint.fmpz_mat([[ev(k, e) for e in mons if sum(e) < 9] for k in A])
    r8, r9 = M8.rank(), M9.rank()
    X, nul = M8.transpose().nullspace()
    c = [int(X[i, 0]) for i in range(len(A))]
    print(Y, "|A_Y|", len(A), "rank<8", r8, "rank<9", r9, "dim V(A_Y,8)", nul, flush=True)
    # exact sign test at random integer points
    def Phi(u): return sum(ci * sum(a*b for a, b in zip(k, u))**8 for ci, k in zip(c, A) if ci)
    pos = neg = 0; wpos = wneg = None
    for _ in range(3000):
        u = [random.randint(-30, 30) for _ in range(7)]
        v = Phi(u)
        if v > 0: pos += 1; wpos = wpos or u
        elif v < 0: neg += 1; wneg = wneg or u
    print("  exact signs at 3000 random integer points: pos", pos, "neg", neg, "witness+", wpos, "witness-", wneg, flush=True)
    if "--factor" in sys.argv:
        us = sp.symbols("u0:7")
        P = sp.Poly(sum(ci * sum(a*b for a, b in zip(k, us))**8 for ci, k in zip(c, A) if ci), *us)
        fl = sp.factor_list(P)
        print("  factor_list: content", fl[0], "factors (deg, mult):", [(f.total_degree(), m) for f, m in fl[1]], f"t={time.time()-t0:.0f}s", flush=True)
