# all 30 translation classes of 8-point sections (Y = 7 finite points + infinity) of Paley-12, L=2, m=8:
# exact dim V(A_Y,8), exact sign change of the top form, irreducibility over Q.
import itertools, random, flint, sympy as sp, time
q = 11; D0 = [0,1,3,4,5,9]; L = 2
cols = [sorted((i + x) % q for x in D0) for i in range(q)]
reps = sorted({min(tuple(sorted((y + t) % q for y in Y)) for t in range(q)) for Y in itertools.combinations(range(q), 7)})
print("classes", len(reps), flush=True)
mons = sorted([e for e in itertools.product((0,1,2), repeat=7)], key=sum)
def ev(k, e):
    v = 1
    for a, b in zip(k, e):
        if b: v *= a**b
    return v
random.seed(5); us = sp.symbols("u0:7"); summary = []
for Y in reps:
    t0 = time.time()
    A = [k for k in itertools.product((-1,0,1), repeat=7)
         if all(abs(sum(dict(zip(Y, k)).get(x, 0) for x in c)) <= L for c in cols)]
    M8 = flint.fmpz_mat([[ev(k, e) for e in mons if sum(e) < 8] for k in A])
    M9r = flint.fmpz_mat([[ev(k, e) for e in mons if sum(e) < 9] for k in A]).rank()
    X, nul = M8.transpose().nullspace()
    c = [int(X[i, 0]) for i in range(len(A))]
    def Phi(u): return sum(ci * sum(a*b for a, b in zip(k, u))**8 for ci, k in zip(c, A) if ci)
    pos = neg = 0
    for _ in range(400):
        v = Phi([random.randint(-30, 30) for _ in range(7)]); pos += v > 0; neg += v < 0
        if pos and neg: break
    Pp = sp.Poly(sum(ci * sum(a*b for a, b in zip(k, us))**8 for ci, k in zip(c, A) if ci), *us)
    fl = sp.factor_list(Pp)
    irr = len(fl[1]) == 1 and fl[1][0][1] == 1 and fl[1][0][0].total_degree() == 8
    summary.append((Y, len(A), nul, M9r == len(A), bool(pos and neg), irr))
    print(Y, "|A_Y|", len(A), "dimV8", nul, "rho<=8", M9r == len(A), "indefinite(exact)", bool(pos and neg), "irreducible/Q", irr, f"{time.time()-t0:.0f}s", flush=True)
print("all classes: dimV8=1:", all(s[2] == 1 for s in summary), "rho=8:", all(s[3] for s in summary),
      "indefinite:", sum(s[4] for s in summary), "irreducible:", sum(s[5] for s in summary), "of", len(summary))
