# (a),(b),(c): symbolic identities over F_2[x1..x5,Y]; Jacobian ranks at random points of F_{2^20}.
import itertools, random, sys
from gf import *
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
n = 6  # variables x0..x4 (indices 0..4), Y (index 5)
X = [mp_var(i, n) for i in range(5)]; Yv = mp_var(5, n)

def span_nonzero(idx):
    out = []
    for coeffs in itertools.product([0, 1], repeat=len(idx)):
        if any(coeffs):
            s = set()
            for c, i in zip(coeffs, idx):
                if c: s = mp_add(s, X[i])
            out.append(s)
    return out

def prod(ps):
    r = mp_const1(n)
    for p in ps: r = mp_mul(r, p)
    return r

# y_j = product of the 15 nonzero elements of H_j = span(x_k, k != j)
ys = []
for j in range(5):
    idx = [k for k in range(5) if k != j]
    ys.append(prod(span_nonzero(idx)))
# Moore determinant check: det (x_k^(2^a)) a=0..3, k in idx, Leibniz over F_2
def moore(idx):
    s = set()
    for perm in itertools.permutations(range(4)):
        term = mp_const1(n)
        for a, kpos in enumerate(perm):
            term = mp_mul(term, mp_pow2(X[idx[kpos]], a))
        s = mp_add(s, term)
    return s
for j in range(5):
    idx = [k for k in range(5) if k != j]
    assert moore(idx) == ys[j], "Moore det != product of 15 forms"
print("(b) y_j = 4x4 Moore minor deleting column j: verified symbolically for all j")
# kernel identities sum_j y_j x_j^(2^a) = 0, a = 0..3
for a in range(4):
    s = set()
    for j in range(5): s = mp_add(s, mp_mul(ys[j], mp_pow2(X[j], a)))
    assert not s
print("(b) sum_j y_j x_j^(2^a) = 0 for a=0,1,2,3: verified symbolically")
# a=4 is NOT zero (sanity)
s = set()
for j in range(5): s = mp_add(s, mp_mul(ys[j], mp_pow2(X[j], 4)))
print("    sanity a=4: nonzero ->", bool(s), " (it is the 5x5 Moore det = prod of 31 forms):",
      s == prod(span_nonzero(range(5))))
# (c) derivative identities sum_j x_j^(2^a) dy_j/dx_k = 0 for a=1,2,3
J = [[mp_diff(ys[j], k) for k in range(5)] for j in range(5)]
for a in (1, 2, 3):
    for k in range(5):
        s = set()
        for j in range(5): s = mp_add(s, mp_mul(mp_pow2(X[j], a), J[j][k]))
        assert not s
print("(c) sum_j x_j^(2^a) d y_j/d x_k = 0 for a=1,2,3, all k: verified symbolically")

# (a) F_j(Y) = P_j(Y)/Y = prod_{v in H_j \ 0} (Y+v): linearized support check, deg_Y(F_i - F_j) <= 7
Fs = []
for j in range(5):
    idx = [k for k in range(5) if k != j]
    Fs.append(prod([mp_add(Yv, v) for v in span_nonzero(idx)]))
for j in range(5):
    degsY = sorted({e[5] for e in Fs[j]})
    assert set(degsY) <= {0, 1, 3, 7, 15}, degsY
print("(a) each F_j has Y-exponents in {0,1,3,7,15} (additivity of P_j): verified")
mx = 0
for i in range(5):
    for j in range(i + 1, 5):
        d = max(e[5] for e in mp_add(Fs[i], Fs[j]))
        mx = max(mx, d)
print("(a) max_{i<j} deg_Y(F_i - F_j) =", mx, " against deg_Y F_j = 15")
# term count of F_j in the split model (Y=t, x constants): number of distinct Y exponents
print("    split model monomials in t per mark:", [len({e[5] for e in F}) for F in Fs])

# Jacobian ranks at random points over F_{2^20}
F = GF2k(20)
def proj_rank(cols_fn, vals):
    pass
dF = [[mp_diff(Fs[j], k) for k in range(6)] for j in range(5)]
for trial in range(4):
    pt = [F.rand() for _ in range(6)]
    # x-only Jacobian of y at Y=0: projective rank = rank [J | y] - 1
    pt0 = pt[:5] + [0]
    Jy = [[mp_eval(F, J[j][k], pt0) for k in range(5)] + [mp_eval(F, ys[j], pt0)] for j in range(5)]
    r1 = rank_gf(F, Jy) - 1
    JyonlyX = rank_gf(F, [row[:5] for row in Jy])
    # full family (x, Y) -> [F_1(Y):...:F_5(Y)]
    Jf = [[mp_eval(F, dF[j][k], pt) for k in range(6)] + [mp_eval(F, Fs[j], pt)] for j in range(5)]
    r2 = rank_gf(F, Jf) - 1
    rX = rank_gf(F, [row[:5] + [row[6]] for row in Jf]) - 1
    rY = rank_gf(F, [[row[5], row[6]] for row in Jf]) - 1
    print(f"trial {trial}: Y=0 slice: proj. rank d[y] = {r1} (affine rank of J = {JyonlyX});"
          f" full family (x,Y): proj. rank = {r2}; x-directions only: {rX}; Y-direction only: {rY}")
