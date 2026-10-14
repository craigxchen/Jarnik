# Judge check J8: BLOCK(t0=5) is false for every j0 and every m0 > 8 j0;
# BLOCK(t0=6, j0=2) is false for every m0 > 16.
# rho upper bounds are rigorous: mod-p rank of the evaluation matrix equal to |Y| implies rational rank |Y|.
import itertools, random, sys
from flint import nmod_mat
P = 2147483629

def monos(nv, deg):
    out = []
    def rec(i, left, cur):
        if i == nv - 1:
            out.append(tuple(cur + [left])); return
        for e in range(left + 1):
            rec(i + 1, left - e, cur + [e])
    for dd in range(deg + 1):
        rec(0, dd, [])
    return out

def rank_eval(pts, deg):
    nv = len(pts[0]); ms = monos(nv, deg)
    rows = []
    for p in pts:
        row = []
        for m in ms:
            v = 1
            for x, e in zip(p, m):
                v = v * pow(x, e, P) % P
            row.append(v % P)
        rows.append(row)
    M = nmod_mat(len(pts), len(ms), [x for r in rows for x in r], P)
    return M.rank(), len(ms)

print("(1) parity box P_j = {|y_i|<=j, y1+y2+y3 even}: rho(P_j) <= 3j")
for j in range(1, 6):
    pts = [y for y in itertools.product(range(-j, j + 1), repeat=3) if sum(y) % 2 == 0]
    r, nm = rank_eval(pts, 3 * j)
    print(f"  j={j} |P_j|={len(pts)} rank_p(deg<={3*j})={r}  full={r==len(pts)}")

print("(2) identity k(S_a)+k(S_b)+k(S_a+b) = 2 k(S_a & S_b) for sum k = 0 (random exact tests, F_2^5)")
pts5 = list(itertools.product([0, 1], repeat=5))
dot = lambda a, x: sum(ai * xi for ai, xi in zip(a, x)) % 2
ok = True
for _ in range(2000):
    k = [random.choice([-1, 0, 1]) for _ in pts5]
    s = sum(k); k[0] -= s  # force sum 0 (k[0] may leave {-1,0,1}; identity is linear, fine)
    a = tuple(random.randint(0, 1) for _ in range(5)); b = tuple(random.randint(0, 1) for _ in range(5))
    if not any(a) or not any(b) or a == b: continue
    ab = tuple((x + y) % 2 for x, y in zip(a, b))
    kS = lambda c: sum(kk for kk, x in zip(k, pts5) if dot(c, x) == 0)
    lhs = kS(a) + kS(b) + kS(ab)
    rhs = 2 * sum(kk for kk, x in zip(k, pts5) if dot(a, x) == 0 and dot(b, x) == 0)
    ok &= (lhs == rhs)
print("  identity holds on all samples:", ok)

print("(3) collision-freeness of span{1_S_c : c in F} + constants, F = {e1,e2,e1+e2, e3,e4,e3+e4, e5} on F_2^5")
from fractions import Fraction
import sympy
E = [tuple(1 if i == t else 0 for i in range(5)) for t in range(5)]
add = lambda a, b: tuple((x + y) % 2 for x, y in zip(a, b))
F = [E[0], E[1], add(E[0], E[1]), E[2], E[3], add(E[2], E[3]), E[4]]
cols = [[1 if dot(c, x) == 0 else 0 for x in pts5] for c in F] + [[1] * 32]
target = [x[0] + 2 * x[1] + 4 * x[2] + 8 * x[3] + 16 * x[4] for x in pts5]
Mx = sympy.Matrix(cols).T
sol = Mx.gauss_jordan_solve(sympy.Matrix(target))[0]
print("  binary code (32 distinct values) lies in span exactly:", Mx * sol == sympy.Matrix(target),
      " distinct:", len(set(target)) == 32)

print("(4) Y_3(j) = {(sum_c K_c [a.c=0])_{a!=0} : K in Z^8, sum K = 0} cap [-j,j]^7; rho(Y_3(2)) <= 8")
cls = list(itertools.product([0, 1], repeat=3)); A7 = [a for a in cls if any(a)]
Mmap = sympy.Matrix([[1 if dot(a, c) == 0 else 0 for c in cls] for a in A7])  # 7 x 8
# lattice = image of {K: sum K = 0}; basis K = e_c - e_0
B = sympy.Matrix.hstack(*[Mmap[:, i] - Mmap[:, 0] for i in range(1, 8)])  # 7x7
Binv = B.inv(); print("  det of lattice basis:", B.det())
def inlat(y):
    v = Binv * sympy.Matrix(y)
    return all(x.q == 1 for x in v)
# faster lattice test: integer combination via adjugate
adj = B.adjugate(); D = B.det()
def inlat2(y):
    v = adj * sympy.Matrix(y)
    return all(x % D == 0 for x in v)
for j in (1, 2):
    Y = [y for y in itertools.product(range(-j, j + 1), repeat=7) if inlat2(y)]
    r, nm = rank_eval(Y, 4 * j)
    print(f"  j={j} |Y_3(j)|={len(Y)} rank_p(deg<={4*j})={r} full={r==len(Y)}")
    if j == 2:
        r7, _ = rank_eval(Y, 7)
        print(f"  j=2 rank_p(deg<=7)={r7} (only an upper bound on the rational rank)")
print("(5) two complementary 3-dim blocks on F_2^6: span of 14 indicators + const contains an injective function")
pts6 = list(itertools.product([0, 1], repeat=6))
V1 = [a + (0, 0, 0) for a in A7]; V2 = [(0, 0, 0) + a for a in A7]
cols = [[1 if dot(c, x) == 0 else 0 for x in pts6] for c in V1 + V2] + [[1] * 64]
target = [sum(x[i] << i for i in range(6)) for x in pts6]
Mx = sympy.Matrix(cols).T
sol = Mx.gauss_jordan_solve(sympy.Matrix(target))[0]
print("  6-bit code in span exactly:", Mx * sol == sympy.Matrix(target))
