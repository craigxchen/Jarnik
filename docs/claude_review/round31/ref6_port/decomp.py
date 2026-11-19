"""(1) The multiple-decomposition identity in characteristic 3 versus characteristic 0.

char 3: for V = i<1,t>_{F_27} and the four kernel lines H of x0, x1, x0+x1, x0-x1, check that
A_H^q - c_H^{q-1} A_H (A_H = X^q - w^{q-1} X, w spanning H, c_H = A_H(v0)) is one and the same polynomial
P_V(X) in F_729[t][X], i.e. P_V has four different right composition factors of the same degree q.
char 0: sympy check that z_g = A (A^s - g)^2, conj(A) = -A, has norm independent of g only for g = +-g0,
and of the decomposition lemma on an explicit example.
"""
import sys
sys.path.insert(0, '.')
from tern4 import f729_mul, f729_add, f729_neg, f729_pow, enc, I, polmul729, p27, M27, A27, N27

q = 27
def padd(a, b):
    n = max(len(a), len(b)); r = [0]*n
    for k in range(n):
        x = a[k] if k < len(a) else 0; y = b[k] if k < len(b) else 0
        r[k] = f729_add(x, y)
    while len(r) > 1 and r[-1] == 0: r.pop()
    return r
def pneg(a): return [f729_neg(x) for x in a]
def pfrob(a, e=q):  # (sum c_k t^k)^q = sum c_k^q t^{qk}, char 3, q = 3^3
    r = [0]*((len(a)-1)*e + 1)
    for k, c in enumerate(a):
        r[k*e] = f729_pow(c, e)
    return r
def ppow(a, e):
    r = [1]
    for _ in range(e): r = polmul729(r, a)
    return r

res = []
for (a, b) in [(1, 0), (0, 1), (1, 1), (1, 2)]:
    # kernel of a x0 + b x1 is spanned by (x0,x1) = (-b, a); w = i(-b + a t)
    w = [f729_mul(I, enc(N27[b] if b else 0, 0)), f729_mul(I, enc(a, 0))]
    v0x = (1, 0) if a else (0, 1)   # a vector with form value 1 (coefficient 1)
    v0 = [f729_mul(I, enc(v0x[0], 0)), f729_mul(I, enc(v0x[1], 0))]
    wq1 = ppow(w, q - 1)
    # c = v0^q - w^{q-1} v0
    c = padd(pfrob(v0), pneg(polmul729(wq1, v0)))
    cq1 = ppow(c, q - 1)
    # A^q - c^{q-1} A with A = X^q - wq1 X:  X^{q^2} + (-(wq1)^q - cq1) X^q + (cq1*wq1) X
    coefXq = padd(pneg(pfrob(wq1)), pneg(cq1))
    coefX = polmul729(cq1, wq1)
    tr=lambda p: (lambda r: r)(p[:max([k+1 for k,x in enumerate(p) if x] or [1])])
    coefXq, coefX = tr(coefXq), tr(coefX)
    res.append((coefXq, coefX))
    print("form", (a, b), "deg_t coeff X^q:", len(coefXq) - 1, " deg_t coeff X:", len(coefX) - 1)
print("char 3: the four decompositions give the same P_V:", all(r == res[0] for r in res))

# Dickson-invariant check: P_V(v) = 0 for all v in V
coefXq, coefX = res[0]
def peval_lin(v):  # v in F_729[t] (list); returns P_V(v) = v^{q^2} + coefXq v^q + coefX v
    return padd(padd(pfrob(pfrob(v)), polmul729(coefXq, pfrob(v))), polmul729(coefX, v))
import itertools, random
random.seed(1)
ok = True
for _ in range(30):
    x0, x1 = random.randrange(27), random.randrange(27)
    v = [f729_mul(I, enc(x0, 0)), f729_mul(I, enc(x1, 0))]
    ok &= (peval_lin(v) == [0])
print("P_V vanishes on 30 random elements of V:", ok)

# ---------- characteristic 0 ----------
import sympy as sp
a, g, g0, X = sp.symbols('a g g0 X')
s = 3
A = sp.I * a
z = A * (A**s - g)**2
zb = (-A) * ((-A)**s - sp.conjugate(g))**2
# with g real: norm
nz = sp.expand(z * zb.subs(sp.conjugate(g), g))
print("char 0, s=3, g real: norm =", sp.factor(nz))
# norm as a function of complex g = g0*zeta: equal to the g0 norm iff zeta real
zeta = sp.symbols('zeta')
zz = A * (A**s - g0*zeta)**2
zzb = (-A) * ((-A)**s - g0*sp.conjugate(zeta))**2
diff = sp.expand(zz * zzb - nz.subs(g, g0))
print("norm(g0*zeta) - norm(g0) =", sp.factor(diff))
# decomposition lemma example: f(Y)=Y^7 - C Y, A = X^3 + X + 1; search A' = X^3 + p X^2 + r X + u, f'(Y) = Y^7 + ...
# (generic statement is proved in the report; here we only check that f(alpha A + beta) = f(A) forces beta=0, alpha^7=1, alpha*C'=C)
al, be, C, Cp, Y = sp.symbols('alpha beta C Cp Y')
expr = sp.expand((al*Y + be)**7 - Cp*(al*Y + be) - (Y**7 - C*Y))
eqs = [sp.Eq(expr.coeff(Y, k), 0) for k in range(8)]
print("f(alpha Y + beta) = f(Y) solutions:", sp.solve(eqs, [al, be, Cp], dict=True)[:3], "...")
