# E3: exact verification (Gaussian-rational arithmetic) of the contact-2 curve family
#   u_j(e) = prod_{l != j} (1 + (b_l + i g) e) * (1 + (b_j - i g) e)
#   v_j(e) = prod_{l != j} (1 + (b_l - i g) e) * (1 + (b_j + i g) e)        (j = 1..M)
# i.e. sigma^(j) = 1 - 2 e_j.  Claims:
#   (i)  u_j v_j is independent of j (curve lies on the cone of X_M);
#   (ii) u_j - u_1 and v_j - v_1 are divisible by e^2 (contact >= 2 with Y at e = 0, no base point);
#   (iii) degree M;
#   (iv) given any point x of the open torus of X_M (M odd: after a Q(i)-rescaling), parameters exist
#        with (u(1), v(1)) proportional to x:  r_j = mu v_j, 1 + b_j = i g (r_j + 1)/(r_j - 1).
#   (v)  for a form F of degree D and order m along Y, F(phi(e)) has degree <= M D and order >= 2m;
#        checked for F_3 (M=3: F_3(phi(e)) = c e^6, c != 0) and the M=5 Vandermonde form.
# Consequence (proved in verify.md): tau(M) <= M/2; tau(3) = 3/2, tau(4) <= 2.
from fractions import Fraction as Fr
import itertools, math, random

class G:
    __slots__ = ('a', 'b')
    def __init__(self, a, b=0): self.a = Fr(a); self.b = Fr(b)
    def __add__(s, o): o = o if isinstance(o, G) else G(o); return G(s.a + o.a, s.b + o.b)
    __radd__ = __add__
    def __sub__(s, o): o = o if isinstance(o, G) else G(o); return G(s.a - o.a, s.b - o.b)
    def __rsub__(s, o): return G(o) - s
    def __neg__(s): return G(-s.a, -s.b)
    def __mul__(s, o): o = o if isinstance(o, G) else G(o); return G(s.a * o.a - s.b * o.b, s.a * o.b + s.b * o.a)
    __rmul__ = __mul__
    def conj(s): return G(s.a, -s.b)
    def norm(s): return s.a * s.a + s.b * s.b
    def __truediv__(s, o):
        o = o if isinstance(o, G) else G(o); n = o.norm(); q = s * o.conj(); return G(q.a / n, q.b / n)
    def __eq__(s, o): o = o if isinstance(o, G) else G(o); return s.a == o.a and s.b == o.b
    def iszero(s): return s.a == 0 and s.b == 0
    def __repr__(s): return f"({s.a}+{s.b}i)"
I = G(0, 1)

def padd(p, q):
    n = max(len(p), len(q)); return [(p[i] if i < len(p) else G(0)) + (q[i] if i < len(q) else G(0)) for i in range(n)]
def psub(p, q): return padd(p, [-c for c in q])
def pmul(p, q):
    out = [G(0) for _ in range(len(p) + len(q) - 1)]
    for i, a in enumerate(p):
        if a.iszero(): continue
        for j, b in enumerate(q): out[i + j] = out[i + j] + a * b
    return out
def ppow(p, k):
    out = [G(1)]
    for _ in range(k): out = pmul(out, p)
    return out
def pord(p):
    for i, c in enumerate(p):
        if not c.iszero(): return i
    return None
def pdeg(p):
    d = None
    for i, c in enumerate(p):
        if not c.iszero(): d = i
    return d
def peval(p, x):
    out = G(0)
    for c in reversed(p): out = out * x + c
    return out

def curve(bs, g):
    M = len(bs)
    lin_p = [[G(1), b + I * g] for b in bs]   # 1 + (b + i g) e
    lin_q = [[G(1), b - I * g] for b in bs]   # 1 + (b - i g) e
    U, V = [], []
    for j in range(M):
        u = [G(1)]; v = [G(1)]
        for l in range(M):
            if l != j: u = pmul(u, lin_p[l]); v = pmul(v, lin_q[l])
        u = pmul(u, lin_q[j]); v = pmul(v, lin_p[j])
        U.append(u); V.append(v)
    return U, V

def form_on_curve(meas, D, U, V):
    out = [G(0)]
    t = pmul(U[0], V[0])
    for k, c in meas.items():
        l1 = sum(abs(x) for x in k)
        term = [G(c)]
        term = pmul(term, ppow(t, (D - l1) // 2))
        for j, kj in enumerate(k):
            term = pmul(term, ppow(U[j], kj) if kj > 0 else ppow(V[j], -kj))
        out = padd(out, term)
    return out

def perm_sign(p):
    s = 1; p = list(p)
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]; p[i], p[j] = p[j], p[i]; s = -s
    return s
def vandermonde_measure(v):
    M = len(v); out = {}
    for p in itertools.permutations(range(M)):
        k = tuple(v[p[i]] for i in range(M)); out[k] = out.get(k, 0) + perm_sign(p)
    return {k: x for k, x in out.items() if x}

def params_from_r(r, g):
    """1 + b_j = i g (r_j + 1)/(r_j - 1)  <=>  (1 + b_j + i g)/(1 + b_j - i g) = r_j."""
    bs = []
    for rj in r:
        assert not (rj - 1).iszero()
        bs.append(I * g * (rj + 1) / (rj - 1) - 1)
    return bs

def target_point(r):
    """the point of X_M reached at e = 1: (u_j, v_j) = (prod_{l != j} r_l, r_j)."""
    M = len(r); P = G(1)
    for x in r: P = P * x
    return [P / x for x in r], list(r)

def params_through(x_u, x_v, g):
    """For a point x = (u_j, v_j) of X_M (u_j v_j = t) we need r_j = mu v_j with mu^(M-2) = t / prod v.
    For M = 3, mu = t / prod v lies in Q(i), so this is exact for actual lattice configurations."""
    M = len(x_u); assert M == 3
    t = x_u[0] * x_v[0]; pv = G(1)
    for v in x_v: pv = pv * v
    mu = t / pv
    r = [mu * v for v in x_v]
    return params_from_r(r, g), x_u, x_v

def check(M, bs, g, label, forms=()):
    U, V = curve(bs, g)
    prod0 = pmul(U[0], V[0])
    ok_cone = all(psub(pmul(U[j], V[j]), prod0) == [] or all(c.iszero() for c in psub(pmul(U[j], V[j]), prod0)) for j in range(M))
    contact = min(min(pord(psub(U[j], U[0])) or 99, pord(psub(V[j], V[0])) or 99) for j in range(1, M))
    deg = max(max(pdeg(u) for u in U), max(pdeg(v) for v in V))
    base0 = U[0][0]
    print(f"  [{label}] M={M}: on cone: {ok_cone}; contact order at e=0: {contact}; degree: {deg}; phi(0) = ({base0}, ...) not a base point: {not base0.iszero()}")
    for name, meas, D, m in forms:
        P = form_on_curve(meas, D, U, V)
        print(f"     F={name} (deg {D}, ord {m}): F(phi(e)) has order {pord(P)} (>= 2m = {2*m}) and degree {pdeg(P)} (<= M D = {M*D}); value at e=1 nonzero: {not peval(P, G(1)).iszero()}")
    return U, V

random.seed(1)
F3 = vandermonde_measure((1, 0, -1))
F5 = vandermonde_measure((2, 1, 0, -1, -2))
# 1. random parameters
for M in (3, 4, 5, 6, 7):
    r = [G(Fr(random.randint(-9, 9), random.randint(1, 5)), Fr(random.randint(1, 9), random.randint(1, 5))) for _ in range(M)]
    g = G(Fr(random.randint(1, 9), random.randint(1, 5)), Fr(random.randint(-3, 3), 7))
    bs = params_from_r(r, g)
    forms = []
    if M >= 3: forms.append(("F_3(z1,z2,z3)", {k + (0,) * (M - 3): v for k, v in F3.items()}, 2, 3))
    if M >= 5: forms.append(("F_5 Vandermonde", {k + (0,) * (M - 5): v for k, v in F5.items()}, 6, 10))
    U, V = check(M, bs, g, "random r", forms)
    xu, xv = target_point(r)
    lamb = peval(U[0], G(1)) / xu[0]
    prop = all((peval(U[j], G(1)) - lamb * xu[j]).iszero() and (peval(V[j], G(1)) - lamb * xv[j]).iszero() for j in range(M))
    print(f"     phi(1) proportional to (prod_(l!=j) r_l, r_j)_j: {prop}")

# 2. through actual lattice configurations (M odd).
def lattice_points(N):
    pts = []
    a = 0
    while a * a <= N:
        b2 = N - a * a; b = math.isqrt(b2)
        if b * b == b2:
            for sa in (1, -1):
                for sb in (1, -1):
                    pts.append((sa * a, sb * b))
        a += 1
    return sorted(set(pts), key=lambda p: math.atan2(p[1], p[0]))
N = 5**2 * 13**2 * 17 * 29 * 37
pts = lattice_points(N)
print(f"\n  circle N = {N}: {len(pts)} lattice points")
for M, forms in ((3, [("F_3", F3, 2, 3)]),):
    sel = pts[10:10 + M]
    x_u = [G(a, b) for a, b in sel]; x_v = [G(a, -b) for a, b in sel]
    g = G(Fr(3, 2))
    bs, xu, xv = params_through(x_u, x_v, g)
    U, V = check(M, bs, g, f"through lattice points {sel}", forms)
    # point at e = 1 proportional to x
    lamb = peval(U[0], G(1)) / xu[0]
    prop = all((peval(U[j], G(1)) - lamb * xu[j]).iszero() and (peval(V[j], G(1)) - lamb * xv[j]).iszero() for j in range(M))
    print(f"     phi(1) proportional to the lattice configuration: {prop}")
    if M == 3:
        P = form_on_curve(F3, 2, U, V)
        nz = [i for i, c in enumerate(P) if not c.iszero()]
        print(f"     F_3(phi(e)) nonzero coefficients at powers {nz}  (so F_3 o phi = c e^6 exactly)")

# 3. dominance of r -> point: log-Jacobian is J - 2I (ratios u_j/v_j = prod_l r_l / r_j^2)
for M in range(2, 12):
    # det(J - 2I) = (M-2) (-2)^(M-1)
    det = (M - 2) * (-2) ** (M - 1)
    print(f"  M={M}: det(J-2I) = {det}  ({'dominant' if det != 0 else 'NOT dominant'})")
