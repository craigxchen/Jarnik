"""Referee check of p3Dlower Def. 4.1 / Lemma 4.2 / Thm 4.3 at ALL levels a <= A_q (author checked level 1 only).
Exact arithmetic in T_Q = norm-one units of Z[i]/Q, Q = q^A_q <= X, with discrete logs by table (independent code).
Twisted labels are defined at the top level by Def. 4.1; for every v in [-2,2]^5 minus 0, every q, every level a, every s:
    <v,ell(q)> = s*iota mod m_{q^a}   <=>   q^a | D_s(Z_v rho_t(q)^(d_q(v))),
and m_(v,s) is compared with the cap (4.2)."""
import itertools, math
from sympy import primerange, legendre_symbol, factorint
from tq import two_sq, D
X = 200
def gmul(x, y, Q): return ((x[0]*y[0] - x[1]*y[1]) % Q, (x[0]*y[1] + x[1]*y[0]) % Q)
def gpow(x, e, Q):
    r = (1, 0)
    while e:
        if e & 1: r = gmul(r, x, Q)
        x = gmul(x, x, Q); e >>= 1
    return r
def ginv(x, Q):
    n = (x[0]*x[0] + x[1]*x[1]) % Q; ni = pow(n, -1, Q); return (x[0]*ni % Q, (-x[1])*ni % Q)
class TQ:
    def __init__(self, q, A):
        self.q, self.Q = q, q**A; Q = self.Q
        self.m = q**(A-1) * (q - 1 if q % 4 == 1 else q + 1)
        primes = list(factorint(self.m)); gen = None
        for a in range(Q):
            for b in range(1, Q):
                if (a*a + b*b) % q == 0: continue
                t = gmul((a, b), ginv((a, (-b) % Q), Q), Q)
                if all(gpow(t, self.m // p, Q) != (1, 0) for p in primes): gen = t; break
            if gen: break
        self.tab = {}; x = (1, 0)
        for e in range(self.m): self.tab[x] = e; x = gmul(x, gen, Q)
        assert len(self.tab) == self.m
        self.iota = self.tab[(0, 1)]
    def L(self, a, b):
        Q = self.Q; return self.tab[gmul((a % Q, b % Q), ginv((a % Q, (-b) % Q), Q), Q)]
split = [p for p in primerange(5, 10**5) if p % 4 == 1]
zp = [p for p in split if p > X][:5]; zs = [two_sq(p) for p in zp]           # zeta_j over p'_j > X
genp = [p for p in primerange(10**6, 10**6 + 2000) if p % 4 == 1][:5]                                    # genuine p_j ~ 1e6
ok = True; ncoll = 0; worst = 0.0; capviol = 0
qs = list(primerange(3, X + 1))
data = {}
for q in qs:
    A = int(math.floor(math.log(X) / math.log(q) + 1e-12))
    while q**(A+1) <= X: A += 1
    t = TQ(q, A)
    nu = [1 if legendre_symbol(p % q, q) == -1 else 0 for p in genp]
    coh = [t.L(a, b) for a, b in zs]
    tq_ = next(i for i, l in enumerate([2] + split) if legendre_symbol(l % q, q) == -1)
    ell_t = ([2] + split)[tq_]; rho = (1, 1) if tq_ == 0 else two_sq(ell_t)
    tw = t.L(*rho)
    delta = [(coh[j] % 2) ^ nu[j] for j in range(5)]
    lab = [(coh[j] + delta[j] * tw) % t.m for j in range(5)]
    ok &= all(lab[j] % 2 == nu[j] for j in range(5))
    data[q] = (t, A, lab, delta, tq_, rho, ell_t)
def gpow_int(z, e):           # exact Gaussian integer power (e >= 0)
    r = (1, 0)
    for _ in range(e): r = (r[0]*z[0] - r[1]*z[1], r[0]*z[1] + r[1]*z[0])
    return r
def conj(z): return (z[0], -z[1])
def gm(x, y): return (x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0])
T_ = max(d[4] for d in data.values()); lam = math.log(max(d[6] for d in data.values()))
for v in itertools.product(range(-2, 3), repeat=5):
    if not any(v): continue
    k = sum(abs(x) for x in v)
    Z = (1, 0)
    for j, x in enumerate(v): Z = gm(Z, gpow_int(zs[j] if x > 0 else conj(zs[j]), abs(x)))
    logm = [0.0] * 4
    for q, (t, A, lab, delta, tq_, rho, ell_t) in data.items():
        d = sum(x * dl for x, dl in zip(v, delta))
        W = gm(Z, gpow_int(rho if d >= 0 else conj(rho), abs(d)))
        val = sum(x * l for x, l in zip(v, lab)) % t.m
        for a in range(1, A + 1):
            ma = q**(a-1) * (q - 1 if q % 4 == 1 else q + 1)
            for s in range(4):
                lhs = (val - s * t.iota) % ma == 0
                Ds = D(W[0], W[1], s); ok &= (Ds != 0)
                rhs = Ds % (q**a) == 0
                ok &= (lhs == rhs)
                if lhs: ncoll += 1; logm[s] += math.log(q) if True else 0
    absZ = math.sqrt(Z[0]**2 + Z[1]**2)
    cap = (T_ * (2*k + 1) + 4) * math.log(math.sqrt(2) * absZ) + T_ * k * (k + 1) / 2 * lam
    for s in range(4):
        if logm[s] > cap + 1e-9: capviol += 1
        worst = max(worst, logm[s] / cap)
print(f"X={X}, T(X)={T_}, levels up to A_q (3^4,5^3,7^2,...,13^2 included)")
print("Lemma 4.2 at all levels (dlog collision <=> q^a | D_s(W)), D_s != 0, parity of twisted labels:", ok, " level-collisions:", ncoll)
print("Thm 4.3 cap violations:", capviol, "  max m_(v,s)/cap in log scale:", round(worst, 3))
print("ALL OK" if ok and capviol == 0 else "FAIL")
