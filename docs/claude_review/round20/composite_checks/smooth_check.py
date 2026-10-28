"""Check 2 (double precision; evidence for Prop E, not used in any proof):
smooth all-ideal Hecke sums F_k(X) = sum_alpha f(N alpha / X) e^{4 i k theta(alpha)}, f a C^infty bump on (1/2, 2).
Genuine: O(X^-B) (Poisson / entire L). Multiplicative graft and minimal graft: O(1)-size oscillation (poles on Re s = 0).
Near-genuine grafts phi = theta + eps: size proportional to eps (quantitative pinning)."""
import math, numpy as np
from common import ideals, grafts, thetas
from sympy import primerange

def f(u):
    out = np.zeros_like(u)
    m = (u > 0.5) & (u < 2.0)
    out[m] = np.exp(-1.0 / ((u[m] - 0.5) * (2.0 - u[m])))
    return out

Xmax = 2_000_000
x, y, n = ideals(2 * Xmax)
S = [p for p in primerange(101, 2000) if p % 4 == 1][:6]
rng = np.random.default_rng(7)
phi = rng.uniform(0, math.pi / 2, size=len(S))
th, mult, mini, smooth, spart = grafts(x, y, n, S, phi)
print("S =", S, " #ideals", len(n), " #S-smooth nontrivial", int((smooth & (n > 1)).sum()))
nf = n.astype(np.float64)
for k in (1, 2, 5):
    print(f"k={k}")
    for X in (20000, 100000, 500000, 1000000, 2000000):
        w = f(nf / X)
        g = np.sum(w * np.exp(4j * k * th)); m_ = np.sum(w * np.exp(4j * k * mult)); mi = np.sum(w * np.exp(4j * k * mini))
        print(f"  X={X:8d} |genuine|={abs(g):.3e} |mult graft|={abs(m_):.3e} |minimal graft|={abs(mi):.3e}  (sum of weights {w.sum():.3e})")
# pinning: phi = theta + eps * (random signs), multiplicative graft
print("pinning test, k=1, X=1e6: phi_p = theta_p + eps*u_p, |u_p|<=1")
u = rng.uniform(-1, 1, size=len(S))
X = 1000000; w = f(nf / X)
for eps in (1e-1, 1e-2, 1e-3, 1e-4):
    _, mlt, _, _, _ = grafts(x, y, n, S, thetas(S) + eps * u)
    print(f"  eps={eps:.0e}  |F_1 graft| = {abs(np.sum(w*np.exp(4j*mlt))):.3e}")
