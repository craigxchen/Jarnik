"""Check 1 (numerical linear algebra, double precision; illustration/evidence):
all-ideal Hecke large-sieve optimal constants Lambda(X,K) = lambda_max of the Toeplitz matrix T(k-l),
T(m) = sum_{N alpha <= X} e^{4 i m theta(alpha)}, for genuine angles, multiplicative graft, minimal graft.
Also checks the proved lower bound Lambda >= max(#ideals, (2K+1) floor(sqrt X)) and the proved
upper bound Lambda_all(Y,K) <= 56 (Y + K Y^(1/2))."""
import math, sys, numpy as np
from scipy.linalg import toeplitz
from common import ideals, grafts
from sympy import primerange

def lam(theta, K):
    m = np.arange(0, 2 * K + 1)
    T = np.array([np.exp(4j * mm * theta).sum() for mm in m])
    H = toeplitz(T, np.conj(T))   # H[k,l] = T(k-l) with T(-m) = conj T(m)
    return np.linalg.eigvalsh(H).max()

rng = np.random.default_rng(20)
for (X, P0, s) in [(20000, 101, 4), (100000, 211, 6), (300000, 401, 8)]:
    S = [p for p in primerange(P0, 10 * P0) if p % 4 == 1][:s]
    x, y, n = ideals(X)
    phi = rng.uniform(0, math.pi / 2, size=s)
    th, mult, mini, smooth, spart = grafts(x, y, n, S, phi)
    nS = int((spart > 1).sum()); nsm = int((smooth & (n > 1)).sum())
    tau_half = np.prod([(1 - p ** -0.5) ** -2 for p in S]) - 1
    print(f"X={X} S={S} #ideals={len(n)} #S-divisible={nS} #S-smooth(nontrivial)={nsm} tau_1/2={tau_half:.4f}")
    for K in (5, 20, 60, 150):
        lg, lm, ln = lam(th, K), lam(mult, K), lam(mini, K)
        lb = max(len(n), (2 * K + 1) * math.isqrt(X))
        ub = 56 * (X + K * math.sqrt(X))
        print(f"  K={K:4d} genuine={lg:10.1f} mult={lm:10.1f} ({lm/lg-1:+.4f}) minimal={ln:10.1f} ({ln/lg-1:+.4f}) "
              f"proved_lb={lb} lb_ok={lg >= lb - 1e-6} proved_ub_ok={lg <= ub} X+K*sqrtX={X + K*math.sqrt(X):.0f}")
    sys.stdout.flush()
