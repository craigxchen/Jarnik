# Li-Liu residual kernel E_N(s,sbar) at real s=sigma for N=2p (p odd prime).
# E = P - D - K ; P = (2/phi(N)) sum_{chi odd} |L'/L(1/2+sigma,chi)|^2
# D(w) = sum_{(n,N)=1} Lambda(n)^2 n^{-1-w}, w=2 sigma ; K = S(N)(pi/N)^w P_N(w) zeta(1-w) X1(sigma)^2
# Lower bound (5.8): ||T_N|| >= 8 d(s) |E(s,sbar)| / |1+s|^2, d(s)=min(sigma-alpha_N, 1/8-sigma).
import sys, math, numpy as np
import flint
from flint import arb, acb
import mpmath as mp
from sympy import primitive_root, primerange, isprime
flint.ctx.prec = 90
mp.mp.dps = 30

def hurwitz_and_deriv(s, a):
    h = arb("1e-8")
    f0 = acb(s).zeta(acb(a))
    fp = acb(s + h).zeta(acb(a)); fm = acb(s - h).zeta(acb(a))
    return float(f0.real), float(((fp - fm) / (2 * h)).real)

C2 = None
def twin_const():
    # C2 = prod_{p>2} (1-1/(p-1)^2)
    return float(mp.mpf('0.66016181584686957392781211001455577843262336028473'))

def prime_sum_log2(w, exclude):
    # sum_{p not in exclude} log^2 p /(p^{1+w}-1) = sum_k sum_p log^2p p^{-k(1+w)}
    tot = mp.mpf(0)
    k = 1
    while True:
        s = k * (1 + w)
        if s > 60: break
        term = mp.diff(lambda x: mp.primezeta(x), s, 2)
        tot += term
        k += 1
    for q in exclude:
        tot -= mp.log(q) ** 2 / (mp.mpf(q) ** (1 + w) - 1)
    return float(tot)

def logPN(w, exclude, B=2*10**6):
    # log prod_{p>2, p not in exclude} (1 + p^{-w}/(p-2)); direct to B, tail via p^{-1-w}+2p^{-2-w} using primezeta
    tot = 0.0
    ps = np.array(list(primerange(3, B)), dtype=float)
    tot = float(np.sum(np.log1p(ps ** (-w) / (ps - 2))))
    # tail: sum_{p>=B} [p^{-1-w} + 2 p^{-2-w} - p^{-2-2w}/2 ...] approx by primezeta differences
    P1 = mp.primezeta(1 + w) - sum(mp.mpf(q) ** (-1 - w) for q in [2]) 
    head = float(np.sum(ps ** (-1 - w))) 
    tail = float(P1) - head  # sum_{p>=B} p^{-1-w}
    tot += tail
    for q in exclude:
        if q > 2:
            tot -= math.log1p(q ** (-w) / (q - 2))
    return tot

def E_at(p, sigma):
    N = 2 * p
    g = primitive_root(p)
    # discrete logs mod p
    dlog = [0] * p
    x = 1
    for m in range(p - 1):
        dlog[x] = m; x = x * g % p
    s = 0.5 + sigma
    # F(a) for a odd in 1..N-1, a != p : contributes chi_p(a mod p)
    A = np.zeros(p - 1, dtype=complex); Ad = np.zeros(p - 1, dtype=complex)
    for a in range(1, N, 2):
        if a == p: continue
        z0, z1 = hurwitz_and_deriv(s, a / N)
        A[dlog[a % p]] += z0; Ad[dlog[a % p]] += z1
    # chi_j(g^m)=e(jm/(p-1)); sum_m chi_j(g^m) A[m] = fft with conj sign
    FA = np.fft.ifft(A) * (p - 1); FAd = np.fft.ifft(Ad) * (p - 1)
    j = np.arange(p - 1)
    odd = (j % 2 == 1)
    L = N ** (-s) * FA; Ld = N ** (-s) * (FAd - math.log(N) * FA)
    gvals = Ld[odd] / L[odd]
    phi = p - 1
    P = 2.0 / phi * float(np.sum(np.abs(gvals) ** 2))
    w = 2 * sigma
    D = prime_sum_log2(w, [2, p])
    S = 2 * twin_const() * (p - 1) / (p - 2)
    lP = logPN(w, [2, p])
    X1 = float(mp.gamma((1.5 - sigma) / 2) / mp.gamma((1.5 + sigma) / 2))
    K = S * (math.pi / N) ** w * math.exp(lP) * float(mp.zeta(1 - w)) * X1 ** 2
    E = P - D - K
    alpha = 1 / math.log(N)
    d = min(sigma - alpha, 0.125 - sigma) if sigma > alpha else float("nan")
    lb = 8 * d * abs(E) / (1 + sigma) ** 2
    delta = math.sqrt(S / N) * math.log(N) ** 1.5
    return dict(N=N, sigma=sigma, P=P, D=D, K=K, E=E, lb=lb, delta=delta)

if __name__ == "__main__":
    for p in [int(x) for x in sys.argv[1:]]:
        N = 2 * p
        for sigma in [1 / 9, 0.05, 1.5 / math.log(N)]:
            if sigma >= 0.125: continue
            r = E_at(p, sigma)
            print(" ".join(f"{k}={v:.6g}" for k, v in r.items()), flush=True)
