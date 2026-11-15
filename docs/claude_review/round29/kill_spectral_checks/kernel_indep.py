# Independent recomputation of Li-Liu E_N(sigma,sigma)=P-D-K for N=2p.
# P via mpmath Hurwitz zeta with exact derivative (derivative=1), char sums by explicit DFT over discrete logs.
# D via prime zeta from Moebius-log-zeta (P(s)=sum mu(k)/k log zeta(ks)), differentiated analytically by mpmath diff on log zeta.
# K via (2.12) with P_N computed as primezeta part + direct correction.
import sys, math, numpy as np, mpmath as mp
from sympy import primitive_root, primerange, mobius
mp.mp.dps = 25
def Pdd(s):  # sum_p log^2 p p^{-s} = sum_k mu(k) k (log zeta)''(ks)
    t = mp.mpf(0); k = 1
    while k*s < 80:
        mu = mobius(k)
        if mu: t += mu*k*mp.diff(lambda x: mp.log(mp.zeta(x)), k*s, 2)
        k += 1
    return t
def run(p, sigma):
    N = 2*p; s = mp.mpf(0.5)+sigma
    g = primitive_root(p); dl = {}; x = 1
    for m in range(p-1): dl[x] = m; x = x*g % p
    A = np.zeros(p-1, complex); Ad = np.zeros(p-1, complex)
    for a in range(1, N, 2):
        if a == p: continue
        A[dl[a % p]] += float(mp.zeta(s, mp.mpf(a)/N)); Ad[dl[a % p]] += float(mp.zeta(s, mp.mpf(a)/N, 1))
    m = np.arange(p-1); P = 0.0
    for j in range(1, p-1, 2):   # odd characters chi_j(g^m)=exp(2 pi i j m/(p-1)); chi(-1)=(-1)^j
        e = np.exp(2j*np.pi*j*m/(p-1))
        L = np.dot(e, A); Ld = np.dot(e, Ad) - math.log(N)*L   # common factor N^{-s} cancels in L'/L
        P += abs(Ld/L)**2
    P *= 2.0/(p-1)
    w = 2*sigma
    D = Pdd(1+w)
    for q in (2, p): D -= mp.log(q)**2/(mp.mpf(q)**(1+w)-1)
    C2 = mp.mpf('0.6601618158468695739278121100145557784326')
    S = 2*C2*(p-1)/mp.mpf(p-2)
    lp = mp.primezeta(1+w) - mp.mpf(2)**(-1-w) - mp.mpf(p)**(-1-w)
    for q in primerange(3, 300000):
        if q == p: continue
        lp += math.log1p(q**(-w)/(q-2)) - q**(-1-w)
    PN = mp.e**lp
    X1 = mp.gamma((1.5-sigma)/2)/mp.gamma((1.5+sigma)/2)
    Km = S*(mp.pi/N)**w*PN*mp.zeta(1-w)*X1**2
    E = P - D - Km
    delta = mp.sqrt(S/N)*mp.log(N)**1.5
    print(f"N={N} sigma={sigma:.6f} P={P:.6f} D={float(D):.6f} K={float(Km):.6f} E={float(E):.6f} delta={float(delta):.4f}", flush=True)
for p in map(int, sys.argv[1:]):
    for sg in (1/9, 0.05): run(p, sg)
