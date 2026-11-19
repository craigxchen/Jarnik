"""Fixed integer T, fixed seeds, n -> infinity: exact Gaussian integers, primitive height, Delta*H^(1/4)."""
from fractions import Fraction as Fr
from math import log, exp
import sys

def matpow(T, n):
    A = ((T*T+1, T), (T, 1))
    R = ((1, 0), (0, 1))
    def mm(X, Y):
        return ((X[0][0]*Y[0][0]+X[0][1]*Y[1][0], X[0][0]*Y[0][1]+X[0][1]*Y[1][1]),
                (X[1][0]*Y[0][0]+X[1][1]*Y[1][0], X[1][0]*Y[0][1]+X[1][1]*Y[1][1]))
    for _ in range(n):
        R = mm(R, A)
    return R

def mul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def cj(a): return (a[0], -a[1])

def ggcd(a, b):
    def mod(x, y):
        nr = x[0]*y[0]+x[1]*y[1]; ni = x[1]*y[0]-x[0]*y[1]; d = y[0]**2+y[1]**2
        q = ((2*nr+d)//(2*d), (2*ni+d)//(2*d))
        return (x[0]-(q[0]*y[0]-q[1]*y[1]), x[1]-(q[0]*y[1]+q[1]*y[0]))
    while b != (0, 0):
        a, b = b, mod(a, b)
    return a

def marks(T, n, seeds):
    # seeds integers (pass already-cleared integer seed vectors (p,q) meaning (q, p) start vector q*(1,p/q))
    M = matpow(T, n)
    G = []
    for (p, q) in seeds:
        x = M[0][0]*q + M[0][1]*p; y = M[1][0]*q + M[1][1]*p
        G.append((x, y))
    beta = (T, 2)
    def prodc(idx):
        r = (1, 0)
        for l in idx: r = mul(r, cj(G[l]))
        return r
    z = [mul(beta, prodc(range(4)))]
    for (j, k) in [(0, 1), (0, 2), (1, 2), (0, 3)]:
        z.append(mul(mul(cj(beta), mul(G[j], G[k])), prodc([l for l in range(4) if l not in (j, k)])))
    return z

def stats(z):
    g = z[0]
    for v in z[1:]: g = ggcd(g, v)
    gn = g[0]**2+g[1]**2
    R2 = z[0][0]**2+z[0][1]**2
    assert all(v[0]**2+v[1]**2 == R2 for v in z)
    H = R2//gn
    d2 = [(z[a][0]-z[b][0])**2+(z[a][1]-z[b][1])**2 for a in range(5) for b in range(a+1, 5)]
    assert min(d2) > 0
    # Delta = diam / R ; log(Delta*H^(1/4))
    lDH = 0.5*(log(max(d2))-log(R2)) + 0.25*log(H)
    gam = log(H)/(-0.5*(log(max(d2))-log(R2)))
    return gn, H, exp(lDH), gam

if __name__ == '__main__':
    for T in (1, 2, 3):
        for seeds in ([(1, 1), (2, 1), (3, 1), (5, 1)], [(1, 2), (-2, 1), (7, 3), (4, 1)]):
            print(f'T={T} seeds={seeds}')
            for n in (2, 5, 10, 20, 40, 80, 160):
                gn, H, DH, gam = stats(marks(T, n, seeds))
                print(f'   n={n:4d} N(gcd)={gn} log10H={log(H)/log(10):8.1f} Delta*H^(1/4)={DH:.8f} gamma={gam:.6f}')
