"""Exact Gaussian-integer helpers (pairs of Python ints)."""
from itertools import combinations
import math

EDGES = list(combinations(range(5), 2))

def mul(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])

def conj(a):
    return (a[0], -a[1])

def norm(a):
    return a[0]*a[0]+a[1]*a[1]

def divround(a, b):
    # nearest Gaussian integer to a/b
    n = norm(b)
    num = mul(a, conj(b))
    def rnd(p, q):
        return (2*p+q)//(2*q)
    return (rnd(num[0], n), rnd(num[1], n))

def sub(a, b):
    return (a[0]-b[0], a[1]-b[1])

def ggcd(a, b):
    while b != (0, 0):
        q = divround(a, b)
        a, b = b, sub(a, mul(q, b))
    return a

def exactdiv(a, b):
    n = norm(b)
    num = mul(a, conj(b))
    assert num[0] % n == 0 and num[1] % n == 0, (a, b)
    return (num[0]//n, num[1]//n)

def prod(lst):
    r = (1, 0)
    for x in lst:
        r = mul(r, x)
    return r

def tuple_from_blocks(G, U=None):
    """z_k = U_k * prod_e (conj G_e if k in e else G_e); returns primitive tuple."""
    z = []
    for k in range(5):
        f = [conj(g) if k in e else g for e, g in zip(EDGES, G)]
        zk = prod(f)
        if U is not None:
            zk = mul(zk, U[k])
        z.append(zk)
    g = z[0]
    for w in z[1:]:
        g = ggcd(g, w)
    z = [exactdiv(w, g) for w in z]
    return z

def width_and_H(z):
    H = norm(z[0])
    assert all(norm(w) == H for w in z)
    m = max(norm(sub(a, b)) for a, b in combinations(z, 2))
    # Delta = sqrt(m/H)
    logD = 0.5*(math.log(m) - math.log(H)) if m > 0 else float('-inf')
    return H, logD

def cofactors(G, U=None):
    """t_ij = [U_i P_i, U_j P_j]/|G_ij|^2, exact integers."""
    if U is None:
        U = [(1, 0)]*5
    P = [mul(U[i], prod([g for e, g in zip(EDGES, G) if i in e])) for i in range(5)]
    t = {}
    for (i, j), g in zip(EDGES, G):
        x = mul(conj(P[i]), P[j])[1]   # Im(conj(P_i) P_j) = [P_i,P_j]
        n = norm(g)
        assert x % n == 0
        t[(i, j)] = x // n
    return t
