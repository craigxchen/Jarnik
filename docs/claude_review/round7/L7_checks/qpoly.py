"""Exact univariate polynomials over Q (coefficient lists, low degree first, Fractions).
Small, dependency-free helper used by the L7 round-7 scripts."""
from fractions import Fraction as Fr


def norm(p):
    p = [Fr(c) for c in p]
    while p and p[-1] == 0:
        p.pop()
    return p


def deg(p):
    p = norm(p)
    return len(p) - 1  # -1 for the zero polynomial


def add(p, q):
    n = max(len(p), len(q))
    return norm([(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)])


def sub(p, q):
    return add(p, [-c for c in q])


def scal(c, p):
    return norm([Fr(c) * x for x in p])


def mul(p, q):
    p, q = norm(p), norm(q)
    if not p or not q:
        return []
    r = [Fr(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a == 0:
            continue
        for j, b in enumerate(q):
            r[i + j] += a * b
    return norm(r)


def prod(ps):
    r = [Fr(1)]
    for p in ps:
        r = mul(r, p)
    return r


def lin(root):
    """s - root"""
    return [Fr(-root), Fr(1)]


def ev(p, x):
    r = Fr(0)
    for c in reversed(p):
        r = r * x + c
    return r


def divmod_(p, q):
    p, q = norm(p), norm(q)
    if not q:
        raise ZeroDivisionError
    p = p[:]
    out = [Fr(0)] * max(len(p) - len(q) + 1, 1)
    while len(p) >= len(q) and p:
        c = p[-1] / q[-1]
        k = len(p) - len(q)
        out[k] = c
        for i, b in enumerate(q):
            p[i + k] -= c * b
        p = norm(p)
    return norm(out), p


def gcd(p, q):
    p, q = norm(p), norm(q)
    while q:
        _, r = divmod_(p, q)
        p, q = q, r
    if not p:
        return []
    return scal(1 / p[-1], p)


def deriv(p):
    return norm([i * p[i] for i in range(1, len(p))])


def squarefree(p):
    return deg(gcd(p, deriv(p))) == 0


def resultant(p, q):
    """Resultant of p, q (as polynomials of their actual degrees) via Euclid."""
    p, q = norm(p), norm(q)
    if not p or not q:
        return Fr(0)
    dp, dq = len(p) - 1, len(q) - 1
    if dq == 0:
        return q[0] ** dp
    if dp == 0:
        return p[0] ** dq
    # Res(p,q) = (-1)^{dp dq} Res(q,p);  Res(q, p) with p = a q + r : lc(q)^{dp - dr} Res(q, r)
    _, r = divmod_(p, q)
    if not r:
        return Fr(0)
    dr = len(r) - 1
    res_qr = resultant(q, r)
    return (-1) ** (dp * dq) * q[-1] ** (dp - dr) * res_qr


def interp(xs, ys):
    """Lagrange interpolation, exact."""
    n = len(xs)
    res = []
    for i in range(n):
        num = [Fr(1)]
        den = Fr(1)
        for j in range(n):
            if j != i:
                num = mul(num, lin(xs[j]))
                den *= (xs[i] - xs[j])
        res = add(res, scal(Fr(ys[i]) / den, num))
    return res


def rational_roots(p):
    """All rational roots of p (with multiplicity ignored), via rational root theorem on the
    integer-scaled polynomial. Only used for small-degree checks."""
    from math import gcd as igcd
    p = norm(p)
    if not p:
        raise ValueError
    # scale to integers
    L = 1
    for c in p:
        L = L * c.denominator // igcd(L, c.denominator)
    q = [int(c * L) for c in p]
    # strip zero roots
    roots = set()
    while q and q[0] == 0:
        roots.add(Fr(0))
        q = q[1:]
    if len(q) <= 1:
        return roots
    a0, an = abs(q[0]), abs(q[-1])

    def divisors(n):
        ds = set()
        i = 1
        while i * i <= n:
            if n % i == 0:
                ds.add(i)
                ds.add(n // i)
            i += 1
        return ds
    for a in divisors(a0):
        for b in divisors(an):
            for sgn in (1, -1):
                x = Fr(sgn * a, b)
                if ev([Fr(c) for c in q], x) == 0:
                    roots.add(x)
    return roots
