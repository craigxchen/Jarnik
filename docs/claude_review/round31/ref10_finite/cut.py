# Referee's independent check of the X.5 least-cutoff table (exact integer arithmetic).
import sys
from math import gcd, isqrt
from sympy import primerange, isprime

def gprime(p):
    for a in range(1, isqrt(p) + 1):
        b2 = p - a * a
        b = isqrt(b2)
        if b * b == b2:
            return (a, b)

def mul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])

def conj(z):
    return (z[0], -z[1])

def canon(z):
    # representative of z up to the 4 units and conjugation
    x, y = z
    imgs = []
    for w in [(x, y), (-y, x), (-x, -y), (y, -x)]:
        imgs.append(w); imgs.append(conj(w))
    return min(imgs)

def labels(P, five):
    return [p for p in primerange(P, 2 * P + 1) if (p % 8 == 5 if five else p % 4 == 1)]

def Lcm(X, S):
    L = 1
    for q in primerange(3, X + 1):
        if q in S:
            continue
        qa = q
        while qa * q <= X:
            qa *= q
        L *= qa
    return L

def enum_lifts(S, H):
    """Independent enumeration: all r = pi*gamma/conj(gamma), Norm pi = p in S, pi either prime
    of norm p, gamma in Z[i] nonzero with Norm gamma <= H.  Reduce r = alpha/b, b>0 minimal.
    Keep b>1 with every prime of b in S (local-unit at all exterior primes <= X when X >= H).
    Return set of (p, b, canon(alpha))."""
    out = set()
    allb = set()
    Sset = set(S)
    SPF = list(range(H + 1))
    for q in range(2, isqrt(H) + 1):
        if SPF[q] == q:
            for k in range(q * q, H + 1, q):
                if SPF[k] == k:
                    SPF[k] = q
    R = isqrt(H)
    for p in S:
        a, bb = gprime(p)
        for pi in [(a, bb), (a, -bb)]:
            for u in range(-R, R + 1):
                for v in range(-R, R + 1):
                    Ng = u * u + v * v
                    if Ng == 0 or Ng > H:
                        continue
                    g2 = mul((u, v), (u, v))
                    num = mul(pi, g2)
                    d = gcd(gcd(num[0], num[1]), Ng)
                    alpha = (num[0] // d, num[1] // d)
                    b = Ng // d
                    assert alpha[0] ** 2 + alpha[1] ** 2 == p * b * b
                    assert gcd(alpha[0], alpha[1]) == 1
                    allb.add(b)
                    if b == 1:
                        assert alpha[0]**2 + alpha[1]**2 == p
                        continue
                    # primes of b
                    bad = False
                    t = b
                    while t > 1:
                        q = SPF[t]
                        if q not in Sset:
                            bad = True
                        while t % q == 0:
                            t //= q
                    if not bad:
                        out.add((p, b, canon(alpha)))
    return out, allb

def feasible_n(alpha, p, n_pow, L):
    """alpha^n given as n_pow=(x,y); p^n. Exists phi with all four inequalities at this n?"""
    pass

def single_n_kill(z, pn, L):
    x, y = z
    m0, m1, m2, m3 = gcd(y, L), gcd(x - y, L), gcd(x, L), gcd(x + y, L)
    # c = p^n cos(2n phi), d = p^n sin(2n phi)
    cl, ch = 2 * m2 * m2 - pn, pn - 2 * m0 * m0
    dl, dh = m3 * m3 - pn, pn - m1 * m1
    if cl > ch or dl > dh:
        return True
    near = (max(0, cl, -ch)) ** 2 + (max(0, dl, -dh)) ** 2
    far = max(cl * cl, ch * ch) + max(dl * dl, dh * dh)
    return not (near <= pn * pn <= far)

def killed(alpha, p, L, nmax):
    z = (1, 0); pn = 1
    for n in range(1, nmax + 1):
        z = mul(z, alpha); pn *= p
        if single_n_kill(z, pn, L):
            return n
    return 0

def survives_with(alpha, p, L, w):
    """phi = arg w, Norm w = 5 exterior.  Check (1.2) at all four units for all n
    exactly; beyond n0 with (p/5)^n >= L^2 it holds automatically."""
    z = (1, 0); t = (1, 0); pn = 1; fn = 1; n = 0
    while pn < fn * L * L:
        n += 1
        z = mul(z, alpha); t = mul(t, w); pn *= p; fn *= 5
        x, y = z; a, b = t
        for D, E in [(y, b), (x - y, a - b), (x, a), (x + y, a + b)]:
            # nu_s^2 E^2/5^n >= nu_s^2 m^2/p^n  <=>  E^2 p^n >= m^2 5^n
            m = gcd(D, L)
            if E * E * pn < m * m * fn:
                return None
    return n

if __name__ == "__main__":
    P = int(sys.argv[1]); five = sys.argv[2] == "5"; X = int(sys.argv[3]); nmax = int(sys.argv[4])
    S = labels(P, five)
    s = len(S)
    H = 2 * P
    lifts, allb = enum_lifts(S, H)
    print("P", P, "five", five, "labels", s, "nonintegral local-unit lifts up to sym:", len(lifts), "s(2s-1)=", s * (2 * s - 1))
    print("denominators b:", sorted(set(b for _, b, _ in lifts)) == sorted(S))
    T = max(((l - 1) & -(l - 1)).bit_length() - 1 for l in S)
    print("T =", T)
    LX = Lcm(X, set(S)); LX1 = Lcm(X - 1, set(S))
    print("L(X) > L(X-1):", LX > LX1)
    # genuine primes pass at X (check up to nmax only via single-n)
    for p in S:
        assert killed(gprime(p), p, LX, nmax) == 0
    notkilled = []
    from collections import Counter
    cnt = Counter()
    for (p, b, al) in lifts:
        k = killed(al, p, LX, nmax)
        cnt[k] += 1
        if k == 0:
            notkilled.append((p, b, al))
    print("at X =", X, "kill exponent distribution:", sorted(cnt.items()))
    # survivors at X-1
    surv = []
    for (p, b, al) in lifts:
        if killed(al, p, LX1, nmax) == 0:
            for w in [(2, 1), (2, -1), (1, 2), (1, -2)]:
                for a2 in [w8 for u8 in [al,(-al[1],al[0]),(-al[0],-al[1]),(al[1],-al[0])] for w8 in (u8,conj(u8))]:
                    r = survives_with(a2, p, LX1, w)
                    if r is not None:
                        surv.append((p, b, a2, w, r))
                        break
                else:
                    continue
                break
    print("at X-1 =", X - 1, "survivors with phi = arg(2+-i) etc:", len(surv), surv[:5])
