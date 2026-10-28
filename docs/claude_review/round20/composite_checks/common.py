"""Shared helpers: ideals of Z[i] (first-quadrant reps), genuine angles, S-valuations, grafts."""
import math, numpy as np

def ideals(X):
    """first-quadrant representatives x>0, y>=0 of nonzero ideals with norm <= X."""
    r = math.isqrt(X)
    xs, ys = np.meshgrid(np.arange(1, r + 1, dtype=np.int64), np.arange(0, r + 1, dtype=np.int64), indexing='ij')
    x, y = xs.ravel(), ys.ravel()
    n = x * x + y * y
    keep = n <= X
    return x[keep], y[keep], n[keep]

def gauss_prime(p):
    """a+bi with a>b>0, a^2+b^2=p (p = 1 mod 4)."""
    for b in range(1, math.isqrt(p) + 1):
        a2 = p - b * b
        a = math.isqrt(a2)
        if a * a == a2 and a > b:
            return a, b
    raise ValueError(p)

def val(x, y, a, b, p):
    """valuation of x+iy at a+bi (vectorised)."""
    x = x.copy(); y = y.copy(); v = np.zeros(x.shape, dtype=np.int64)
    active = np.ones(x.shape, dtype=bool)
    while active.any():
        re = x * a + y * b; im = y * a - x * b
        ok = active & (re % p == 0) & (im % p == 0)
        if not ok.any(): break
        x = np.where(ok, re // p, x); y = np.where(ok, im // p, y)
        v += ok; active = ok
    return v

def s_data(x, y, n, S):
    """for each p in S: exponent difference v_pi - v_pibar; and S-part norm."""
    diffs = []; spart = np.ones(n.shape, dtype=np.float64)
    for p in S:
        a, b = gauss_prime(p)
        v1 = val(x, y, a, b, p); v2 = val(x, y, a, -b, p)
        diffs.append(v1 - v2); spart *= float(p) ** (v1 + v2)
    return np.array(diffs), spart

def thetas(S):
    return np.array([math.atan2(gauss_prime(p)[1], gauss_prime(p)[0]) for p in S])

def grafts(x, y, n, S, phi):
    """returns genuine angle, multiplicative graft, minimal graft (all mod pi/2)."""
    th = np.arctan2(y, x).astype(np.float64)
    D, spart = s_data(x, y, n, S)
    shift = (D * (np.asarray(phi) - thetas(S))[:, None]).sum(axis=0)
    mult = np.mod(th + shift, math.pi / 2)
    smooth = np.isclose(spart, n.astype(np.float64))
    minimal = np.where(smooth, mult, th)
    return th, mult, minimal, smooth, spart
