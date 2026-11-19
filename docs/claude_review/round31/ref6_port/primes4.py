"""Prime degrees (over F_9) of the blocks Q_O of the four-mark ternary example, by distinct-degree
factorisation, and the induced primes of F_3[t] (all split)."""
import sys, itertools
sys.path.insert(0, '.')
import numpy as np
from collections import Counter
import tern4 as T

def ddf(Q):
    Q = T.trim(Q)
    out = Counter()
    tpol = (np.array([0, 1]), np.array([0, 0]))
    h = tpol
    k = 0
    f = Q
    while T.deg9(f) >= 2 * (k + 1) or (T.deg9(f) > 0 and k < 40):
        k += 1
        # h = h^9 mod Q
        for _ in range(2):
            hh = h
            h3 = T.mul9(T.mul9(hh, hh), hh)
            h = T.mod9(h3, Q)
        diff = T.sub9(h, tpol)
        gg = T.gcd9(f, diff) if T.deg9(diff) >= 0 else f
        dg = T.deg9(gg)
        if dg > 0:
            out[k] += dg // k
            # remove
            while True:
                gg2 = T.gcd9(f, gg)
                if T.deg9(gg2) <= 0: break
                # exact division f / gg2 via long division
                f = divexact(f, gg2)
        if T.deg9(f) <= 0: break
    return out

def divexact(a, b):
    a = [a[0].copy(), a[1].copy()]; b = T.trim(b)
    db = len(b[0]) - 1; da = len(a[0]) - 1
    inv = T.F9INV[(int(b[0][db]), int(b[1][db]))]
    qre = np.zeros(da - db + 1, dtype=np.int64); qim = np.zeros(da - db + 1, dtype=np.int64)
    while da >= db:
        lr, li = int(a[0][da]), int(a[1][da])
        if lr or li:
            cr = (lr * inv[0] - li * inv[1]) % 3; ci = (lr * inv[1] + li * inv[0]) % 3
            s = da - db; qre[s], qim[s] = cr, ci
            a[0][s:da + 1] = (a[0][s:da + 1] - (cr * b[0] - ci * b[1])) % 3
            a[1][s:da + 1] = (a[1][s:da + 1] - (cr * b[1] + ci * b[0])) % 3
        da -= 1
    assert not a[0][:db].any() and not a[1][:db].any()
    return T.trim((qre, qim))

def main():
    q, L = 27, 11
    vecs = list(itertools.product(range(27), repeat=2))
    seen, orbits = set(), []
    for x in vecs:
        if x in seen or x == (0, 0): continue
        orb = [x]; y = x
        while True:
            y = (T.p27(y[0], 9), T.p27(y[1], 9))
            if y == x: break
            orb.append(y)
        seen.update(orb); orbits.append(orb)
    tot = Counter()
    for orb in orbits:
        poly = [1]
        for (x0, x1) in orb:
            lin = [0] * (L + 1); lin[L] = 1
            lin[0] = T.f729_neg(T.f729_mul(T.I, T.enc(x0, 0)))
            lin[1] = T.f729_add(lin[1], T.f729_neg(T.f729_mul(T.I, T.enc(x1, 0))))
            poly = T.polmul729(poly, lin)
        Q = T.to_f9(poly)
        if orb[0][0] == 0:  # strip the factor t^{|O|}
            Q = (Q[0][len(orb):], Q[1][len(orb):])
        c = ddf(Q)
        tot.update(c)
    print("F_9-prime degree histogram over all nonzero labels (both members of each conjugate pair):", sorted(tot.items()))
    print("smallest prime degree:", min(tot))

if __name__ == "__main__":
    main()
