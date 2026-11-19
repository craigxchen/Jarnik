"""Residual gcd enumeration (check of (V.4) and of Astra's d=4 count).

Method E (Euclid): DFS over forms, choice skip/kernel/complement, running Euclidean gcd in F_4[t].
Method P (valuations): factor every Y+v over F_4 by trial division; gcd degree =
   sum_f deg f * min_i sum_{v in S_i} ord_f(Y+v).
residual = deg gcd - sum_{v in C} deg(Y+v),  C = intersection of the chosen sides.
"""
import itertools, random, sys, time
import numpy as np
from f4poly import *

random.seed(5)


def monics(k):
    for coeffs in itertools.product(range(4), repeat=k):
        yield tuple(coeffs) + (1,)


def factor(p):
    """trial division; returns dict monic irreducible -> exponent"""
    out = {}
    p = trim(p)
    if deg(p) <= 0:
        return out
    k = 1
    while 2 * k <= deg(p):
        for m in monics(k):
            while True:
                q, r = divmod_(p, m)
                if r:
                    break
                out[m] = out.get(m, 0) + 1
                p = q
            if 2 * k > deg(p):
                break
        k += 1
    if deg(p) >= 1:
        c = INV[p[-1]]
        p = trim([MUL[c][x] for x in p])
        out[p] = out.get(p, 0) + 1
    return out


def prime_data(d, T):
    Q = 1 << d
    facs = [factor(linfac(T, v)) for v in range(Q)]
    for v in range(Q):  # self-check of factorization
        lf = linfac(T, v)
        pr = prod(f for f, e in facs[v].items() for _ in range(e))
        assert deg(pr) == max(0, deg(lf)), (v, lf, pr)
    primes = sorted({f for F in facs for f in F})
    ordm = np.array([[facs[v].get(f, 0) for v in range(Q)] for f in primes], dtype=np.int32)
    degs = np.array([deg(f) for f in primes], dtype=np.int64)
    return primes, ordm, degs, [max(0, deg(linfac(T, v))) for v in range(Q)]


def bound_unsimplified(d):
    """(Q/2) sum_{e even, 2<=e<d} e*pi_2(e) * sum_{j>=1, ej<d} 2^{-ej}  with exact pi_2"""
    Q = 1 << d
    def pi2(e):
        # Moebius
        def mu(n):
            r, x, p = 1, n, 2
            while p * p <= x:
                if x % p == 0:
                    x //= p
                    if x % p == 0:
                        return 0
                    r = -r
                p += 1
            return -r if x > 1 else r
        return sum(mu(e // k) * 2 ** k for k in range(1, e + 1) if e % k == 0) // e
    tot = 0.0
    for e in range(2, d, 2):
        tot += e * pi2(e) * sum(2.0 ** (-e * j) for j in range(1, d) if e * j < d)
    return Q / 2 * tot


def method_E(d, T):
    Q = 1 << d
    fs = forms_and_sides(d)
    fdeg = [max(0, deg(linfac(T, v))) for v in range(Q)]
    sides = []
    for (w, k, c) in fs:
        sides.append([(S, prod(linfac(T, v) for v in sorted(S))) for S in (k, c)])
    full = frozenset(range(Q))
    stats = {"cases": 0, "peak": 0, "Vgap": 10**9}

    def visit(i, cnt, C, G):
        if i == len(sides):
            if cnt >= 2:
                stats["cases"] += 1
                res = deg(G) - sum(fdeg[v] for v in C)
                stats["peak"] = max(stats["peak"], res)
                B = max(d, f2deg(T) if T else -1)
                stats["Vgap"] = min(stats["Vgap"], len(C) * B + Q * d / 4 - deg(G))
            return
        visit(i + 1, cnt, C, G)
        for S, F in sides[i]:
            visit(i + 1, cnt + 1, C & S, gcd(F, G) if G else F)

    visit(0, 0, full, ())
    return stats


def method_P(d, T, sample=None):
    """exhaustive (numpy meet-in-the-middle) or sampled valuation computation"""
    Q = 1 << d
    fs = forms_and_sides(d)
    primes, ordm, degs, fdeg = prime_data(d, T)
    B = max(d, f2deg(T) if T else -1)
    nf = len(fs)
    # per side: valuation vector and mask
    side_ord, side_mask = [], []
    for (w, k, c) in fs:
        row_o, row_m = [], []
        for S in (k, c):
            idx = sorted(S)
            row_o.append(ordm[:, idx].sum(axis=1) if len(primes) else np.zeros(0, np.int32))
            row_m.append(sum(1 << v for v in S))
        side_ord.append(row_o)
        side_mask.append(row_m)
    fdeg_arr = np.array(fdeg, dtype=np.int64)
    BIG = 10**6
    kP = len(primes)

    def half_table(forms):
        cfgs = list(itertools.product(range(3), repeat=len(forms)))
        mins = np.full((len(cfgs), kP), BIG, dtype=np.int64)
        masks = np.zeros(len(cfgs), dtype=np.int64)
        cnts = np.zeros(len(cfgs), dtype=np.int64)
        for n, cfg in enumerate(cfgs):
            m = (1 << Q) - 1
            mv = np.full(kP, BIG, dtype=np.int64)
            c = 0
            for f, ch in zip(forms, cfg):
                if ch:
                    mv = np.minimum(mv, side_ord[f][ch - 1])
                    m &= side_mask[f][ch - 1]
                    c += 1
            mins[n], masks[n], cnts[n] = mv, m, c
        return mins, masks, cnts

    def cdeg(masks):
        # sum of fdeg over bits of mask, via byte tables
        out = np.zeros(masks.shape, dtype=np.int64)
        for byte in range((Q + 7) // 8):
            tab = np.array([sum(fdeg[8 * byte + b] for b in range(8) if (x >> b) & 1 and 8 * byte + b < Q)
                            for x in range(256)], dtype=np.int64)
            out += tab[(masks >> (8 * byte)) & 255]
        return out

    def popc(masks):
        out = np.zeros(masks.shape, dtype=np.int64)
        tab = np.array([bin(x).count("1") for x in range(256)], dtype=np.int64)
        for byte in range((Q + 7) // 8):
            out += tab[(masks >> (8 * byte)) & 255]
        return out

    assert Q <= 62  # masks fit in int64
    if sample is None:
        h1 = list(range(nf // 2))
        h2 = list(range(nf // 2, nf))
        m1, k1, c1 = half_table(h1)
        m2, k2, c2 = half_table(h2)
        peak, cases, gap = 0, 0, 10**9
        hist = {}
        for a in range(len(c1)):
            mv = np.minimum(m1[a], m2)
            cnt = c1[a] + c2
            ok = cnt >= 2
            if not ok.any():
                continue
            mask = k1[a] & k2
            g = (mv * degs).sum(axis=1) if kP else np.zeros(len(cnt), np.int64)
            res = g - cdeg(mask)
            res, g, mask = res[ok], g[ok], mask[ok]
            cases += int(ok.sum())
            peak = max(peak, int(res.max()))
            gap = min(gap, float((popc(mask) * B + Q * d / 4 - g).min()))
            for r in np.unique(res):
                hist[int(r)] = hist.get(int(r), 0) + int((res == r).sum())
        return {"cases": cases, "peak": peak, "Vgap": gap, "hist": hist, "nprimes": kP}
    else:
        peak = 0
        recs = []
        for _ in range(sample):
            s = random.randint(2, nf)
            sel = random.sample(range(nf), s)
            ch = [random.randint(0, 1) for _ in sel]
            mv = np.full(kP, BIG, dtype=np.int64)
            m = (1 << Q) - 1
            for f, c in zip(sel, ch):
                mv = np.minimum(mv, side_ord[f][c])
                m &= side_mask[f][c]
            g = int((mv * degs).sum())
            res = g - sum(fdeg[v] for v in range(Q) if (m >> v) & 1)
            peak = max(peak, res)
            recs.append((sel, ch, g, res))
        return {"peak": peak, "recs": recs}


def euclid_selection(d, T, sel, ch):
    fs = forms_and_sides(d)
    G = ()
    for f, c in zip(sel, ch):
        S = fs[f][1 + c]
        F = prod(linfac(T, v) for v in sorted(S))
        G = gcd(F, G) if G else F
    return deg(G)


if __name__ == "__main__":
    what = sys.argv[1]
    if what == "d4":
        d = 4
        print("unsimplified prime bound d=4:", bound_unsimplified(4), " Qd/4 =", 16 * 4 / 4)
        tot = 0
        pk = 0
        for T in (0, 16, 32, 48):
            t0 = time.time()
            sE = method_E(d, T)
            sP = method_P(d, T)
            assert sE["cases"] == sP["cases"] and sE["peak"] == sP["peak"], (sE, sP)
            tot += sE["cases"]
            pk = max(pk, sE["peak"])
            print(f"d=4 T={T}: cases {sE['cases']}, peak residual E={sE['peak']} P={sP['peak']}, "
                  f"min slack in g<=hB+Qd/4: {sE['Vgap']}, hist {dict(sorted(sP['hist'].items()))}  ({time.time()-t0:.1f}s)")
        print(f"d=4 TOTAL cases {tot}, max residual {pk}")
        # extension: all T of degree < 10 modulo V (64 classes), valuation method only
        pk2 = 0
        for k in range(64):
            T = 16 * k
            sP = method_P(d, T)
            pk2 = max(pk2, sP["peak"])
            assert sP["Vgap"] >= 0
        print(f"d=4, all 64 classes T mod V with deg T < 10: max residual {pk2}")
        # Euclid cross-check on a few extra classes
        for T in (16 * 13, 16 * 37, 16 * 63):
            assert method_E(d, T)["peak"] == method_P(d, T)["peak"]
        print("Euclid = valuation on T classes 13,37,63 (x16)")
    elif what == "d5":
        d = 5
        print("unsimplified prime bound d=5:", bound_unsimplified(5), " Qd/4 =", 32 * 5 / 4)
        pk = 0
        for k in range(int(sys.argv[2])):
            T = 32 * k
            t0 = time.time()
            sP = method_P(d, T)
            pk = max(pk, sP["peak"])
            print(f"d=5 T={T} (deg {f2deg(T) if T else '-inf'}): {sP['nprimes']} primes, cases {sP['cases']}, "
                  f"peak residual {sP['peak']}, min slack g<=hB+Qd/4 {sP['Vgap']}, "
                  f"hist {dict(sorted(sP['hist'].items()))} ({time.time()-t0:.0f}s)", flush=True)
        print("d=5 max residual over classes:", pk)
    elif what == "d5euclid":
        d = 5
        n = int(sys.argv[2])
        for T in (0, 32, 64, 96):
            r = method_P(d, T, sample=n)
            for sel, ch, g, res in r["recs"]:
                assert euclid_selection(d, T, sel, ch) == g
            print(f"d=5 T={T}: {n} random selections, Euclid gcd degree == valuation gcd degree; sample peak {r['peak']}")
    elif what == "d6sample":
        d = 6
        n = int(sys.argv[2])
        print("unsimplified prime bound d=6:", bound_unsimplified(6), " Qd/4 =", 64 * 6 / 4)
        for T in (0, 64, 128, 192, 64 * 5):
            r = method_P(d, T, sample=n)
            print(f"d=6 T={T}: sample {n}, peak residual {r['peak']}")
