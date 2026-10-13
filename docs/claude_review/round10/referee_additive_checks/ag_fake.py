"""AG test on explicit skew double core fakes (round5/exact.md Thm 5.1) and on generic arcs.

For a configuration w_x = R exp(i delta_x), delta in [0, Delta]^M, R = e^(W/2), Delta = C e^(-W/4):
  * search, by exact integer LLL (python-flint) on a sweep of weights lambda, for zero-sum integer c with
    0 < |sum c_x w_x| < 1, recording the least ||c||_1 found  (h_found: an UPPER bound for the true least
    AG-violating height h_AG);
  * re-evaluate every candidate with mpmath at 200 digits;
  * for the fake: v(c) = sum c_x a_x, w(v(c)) = sum |v_j| log p_j, and test referee Prop R1(i):
    every c with |sum c w| < 1 and v(c) != 0 has w(v(c)) > W + log 2 - 2 log(1 + n C^3 N^(-1/4)/48);
  * check the sampled delta against the Thm 5.1 requirement for all characters with ||c||_1 <= NMAX
    (dist(c.delta/2, (pi/4)Z) >= arcsin(L^(n/2) e^(-w/2))); all other characters are covered only by
    the theorem's union bound (each random delta is a member with probability >= 0.586);
  * compare h_found with: (a) the author's disc pigeonhole l1 height 2(M-1)H_disc,
    (H+1)^(M-1) > (3(M-1)HL+1)^2; (b) the referee's rectangle pigeonhole 2(M-1)H_rect,
    (H+1)^(M-1) > (3(M-1)H C^2/8 + 1)(3(M-1)HL + 1); (c) the heuristic typical height.
usage: python3 ag_fake.py
"""
import sys, math, random, itertools
import numpy as np
import mpmath as mp
from flint import fmpz_mat
from sympy import isprime

mp.mp.dps = 200


def paley_skew(q):
    chi = lambda a: 0 if a % q == 0 else (1 if pow(a % q, (q - 1) // 2, q) == 1 else -1)
    M = q + 1
    S = np.zeros((M, M), dtype=np.int64)
    for i in range(M):
        for j in range(M):
            if i == j:
                continue
            if i == 0:
                S[i, j] = -1
            elif j == 0:
                S[i, j] = 1
            else:
                S[i, j] = chi((j - 1) - (i - 1))
    K = np.eye(M, dtype=np.int64) + S
    d = K[:, 0] * K[0, 0]
    K = (d[:, None] * K) * d[None, :]
    assert np.all(K + K.T == 2 * np.eye(M, dtype=np.int64)) and np.all(K @ K.T == M * np.eye(M, dtype=np.int64))
    return K


def skew_profile(q):
    K = paley_skew(q)
    M = K.shape[0]
    S = [frozenset(x for x in range(M) if K[x, a] == -1) for a in range(1, M)]
    S2 = [s ^ {0, a} for a, s in zip(range(1, M), S)]
    cutsl = S + S2
    a = [[1 if x in s else 0 for s in cutsl] for x in range(M)]
    return M, a


def primes_1mod4_from(start, count):
    out = []
    n = start + (1 - start) % 4
    while len(out) < count:
        if isprime(n):
            out.append(n)
        n += 4
    return out


def lll_search(ws, lam_list, extra_combos=True):
    """Return list of (l1, c, |sum c w|) for candidates with |sum c w| < 1 found by LLL."""
    M = len(ws)
    V = [ws[x] - ws[0] for x in range(1, M)]
    found = {}
    SIG = mp.mpf(10) ** 30
    for lam in lam_list:
        rows = []
        for i in range(M - 1):
            row = [0] * (M - 1)
            row[i] = int(SIG)
            row += [int(mp.nint(SIG * lam * mp.re(V[i]))), int(mp.nint(SIG * lam * mp.im(V[i])))]
            rows.append(row)
        L = fmpz_mat(rows).lll()
        basis = []
        for i in range(L.nrows()):
            cc = [int(L[i, k]) for k in range(M - 1)]
            # coefficients are SIG * c'
            cc = [int(mp.nint(mp.mpf(t) / SIG)) for t in cc]
            basis.append(cc)
        cands = list(basis)
        if extra_combos:
            for i, j in itertools.combinations(range(len(basis)), 2):
                for s in (1, -1):
                    cands.append([basis[i][k] + s * basis[j][k] for k in range(M - 1)])
        for cp in cands:
            if not any(cp):
                continue
            c = [-sum(cp)] + cp
            val = abs(sum(c[x] * ws[x] for x in range(M)))
            if val < 1:
                key = tuple(c) if c[next(k for k in range(M) if c[k])] > 0 else tuple(-t for t in c)
                found[key] = (sum(abs(t) for t in c), val)
    return sorted(((v[0], list(k), v[1]) for k, v in found.items()), key=lambda t: t[0])


def pigeonhole_heights(M, L, C):
    def least(f):
        lo, hi = 1, 2
        while not f(hi):
            hi *= 2
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if f(mid):
                hi = mid
            else:
                lo = mid
        return hi
    disc = least(lambda H: (H + 1) ** (M - 1) > (3 * (M - 1) * H * L + 1) ** 2)
    rect = least(lambda H: (H + 1) ** (M - 1) > (3 * (M - 1) * H * C * C / 8 + 1) * (3 * (M - 1) * H * L + 1))
    return 2 * (M - 1) * disc, 2 * (M - 1) * rect


def heuristic_typical(M, L, C):
    # expected number of c in Lambda_0 with ||c||_1 <= h and |sum c w| < 1:
    #   ~ #{||c||_1 <= h} * E[1/(|Im-spread| |Re-spread|)];  #{} ~ (2h)^(M-1)/(M-1)!,
    #   Im part ~ R c.delta (spread ~ |c|_2 L / sqrt(12)), Re part ~ R c.delta^2/2 (spread ~ |c|_2 C^2 / sqrt(45)/2);
    #   P(both < 1) ~ 4 / (2pi sigma1 sigma2) * ... use Gaussian approx; |c|_2^2 ~ h^2 * 2/(M(M+1)) for uniform l1-ball.
    def expected(h):
        cnt = (2 * h) ** (M - 1) / math.factorial(M - 1)
        c2 = h * math.sqrt(2.0 / (M * (M + 1)))
        s1 = c2 * float(L) / math.sqrt(12)
        s2 = c2 * float(C) ** 2 / (2 * math.sqrt(45))
        p = 4 / (2 * math.pi * s1 * s2)
        return cnt * p
    lo, hi = 1.0, 2.0
    while expected(hi) < 1:
        hi *= 2
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if expected(mid) < 1:
            lo = mid
        else:
            hi = mid
    return hi


def run_fake(q, C, mode, nsamp, rng, NMAX=6):
    M, a = skew_profile(q)
    r = 2 * (M - 1)
    if mode == 'A0':
        Lres = 1
    else:
        Q0 = int(mode) * M
        Lres = primes_1mod4_from(Q0 + 1, M)[-1]
    Pst = max(math.ceil(160000 * M ** 4 * Lres ** 2 / C ** 2), 21 * r * r)
    ps = primes_1mod4_from(Pst, r)
    W = sum(mp.log(p) for p in ps)
    R = mp.exp(W / 2)
    N4 = mp.exp(W / 4)
    Delta = mp.mpf(C) / N4
    Larc = mp.mpf(C) * N4
    print(f'== skew double core q={q} M={M} r={r} C={C} mode={mode} L_res={Lres}: primes {ps[0]}..{ps[-1]}, '
          f'W={float(W):.2f} (W/(M log M)={float(W)/(M*math.log(M)):.2f}), log L_arc={float(mp.log(Larc)):.2f}')
    hd, hr = pigeonhole_heights(M, int(mp.ceil(Larc)), C)
    ht = heuristic_typical(M, Larc, C)
    print(f'   log heights: author disc pigeonhole 2(M-1)H_disc: {math.log(hd):.2f};  referee rectangle '
          f'2(M-1)H_rect: {math.log(hr):.2f};  heuristic typical: {math.log(ht):.2f}')
    # characters with ||c||_1 <= NMAX for the gap pre-check
    chars = []
    def rec(i, cur, l1, s):
        if i == M - 1:
            v = -s
            if l1 + abs(v) <= NMAX:
                c = cur + [v]
                if any(c) and c[next(k for k in range(M) if c[k])] > 0:
                    chars.append(c)
            return
        for t in range(-(NMAX - l1), NMAX - l1 + 1):
            rec(i + 1, cur + [t], l1 + abs(t), s + t)
    rec(0, [], 0, 0)
    logp = [mp.log(p) for p in ps]
    def wv(c):
        v = [sum(c[x] * a[x][j] for x in range(M)) for j in range(r)]
        return sum(abs(v[j]) * logp[j] for j in range(r)), any(v)
    chars_w = [(c, wv(c)[0], sum(abs(t) for t in c)) for c in chars]
    results = []
    kept = 0
    tries = 0
    while kept < nsamp:
        tries += 1
        u = [(mp.mpf(rng.getrandbits(900)) / mp.mpf(2) ** 900) for _ in range(M)]
        delta = [Delta * t for t in u]
        ok = True
        for c, w, n in chars_w:
            val = abs(sum(c[x] * delta[x] for x in range(M))) / 2
            req = mp.asin(min(mp.mpf(1), mp.mpf(Lres) ** (mp.mpf(n) / 2) * mp.exp(-w / 2)))
            # dist to (pi/4)Z equals val since val <= n Delta/4 < pi/8
            if val < req:
                ok = False
                break
        if not ok:
            continue
        kept += 1
        ws = [R * mp.expj(d) for d in delta]
        lams = [mp.exp(t) for t in np.arange(10, 40, 1.0)]
        cand = lll_search(ws, lams)
        if not cand:
            print('   sample', kept, 'no candidate found')
            continue
        best = cand[0]
        # Prop R1(i) on every candidate
        bad = 0
        minw = None
        for l1, c, val in cand:
            w, nz = wv(c)
            bound = W + mp.log(2) - 2 * mp.log(1 + l1 * mp.mpf(C) ** 3 / N4 / 48)
            if nz and w <= bound:
                bad += 1
            mw = w - W
            minw = mw if minw is None or mw < minw else minw
        exact_rel = any(val < mp.mpf(10) ** -100 for _, _, val in cand)
        results.append(float(mp.log(best[0])))
        print(f'   sample {kept} (accepted {kept}/{tries}): least l1 found log h = {float(mp.log(best[0])):.2f} '
              f'(||c||_1 = {best[0]}, |sum c w| = {float(best[2]):.3e}); #cands {len(cand)}; '
              f'exact relation among candidates: {exact_rel}; R1(i) violations: {bad}; '
              f'min over cands of w(v(c)) - W = {float(minw):.2f}')
    return results


def run_generic(M, logL, C, nsamp, rng):
    """Generic arc: R chosen so that L = C N^(1/4) = e^logL."""
    N4 = mp.exp(logL) / C
    R = N4 ** 2
    Delta = mp.mpf(C) / N4
    hd, hr = pigeonhole_heights(M, int(mp.ceil(mp.exp(logL))), C)
    ht = heuristic_typical(M, mp.exp(logL), C)
    out = []
    for s in range(nsamp):
        delta = [Delta * (mp.mpf(rng.getrandbits(900)) / mp.mpf(2) ** 900) for _ in range(M)]
        ws = [R * mp.expj(d) for d in delta]
        lo = max(1.0, (logL - 10) / (M - 3) - 3)
        lams = [mp.exp(t) for t in np.arange(lo, lo + 2 * logL / (M - 3) + 12, 0.5)]
        cand = lll_search(ws, lams)
        out.append(float(mp.log(cand[0][0])) if cand else float('nan'))
    return math.log(hd), math.log(hr), math.log(ht), out


if __name__ == '__main__':
    rng = random.Random(7)
    print('### Part A: explicit skew double core fakes (Thm 5.1 of round5/exact.md)')
    for (q, C, mode) in [(7, 1, 'A0'), (7, 1, '2'), (11, 1, 'A0')]:
        run_fake(q, C, mode, 3, rng, NMAX=(6 if q == 7 else 4))
    print()
    print('### Part B: generic random arcs, least AG-violating l1 height vs log L')
    print('M  logL  | author-disc  referee-rect  heuristic | found (min, median, max over samples)')
    for M in (6, 8):
        for logL in (20, 40, 80, 160):
            hd, hr, ht, out = run_generic(M, logL, 1, 4, rng)
            o = sorted(out)
            print(f'{M}  {logL:4d}  | {hd:8.2f}  {hr:8.2f}  {ht:8.2f} | {o[0]:.2f} {o[len(o)//2]:.2f} {o[-1]:.2f}')
