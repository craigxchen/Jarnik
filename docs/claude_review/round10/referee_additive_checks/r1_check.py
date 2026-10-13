"""Referee Prop R1(i) on explicit skew double core fakes: for every zero-sum c with ||c||_1 <= NMAX,
|sum c_x w_x| >= sqrt2 e^((W - w(v(c)))/2) - n C^3 N^(-1/4)/48 (positions measured from the arc centre).
Also records min over c of |sum c w| and of the bound, and how many c have w(v(c)) <= W (where R1(i)
forces |sum c w| >= 1)."""
import random, math
import mpmath as mp
from ag_fake import skew_profile, primes_1mod4_from
mp.mp.dps = 200
rng = random.Random(11)
for (q, NMAX) in [(7, 6), (11, 4)]:
    M, a = skew_profile(q)
    r = 2 * (M - 1)
    C = 1
    Pst = max(math.ceil(160000 * M ** 4 / C ** 2), 21 * r * r)
    ps = primes_1mod4_from(Pst, r)
    logp = [mp.log(p) for p in ps]
    W = sum(logp); R = mp.exp(W / 2); N4 = mp.exp(W / 4); Delta = mp.mpf(C) / N4
    for s in range(2):
        u = [mp.mpf(rng.getrandbits(900)) / mp.mpf(2) ** 900 for _ in range(M)]
        delta = [Delta * (t - mp.mpf(1) / 2) for t in u]          # centred arc
        ws = [R * mp.expj(d) for d in delta]
        chars = []
        def rec(i, cur, l1, sm):
            if i == M - 1:
                v = -sm
                if l1 + abs(v) <= NMAX:
                    c = cur + [v]
                    if any(c) and c[next(k for k in range(M) if c[k])] > 0:
                        chars.append(c)
                return
            for t in range(-(NMAX - l1), NMAX - l1 + 1):
                rec(i + 1, cur + [t], l1 + abs(t), sm + t)
        rec(0, [], 0, 0)
        viol = 0; low = 0; minval = None; minratio = None
        for c in chars:
            n = sum(abs(t) for t in c)
            v = [sum(c[x] * a[x][j] for x in range(M)) for j in range(r)]
            w = sum(abs(v[j]) * logp[j] for j in range(r))
            val = abs(sum(c[x] * ws[x] for x in range(M)))
            bound = mp.sqrt(2) * mp.exp((W - w) / 2) - n * mp.mpf(C) ** 3 / N4 / 48
            if val < bound:
                viol += 1
            if w <= W:
                low += 1
            ratio = val / mp.sqrt(2) / mp.exp((W - w) / 2)
            minratio = ratio if minratio is None or ratio < minratio else minratio
            minval = val if minval is None or val < minval else minval
        print(f'q={q} M={M} W={float(W):.1f} sample {s}: {len(chars)} chars (||c||_1<={NMAX}), '
              f'{low} with w(v(c)) <= W; R1(i) violations {viol}; min |sum c w| = {mp.nstr(minval, 5)}; '
              f'min |sum c w| / (sqrt2 e^((W-w)/2)) = {mp.nstr(minratio, 6)}')
