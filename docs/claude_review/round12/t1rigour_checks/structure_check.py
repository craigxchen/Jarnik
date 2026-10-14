"""Sanity checks (exact integer arithmetic on enumerated boxes) for the cheap-family structure theorem (t1rigour.md §6).
For small random profiles: A = {k in Z^n_0 : k' in {-1,0,1}^(n-1), |<k,beta_p>| <= w}.
 (i)  h_A(v) >= 2 tau_j(v) and tau_j(v) >= j*ell(v) + e(v) for random integer directions (v_n = 0);
 (ii) for random cheap families (total width < m'), the classes of core points number <= prod(ell_i+1) <= 2^b,
      and the exceptional set has size <= (cost-1)/2 - j*sum(ell)  [cost = total width].
"""
import itertools, random
import numpy as np

def tau(v, j):
    s = sorted(v)
    n = len(s)
    return sum(s[n - 1 - i] - s[i] for i in range(j))

def window(v, j):
    s = sorted(v); n = len(s)
    lo, hi = s[j - 1], s[n - j]
    e = sum(1 for x in v if x < lo or x > hi)
    return lo, hi, hi - lo, e

random.seed(12)
tot_dirs = tot_fams = 0
for trial in range(6):
    n = random.choice([12, 13, 14])
    d = n - 1
    r = random.randint(2, 6)
    betas = [np.array([random.choice((-1, 1)) for _ in range(n)]) for _ in range(r)]
    w = random.choice([4.5, 5.0, 6.2, 7.0])
    j = int(w // 2)
    assert 2 * j <= n
    K = np.array(list(itertools.product((-1, 0, 1), repeat=d)), dtype=np.int64)
    kn = -K.sum(axis=1)
    keep = np.ones(len(K), bool)
    for b in betas:
        val = K @ b[:d] + kn * b[d]
        keep &= np.abs(val) <= w
    A = K[keep]
    def width(v):  # v in Z^n, functional <k,v> on Z^n_0; use k' coords: <k',v'-v_n>
        vv = np.array(v[:d]) - v[d]
        y = A @ vv
        return int(y.max() - y.min())
    # (i)
    for _ in range(400):
        v = [random.randint(-3, 3) if random.random() < 0.5 else 0 for _ in range(n)]
        h = width(v)
        t = tau(v, j)
        lo, hi, ell, e = window(v, j)
        assert h >= 2 * t, (v, h, t)
        assert t >= j * ell + e, (v, t, ell, e)
        tot_dirs += 1
    # (ii) candidate directions: columns (0/1), coordinates, pair differences, sums of columns, random sparse
    cands = []
    for b in betas:
        cands.append([1 if x > 0 else 0 for x in b])
    for x in range(n):
        cands.append([1 if y == x else 0 for y in range(n)])
    for _ in range(30):
        v = [0] * n
        for x in random.sample(range(n), random.randint(1, 3)):
            v[x] = random.choice((-1, 1, 2))
        cands.append(v)
    for _ in range(10):
        a, bb = random.sample(range(len(betas)), 2) if len(betas) >= 2 else (0, 0)
        cands.append([cands[a][x] + cands[bb][x] for x in range(n)])
    widths = [width(v) for v in cands]
    mprime = int(np.ceil(4 * w))
    for _ in range(3000):
        fam = [i for i in range(len(cands)) if random.random() < 0.15]
        cost = sum(widths[i] for i in fam)
        if not fam or cost >= mprime:
            continue
        tot_fams += 1
        E = set(); ells = []
        for i in fam:
            lo, hi, ell, e = window(cands[i], j)
            ells.append(ell)
            E |= {x for x in range(n) if cands[i][x] < lo or cands[i][x] > hi}
        core = [x for x in range(n) if x not in E]
        classes = {tuple(cands[i][x] for i in fam) for x in core}
        b_ = (mprime - 1) // (2 * j)
        assert len(classes) <= np.prod([l + 1 for l in ells]) <= 2 ** b_
        assert sum(ells) <= b_
        assert len(E) <= (cost - 2 * j * sum(ells)) // 2 if cost > 0 else len(E) == 0
print(f"structure checks passed: {tot_dirs} directions (h_A >= 2 tau_j, tau_j >= j ell + e), {tot_fams} cheap families")
