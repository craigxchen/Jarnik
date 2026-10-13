#!/usr/bin/env python3
"""Exact checks for L5.md sections 2-3 (label arithmetic in Z/m_Q, all levels).

Checks:
  (A) Prop 2.1: (w, eps e_y) collides at q^a with unit i^s  <=>  kappa_y = <phi_y - eps w, f> + eps s iota.
  (B) Prop 2.3: if kappa_x = kappa_y mod m_(q^a), the twin monomials of y and the translated twin
      monomials of x have identical collision data at that level.
  (C) Lemma 3.1: greedy per-modulus avoidance with injectivity (constructive).
  (D) Theorem 3.2: switching construction; injectivity and the domination inequality
      a_(q,s)(w, eps e_y) <= sum_i a^lambda_(q,s)(w - 2 eps w_i, eps e_y), for all small w.
Run: python3 check_identities.py   (output: check_identities_out.txt)
"""
import random, itertools, sys
random.seed(20261013)

def chi4(q):  # (-1/q)
    return 1 if q % 4 == 1 else -1

def primes(lo, hi):
    out = []
    for p in range(max(lo, 3), hi + 1):
        if all(p % d for d in range(2, int(p ** 0.5) + 1)):
            out.append(p)
    return out

def twin_profile(M, nJ):
    """Random 0/1 twin profile: columns J (random nonconstant), f(x), s(x) with a_f - a_s = eps_x e_x."""
    r = nJ + 2 * M
    cols = []
    for _ in range(nJ):
        while True:
            c = [random.randint(0, 1) for _ in range(M)]
            if 0 < sum(c) < M:
                break
        cols.append(c)
    eps = []
    fcol, scol = [], []
    for x in range(M):
        e = random.choice([1, -1])
        eps.append(e)
        base = [random.randint(0, 1) for _ in range(M)]
        if e == 1:
            base[x] = 0
            sc = base[:]
            fc = base[:]; fc[x] = 1
        else:
            base[x] = 1
            sc = base[:]
            fc = base[:]; fc[x] = 0
        fcol.append(len(cols)); cols.append(fc)
        scol.append(len(cols)); cols.append(sc)
    A = [[cols[j][x] for j in range(r)] for x in range(M)]
    J = list(range(nJ))
    return A, J, fcol, scol, eps

def setup(M=6, nJ=8):
    A, J, fcol, scol, eps = twin_profile(M, nJ)
    r = len(A[0])
    # phi_y = (a_yJ, (a_(y, s(x)))_x)
    phi = [[A[y][j] for j in J] + [A[y][scol[x]] for x in range(M)] for y in range(M)]
    return A, J, fcol, scol, eps, phi, r

def labels_random(m, nu):
    return [nu_j + 2 * random.randrange(m // 2) for nu_j in nu]

def dot(u, v, m):
    return sum(a * b for a, b in zip(u, v)) % m

def run():
    out = []
    ok = True
    M, nJ = 6, 8
    A, J, fcol, scol, eps, phi, r = setup(M, nJ)
    n = r - M
    qs = [q for q in primes(3, 60)]
    nA = 0; nB = 0
    for q in qs:
        Aq = 2 if q < 8 else 1
        mq = q - chi4(q)
        mQ = mq * q ** (Aq - 1)
        nu = [random.randint(0, 1) for _ in range(r)]
        ell = labels_random(mQ, nu)
        iota = mQ // 4  # a log of i (either generator choice gives m/4 or 3m/4; both span T)
        f = [ell[j] for j in J] + [(ell[fcol[x]] + ell[scol[x]]) % mQ for x in range(M)]
        Lam = [(eps[x] * ell[fcol[x]]) % mQ for x in range(M)]
        F = [dot(phi[y], f, mQ) for y in range(M)]
        kappa = [dot(A[y], ell, mQ) for y in range(M)]
        for y in range(M):
            if (F[y] + Lam[y] - kappa[y]) % mQ:
                ok = False; out.append('kappa != F + Lambda')
        for a in range(1, Aq + 1):
            ma = mq * q ** (a - 1)
            for _ in range(300):
                y = random.randrange(M); e = random.choice([1, -1])
                w = [0] * n
                for _ in range(random.randint(0, 3)):
                    w[random.randrange(n)] += random.choice([1, -1])
                # original-coordinate monomial v with twin coordinates (w, e e_y)
                v = [0] * r
                for idx, j in enumerate(J):
                    v[j] = w[idx]
                for x in range(M):
                    v[scol[x]] = w[nJ + x]
                    v[fcol[x]] = w[nJ + x] + (eps[x] * e if x == y else 0)
                val_v = dot(v, ell, ma)
                val_tw = (dot(w, f, ma) + e * Lam[y]) % ma
                if val_v != val_tw:
                    ok = False; out.append('(2.2) fails')
                for s in range(4):
                    lhs = (val_v - s * (ma // 4)) % ma == 0
                    target = (dot([p - e * ww for p, ww in zip(phi[y], w)], f, ma) + e * s * (ma // 4)) % ma
                    rhs = (kappa[y] - target) % ma == 0
                    nA += 1
                    if lhs != rhs:
                        ok = False; out.append('Prop 2.1 fails q=%d a=%d' % (q, a))
            # (B) class-mate transfer: force a collision of rows 0,1 by adjusting Lambda_1 at this level
            kap = kappa[:]
            Lam2 = Lam[:]
            Lam2[1] = (kappa[0] - F[1]) % mQ  # now kappa_1 = kappa_0 at every level
            if (Lam2[1] - nu[fcol[1]] * eps[1]) % 2:  # parity of Lambda_1 must equal nu_f(1) (eps odd)
                Lam2[1] = (Lam2[1]) % mQ  # parity mismatch would mean different sigma classes; skip
                continue
            for _ in range(200):
                e = random.choice([1, -1])
                w = [0] * n
                for _ in range(random.randint(0, 3)):
                    w[random.randrange(n)] += random.choice([1, -1])
                lhs = (dot(w, f, ma) + e * Lam2[1]) % ma
                w2 = [ww + e * (p0 - p1) for ww, p0, p1 in zip(w, phi[0], phi[1])]
                rhs = (dot(w2, f, ma) + e * Lam2[0]) % ma
                nB += 1
                if lhs != rhs:
                    ok = False; out.append('Prop 2.3 fails')
    out.append('(A) Prop 2.1 instances: %d; (B) Prop 2.3 instances: %d; ok=%s' % (nA, nB, ok))

    # (C) Lemma 3.1, constructive, at moduli with n_q > M - 1 + 4V
    V = 5; okC = True; cnt = 0
    for q in primes(80, 200):
        mq = q - chi4(q); nq = mq // 2
        if nq <= M - 1 + 4 * V:
            continue
        nuf = [random.randint(0, 1) for _ in range(n)]
        f = labels_random(mq, nuf)
        F = [dot(phi[y], f, mq) for y in range(M)]
        nu_fy = [random.randint(0, 1) for _ in range(M)]
        fam = {y: [] for y in range(M)}
        for y in range(M):
            for _ in range(V):
                w = [0] * n
                for _ in range(random.randint(1, 4)):
                    w[random.randrange(n)] += random.choice([1, -1])
                fam[y].append((w, random.choice([1, -1])))
        used = set(); kappa = []
        T = [k * (mq // 4) for k in range(4)]
        for y in range(M):
            par = (F[y] + nu_fy[y]) % 2
            bad = set(used)
            for (w, e) in fam[y]:
                base = dot([p - e * ww for p, ww in zip(phi[y], w)], f, mq)
                for t in T:
                    bad.add((base + t) % mq)
            choice = next(k for k in range(par, mq, 2) if k not in bad)
            kappa.append(choice); used.add(choice)
        # verify
        if len(set(kappa)) != M:
            okC = False
        for y in range(M):
            Lam_y = (kappa[y] - F[y]) % mq
            if Lam_y % 2 != nu_fy[y]:
                okC = False
            for (w, e) in fam[y]:
                if (dot(w, f, mq) + e * Lam_y) % (mq // 4) == 0:
                    okC = False
        cnt += 1
    out.append('(C) Lemma 3.1 constructive: %d moduli, ok=%s' % (cnt, okC))

    # (D) Theorem 3.2 switching with t candidates, Hall via augmenting-path matching
    t, c0 = 5, 2
    S = {y: [[0] * n] for y in range(M)}
    for y in range(M):
        while len(S[y]) < t:
            w = [0] * n
            for _ in range(random.randint(1, c0)):
                w[random.randrange(n)] += random.choice([1, -1])
            if sum(map(abs, w)) <= c0 and w not in S[y]:
                S[y].append(w)
    okD = True; cntD = 0; dom = 0; hallfail = 0
    for q in primes(2 * 2 * M + 1, 400):
        mq = q - chi4(q); nq = mq // 2
        nuf = [random.randint(0, 1) for _ in range(n)]
        f = labels_random(mq, nuf)
        nu_fy = [random.randint(0, 1) for _ in range(M)]
        lam = [nu_fy[y] + 2 * random.randrange(nq) for y in range(M)]
        F = [dot(phi[y], f, mq) for y in range(M)]
        cand = {y: [(F[y] + lam[y] - 2 * dot(w, f, mq)) % mq for w in S[y]] for y in range(M)}
        # bipartite matching rows -> values
        match = {}
        def aug(y, seen):
            for i, val in enumerate(cand[y]):
                if val in seen:
                    continue
                seen.add(val)
                if val not in match or aug(match[val][0], seen):
                    match[val] = (y, i)
                    return True
            return False
        good = all(aug(y, set()) for y in range(M))
        if not good:
            hallfail += 1
            continue
        idx = {y: i for val, (y, i) in match.items()}
        Lam = [(lam[y] - 2 * dot(S[y][idx[y]], f, mq)) % mq for y in range(M)]
        kappa = [(F[y] + Lam[y]) % mq for y in range(M)]
        if len(set(kappa)) != M:
            okD = False
        for y in range(M):
            if Lam[y] % 2 != nu_fy[y]:
                okD = False
        # domination: collision of (w, e e_y) under Lam  =>  collision of some (w - 2 e w_i, e e_y) under lam
        for _ in range(400):
            y = random.randrange(M); e = random.choice([1, -1])
            w = [0] * n
            for _ in range(random.randint(0, 3)):
                w[random.randrange(n)] += random.choice([1, -1])
            for s in range(4):
                hit = (dot(w, f, mq) + e * Lam[y] - s * (mq // 4)) % mq == 0
                if hit:
                    i = idx[y]
                    w2 = [ww - 2 * e * wi for ww, wi in zip(w, S[y][i])]
                    hit2 = (dot(w2, f, mq) + e * lam[y] - s * (mq // 4)) % mq == 0
                    dom += 1
                    if not hit2:
                        okD = False
        cntD += 1
    out.append('(D) Theorem 3.2 switching: %d moduli matched, %d Hall failures (random model), '
               '%d domination instances, ok=%s' % (cntD, hallfail, dom, okD))
    allok = ok and okC and okD
    out.append('ALL OK' if allok else 'FAILURE')
    return out

if __name__ == '__main__':
    res = run()
    txt = '\n'.join(res)
    print(txt)
    open(__file__.replace('.py', '_out.txt'), 'w').write(txt + '\n')
