#!/usr/bin/env python3
"""Independent check of stab.md Cor 3.4 on the Losev-Manin 4-fold (permutohedral fan in Z^5/Z1, 30 rays 1_A),
and a search for pencils of small defect s = D.beta - 8t with h^0(pi^*D - tE) >= 2.
Lattice M = {m in Z^5: sum m = 0}, coordinates m_1..m_4 (m_5 = -sum).  Sections of pi^*D - tE =
Laurent polynomials supported on P_D cap M with order >= t at e = measures annihilating polys of degree < t.
"""
import itertools, random, sys, time
from fractions import Fraction
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/kill_stab_checks')
from thm33 import dimA, MW

RAYS = [frozenset(A) for r in range(1, 5) for A in itertools.combinations(range(5), r)]  # 30 proper nonempty
def pair(m4, A):
    m = list(m4) + [-sum(m4)]
    return sum(m[j] for j in A)

def PD_points(d, R=12):
    """lattice points of P_D = {m: <m,1_A> >= -d_A}; box search |m_j| <= R"""
    pts = []
    # bound: <m, 1_{j}> >= -d_{j} and <m,1_{[5]-j}> = -m_j >= -d_{[5]-j}  => -d_j <= m_j <= d_{comp}
    lo = [-d[frozenset({j})] for j in range(4)]
    hi = [d[frozenset(range(5)) - frozenset({j})] for j in range(4)]
    for m in itertools.product(*[range(l, h + 1) for l, h in zip(lo, hi)]):
        if all(pair(m, A) >= -d[A] for A in RAYS):
            pts.append(m)
    return pts

def cutwidth_sum(S):
    tot = 0
    for A in RAYS:
        if 4 in A: continue  # cut named by side without point 5 (index 4)
        v = [pair(m, A) for m in S]
        tot += max(v) - min(v)
    return tot

def minimal_D(S):
    return {A: -min(pair(m, A) for m in S) for A in RAYS}

def zonotope(edges):
    Z = {(0, 0, 0, 0)}
    for (i, j) in edges:
        e = [0] * 5; e[i] += 1; e[j] -= 1
        seg = [(0, 0, 0, 0), tuple(e[:4])]
        Z = {tuple(a + b for a, b in zip(x, y)) for x in Z for y in seg}
    return sorted(Z)

def main():
    rnd = random.Random(271)
    pairs = [(i, j) for i in range(5) for j in range(i + 1, 5)]
    # ---- (A) root-product classes: D = minimal divisor of a root-binomial zonotope, t = #edges ----
    fails = 0; tested = 0; t0 = time.time()
    for trial in range(120):
        k = rnd.randint(1, 4)
        edges = [rnd.choice(pairs) for _ in range(k)]
        Z = zonotope(edges)
        d = minimal_D(Z)
        Dbeta = sum(d.values())
        S = PD_points(d)
        assert Dbeta == cutwidth_sum(Z) == 8 * k, (Dbeta, cutwidth_sum(Z), k)
        dim = dimA(S, k)
        tested += 1
        if dim > 1: fails += 1; print('Cor3.4 FAIL', edges, len(S), dim)
        if dim != 1: print('  note: dim', dim, 'for root product', edges, '(expected 1)')
    print(f'(A) root-product classes D (D.beta = 8t): {tested} tested, {fails} with h^0(pi*D - tE) >= 2  [{time.time()-t0:.1f}s]', flush=True)

    # ---- (B) perturbed classes: D' = D + random nonneg increments, keep D'.beta = 8t' ----
    fails = 0; tested = 0; nonzero = 0; t0 = time.time()
    for trial in range(200):
        k = rnd.randint(1, 3)
        edges = [rnd.choice(pairs) for _ in range(k)]
        d = minimal_D(zonotope(edges))
        extra = rnd.randint(1, 3) * 8
        for _ in range(extra):
            A = rnd.choice(RAYS); d[A] += 1
        Dbeta = sum(d.values()); t = Dbeta // 8
        S = PD_points(d)
        if len(S) > 400: continue
        dim = dimA(S, t)
        tested += 1; nonzero += dim >= 1
        if dim > 1: fails += 1; print('Cor3.4 FAIL', d, len(S), dim)
    print(f'(B) random classes with D.beta = 8t: {tested} tested, {nonzero} effective, {fails} with h^0 >= 2  [{time.time()-t0:.1f}s]', flush=True)

    # ---- (C) defect search: supports S in M; t2 = max{t: dim A_t(S) >= 2}; s = sum cut widths - 8 t2 ----
    best = {}; t0 = time.time()
    small = [m for m in itertools.product(range(-2, 3), repeat=4) if abs(sum(m)) <= 2]
    roots = [tuple((1 if q == i else 0) - (1 if q == j else 0) for q in range(5))[:4] for i in range(5) for j in range(5) if i != j]
    cands = []
    for trial in range(3000):
        mode = rnd.randint(0, 2)
        if mode == 0:
            S = set(rnd.sample(small, rnd.randint(3, 9)))
        elif mode == 1:  # zonotope + extra point(s)
            Z = set(zonotope([rnd.choice(pairs) for _ in range(rnd.randint(1, 3))]))
            for _ in range(rnd.randint(1, 2)):
                base = rnd.choice(sorted(Z)); r = rnd.choice(roots)
                Z.add(tuple(a + b for a, b in zip(base, r)))
            S = Z
        else:  # union of two zonotopes with equal edge count
            k = rnd.randint(1, 3)
            S = set(zonotope([rnd.choice(pairs) for _ in range(k)])) | set(zonotope([rnd.choice(pairs) for _ in range(k)]))
        S = sorted(S)
        if len(S) < 3: continue
        W = cutwidth_sum(S)
        t2 = None
        for t in range(W // 8 + 1, 0, -1):
            if dimA(S, t, exact_if_ge=3) >= 2:
                t2 = t; break
        if t2 is None: continue
        s = W - 8 * t2
        if s not in best: best[s] = S
    print(f'(C) defect search over 3000 supports [{time.time()-t0:.1f}s]: defects found', sorted(best)[:8])
    for s in sorted(best)[:3]: print('   s =', s, 'example support', best[s])

    # ---- (D) exhaustive small-class search: D with small coefficients, s = D.beta - 8t in {1,2,3} ----
    t0 = time.time(); found = {}
    cnt = 0
    for trial in range(4000):
        d = {A: 0 for A in RAYS}
        for _ in range(rnd.randint(8, 27)):
            d[rnd.choice(RAYS)] += 1
        Dbeta = sum(d.values())
        for s in (1, 2, 3):
            if (Dbeta - s) % 8 == 0 and Dbeta - s > 0:
                t = (Dbeta - s) // 8
                S = PD_points(d)
                if not S or len(S) > 300: continue
                cnt += 1
                if dimA(S, t, exact_if_ge=3) >= 2:
                    found.setdefault(s, d)
    print(f'(D) {cnt} random classes with D.beta - 8t in {{1,2,3}}: pencils found for s =', sorted(found), f'[{time.time()-t0:.1f}s]')

if __name__ == '__main__':
    main()
