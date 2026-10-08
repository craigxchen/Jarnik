"""Referee check of Theorem B(i),(iv), (3.2) and the algebraic core of Lemma 6.2, using a DIFFERENT
algorithm from the author's (p-adic cluster discs with the anchor sent to infinity, not the
quartet min-formula), on actual common-unit lattice-point tuples (k = 5..8), including 8 consecutive
points in angular order of one unit class on circles with many points.

Item-532 blocks: layers (p,t), t=1..e_p; T(p,t) = {i in [m] : [a_i>=t] != [a_0>=t]}.
Checks:
 (B1) det(P_i,P_j) = t_ij * prod_{T contains i,j} n_T, t_ij integer != 0;  det(P_0,P_j) = Y_j != 0
 (B1') Norm gcd_{Z[i]}(P_i : i in S) = prod_{T >= S} n_T  (item-532 (7)) for |S| = 1,2,3
 (B2) z_i - z_j = -2i z_0 det(P_i,P_j)/(conj P_i conj P_j)
 (B3) cluster edge length l_T(p) vs v_p(n_T): equal if p divides no t, else |diff| <= 2 max v_p(t)
 (B4) aggregate: |sum_p l_T(p) log p - log n_T| <= 2 H_Sigma for every boundary T
 (B5) cross-ratio factorization used in Lemma 6.2:
      det(a,b)det(c,d)/(det(a,c)det(b,d)) = prod_T n_T^{eps_T} * t_ab t_cd/(t_ac t_bd),
      with #{eps=+1} = #{eps=-1} = 2^(m-3) for every quartet crossing a boundary split.
"""
import sys, math, itertools, random
from fractions import Fraction
from ref_common import *

random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 5)

def t_of(i, j, t):
    return t[(i, j)] if (i, j) in t else -t[(j, i)]   # antisymmetric, like det

stats = dict(tuples=0, consecutive=0, primechecks=0, resid_primes=0, maxdev=0, B5=0, maxH=0.0)
for trial in range(120):
    r = random.randint(2, 4)
    primes = random.sample(SPLIT[:10], r)
    exps = [random.randint(1, 4) for _ in primes]
    pis = {p: two_squares(p) for p in primes}
    allocs = list(itertools.product(*[range(e + 1) for e in exps]))
    if len(allocs) < 5:
        continue
    def zof(al):
        z = (1, 0)
        for p, e, a in zip(primes, exps, al):
            z = gmul(z, gmul(gpow(pis[p], a), gpow(gconj(pis[p]), e - a)))
        return z
    pts = sorted(((zof(al), al) for al in allocs), key=lambda q: math.atan2(q[0][1], q[0][0]))
    n = len(pts)
    choices = []
    k = min(8, n)
    if n >= 5:
        s0 = random.randrange(n)
        choices.append(('consec', [pts[(s0 + j) % n] for j in range(random.randint(5, k))]))
        choices.append(('random', random.sample(pts, random.randint(5, k))))
    for kind, tup in choices:
        k = len(tup); m = k - 1
        random.shuffle(tup)            # arbitrary anchor and labelling
        z = [q[0] for q in tup]; al = [q[1] for q in tup]
        a0 = al[0]
        P = []
        for i in range(k):
            Pi = (1, 0)
            for idx, p in enumerate(primes):
                d = al[i][idx] - a0[idx]
                Pi = gmul(Pi, gpow(pis[p], d) if d > 0 else gpow(gconj(pis[p]), -d))
            P.append(Pi)
        assert P[0] == (1, 0)
        # (B2)
        for i, j in itertools.permutations(range(k), 2):
            num = gmul((0, -2), gmul(z[0], (det2(P[i], P[j]), 0)))
            den = gmul(gconj(P[i]), gconj(P[j]))
            assert gsub(z[i], z[j]) == gdivexact(num, den)
        # blocks
        nT = {}
        for idx, (p, e) in enumerate(zip(primes, exps)):
            for lay in range(1, e + 1):
                T = frozenset(i for i in range(1, k) if (al[i][idx] >= lay) != (a0[idx] >= lay))
                if T:
                    nT[T] = nT.get(T, 1) * p
        def Gset(S):
            return math.prod(v for T, v in nT.items() if S <= T)
        # (B1') gcd identity via rational content trick: Norm gcd = gcd over pairwise? use direct Gaussian gcd
        def ggcd(a, b):
            while b != (0, 0):
                # Euclid in Z[i] with rounding
                nb = gnorm(b); num = gmul(a, gconj(b))
                q = (round(Fraction(num[0], nb)), round(Fraction(num[1], nb)))
                a, b = b, gsub(a, gmul(q, b))
            return a
        for sz in (1, 2, 3):
            for S in itertools.combinations(range(1, k), sz):
                g = P[S[0]]
                for i in S[1:]:
                    g = ggcd(g, P[i])
                assert gnorm(g) == Gset(frozenset(S)), "item-532 (7)"
        # residues
        t = {}
        for i, j in itertools.combinations(range(k), 2):
            d = det2(P[i], P[j])
            assert d != 0
            G = 1 if i == 0 else Gset(frozenset((i, j)))
            assert d % G == 0
            t[(i, j)] = d // G
        H = sum(math.log(abs(v)) for v in t.values())
        stats['maxH'] = max(stats['maxH'], H)
        # all primes involved
        plist = set(primes)
        for v in t.values():
            plist |= set(factor(v))
        # cluster algorithm with anchor at infinity: x_i = X_i/Y_i (Y_i = det(P_0,P_i) != 0)
        bsplits = [frozenset(T) for s in range(2, m) for T in itertools.combinations(range(1, k), s)]
        L = {T: 0.0 for T in bsplits}
        for p in plist:
            def vdist(i, j):  # v_p(x_i - x_j)
                return vp(det2(P[i], P[j]), p) - vp(P[i][1], p) - vp(P[j][1], p)
            mu = max(vp(v, p) for v in t.values())
            for T in bsplits:
                rT = min(vdist(i, j) for i, j in itertools.combinations(sorted(T), 2))
                sT = max(vdist(i, j) for i in T for j in range(1, k) if j not in T)
                ell = rT - sT if rT > sT else 0
                nu = vp(nT.get(T, 1), p)
                if mu == 0:
                    assert ell == nu, ("B3 equality", p, T, ell, nu)
                assert abs(ell - nu) <= 2 * mu, ("B3", p, T, ell, nu, mu)
                stats['maxdev'] = max(stats['maxdev'], abs(ell - nu))
                L[T] += ell * math.log(p)
            stats['primechecks'] += 1
            stats['resid_primes'] += mu > 0
        for T in bsplits:
            assert abs(L[T] - math.log(nT.get(T, 1))) <= 2 * H + 1e-9, ("B4", T)
        # (B5) for every boundary T and every crossing quartet
        if k == 8 or random.random() < 0.3:
            for T in bsplits:
                Tc = [x for x in range(k) if x not in T]
                for a_, b_ in itertools.combinations(sorted(T), 2):
                    for c_, d_ in itertools.combinations(Tc, 2):
                        def D(i, j):
                            return det2(P[i], P[j])
                        chi = Fraction(D(a_, b_) * D(c_, d_), D(a_, c_) * D(b_, d_))
                        prodn = Fraction(1); plus = minus = 0
                        for TT, v in nT.items():
                            if len(TT) < 2 or len(TT) > m - 1:
                                # invisible blocks must cancel
                                pass
                            def G2(i, j):
                                return 0 if (i == 0 or j == 0) else int(i in TT and j in TT)
                            e = G2(a_, b_) + G2(c_, d_) - G2(a_, c_) - G2(b_, d_)
                            if e:
                                assert 2 <= len(TT) <= m - 1
                            prodn *= Fraction(v) ** e
                        tr = Fraction(t_of(a_, b_, t) * t_of(c_, d_, t), t_of(a_, c_, t) * t_of(b_, d_, t))
                        assert chi == prodn * tr, "B5 factorization"
                        # count of +-1 blocks among ALL 2^m-1 subsets (not just nonempty layers)
                        for TT in [frozenset(S) for s in range(1, m + 1) for S in itertools.combinations(range(1, k), s)]:
                            def G2(i, j):
                                return 0 if (i == 0 or j == 0) else int(i in TT and j in TT)
                            e = G2(a_, b_) + G2(c_, d_) - G2(a_, c_) - G2(b_, d_)
                            plus += e == 1; minus += e == -1
                            assert e in (-1, 0, 1)
                        assert plus == minus == 2 ** (m - 3), (plus, minus, m)
                        stats['B5'] += 1
        stats['tuples'] += 1
        stats['consecutive'] += kind == 'consec'
print("Theorem B / (3.2) / Lemma 6.2 core: independent check passed:", stats)
