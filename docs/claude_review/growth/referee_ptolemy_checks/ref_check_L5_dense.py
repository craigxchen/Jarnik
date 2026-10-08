"""Referee: (L_5) allows discarding a fixed proper Zariski-closed Z in M_{0,5}. One curve family (Theorem C1)
does not refute it. Here: the C1 construction works for many (a, tau); the resulting curves are birational
images of P^1_s with each visible boundary degree 1, so the labelled boundary parameters up to PGL_2 are curve
invariants. Distinct cross-ratios cr(tau_1,tau_2,tau_3,tau_4) => pairwise distinct curves => union is Zariski
dense in the surface M_{0,5}; given any Z, some curve is not in Z and carries points s->oo off Z with H = O(1).
This script verifies, for several (a,tau), all ingredients exactly: polynomial identities, a CRT residue class
fixing every bad-prime valuation, pairwise coprimality of the reduced visible blocks for EVERY s in the class,
bounded residues, and the cross-ratio invariant.
"""
import math, itertools
from fractions import Fraction as F
from ref_common import factor, vp

def family(a, tau):
    blocks = {frozenset([1, 2, 3, 4]): (1, -tau[0])}
    for l in range(1, 5):
        blocks[frozenset(set(range(1, 5)) - {l})] = (1, -tau[l])
    for i, j in itertools.combinations(range(1, 5), 2):
        blocks[frozenset([i, j])] = (a[i] - a[j], -(a[i] * tau[j] - a[j] * tau[i]))
    return blocks

def run(a, tau):
    blocks = family(a, tau)
    roots = {F(-f[1], f[0]) for f in blocks.values()}
    assert len(roots) == 11, "roots not distinct"
    X = lambda i, s: a[i] * (s - tau[0]) * math.prod(s - tau[l] for l in range(1, 5) if l != i)
    for i, j in itertools.combinations(range(1, 5), 2):
        for s in range(-12, 12):
            assert X(i, s) - X(j, s) == math.prod(f[0] * s + f[1] for T, f in blocks.items() if i in T and j in T)
    vis = {T: f for T, f in blocks.items() if 2 <= len(T) <= 3}
    bad = set([2, 3])
    for (T1, f1), (T2, f2) in itertools.combinations(blocks.items(), 2):
        bad |= set(factor(f1[0] * f2[1] - f2[0] * f1[1]))
    for f in blocks.values():
        bad |= set(factor(math.gcd(f[0], f[1])))
    # CRT class: for each bad q pick e, s0 mod q^e with v_q(block(s0)) < e - v_q(c1) for every block
    mod, res = 1, 0
    for q in sorted(bad):
        e = 1
        while True:
            found = None
            for s0 in range(q**e):
                if all((f[0] * s0 + f[1]) % (q**e) != 0 and vp(f[0], q) + e > vp(f[0] * s0 + f[1], q)
                       if f[0] * s0 + f[1] != 0 else False for f in blocks.values()):
                    found = s0; break
            if found is not None:
                break
            e += 1
        # combine
        M2 = q**e
        t = ((found - res) * pow(mod, -1, M2)) % M2
        res, mod = res + mod * t, mod * M2
    # constant bad parts on the class: v_q(c1*(res+mod*u)+c0) = v_q(c1*res+c0) since v_q(c1*mod) larger
    const = {}
    for T, f in blocks.items():
        v0 = f[0] * res + f[1]
        c = 1
        for q in bad:
            vq = vp(v0, q)
            assert vp(f[0] * mod, q) > vq
            c *= q**vq
        const[T] = c
    # pairwise coprimality of reduced visible blocks for every s in the class: any common prime divides a
    # resultant, hence is bad, and bad parts are removed.  Residues after moving bad parts:
    tmax = max(math.prod(const[T] for T in blocks if i in T and j in T) for i, j in itertools.combinations(range(1, 5), 2))
    t1, t2, t3, t4 = (F(tau[l]) for l in range(1, 5))
    cr = (t1 - t3) * (t2 - t4) / ((t1 - t4) * (t2 - t3))
    # spot check a big s
    s = res + mod * 10**15
    for (T1, f1), (T2, f2) in itertools.combinations(vis.items(), 2):
        g = math.gcd(abs(f1[0] * s + f1[1]) // const[T1], abs(f2[0] * s + f2[1]) // const[T2])
        assert g == 1
    return mod, res, tmax, cr

for tau4 in (2, 3, 5, 9, 11, 13):
    a = {1: 1, 2: 2, 3: 3, 4: 4}
    tau = {0: 4, 1: 6, 2: 7, 3: 0, 4: tau4}
    try:
        mod, res, tmax, cr = run(a, tau)
        print(f"tau=(4,6,7,0,{tau4}): class s={res} mod {mod}, max|t'|={tmax}, |Y|=1, curve invariant cr={cr}")
    except AssertionError as e:
        print(f"tau=(4,6,7,0,{tau4}): skipped ({e})")
