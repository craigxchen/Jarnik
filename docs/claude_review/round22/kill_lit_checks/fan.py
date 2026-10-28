# Independent check of the toric geometry of Lbar_k^tw (permutohedral fan in Z^k/Z1).
# (1) Is -K ample (strictly convex support function)?  (2) nef?  (3) Galois-invariant Picard rank.
# (4) polytope of -K: vertices; widths in directions 1_A; Thm-B mean-width bound.
# (5) -K-degree of the 1-PS closure C_v = width of P in direction v.
import itertools
from fractions import Fraction

def rays(k):
    out = []
    for r in range(1, k):
        for A in itertools.combinations(range(k), r):
            out.append(frozenset(A))
    return out

def ind(A, k):
    return tuple(1 if i in A else 0 for i in range(k))

def dot(u, v):
    return sum(a * b for a, b in zip(u, v))

def check(k):
    R = rays(k)
    # maximal cones <-> permutations s: chain {s0} < {s0,s1} < ... < {s0..s_{k-2}}
    # the linear function m with <m,1_A> = 1 on the chain rays, sum m = 0 (well defined mod 1 on N=Z^k/Z1:
    # m in M = {sum = 0}).  Solve: m_{s0}=1, m_{s1}=...=m_{s_{k-2}}=0, m_{s_{k-1}} = -1.
    n_cones = 0; n_fail_strict = 0; n_fail_nef = 0; ms = set(); example = None
    for s in itertools.permutations(range(k)):
        n_cones += 1
        m = [0] * k; m[s[0]] = 1; m[s[-1]] = -1
        chain = [frozenset(s[:j]) for j in range(1, k)]
        for c in chain:
            assert dot(m, ind(c, k)) == 1
        ms.add(tuple(m))
        for A in R:
            val = dot(m, ind(A, k))
            if val > 1:
                n_fail_nef += 1
            if A not in chain and val >= 1:
                n_fail_strict += 1
                if example is None:
                    example = (s, sorted(A), val)
    print(f"k={k}: max cones {n_cones}, distinct m_sigma {len(ms)}, nef violations {n_fail_nef}, "
          f"strict-convexity violations {n_fail_strict}; ample={n_fail_strict==0}; example {example}")
    # Picard rank over Q: Galois orbits of rays (A <-> A^c) minus rank M^Gal (=0, Galois acts by -1)
    orbits = {frozenset([A, frozenset(range(k)) - A]) for A in R}
    print(f"   rays {len(R)}, Galois orbits {len(orbits)} -> rank Pic(X)^Gal = {len(orbits)}; "
          f"geometric Pic rank = {len(R) - (k - 1)}")
    # polytope P of -K = conv(m_sigma); width in direction 1_A
    V = list(ms)
    widths = {}
    for A in [frozenset(B) for r in range(0, k + 1) for B in itertools.combinations(range(k), r)]:
        vals = [dot(u, ind(A, k)) for u in V]
        widths[A] = max(vals) - min(vals)
    MW = Fraction(sum(widths.values()), 2 ** (k - 1))
    print(f"   widths of P(-K) over proper A: {sorted(set(w for A, w in widths.items() if 0 < len(A) < k))}; "
          f"Thm B mean-width bound per m = {MW} (claimed 4-2^(3-k) = {4 - Fraction(1, 2 ** (k - 3)) if k>=3 else None})")
    return V

def deg_on_1ps(V, v):
    vals = [dot(u, v) for u in V]
    return max(vals) - min(vals)

if __name__ == '__main__':
    for k in (3, 4, 5, 6):
        V = check(k)
        if k == 5:
            import random
            random.seed(1)
            ok = True
            for _ in range(2000):
                v = [random.randint(-6, 6) for _ in range(5)]
                if deg_on_1ps(V, v) != 2 * (max(v) - min(v)):
                    ok = False; print("deg mismatch", v)
            print("   -K.C_v == 2(max v - min v) on 2000 random v:", ok)
            # vertices of P: are they exactly e_i - e_j?
            ee = {tuple((1 if t == i else 0) - (1 if t == j else 0) for t in range(5)) for i in range(5) for j in range(5) if i != j}
            print("   vertex set of P == {e_i - e_j}:", set(V) == ee)
