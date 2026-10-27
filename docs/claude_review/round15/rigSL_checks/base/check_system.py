# check_system.py -- independent exact checker for an explicit joint witness (round 14).
# Input: a jw_*.pkl (labels at every forced modulus: the deepest stored lift ell_top mod m q^(A-1)).
# Recomputes from scratch, for every v with 1 <= ||v||_1 <= K (both signs) and every unit s:
#   m_(v,s) = prod_q q^(b_q(v,s)), b_q = number of levels at which v collides with unit s,
# requiring that no collision persists at the deepest stored level A (so b_q is exact), and checks
#   m_(v,s)^2 <= Norm(Z_v) * (2 if s odd).  Also: parity, coverage of the exact forced set, and the class status
# (genuine class / some column-sign class / none) of the label at every forced modulus.
import sys, pickle, json, math, itertools
from rl import mq, nu, gauss_prime, Torus, cap, units_mod, split_primes

def all_monos(r, K):
    out = []
    for v in itertools.product(range(-K, K + 1), repeat=r):
        if 0 < sum(abs(t) for t in v) <= K: out.append(v)
    return out

def eps_set(L, n, r, c):
    out = set()
    for bits in range(2 ** r):
        eps = tuple((-1 if (bits >> j) & 1 else 1) for j in range(r))
        if any(tuple((u * e * x) % n for e, x in zip(eps, L)) == c for u in units_mod(n)): out.add(eps)
    return out

def classes(L, n, r):
    gen = set(); anyeps = set()
    for bits in range(2 ** r):
        eps = [(-1 if (bits >> j) & 1 else 1) for j in range(r)]
        for u in units_mod(n):
            c = tuple((u * e * x) % n for e, x in zip(eps, L))
            anyeps.add(c)
            if bits == 0: gen.add(c)
    return gen, anyeps

def check(fn, thrfn):
    W = pickle.load(open(fn, "rb")); P, K = W["P"], W["K"]; r = len(P)
    F = json.load(open(thrfn))["forced"]
    assert sorted(W["labels"]) == sorted(F), "witness does not cover the forced set exactly"
    V = all_monos(r, K); pis = [gauss_prime(p) for p in P]
    m_vs = {}
    status = {}; epsq = {}
    for q in F:
        lab = W["labels"][q]; top = lab["lifts"][-1]
        m = mq(q); n = m // 4; T = Torus(q); iota = T.iota
        # depth A of the stored top label: modulus m q^(A-1) with all entries < it (lifts are exact integers)
        A = len(lab["lifts"]) + (W.get("deep", {}).get(q, 1) - 1)
        Mtop = m * q ** (A - 1)
        assert all(0 <= x < Mtop for x in top), (q, A)
        assert all(x % 2 == p for x, p in zip(top, nu(P, q))), "parity"
        for v in V:
            d = sum(a * b for a, b in zip(v, top)) % Mtop
            if d % n: continue
            s = [u for u in range(4) if (u * iota) % m == d % m][0]
            b = 1
            while b < A and d % (n * q ** b) == 0: b += 1
            assert not (b == A and d % (n * q ** (A - 1)) == 0 and A > 0 and False)
            # no collision at the deepest stored level unless it is level A and the next level is avoided:
            if d % (n * q ** (A - 1)) == 0 and b == A:
                # collision at level A itself: deeper levels undetermined -> must not happen for a complete witness
                raise AssertionError(f"collision persists at deepest stored level {A} at q={q}, v={v}")
            m_vs[(v, s)] = m_vs.get((v, s), 1) * q ** b
        L = [T.L(pi) for pi in pis]
        gen, anyeps = classes([x % n for x in L], n, r)
        c = tuple(x % n for x in top)
        status[q] = "gen" if c in gen else ("eps" if c in anyeps else "NONE")
        epsq[q] = eps_set([x % n for x in L], n, r, c)
    for (v, s), mm in m_vs.items():
        assert mm * mm <= cap(P, v, s), ("budget", v, s, mm)
    return status, epsq

if __name__ == '__main__':
    st, epsq = check(sys.argv[1], sys.argv[2])
    W = pickle.load(open(sys.argv[1], "rb"))
    print(f"OK: witness {sys.argv[1]} is parity-respecting and good for every ||v||_1 <= {W['K']} at every level on the "
          f"exact forced set ({len(st)} moduli); q0={W['q0']} status {st[W['q0']]}")
    print("  non-genuine-class moduli:", {q: s for q, s in st.items() if s != "gen"})
    # longest initial segment of forced moduli on which ONE column-sign vector (up to global sign) fits
    F = sorted(st); common = None; seg = None
    for q in F:
        c2 = epsq[q] if common is None else common & epsq[q]
        if not c2: break
        common = c2; seg = q
    print(f"  system-level coherent (one eps for all) on forced moduli <= {seg}; eps = {sorted(common)[:2] if common else None}")
