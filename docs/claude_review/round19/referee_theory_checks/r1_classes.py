"""Referee check r1 (independent; imports nothing from the author).
Boundary divisors of Mbar_{0,n} as frozensets S of {0..n-1} with 0 in S, 2<=|S|<=n-2.
Checks:
 (a) Keel relations for beta (all D_T.beta = 1) on Mbar_{0,8};
 (b) K.beta = 10, (K+bd).beta_k for k=4..8, pulled-back degrees, Lambda_k.beta;
 (c) Lambda_k := pi_k^*(K_k+bd_k) - pi_{k-1}^*(K_{k-1}+bd_{k-1}) equals pi_k^* psi_k AS A CLASS:
     compare intersection numbers with all F-curves of Mbar_{0,8} (F-curves span N_1);
 (d) psi_j.beta = 31 for every j, from psi_j = sum_{j in T; a,b notin T} D_T for several (a,b).
Exact arithmetic (Fraction).
"""
from fractions import Fraction as Fr
from itertools import combinations, product

def divisors(n):
    out = []
    rest = list(range(1, n))
    for r in range(1, n - 2):  # |S| = r+1 from 2 to n-2
        for c in combinations(rest, r):
            out.append(frozenset((0,) + c))
    return out

def norm(S, n):
    S = frozenset(S)
    return S if 0 in S else frozenset(range(n)) - S

def K_plus_bd(n):
    return {S: Fr(len(S) * (n - len(S)), n - 1) - 1 for S in divisors(n)}

def K_only(n):
    return {S: Fr(len(S) * (n - len(S)), n - 1) - 2 for S in divisors(n)}

def pullback(div, k, n=8):
    """pull back a divisor on Mbar_{0,k} (points 0..k-1) to Mbar_{0,n} forgetting k..n-1"""
    out = {}
    extra = list(range(k, n))
    for S, c in div.items():
        for bits in product((0, 1), repeat=len(extra)):
            T = set(S) | {extra[i] for i in range(len(extra)) if bits[i]}
            T = norm(T, n)
            out[T] = out.get(T, 0) + c
    return out

def psi(j, a, b, n):
    d = {}
    for S in divisors(n):
        Sc = frozenset(range(n)) - S
        for side in (S, Sc):
            if j in side and a not in side and b not in side:
                d[S] = d.get(S, 0) + 1
    return d

def beta(div):
    return sum(div.values())

def sub(d1, d2):
    out = dict(d1)
    for S, c in d2.items():
        out[S] = out.get(S, 0) - c
    return {S: c for S, c in out.items() if c != 0}

def fcurves(n):
    """all partitions of range(n) into 4 nonempty unlabelled parts"""
    res = []
    def rec(i, parts):
        if i == n:
            if len(parts) == 4:
                res.append([frozenset(p) for p in parts])
            return
        for p in parts:
            p.append(i); rec(i + 1, parts); p.pop()
        if len(parts) < 4:
            parts.append([i]); rec(i + 1, parts); parts.pop()
    rec(0, [])
    return res

def dot_F(div, F, n):
    tot = 0
    full = frozenset(range(n))
    for S, c in div.items():
        Sc = full - S
        v = 0
        for i in range(4):
            if S == F[i] or Sc == F[i]:
                v = -1
        for i, j in combinations(range(4), 2):
            if S == F[i] | F[j]:
                v = 1
        tot += c * v
    return tot

if __name__ == "__main__":
    n = 8
    D = divisors(n)
    print("number of boundary divisors of Mbar_{0,8}:", len(D))
    # (a) Keel relations for beta
    bad = 0; cnt = 0
    for quad in combinations(range(n), 4):
        for (i, j, k, l) in [(quad[0], quad[1], quad[2], quad[3]), (quad[0], quad[2], quad[1], quad[3]), (quad[0], quad[3], quad[1], quad[2])]:
            # relation: sum_{ij|kl} = sum_{ik|jl}; check all three pairings pairwise
            pass
    pairings = []
    for quad in combinations(range(n), 4):
        a, b, c, d = quad
        ps = [((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c))]
        vals = []
        for (p, q) in ps:
            v = 0
            for S in D:
                Sc = frozenset(range(n)) - S
                if (set(p) <= S and not (set(q) & S)) or (set(p) <= Sc and not (set(q) & Sc)):
                    v += 1
            vals.append(v)
        cnt += 1
        if len(set(vals)) != 1:
            bad += 1
    print("(a) quadruples:", cnt, "x 3 pairings (3*70 = 210 quadruple-pairings, giving 2*70*... relations); violations:", bad, "; common value", vals)
    # (b)
    Kb = beta(K_only(8))
    print("(b) K.beta on Mbar_{0,8} =", Kb)
    pulled = {}
    for k in range(4, 9):
        kb = K_plus_bd(k)
        pb = pullback(kb, k)
        pulled[k] = pb
        print(f"    k={k}: (K+bd).beta_k = {beta(kb)}  formula (k-4)2^(k-3)+1 = {(k-4)*2**(k-3)+1};  pulled-back beta-degree = {beta(pb)}  formula 32(k-4)+2^(8-k) = {32*(k-4)+2**(8-k)}")
    lam = {}
    for k in range(4, 9):
        lam[k] = pulled[k] if k == 4 else sub(pulled[k], pulled[k - 1])
        print(f"    Lambda_{k}.beta = {beta(lam[k])}   (32 - 2^(8-k) = {32 - 2**(8-k)})")
    # (c) class identity Lambda_k = pi_k^* psi_k (psi of point k-1, the last point of Mbar_{0,k})
    Fs = fcurves(8)
    print("(c) number of F-curves on Mbar_{0,8}:", len(Fs))
    for k in range(4, 9):
        ps = pullback(psi(k - 1, 0, 1, k), k)
        diff = sub(lam[k], ps)
        mism = sum(1 for F in Fs if dot_F(diff, F, 8) != 0)
        # also check psi with another choice of (a,b) gives the same class
        ps2 = pullback(psi(k - 1, 1, 2, k) if k > 4 else psi(k - 1, 1, 2, k), k)
        mism2 = sum(1 for F in Fs if dot_F(sub(ps, ps2), F, 8) != 0)
        print(f"    k={k}: F-curves with Lambda_k.F != (pi_k^* psi_k).F : {mism};  psi(a,b) choice-independence mismatches: {mism2};  psi_k beta-degree pulled back = {beta(ps)}")
    # sanity of F-curve pairing: boundary divisor self-check, Keel relation via F-curves, K.F values
    KF = {}
    for F in Fs:
        t = tuple(sorted(len(p) for p in F))
        KF.setdefault(t, set()).add(dot_F(K_only(8), F, 8))
    print("    K.F by F-curve type (expect 2222:2, 1223:1, 1133:0, 1124:0, 1115:-1):", {t: sorted(v) for t, v in KF.items()})
    # (d)
    vals = set()
    for j in range(8):
        for a, b in combinations([x for x in range(8) if x != j], 2):
            vals.add(beta(psi(j, a, b, 8)))
    print("(d) psi_j.beta over all j,(a,b):", vals)
    # Lambda_k effective/nef sanity: Lambda_k.F >= 0 for all F (pullback of psi, psi nef on M_0,n)
    for k in range(4, 9):
        mn = min(dot_F(lam[k], F, 8) for F in Fs)
        print(f"    min_F Lambda_{k}.F = {mn}")
