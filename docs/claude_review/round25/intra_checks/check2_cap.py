"""Check 2 (exact): orders and the balanced cap with shared primes (Lemma 3.2, Thm 4.1, Thm 5.1, Prop 6.1).
 (i)  star partition: Phi^*: M_5^q -> M_M is the identity in the bases e_x - e_anchor (so T_M ~= T_5^q), and every
      route-5 slab monomial k (k' in {-T..T}^{M-1}) has a unique preimage with sum|s_j| <= 8T (s_j in 4T R_5);
 (ii) random F on T_5^n (blocks in [M], shared indices): idx(F) := min|J| <= ord_e(G) <= reg(supp G)
      <= 2^{2-M} sum_cuts w_A(supp G)   (Lemma 3.2 + Theorem B on the sum-zero slice), all exact;
      and with a star partition idx(F) == ord_e(G);
 (iii) cube (fully balanced, simple, w=1) profile on T_M: RHS(*) = 1/2 sum_cuts width_{Phi^*Q_F}(1_A)
      >= 2^{M-3} ord_e(G) >= 2^{M-3} idx(F)   (Theorem 5.1 at eta = 0);
 (iv) lacunary inequality of Prop 6.1(a) on the actual 10-point profile: width_{Q_F}(v) >= width_{pr_1 Q_F}(v_1)
      - sum_{i>=2} d_i width_{R_5}(v_i), for random F with Q_F in prod d_i R_5."""
import random, sys, math
from fractions import Fraction as Fr
from itertools import combinations, combinations_with_replacement, product
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from glib import *

def rank_exact(rows):
    rows = [[Fr(x) for x in r] for r in rows]; rk = 0; ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] != 0), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        for i in range(len(rows)):
            if i != rk and rows[i][c] != 0:
                f = rows[i][c]/rows[rk][c]; rows[i] = [x-f*y for x, y in zip(rows[i], rows[rk])]
        rk += 1
    return rk

def monos(nv, t):
    """exponent vectors of total degree < t in nv variables"""
    out = []
    for deg in range(t):
        for c in combinations_with_replacement(range(nv), deg):
            a = [0]*nv
            for i in c: a[i] += 1
            out.append(tuple(a))
    return out

def powprod(k, a):
    r = 1
    for x, e in zip(k, a): r *= x**e
    return r

def reg(S):
    S = list(S); nv = len(S[0]) - 1; t = 0
    while True:
        E = [[powprod(k[:nv], a) for a in monos(nv, t+1)] for k in S]
        if rank_exact(E) == len(S): return t
        t += 1

def ord_e(G):
    """order at e of G = sum c_k x^k (k sum-zero): least t with sum_k c_k k'^a != 0 for some |a| = t"""
    ks = list(G); nv = len(ks[0]) - 1; t = 0
    while True:
        for c in combinations_with_replacement(range(nv), t):
            a = [0]*nv
            for i in c: a[i] += 1
            if sum(G[k]*powprod(k[:nv], a) for k in ks) != 0: return t
        t += 1

def idx(F, n):
    """min |J| with F_J != 0: least total degree of (beta_1..beta_n), beta_i on 4 coords of block i"""
    sv = list(F); t = 0
    while True:
        for c in combinations_with_replacement(range(4*n), t):
            a = [0]*(4*n)
            for i in c: a[i] += 1
            tot = 0
            for s in sv:
                v = F[s]
                for i in range(n): v *= powprod(s[i][:4], a[4*i:4*i+4])
                tot += v
            if tot != 0: return t
        t += 1

def nullvec(rows, ncols, rng):
    """random integer vector in the right kernel of an integer matrix (exact RREF)"""
    A = [[Fr(x) for x in r] for r in rows]; piv = []; rk = 0
    for c in range(ncols):
        p = next((i for i in range(rk, len(A)) if A[i][c] != 0), None)
        if p is None: continue
        A[rk], A[p] = A[p], A[rk]; A[rk] = [x/A[rk][c] for x in A[rk]]
        for i in range(len(A)):
            if i != rk and A[i][c] != 0:
                f = A[i][c]; A[i] = [x-f*y for x, y in zip(A[i], A[rk])]
        piv.append(c); rk += 1
    free = [c for c in range(ncols) if c not in piv]
    if not free: return None
    v = [Fr(0)]*ncols
    for c in free: v[c] = Fr(rng.randint(-2, 2))
    if all(v[c] == 0 for c in free): v[free[0]] = Fr(1)
    for r, c in enumerate(piv): v[c] = -sum(A[r][f]*v[f] for f in free)
    L = 1
    for x in v: L = L*x.denominator//math.gcd(L, x.denominator)
    return [int(x*L) for x in v]

def rand_s(d, rng):
    s = [0]*5
    for _ in range(d):
        if rng.random() < 0.85:
            a, b = rng.sample(range(5), 2); s[a] += 1; s[b] -= 1
    return tuple(s)

def pull(blocks, svec, M):
    k = [0]*M
    for S, s in zip(blocks, svec):
        for j, sj in zip(S, s): k[j] += sj
    return tuple(k)

def cuts(M):
    for r in range(1, M//2+1):
        for A in combinations(range(M), r):
            if 2*r < M or 0 in A: yield set(A)

def wA(S, A):
    vals = [sum(k[j] for j in A) for k in S]; return max(vals)-min(vals)

def part_i(rng):
    print("(i) star partition")
    for M in (9, 13, 17):
        q = (M-1)//4; anchor = M-1; Bs = [list(range(4*j, 4*j+4)) for j in range(q)]
        blocks = [tuple(B+[anchor]) for B in Bs]
        mat = []
        for j, S in enumerate(blocks):
            for x in range(4):
                s = [0]*5; s[x] = 1; s[4] = -1
                svec = [tuple([0]*5)]*q; svec[j] = tuple(s)
                k = pull(blocks, svec, M); mat.append(k[:M-1])
        ident = all(mat[r][c] == (1 if r == c else 0) for r in range(M-1) for c in range(M-1))
        okT = True
        for T in (1, 2):
            for _ in range(2000):
                kp = [rng.randint(-T, T) for _ in range(M-1)]
                svec = []
                for j, B in enumerate(Bs):
                    s = [kp[x] for x in B]; s.append(-sum(s)); svec.append(tuple(s))
                    if sum(abs(x) for x in s) > 8*T: okT = False
                k = pull(blocks, svec, M)
                if list(k[:M-1]) != kp or sum(k) != 0: okT = False
        print(f"  M={M}, q={q}: Phi^* = identity in bases e_x - e_anchor: {ident}; slab monomials pull back with s_j in 4T R_5 (T=1,2): {okT}")

def part_ii_iii(rng, M, trials):
    okA = okB = okC = okD = okS = True; ex = []
    for t in range(trials):
        star = (t % 4 == 0) and M == 9
        if star:
            blocks = [(0, 1, 2, 3, 8), (4, 5, 6, 7, 8)]
        else:
            n = rng.choice([2, 3]); blocks = [tuple(sorted(rng.sample(range(M), 5))) for _ in range(n)]
        n = len(blocks); degs = [rng.choice([1, 2]) for _ in range(n)]
        mode = t % 2   # 0: impose T_5^n index >= tt ; 1: impose T_M order >= tt (route-5 bookkeeping)
        tt = 2 if n == 3 else rng.choice([2, 3])
        nv = (4*n) if mode == 0 else (M-1)
        conds = monos(nv, tt)
        supp = list({tuple(rand_s(d, rng) for d in degs) for _ in range(len(conds)+6)})
        def feat(s):
            if mode == 0:
                flat = [x for i in range(n) for x in s[i][:4]]
            else:
                flat = list(pull(blocks, s, M)[:M-1])
            return flat
        rows = [[powprod(feat(s), a) for s in supp] for a in conds]
        vec = nullvec(rows, len(supp), rng)
        if vec is None: continue
        F = {s: c for s, c in zip(supp, vec) if c}
        G = {}
        for s, c in F.items():
            k = pull(blocks, s, M); G[k] = G.get(k, 0) + c
        G = {k: c for k, c in G.items() if c}
        if len(G) < 2: continue
        i_F = idx(F, n); o = ord_e(G); r = reg(list(G))
        mw = Fr(sum(wA(list(G), A) for A in cuts(M)), 2**(M-2))
        okA &= (i_F <= o); okB &= (o <= r); okC &= (r <= mw)
        if star: okS &= (i_F == o)
        rhs = Fr(sum(wA([pull(blocks, s, M) for s in F], A) for A in cuts(M)), 2)
        okD &= (rhs >= 2**(M-3)*o)
        ex.append((mode, i_F, o, r, float(mw), float(rhs/2**(M-3)))); strict_io = locals().get('strict_io', 0)
    print(f"  M={M}: idx<=ord {okA}; ord<=reg {okB}; reg<=2^(2-M)sum_cuts w_A {okC}; cube RHS >= 2^(M-3) ord {okD};"
          f" star: idx==ord {okS if M == 9 else 'n/a'}; idx<ord in mode0: {sum(1 for e in ex if e[0]==0 and e[1] < e[2])}/{sum(1 for e in ex if e[0]==0)}, mode1: {sum(1 for e in ex if e[0]==1 and e[1] < e[2])}/{sum(1 for e in ex if e[0]==1)}; samples (mode,idx,ord,reg,MW,RHS/2^(M-3)): {ex[:8]}")

def part_iv(rng):
    d, Z, Nn, g0 = normalised_cluster(1176852625, 10); prof = profile(Z, Nn); ok = True; worst = []
    for _ in range(300):
        n = 3; blocks = [tuple(sorted(rng.sample(range(10), 5))) for _ in range(n)]
        degs = [8, 2, 1]
        F = [tuple(rand_s(dd, rng) for dd in degs) for _ in range(6)]
        for p, e, pi, a in prof:
            vy = [[a[j]-min(a[j] for j in S) for j in S] for S in blocks]
            full = [sum(sum(x*y for x, y in zip(s, v)) for s, v in zip(sv, vy)) for sv in F]
            wF = max(full)-min(full)
            w1 = width([sv[0] for sv in F], vy[0])
            corr = sum(degs[i]*2*(max(vy[i])-min(vy[i])) for i in range(1, n))  # width_{R_5}(v) = 2(max-min)
            ok &= (wF >= w1 - corr); worst.append(wF - (w1-corr))
    print(f"(iv) lacunary width inequality on the 10-pt profile, 300 F x primes: {ok}; min slack {min(worst)}")

if __name__ == "__main__":
    rng = random.Random(2525)
    part_i(rng)
    print("(ii)/(iii) orders, Theorem B on the sum-zero slice, cube cap")
    for M, tr in ((6, 40), (7, 30), (9, 24)):
        part_ii_iii(rng, M, tr)
    part_iv(rng)
