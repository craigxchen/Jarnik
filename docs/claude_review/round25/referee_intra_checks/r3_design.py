# Referee check r3: 5-wise uniform (linear) cut designs. Every 5-subtuple is EXACTLY balanced, yet the cluster point Z
# is far from the cube profile required by Theorem 5.1 (most cuts of [M] carry no prime).  Exact (integer counting).
import sys; sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round25/referee_intra_checks')
from itertools import combinations
from math import log
def dot(u, c): return bin(u & c).count('1') & 1
def five_indep_greedy(m, cap=None):
    vecs = []; sums = {0}  # sums of <=4 distinct chosen vectors
    for v in range(1, 2**m):
        if v in sums: continue          # v = sum of <=4 chosen  <=> some <=5-subset incl. v sums to 0
        # update sums of <=4: need sums of <=3 to add v; recompute simply
        vecs.append(v)
        s = {0}
        for k in range(1, 5):
            for T in combinations(vecs, k):
                x = 0
                for t in T: x ^= t
                s.add(x)
        sums = s
        if cap and len(vecs) >= cap: break
    return vecs
def ext_bch(t):
    # columns (1, a, a^3) for a in GF(2^t): parity-check of extended 2-error-correcting BCH, any 5 columns independent
    polys = {3: 0b1011, 4: 0b10011, 5: 0b100101}
    P = polys[t]
    def mul(a, b):
        r = 0
        while b:
            if b & 1: r ^= a
            b >>= 1; a <<= 1
            if a >> t: a ^= P
        return r
    cols = []
    for a in range(2**t):
        a3 = mul(a, mul(a, a))
        cols.append(1 | (a << 1) | (a3 << (t+1)))
    return cols, 2*t+1
def analyse(name, cols, m):
    M = len(cols); us = range(1, 2**m); r = 2**m-1
    # independence check
    ok_ind = True
    for k in range(1, 6):
        for T in combinations(range(M), k):
            x = 0
            for t in T: x ^= cols[t]
            if x == 0: ok_ind = False
    # cut of [M] carried by each prime u (normalised to contain index 0 or its complement)
    cutset = set()
    for u in us:
        pat = tuple(dot(u, c) for c in cols)
        if pat[0] == 1: pat = tuple(1-b for b in pat)
        cutset.add(pat)
    ok_bal = True
    for S in combinations(range(M), 5):
        cnt = {}
        for u in us:
            pat = tuple(dot(u, cols[j]) for j in S)
            if pat[0] == 1: pat = tuple(1-b for b in pat)
            cnt[pat] = cnt.get(pat, 0)+1
        triv = cnt.get((0,)*5, 0); nontriv = [v for p, v in cnt.items() if p != (0,)*5]
        if not (len(nontriv) == 15 and set(nontriv) == {2**(m-4)} and triv == 2**(m-4)-1): ok_bal = False
    ncuts = 2**(M-1)-1
    mu = (2**m-1)/(4*2**(m-4))
    print(f"{name}: M={M}, m={m}, r=2^m-1={r} primes (one per u); any 5 columns independent: {ok_ind}")
    print(f"   every 5-subtuple exactly balanced (15 cuts x 2^(m-4) primes, gcd class 2^(m-4)-1): {ok_bal}")
    print(f"   distinct cuts of [M] carried: {len(cutset)} of {ncuts}  (Thm 5.1 needs all {ncuts}; {ncuts-len(cutset)} empty)")
    print(f"   equal-weight primes, diameter L=C sqrt R, subtuple diam >= L/kappa:  mu(y_S) = 4(1-2^-m) - log(C/kappa)/w_S = {mu:.6f} - log(C/kappa)/w_S  (< 4 iff kappa < C*e^{{w_S 2^(2-m)}})")
cols = five_indep_greedy(6); analyse("greedy F_2^6", cols, 6)
cols = five_indep_greedy(7); analyse("greedy F_2^7", cols, 7)
cols, m = ext_bch(3); analyse("ext-BCH t=3", cols, m)
cols, m = ext_bch(4); analyse("ext-BCH t=4", cols, m)
# route-5 regime arithmetic for ext-BCH, M = 2^t, r = 2^(2t+1)-1 (round10/determinant.md Thm 2.1(d), Thm 4.1)
for t in range(8, 16):
    M = 2**t; r = 2**(2*t+1)-1
    print(f"t={t}: M={M}, r={r}, n-1 >= 300 ln(4r)? {M-1 >= 300*log(4*r)}  ({M-1} vs {300*log(4*r):.0f}); cube would need r={2**(M-1)-1 if M<64 else '2^(M-1)-1'}")
