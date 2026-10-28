# Referee J4, check 3: J4 step 3 on actual clusters (deterministic inequality, exact per layer) and step 4/5 bookkeeping.
# For a set of M >= 9 points of one circle (any arc), per split layer (p,s):
#   P_S[S meets both {a>=s},{a<s}] = 1 - (C(k,5)+C(M-k,5))/C(M,5) <= q(M)   (exact Fractions)
# hence E_S log N''_S <= q(M) * sum_p e_p log p <= q(M) log N.  We also verify the aggregate with exact enumeration.
# Step 4: reduced distances <= chord/|g_S|; the A_5 test chord/|g_S| <= R''^{1/2-delta} is recorded (evidence only:
# these radii are far below R_1(C)).  Step 5 via Lemma 1.2: R'' >= R/min chord.
import sys
from itertools import combinations
from math import comb, log, sqrt
from fractions import Fraction as Fr
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from gz import *
from c2_identity import decompose
def qM(M): return 1 - min(Fr(comb(k, 5)+comb(M-k, 5), comb(M, 5)) for k in range(1, M))

def cluster_check(cl, N, delta=Fr(1, 31), eta=Fr(1, 1000)):
    M = len(cl); e0, inert, split = decompose(N)
    pis = {p: gauss_prime_above(p) for p in split}
    a = {p: [val(z, pis[p]) for z in cl] for p in split}
    q = qM(M); ok = True
    for p, e in split.items():
        for s in range(1, e+1):
            k = sum(1 for x in a[p] if x >= s)
            prob = 1 - Fr(comb(k, 5)+comb(M-k, 5), comb(M, 5))
            ok &= prob <= q
    # aggregate, exact: sum_S log2-free form -> compare exponents layer-weighted via exact enumeration of N''
    tot = 0.0; good = 0; nS = 0; lemma12 = True; a5 = 0
    R = sqrt(N)
    for S in combinations(range(M), 5):
        g = ggcd_list([cl[j] for j in S]); N2 = N // norm(g)
        prod = 1
        for p in split: prod *= p ** (max(a[p][j] for j in S) - min(a[p][j] for j in S))
        ok &= (prod == N2)
        tot += log(N2); nS += 1
        if log(N2) <= float(Fr(15, 16)+eta)*log(N): good += 1
        chords = [sqrt(norm((cl[i][0]-cl[j][0], cl[i][1]-cl[j][1]))) for i, j in combinations(S, 2)]
        lemma12 &= sqrt(N2) >= R/min(chords) - 1e-9
        if max(chords)/sqrt(norm(g)) <= sqrt(N2)**float(Fr(1, 2)-delta): a5 += 1
    mean = tot/nS/log(N)
    ok &= mean <= float(q) + 1e-12
    return M, ok, mean, float(q), good/nS, lemma12, a5, nS

if __name__ == "__main__":
    n = 1176852625; pts = lattice_points(n); L = len(pts)
    best = min(((sqrt(norm((pts[i][0]-pts[(i+9) % L][0], pts[i][1]-pts[(i+9) % L][1]))), i) for i in range(L)))
    cl = [pts[(best[1]+j) % L] for j in range(10)]
    print("record 10-cluster: chord span/sqrt(R)=%.4f" % (best[0]/n**0.25))
    print("  M=%d ok=%s mean logN''/logN=%.4f <= q(M)=%.4f; frac(logN''<=(15/16+.001)logN)=%.3f; Lemma1.2 ok=%s; "
          "subsets already at A_5(1/31) scale: %d/%d" % cluster_check(cl, n))
    # whole lattice-point sets of circles with many points (M up to ~ 48): the inequality is arc-free
    allok = True; cnt = 0; worst = 0
    for N in [5**2*13**2*17*29, 5**3*13*17*29*37, 2*5**2*13*17*29*41, 9*5*13*17*29*37, 5**4*13**2*17**2,
              5*13*17*29*37*41, 5**2*13**2*17**2*29, 65**3, 1105**2, 5**3*13*17*29]:
        pts = lattice_points(N)
        for start in range(0, len(pts), max(1, len(pts)//4)):
            for M in (9, 12, 16):
                if M > len(pts): continue
                cl = [pts[(start+j) % len(pts)] for j in range(M)]
                r = cluster_check(cl, N)
                allok &= r[1] and r[5]; cnt += 1; worst = max(worst, r[2]/r[3])
    print("consecutive-point clusters tested: %d; all exact checks ok: %s; max mean/q(M) = %.4f" % (cnt, allok, worst))
