import random, sys
from fractions import Fraction
from plib import *
import sympy
random.seed(17)
out = []
for q in [7, 11, 19, 23, 31, 43, 47, 59]:
    P1 = patterns(q, +1); P2 = patterns(q, -1)
    # basic structure
    sizes = sorted(set(len(S) for S in P1)); infin = sum(1 for S in P1 if q in S)
    mism = 0; mism2 = 0
    for trial in range(300):
        B = random.choice([1, 1, 2, 3, 5])
        lam = [random.randint(-B, B) if random.random() < random.choice([0.1, 0.3, 1.0]) else 0 for _ in range(q)]
        lam.append(-sum(lam) + random.choice([0, 0, 1, -1, 2]))  # adjust inf so total is small; then fix sum 0
        s = sum(lam); lam[q] -= s
        nd = N_def(lam, P1); nf = N_formula(lam, q); nd2 = N_def(lam, P2)
        mism += (nd != nf); mism2 += (nd != nd2)
    # pairs
    pairs = set()
    for x in range(q+1):
        for y in range(x+1, q+1):
            lam = [0]*(q+1); lam[x] = 1; lam[y] = -1
            pairs.add(N_def(lam, P1) - q)
    out.append(f"q={q}: |S| values {sizes}, #patterns containing inf {infin}/{len(P1)}, distinct {len(set(P1))}; "
               f"formula mismatches {mism}/300, conv(+/-) mismatches {mism2}/300, pair sigmas {sorted(pairs)}")
    print(out[-1], flush=True)
# lattice L_T = {lam in Q^n: sum 0, lam(S) in Z}: compute L_T / Z^n_0 via Smith normal form of the constraint matrix
for q in [7, 11, 19, 23]:
    P = patterns(q); n = q+1
    # basis of sum-zero rationals: v_i = e_i - e_n, i < n. lam = sum c_i v_i. lam(S) = sum c_i (1_S(i) - 1_S(inf)).
    A = sympy.Matrix([[ (1 if i in S else 0) - (1 if q in S else 0) for i in range(n-1)] for S in P])
    # L_T/Z^{n-1} (in v-coords) is dual-lattice quotient; its invariants = nonunit Smith invariants of A (full column rank)
    from sympy.matrices.normalforms import smith_normal_form
    snf = smith_normal_form(A, domain=sympy.ZZ)
    diag = [abs(snf[i, i]) for i in range(min(snf.shape))]
    nonunit = [d for d in diag if d != 1]
    print(f"q={q}: rank {A.rank()}, SNF non-unit invariants {nonunit} (expect [(q+1)/2] = [{(q+1)//2}])", flush=True)
