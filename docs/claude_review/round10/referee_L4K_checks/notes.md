# Referee notes (L4K)

Definitions (round9/L4.md §1): beta-good up to norm K: every v outside the character lattice with
1<=||v||_1=k<=K has log m_(v,s) <= beta k log M for every s. Author's "good" (per-unit height with
B(v,s)=w(v)/2+[s odd]log2/2) is stronger in regime (1.4) when beta >= h_0/2 + log2/(2 log M). Counterexample
direction OK.

Theorem 1(b) arithmetic re-derived by hand: OK (2 Lambda > 2(K+1)KB + K; max a = 2(K+1)^2B+2K+1; K max a = A_K).
Theorem 1(a) pigeonhole OK; for K odd a better set S_h u (e_1+S_h) gives ~2|S_h| (constant only).

Thm 2 (ii) union bound re-derived: OK; only needs 2(K+1)B < n_q for the "one sigma per u" step.
Thm 2 issue candidates:
 - q = p_j > q_T: (2.1) undefined (nu_j(p_j) undefined); statement says "every prime q > q_T". In regime (1.4)
   with K/2 > h_0 no issue; for K small vs h_0 (e.g. K=2, p_j up to M^h0) p_j > q_T possible. Label systems in
   round 9 are defined at q not dividing N (check). Minor.
 - L4^forced with F_K = {n_q <= A_K}: contains free moduli between ~|S_h| and A_K; Thm-2-type non-coherence
   likely there too -> corrected statement must use exact forced set, not {n_q<=A_K}. Test.
 - Remark 2 "not in P3-all*": substantively true (norm 2h+2 monomials with <v,a>=0 as integers collide at all
   levels at a positive fraction of free q). Test concretely.

Computations to do: thm1, free thresholds, genuine check, forced_sig count, instance coherence, L4^forced gap,
Remark 2, K=4 spot check.

RESULTS so far
- ref_genuine.py: 74,080 checks, 0 mismatches (collision iff q^a | D_s), parity 0 mismatches.
- ref_thm1.py: K=1..6, r=2..6, all parity patterns, random sigma: ALL OK, max = A_K.
- ref_free.py: K=2 r=4..8 first/last and K=3 r=4..8 first/last: q* all match the author's table.
- ref_forced_sig.py: 5,066,549,580,791,808 good level-1 systems, 36,864 tuples, 2^37 genuine: exact match; shares match.
- ref_sigclass.py: 16/32, 64/128, 256/1024, 512/4096 match. NOTE q=3 (n=1) and q=7 (n=2): genuine class = ALL parity
  vectors -> rigidity there is vacuous. K=2 "rigid only at {3,7}" = no non-vacuous rigidity; K=3 joint "except {3,7}"
  = every modulus where non-genuine is possible.
- ref_instance.py: K=2 r=4 first H=2000: 0.0827/0.0818 (author 0.076-0.083); last 0.078; r=8 first 0.034/0.026;
  K=3 r=4 Y=1625: 0.065/0.058; H=8000: 0.098/0.100. Restricted to F_K={n<=A_K}: coherent share <= 0.058+0.19 (K=2 r=4),
  <= 0.07+0.09 (r=8), 0.06+0.065 (K=3) -> L4^forced with F_K={n<=A_K} fails at the test; must use exact forced set.
- ref_remark2.py: (2.1) system: norm-4 v with <v,a>=0 in Z at 63% (r=4) / 100% (r=8) of free q -> collide at ALL levels,
  demand unbounded in X. Remark 2 conclusion right, its wording ("about as often as random") wrong (far worse).
- Threshold: L4(K)'s conclusion is Hamming ("all but small part of mass") -> feeds Theorem B (D <= X), not Thm A.
  So K >= 2(h_0+h)(1+o(1)) suffices for the analytic step; "K >= 2L(h_0+h) (twice README)" is overstated and the
  "(proved)" label on "must" is wrong (Linnik bound is an upper bound on the least prime).
- ref_deform (K=4): 47:8/0, 83:36/24, 97:16/0, 101:20/0, 103:48/0, 107:216/144 -- exact match with author.
- ref_core: 1 tuple, 1,966,080 systems, all genuine -- exact match.
- ref_top (K=4, genuine <=103, zero above 199): non-genuine at all 18 top, 0 viol, no zeta pair coherent on all 5 top Q7;
  max pair-coherent share of top Q7 mass 0.41-0.60 -> system coherent share of forced Q7 in [0.5, ~0.71].
- ref_top2: same with genuine at ALL primes in (199,349] and (199,1000]: flexible top still found (2 seeds each).
  => flexible top robust to genuine labels above the forced range.
- ref_deform_mix (norm<=4, genuine on norm-5 forced set <=349 except q0): 107: 36/0 RIGID; 127: 1152/640 flexible.
- ref_deform K=5 (norm<=5, genuine on F_5 except q0): 107: 18/0, 113: 24/0 rigid. (one-modulus; sensitive to reference)
- ref_deform_ext: sign-extended class doesn't change 83 (24) and 107 (144) counts.
- Running: ref_top5 (joint search with norm<=5 on [107,199]).
- K=5 one-modulus (genuine on F_5 ref): 107 18/0, 113 24/0, 127 32/0, 151 36/0 rigid. mix (norm<=4, F_5 ref): 107 36/0,
  127 1152/640, 151 1728/1584 -> norm-5 monomials cause rigidity at 127, 151.
- K=5 joint searches (top5 3 seeds, top5b 2 seeds): fail at 109 after a 13-collision non-genuine at 107.
  topQ7 (non-genuine only at 5 top Q7, genuine elsewhere): 199 ok (10 coll.), 191 fails (2 seeds).
- 211 run OOM-killed (memory); rerunning with ref_deform2.py (chunk by (x1,x2)) for 167, 199, 211.
- K=5 one-modulus: 167 120/96, 199 120/40 (flexible), 211 884/832 (flexible). mix 167: 90816/90048.
