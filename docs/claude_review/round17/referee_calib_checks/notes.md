# Referee notes (calib.md)

## R1 Prop 2.1 re-derivation (done mentally)
- (1) uses only translation invariance + 2(i'-i) not div by q (m<q). Correct. Real projection = complex projection of real fn.
- (2) rationality uses multiplication invariance (pi_t), Pi rational. Correct.
- (4) |u_theta|_inf >= cos(pi/(2q)): union of grid and grid+pi equispaced pi/q (q odd). Correct.
- (2.1) needs Cheb = d^m kappa_V (reduction Thm 2.1(1)), Psi^{1/2} <= Cheb. Correct.
- CAVEAT: "rigorous k_lo whenever Psi exact": the multimodular E values are residuals relative to span(S), S = pivots mod p.
  If rank_Q > rank_p, the true Psi is SMALLER (S misses columns). Rational reconstruction acceptance is heuristic.
  So k_lo is rigorous only modulo the multimodular caveat; the caveat goes in the unfavourable direction for a lower bound.

## R2 kappa_H claims
- kappaH_7_2_6: H=1 lower bound 1.370e-4 vs k_lo 5.006e-4 -> ratio 0.274 (<0.28 claimed). m=5: 6.254e-4/2.248e-3=0.278.
- only two theta values (0, 0.2244); float.
- "kappa_H is kappa_V up to factor <=4 even at H=1" (7.2) not supported: LP kappa_V(theta-grid inf)=1.073e-3 vs kappa_1 lb 1.37e-4 -> ~7.8.
- kappaH_11_2_8 killed; no kappaH_7_2_8 output file for the claimed H_red=42.

## R3 independent exact (exact_full.py, symmetry-free, fmpq; no multimodular) -- ALL MATCH calib rationals
- sdc7 L=1 m=1..4; sdc7 L=2 m=1..9 (m=7: offdiag != 0 -> Psi NOT constant on C_1 at m=q, confirms m<q needed; m=8 Psi=420175/18168=23.1272 matches §6)
- sdc11 L=1 m=2..4; sdc11 L=2 m=2..5 (m=5: 17047.761788 exact); sdc19 L=1 m=2,3; sdc23 L=1 m=2,3.
- full cube q=7 m=2..6 equals cube_q7.txt (Thm 5.1(a) formula).
## R4 cube at regime (cube_regime.py dps 1500/2500 + cube_logprod.py): q=3251 m'=2057:
  ln Psi_cube=27188.36, ln E_cube=28949.37 (Gaussian model 29436.5), decay 0.856/degree (not ln2). q=103 m=65: 0.814.
## R5 LP (lp_check.py): sdc7 L2 m3 6.595e-2 (0.997), m6 1.0726e-3 (0.9946), L1 m3 4.993e-3 (0.856): match.
## R6 enumeration: own DFS enum2.c: 757, 9373, 264291, 818755 (q=23 in 49 s). q=31 count running.
## R7 own float iso engine iso_float.py: sdc7 m4-6 and sdc11 m5,6 match.
