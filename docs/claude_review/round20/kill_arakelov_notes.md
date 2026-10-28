# kill_arakelov working notes
- B(2) normalisation re-derived: sum_{i<j}log|z_i-z_j|^2 <= K log(C^2 R); finite p|N part W(M^2/4-M/2)+Q+2coll;
  W=2log R; difference binom(M,2)logR-(M^2/2-M)logR = (M/2)logR; halve -> (M/4)logR + K log C. Matches (5.1). OK.
- B(1) free levels: unconstrained min automatically has n_0,n_e>=1 (moving a point from lowest level t0 to 0
  changes cut sum by t0(1-M)<0, collisions do not increase). So "free" = Psi^(e). OK.
- A(3) chain Laplacian re-derived (y_tau = x_tau - x_{tau-1} = -S_tau, -Phi^2 = sum y^2). OK.
- multipoint_out.txt pair column sums to 29 (not 28 as stated), expected sums to ~25.66 (not 25.5). check by code.
- Prop B.2 stated for X>=M; depths at q<=X with q^a > X have (q-1)q^{a-1} classes which can be < M when X<2M
  (e.g. q=5, M=24, X=24: depth 2, 20 classes). Fine at X>=2M (witness is at X=2M). Minor.
- Prop 6.1(b): "only Gamma^{>M} not fake data" should be Gamma^{>X}; conservative imprecision.
- kdict.py: Theorem A (1)-(4) all 26 configs pass (independent code). kadelic.py: table 4.2 reproduced to 3 dp.
- kmulti.py: per-row counts reproduced; totals pairs 29 (exp 25.67), triples 1 (exp 0.023). Only triple at q=7 <= 2M=12.
- k2adic.py: norm-N points have classes 1,2,4,4,8,8,... mod (1+i)^a (~2^{a/2}), not 2^{a-1}; genuine 2-adic local min
  sum_a E(M,cl_a) = 85, 2809, 48336 vs K = 45, 1225, 19900 (ratio 1.89, 2.29, 2.43). ADC (ii) includes it;
  construct.md (R) has only one (1+i) per pair; Prop B.2 proof omits (1+i). Repair: lower.md (2-adic at X=2M) or slack. Constant-level.
- (L_8) alone => uniform theorem (overview 4.2) so the "equivalent" needed statement is >= target.
