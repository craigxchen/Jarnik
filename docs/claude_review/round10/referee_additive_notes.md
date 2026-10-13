# Referee notes (additive.md)

## Running issues
- N1. Thm 2.1, 2.2, 2.3, 2.4 proofs: read; 2.4(1) pigeonhole fine (square side 2/3 diam 0.943<1).
- N2. Cor 2.5 hypothesis "admissible positions delta form positive-measure set in R^M": fake positions are
  delta = affine image of phi in R^J (+theta_0). If J+1 < M the image is null in R^M -> corollary vacuous
  as stated. Fix: work in phi-space; F_c(delta(phi)) is an exponential polynomial with distinct
  frequencies 2a_x - e (a_x distinct since positions in arc < pi/2), so zero set is null in R^J. Need
  "positive measure in phi-space". CHECK J vs M in round-5 fakes.
- N3. MAIN: "AG provably false for almost every fake" is a measure statement. Barrier logic needs ONE fake
  satisfying the inputs. Actual circles are themselves angle models (genuine phi) and satisfy AG. So a.e.
  violation does not escape barrier; escaping it is exactly Lemma (AG), open, and Lemma AG for all angle
  models includes actual circles (=> stronger than growth). So "AG is fake-violated" must be downgraded
  to "generic-fake-violated".
- N4. Hybrid: points of norm n in an order O of Q(sqrt -d), embedded in C, satisfy AG at ALL heights
  (min of lattice >= 1). If they can carry a Gaussian profile + angles phi satisfying gaps -> fake
  satisfying AG. Try to construct / test.
- N5. HEIGHT EXPONENTS. Thm 2.4(1): H_1 is the BOX height; AG height is 2(M-1)H_1 = M^{1+c0/2+o(1)}, not
  M^{c0/2+o(1)} as in summary/Cor 2.5/table. Same slip in Prop 2.8 (rho <= 2 + c0/(2(A-1)) M).
- N6. CURVATURE LOSS. Values sum c'(w_x - w_1) lie in a rectangle (radial C^2/2 per unit, tangential L per
  unit), not a disc of radius L. Rectangle pigeonhole: H^{M-3} >~ M^2 L C^2 -> box height M^{c0/4},
  l1 height M^{1+c0/4+o(1)}. Prop 2.8 -> rho <= 2 + (c0/(4(A-1)))M. Heuristic 2.9 count omits (M-1)!
  (l1 ball volume) -> typical threshold M^{1+c0/4} too. Imposed-relation dim count -> 1 + c0/2 (not c0).
  Same curvature loss in round9 Prop 5.2 covol bound: true covol(K) = sqrt(sum_triples A^2)/I
  <= M^{3/2} C^3 N^{1/4}/(2I) => relations at l2 height M^{c0/4+1/2}, not c0/2+1/2 (§4 "geometric mean").
  So stated "AG violated at M^{c0/2+o(1)}" holds (via sharpened proof) only for c0 >= 4; for 3<=c0<4 unproved
  and heuristically false.
- N7. NEW (prop R1): first-order expansion + monomial gap: for any angle model, c in Lambda0, v(c)!=0,
  n=||c||_1, nD<=pi/2: |sum c w| >= sqrt2 e^{(W-w(v(c)))/2} - nC^3N^{-1/4}/48. So the monomial gaps (A_0)
  already imply AG for every c with w(v(c)) <= W + log 2 - o(1); AG-violating c have w > W + log2 - o(1).
  Exact relations (v(c)!=0) have w(v(c)) >= 3W/2 - 6 log C - 2 log n + 8.436. TEST on actual + fakes.
- RESULTS actual_checks (4602 runs, 1790 M>=4, 0 failures): P31 slack -0.0100, T32b slack -0.0248, RW min w/W=1.414,
  RW min slack 3.70 (R1(ii) holds), COV exact identity covol(K)^2 I^2 = sum A^2 holds on all 1790.
- RESULTS ag_fake (after fixing 53-bit artifact!): skew core M=8 A0 W=284: least found AG-viol l1 e^{14.0-14.4};
  M=8 mode2 W=410: e^{20.2-20.6}; M=12 A0 W=482: e^{13.9-14.1}. Predictions (disc/rect/heur): 32.3/17.7/13.6,
  44.9/23.9/19.8, 30.7/17.0/13.6. Generic arcs: found tracks heuristic; slope 1/(M-3) confirmed (M=6: 0.3375; M=8: 0.198).
  First run used rng.random() (53-bit) positions -> rational structure -> spurious low heights (saved as
  ag_fake_out_53bit_artifact.txt). Lesson: positions with small denominators violate AG EARLIER.
- EMPTY INTERIOR: bad set for char c is open and contains hyperplane {c.delta=0}; hyperplanes dense -> admissible
  delta-set of skew core (and any rank-M profile) is closed? (complement of open dense) with empty interior,
  measure >= 0.586 Delta^M. Answers author's open point 7.2(a).
