# kill_interp notes (killer's scratch)

## Re-derivations (done by hand)
- Lemma 3.0: correct. moments of (k+l)^alpha - k^alpha have degree <= |alpha|-1.
- Thm 3.1 (1)-(4): correct. Rationality of the induction OK (decompositions can be chosen over Q).
- Equivalence (3.2) <-> M_e(A) = sum M_e(A_l) <-> J(A)_e = cap J(A_l)_e: correct (V(A_l,e) subset V(A,e),
  kernel of c -> Phi^(e)_c on V(A,e) is V(A,e+1), M_e = J_e^perp).
- Lemma 2.1: correct (|Z^v-1| >= e^{-X_v/2}, Lambda(v) = 2X_v/W).
- (2.2) zonotope cost: correct via Ostrowski (Newton polytope widths add).
- Thm 2.2: correct (Baire + UFD + Ostrowski + (CG_m)). Note: ord_e P <= sum m_j; uses only zeta_j=1 count.
- Prop 3.2: correct (monomial ideals).
- Prop 3.5: correct (faithfulness via Gröbner basis of grid + Comb. Nullstellensatz).
- Cor 4.1: implication correct; r <= W/ln 5 OK.

## Record facts found
- round13 README §3/§7.1: "non-zonotopal nondegenerate certificate on an overdetermined fake box would be the
  first sign that route 5 is CAPPABLE" -> round 13's "certificate" = certificate of a TC/base point (cap side).
  interp §3.1 "Reading" attributes the survival-side phrase to round 13 §3: misattribution.
- round16 README §0/§7.2: "first non-zonotopal LOWER-BOUND certificate at an irrational, collision-free direction
  of an overdetermined REGIME box" (quantitative, for (R_H)/(R_inf)). interp's examples: 1-D and 2-D non-slab sets.
- round26: target ZE^all_H; "ZE^all => ZE" needs V_H = V.
- round13: Loc(b0) is a different statement (generation of J(A) by projections onto few columns).
- t1rigour Lemma 6.10 / Thm 6.11: projection-cheap R^{S'} in coordinate hyperplanes (u_x = 0 off S').

## Ideas for kill
- {+-e_x}-only SHIFT fails at e = 2j on A2(s,j) (sigma^{2j} in all J(A2_{+-e_x}), not in J(A2)_{2j}).
  With L_pair need more. Compute.
- In all slab models SHIFT reach = zonotope order: evidence AGAINST (SHIFT)_A in regime.

## Computations so far (killer's own code)
- k0: cube d=2,3,4 SHIFT holds e<2d (coord and pair) exact; hexagon H=1 2 3 1, pair SHIFT fails at 3, coord-only fails at 2.
- k1: 1-D example reproduces exactly (gcd (y-1)^m m<=4; m=5 cofactor root modulus 2.677, |prod| = 24/11; max|S_l|=3;
  rational cyclotomic products reach order 2).
- k2: 2-D example reproduces exactly (Sidon; H_S = 1 2 3 4 4; SHIFT e<=2 using only +-t; products reach 2).
- k3: sdc7 H = 1 7 28 77 161 251 188 43 1, SHIFT fails at 6 (188 vs 166); syl8 H = 1 7 28 77 161 259 287 175 42,
  fails at 6 (287 vs 286). Two primes.
- k4: single-cut boxes: SHIFT holds up to rho-1 (d5s4L1: rho 7; d6s5L1: rho 8; d6s3L1: rho 10).
- k5: no codim-1 light subtorus {w0 = w1 w6} etc. lies in Bs(V(M8,m)), m<=8 (dimension count).
- k6/k7: **order-7 lattice zonotopes with all generators in L_pair exist in sdc7 and syl8 (L=2)**, verified exactly
  (128 pts in A, c in V(A,7)\V(A,8), P_c = Y^k0 prod(Y^v-1), |c|_inf <= 2). So zeta = 7 in both (LP bound 7 attained).
  Criterion = 6 < zeta = 7 < truth (syl8: 8 by round 11 Thm 5.2; sdc7: 7 since base points at 8).
  => interp's "in every computable model it reaches exactly the zonotope order (6 for sdc7 and syl8)" is FALSE.
- k9/k10/k12: syl8 SHIFT_6 HOLDS with all {-1,0,1}^7 shifts (287=287), fails at 7 (175 vs 141); 98 shifts and +weight-3 fail at 6.
  sdc7 with all cube shifts still fails at 6 (188 vs 166).
- k8/k11: 27 random slab boxes: zeta > reach in 20, equal 7, criterion > zeta 0 (both L_pair and all cube shifts).
- k13/k14: **M8' (syl8 L=3): SHIFT holds e<=10 (2 primes), fails at top e=11 (7 vs 0): criterion = 11 = rho = truth
  (round 11 Thm 5.2(a)); zeta = 10 (exact order-10 zonotope).** Criterion beats zonotopes on this slab box.
  syl8 L=1: H = 1 7 28 21, reach 3 = rho.
- Record: Theorem 3.1(1) with L_pair = round 11 compute §5.3 "ZE-braid", refuted by round 11 referee §2.8 at Sylvester
  orders (T=n), t1rigour Prop 6.5 (T=1, m' > 2jk). interp never cites ZE-braid; Prop 3.4 = that refutation.
