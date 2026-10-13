# referee notes (scratch)
Read: subcluster.md fully; [D]; [RD].
Checked OK so far: Lemma 1.0 (pi/4 gap, odd N), (1.1), Lemma 1.1 chain, Lemma 1.2, Prop 2.1, Prop 2.2,
Thm 3.1, Prop 3.2 (i)-(iii), Thm 4.1(a),(b) (Siegel H=2MB+1 ok; f(3)>W needs W>8lnM+8lnW+4ln2+12lnC', stated
stronger), Cor 4.2, Prop 4.3, Thm 4.4 proof (induction+peeling; base cases), Cor 4.5 chain, Lemma 5.3(ii) logic.
Issues so far:
 A. ZE-tor def "in the ball, Z_F coincides with finite union of rational subtori through e": global subtori can
    re-enter the ball (other sheets); Lemma 5.3(i) proof "near e closed subgroup = identity component" is only
    true in a K-dependent neighbourhood, not at fixed radius 1/(4D). Fix: define ZE-tor in lifted form
    (Z_F cap ball = union (V_i cap ball), V_i rational subspaces containing 1). Then (i),(ii) trivial/valid.
 B. ZE-col(univ) quantified over ALL profiles and all nondegenerate delta' in L: if the univ family restricted
    to L is small (one relation, rho>=3), zero set near 0 has dim >= rho-2 with generic irrational tangent
    directions -> ZE-col false for such profile -> Thm 5.4 hypothesis false -> vacuous. Fix: only need ZE-col
    for profiles with rho>=300ln4r or M>=(rho+1)(rho+2) (and for actual delta). Test with toy example.
 C. Thm 5.4 "R>=R_0(C)": normalised radius may be small; need M <= const(C) for small normalised R. minor.
 D. §6.4 "below the regime honest combined system empty" generalises from n=4,5 tests; 2137 is threshold
    for the specific T=n slab-box only; other supports (Lemma 1.1 univ family, e.g. Thm 4.1(b)) can be
    admissible at small n. Overstatement.
 E. §6.5 100/100 off-locus: zero evidential value for ZE (null set missed a.s. by random points; [D] Cor 4.3).
    "no zeros" of tangent cones: 12 of 280 systems sampled, BFGS 30 starts: not certified. wording.
 F. Prop 3.2 "Reading": window O(|alpha| log) not shown to be harmless; "capped" relies on fakes' margins,
    not proved. Reading is heuristic; summary item says "So the R1 full-box mechanism never appears" -> overclaim?
 G. Prop 5.1 equivalence (i)<=>(ii) is essentially tautological; (i) loosely formalised.
TODO: exact tests: Thm 4.4 random sets (incl. larger, adversarial); Cor 4.5 cross-polytope push-forward;
  honest_small recompute; Thm 4.1(b) on actual cluster (linear relation sum c_x z_x = 0) exact;
  toy for B (one relation, rho=3, irrational zero directions).
## computations done (referee, own code)
- rc_thm44.py: 950 random sets dims 1-4 + hill climbing (50 restarts): reg <= Omega' never violated; equality on
  l1-balls and boxes reproduced; MW values reproduced (6, 9, 7, 10.5).
- rc_honest.py: own pair-separated codes, random samples incl all bases, n=4,5,6: V(A,floor(LP)+1)=0 always (mod-p full rank).
- rc_onerel.py: single relation mu on Z^3, ord 2: leading forms Q (deg2), C (deg3) have common real zero
  directions with x root of irreducible quartic, y in Q(x) with x^2,x^3 coefficients nonzero -> (1,x,y) Q-lin indep;
  Jacobian det = nonzero cubic in x mod quartic -> transversal. IFT -> curve of zeros with irrational tangent ->
  nondegenerate zeros arbitrarily near e. Numerics (60 digits) follow the curve: |ghat| ~ 1e-60 at t = 1e-2..1e-8.
- rc_tc_exact.py: in both 6.5 models, all 280 single systems: no polynomial of degree <= m=6 vanishes on A
  (full column rank mod p) => single tangent cone = {0} EXACTLY. Trivial because m=6 < 2T+1=7. In the regime
  m >= 0.7n(n-1) >> 2T+1 = 2n+1, so prod_j (k_x - j) (deg 2T+1 <= m) vanishes on A and TC always contains
  sparse/column directions. => 6.5 evidence structurally uninformative.
- rc_threshold.py: r=n: slab-box (T=n) LP-supercritical (m>2w) for all n in [440,2400); generators admissible
  (Lemma 1.1, W = 3n ln n - 14n) for all n in [468,2400). r=n ln n: 563; r=n^2: 881. So "below n~2137 the
  honest combined system is empty" is false for the universal family; 2137 is only [D]'s factor-2 margin.
Additional issues:
 H. "Nondegenerate sub-clusters add nothing" (table, summary) contradicts Remark (i) of 1.2 and is unsupported:
    Thm 3.1(c) is about ONE J; full-rank case push-forward = identity; intersection over J's is the question itself.
 I. "reduces exactly to ZE-col": in full-rank case ZE-col(comb) is the question verbatim -> reformulation.
 J. Prop 3.2 reading: products could beat characters inside the O(log W) window; window is the scale of
    log R/loglog R constants; order-level cap plausible but not proved.
 K. lowrank_profiles.py header says "verified with fractions" -- it isn't; float LP with finite test set (s is an
    upper bound so nonexistence conclusions are valid modulo float, margins large).
 L. Summary "fibres are finite": Prop 2.2 allows whole-circle fibres.
