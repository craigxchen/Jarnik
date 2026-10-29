# kill_split working notes
N1. Round23 split.md 2: "section cap tau<=15/4 same for split because widths" -> transfer used widths = symmetric weighting.
    Author's gap claim (asymmetric cut-balanced profiles) is legitimate w.r.t. record. Round23 pair-averaging shows binomials
    ACHIEVE 15/4 at every cut-balanced profile (lower bound on alpha), not that nothing beats it.
N2. Typo: Thm 3.2 proof says cost_{1/2} = (1/8) sum_P w_P; correct is (1/4) sum_P w_P (root binomial: cost 4, sum_P w_P = 16).
N3. Split inversion iota (d -> n/d, u_a -> lcm/u_a) maps split A-points to A-points: admissible Z can be taken iota-stable,
    so Z must ALSO contain iota(cone) = cubic (u1u4u5 - u1u2u3 - u2u4u5 - u3u4u5 + u2u3u4 + u2u3u5 = 0) and iota(W). Author omits.
N4. Near-tight condition X_a equal == balancing condition for curve classes: Pi_nt = cut-balanced curve classes on LM_5.
    WTB at x <=> pi^*gamma_x - 4 ell is in Mov(Bl_e LM_5) = Eff^vee (BDPP). Prop 3.7 = easy direction of BDPP (standard).
N5. M0bar_{0,7} -> LM_5 factors through Bl_e LM_5 (fibre over e is D_{0 infty}); contracts the 15 D_T, T in [5], |T| in {3,4}.
    rank check: 42 - 15 = 27 = 26 + 1. So WTB is a statement about Eff(M0bar_7) (hypertree divisors etc.) - literature link.
N6. Verified exactly: Thm 3.1 (B.C=7, 7 distinct cuts, quadratic place -> e, cone identity, Jac rank 3 at 4 pts), Pell gamma->7/2
    (log(dist^3.5 H) periodic {4.064,4.267,4.960}), C->2sqrt3. Co-divisor (iota) tuples: same lcm, gcd 1, on NONE of 15 cones -> Z must contain iota(cone) cubic.
    m-czech W-curves on none of 15 cones (symbolic, 2 param pts). dim Pi_nt=11, 278 one-sided, 16 classes, code hull 10, +x[5] ->11.
N7. Theorem B (construct.md 2.4): reg <= 2^{1-n} sum_P w_P. So cost_{1/2} = (1/4) sum_P w_P; author wrote 1/8 (false as stated, conclusion right).
N8. Covering-family reformulation: one-sided vertex F, P_ab=4 => condition: for each a<5, c_a A_a(s) - D_a(s) = (c_a-1) q(s), q monic quartic
    (roots = 4 branch params). W-curves = q a perfect square (codim 2). Plan: own complex Newton at all 16 one-sided classes.
N9. k_misc.py: Thm 3.6 identity 200/200; exact V(y)=prod t_ab * prod n_B^{-f(|B|)} 200/200; sum x_B f(|B|) = 40 on cut-balanced (no near-tight needed);
    m-czech phase-0 k=(0,224,225,324,900) at 13 Pell points, C->30. mnn31 k=(0,3,4,9,12) by hand (s-values 0,4,9,3,12).
N10. Castravet-Tevelev (Duke 2015, arXiv 1311.7673, Thm 1.1 via search snippet): morphisms LM~_{n+1} -> M0bar_n -> Bl_e LM_n
    (CT's LM_n has dim n-3). For n=7: M0bar_7 -> Bl_e LM_5(record). Record never makes this link (grep). WTB = movable-curve
    statement on M0bar_7 / Bl_e LM_5; Eff(M0bar_7) (hypertree divisors, Opie, DGJ) is the literature testing ground.
N11. k_cover.py (complex random-start damped Gauss-Newton, Mobius gauge s=0,1,-1,s4): 4-branch dominant families FOUND at 14/16
    one-sided classes in <=40 s each (dominance SVs 4 nonzero, min 3e-4..0.38; 16x20 Jacobian full rank). NOT found at class 4
    (1^5 3^10) and 13 (2^10 4^5 = x_[5], the S_5-invariant even-set vertex) = an iota-pair. Same pair where author failed.
    Author's random-start failure is a solver issue (real?) - complex starts work.
N12. Overclaim: "Subspace with rational forms and no special p-adic behaviour = Liouville" false in general (S-unit eq. uses rational
    forms). Lemma 3.5 covers only single functions with Delta-poles; novelty (b) "fixed-S Siegel/CZ input adds nothing" overstated.
N13. THEOREM (mine, to write up): no nondegenerate curve of class x_[5] (N=1) has mult_e >= 4. Proof: phi_a = c_a Q_a(s)/(s-rho_a);
    e-condition => phi_a - phi_b = (c_a-c_b) q(s)(s-r_ab)/((s-rho_a)(s-rho_b)) with common q (deg 4) => phi_a = Psi + kappa_a q/(s-rho_a),
    Psi poly deg<=3. Shared roots r_ab <=> in plane (y,w=ys) the 10 nodes P_ab of 5 lines L_a: w = rho_a y + kappa_a lie on the
    quartic G = Qtilde + Psitilde (order >=3 at O). Quartics through 10 nodes = span{prod_{b!=c} L_b} (restriction-to-line induction).
    Order 3 at O <=> sum mu_c nu(v_c) = 0 (Veronese, v_c = (rho_c,-1)/kappa_c) <=> 4 of v_c collinear <=> 4 lines concurrent/parallel
    => degenerate. Explains Newton failures at classes 13 and (by iota) 4. Covers cuspidal case too (q = prod (s-sigma)^m).
    Consequence: WTB at x_[5] NOT decided; t*(x_[5]) = inf cost/ord could be < 4 -> section beating 15/4 at even profile. TEST.
N14. 276 one-sided vertices except x_[5] pair: hull affine dim 11, average 1/2 -> (numerical) WTB on full 11-dim nbhd of 1/2.
    60-digit refinement of a W-class solution: |G| 1e-14 -> 5e-27 -> 6e-53; dominance SVs (6.28,1.09,0.24,0.029); 16x20 sigma_min 0.109.
N15. Prop K hypotheses automatic within the class: non-nested merges (r_ab=r_ac, r_ab=r_cd, rho_a=rho_b) change the class
    (vector not in a chain cone); nested merge r_ab = rho_e allowed and harmless. So Prop K: NO rational curve of class x_[5]
    with mult_e>=4 outside boundary/Delta. Even-vertex symmetric polytopes P(alpha,beta): max 4reg/cost = 1 (alpha=3,beta=2:
    Vandermonde-tight), no violation.
N16. k_coverN2.py (analytic Jacobian, poly residuals, k_cover-style init, dominance filter) finds W-class N=1 in 1 s.
    Running N=2 at W-class (control) and at x_[5] (900 s each). Old slow k_coverN.py killed.
N17. Pairs-class (B.C=10) with mult 3 also impossible for rational curves (same Veronese argument on lines y = c_a s - kappa_a):
    so the decomposition gamma_x5 - 4l = (pairs - 3l) + (quads - l) is not realised by rational covering families either.
