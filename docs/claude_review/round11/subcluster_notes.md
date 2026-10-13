# notes (scratch)
Plan:
- Universal admissible system U: every integer c on Z^M_0 with ord>=m, |k|_1<=D on supp, and
  sum_p log p h_{supp c}(a_p) < mW/2 - 2ln|c|_1 - 2m ln(DC/2) + 2 ln m!  => F_c(z)=0. Subcluster systems (any J, any base
  point b in J, any T,w,H) are instances (extension by zero). Embedded V(A_J,m_J) NOT inside V(A_J',m_J') (lower order).
- Prop: dim Z_comb <= n-2 (projection argument; tangent space).
- Fibre lemma: fix J0 (n-1 pts incl base), x outside: slice in zeta_x is finite (<=2T pts) or full.
- Base-point freedom: sparse strata of Prop 5.1 need coincidence with base point; all base points => only delta=const? check.
- Column subtori of rank d with sum_{p in S} e_p h(beta)/2 < m are in every Z_PJ => combined locus contains cheap column
  subtori; at minimal admissible n: m/w ~ 4.5 => rank <=4 (e_p=1). So best possible conclusion: rank <= const => M=O(log W).
- Admissibility threshold: n >= c ln r intrinsic (T cancels in m/(2w)); so rank in (theta log W, n) only seen through
  DEGENERATE J: requires relative ZE (generic point of Col_J).
- TC_J = {u : <.,u>^m in TF_m(A)}  (top forms of deg-m polys vanishing on A).
Computations: n=4..6, M=6..12, all n-subsets, all base points.
## Finding 1 (pushforward)
- Restrict to L = Col (rank rho), t-coords: g_c(delta) = ghat_{pi_* c}(t), pi(k) = (<k,a_p>)_p in Lambda (rank rho).
  pi_*c has order >= m.  Combined system in t-space: Vcomb = sum_J pi_J* V_H(J) (all J containing...); all supported in
  slab region {|y_p| <= e_p w/2}. For full-rank J pi_* iso (no gain).
- Model cube rho=3 (r=3): single pi_*V(A,m) = V(ybox,m), contains Delta^alpha -> but cube profile VIOLATES pair/character
  inequality (||e_i||_omega = 2/3 < 1).  With norm ||y|| = sum_p omega_p (2/e_p)|y_p| >= 1 on Lambda\0 (character ineq),
  any Lambda-parallelepiped inside slab region has |alpha| <= sum alpha_i ||l_i|| = LP(P) <= 2w sum omega <= 2w < m.
  => NO product certificate in any combination (rigorous). Combined system stays slab-cut type in dim rho.
- Packing: S subset Lambda cap {||y||<=w}, 1-separated => |S| <= (2w+1)^rho; V(S,m)!=0 needs (2w+1)^rho > binom(m+rho-1,rho)
  ~ (e m/rho)^rho => rho > e m/(2w+1). At minimal n (m/w ~4.2): rho >= 6 needed for visibility.
- Universal system U (all supports with LP<m) contains all sub-cluster systems; ZE_univ.
Next: honest low-rank models (rho=2,3, e_p>=2, pair-separated via LP weights), compute single vs combined vs universal.
## Finding 2
- cube (non pair-separated): single pi_*V = V(ybox,m) contains Delta^alpha  (pushforward_out.txt)
- cross-polytope rho=2,3 (pair-separated tight, e=2, r=2^(rho-1)): pi_*V = 0 for ALL J, combined 0, U=V(S_slab,m)=0: BLIND (exact)
  (crosspoly2_out.txt, crosspoly3_out.txt). Conjecture: reg(S) <= E_eps h_S(<eps,.>) (cross-polytope criticality);
  Theorem B gives only reg <= MW = Omega' + h(<1,.>)-ish. Evidence only.
- Blindness => M <= binom(rho+2,2) via D=2 support {e_x-e_b} (LP<=2, ord 3), if M >= 2 binom (Siegel height small).
- Cube-like geometries r ~ 2^rho blind => method cannot force rank below ~log2 r ~ 1.44 ln W > 0.587 ln W unless blind
  geometries carry few points (they do: M <= binom(rho+2,2) ~ rho^2/2 -> but that's O(log^2 W) -- fine for growth!)
  => dichotomy: blind => M <= (rho+1)(rho+2)/2 ; visible => ZE in column torus.
Next: experiment F (fakes/full-rank, n=4, M=8): single vs combined TC zero sets, all bases.
## CORRECTION (important)
- Cross-polytope blind profiles have M = 2rho+1 < 300 ln 4r ~ 208 rho: NO admissible sub-cluster -> Z_comb trivial.
  So Prop 5.1 "rank target unattainable" was vacuous. Rewrite: under ZE-tor(single), actual clusters with
  M > 300 ln 4r are Z_comb-blind => rank target <=> combinatorial blindness statement B(theta). Cor 4.5 = sufficient
  criterion for blindness (ratio-2 criticality via Thm 4.4); known examples too small. Prop 4.3 subtori also vacuous
  for actual clusters (rank<=3 => M<=83 < 415).
