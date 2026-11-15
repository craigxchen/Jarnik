# bridge notes
- cube code (coords +r0,-r0,+r1,-r1,...): 0, 1, {±0,±1}, comp, {±0,±2}, comp, {+}, {-}. Each column weight 4 of 8
  => no 5-subset has a constant column => gcd g=0 for every 5-subset. 
- FF exponent for 5-subset at level c, gcd deg g: contact m = D/2-c (=L+1, even D), -K.C = 2(D-g), alpha = 2(D-g)/(D/2-c).
  alpha<4 iff g>2c. Cube: alpha = 16/4 = 4 exactly (if max pair deg = 4). Not A-points; Vojta-critical.
- random 5-subset of large balanced cluster: E[g] >= D/16 -> alpha ~ 15/4: FF version of gauss Lemma 1.2.
- Dominating families (Thm D) rank 4: 4m+1 <= 2D -> alpha >= 4+1/m.
- Template constant (invariant form): for rational t0, primitive normalisation of values z_j(t0) by fractional
  Gaussian gcd g(t0): C(t0) = max|z_i(t0)-z_j(t0)| / sqrt(|z(t0)| |g(t0)|) -> max|gamma_ij| / sqrt|g(t0)|.
  gamma = coeff of t^{D/2} in z_i - z_j. Invariant under rational affine reparametrisation. Only bounded denominators of t0 matter.
- Odd D, c=0: K=(D-1)/2 < D/2 so templates have C -> 0 (Lemma P). M(odd D,0) over Q(i)(t) gives M(C)>=. for all C>0.
- Golden family = Laurent template on Pell torus (units), poles at two conjugate places; condition only at one place:
  contact a/2 for pole order a, total degree 2a: like ff level c = D/4. alpha = BC/(2 mu0) = 4 too.
- Curves are sigma-invariant (Q-curves in twisted T5), so C in cone Q iff C in sigma Q. Test only cones (+W numerically?).
- cube 5-subsets (judge y): all 56 have g=0, mu0=4, alpha=4; 24 have rank 3 (lie in subtorus u0u1=uAuAc type), 32 rank 4; no cone identity.
  Pair (P,Pc) has degree 3, others 4. Judge y: C_inf = 145877.66 (local factors 2: 1.73e-4, 3: 1.52e-4, 5: 25), vs 3.18e9 unreduced.
- S4 CORRECTION: with D = pi^*boundary + E (snc, e off boundary): K_Y+D = r E, truncated FF Vojta: r m <= |S| - 1 + eps h.
  c=0: r=4 allowed for eps>=1/(2D); r=5 needs eps>=1/4+1/(2D) -> forbidden. So (R_4) = Vojta shape, not sharper.
  Theorem D (dominating) = exactly this with eps=0. Q-B's V_E(5) is the r=4 case with E (matches).
- FF-A5: alpha<4 iff primitive 5-pt cluster with K' < D'/2 (gcd removal: c' = c - g/2). Specialises to 5-pt integer clusters with C->0.
  Odd D', K'=(D'-1)/2, rank 4: expected dim 0 (rigid), alpha = 4D/(D+1).
- Large cluster M>=10 => some 5-subset g >= D*2C(M/2,5)/C(M,5) -> alpha <= 4(1-phi(M)) -> 15/4 (FF gauss Lemma 1.2).
- Cube cluster closed under coefficient conjugation (pairs 0/1, A/Ac, B/Bc, P/Pc are z, zbar). So cube templates are AXIS clusters
  (x,±y) near real axis; depths bounded; C^2 -> 8 dmax + 4 y0^2/a (families.md 3.1). Corollary S: C>=4sqrt2.
- Any conj-symmetric 8-pt cluster on squarefree n: every cut separates the 4 conjugate pairs -> cuts among the 8 transversal cuts
  = sub-code of the cube code. Record 5.754 cluster: rank 4, all 8 cuts, matching = conj pairs: EXACTLY the cube code (unequal weights).
  Golden: 4 big cuts (one-e-side) + bounded 5-layers; record 5.674: small primes with multiplicity break it (rank 7).
  10-point: rank 7, 8 cuts.
