# notes (round 20 composite)
Plan of results:
- Lemma 0 (fibre decomposition): every ideal a = d b, d S-smooth, b coprime to S, uniquely. Mult graft: th^S(a)=th^S(d)+th(b) ; rho(d)=th^S(d)-th(d)=<v(d),phi-theta>.
- Lemma 1 (S-smooth counts): S = s split rational primes, all >= P0. tau_sigma := sum_{d S-smooth, d!=1} Nd^{-sigma} <= exp(2s P0^{-sigma}/(1-P0^{-sigma})) - 1.
  Psi_S(X) <= binom(2s+L,L), L=floor(log X/log P0) => Psi_S(X) <= X^{log(3e s)/log P0} for X>=P0 (L<=2s case) ; general.
- Prop 2 (all-ideal LS automatically robust): Lambda_opt(X,K) >= max(#ideals<=X, (2K+1) floor(X^{1/2})) via rational integers (lambda_k(n)=1).
- Thm M (minimal graft): changes only S-smooth ideals; robust inequalities with coefficient growth X^{o(1)} on S-smooth survive if A>1/gamma (A>2 for LS).
- Thm P (mult graft): (a) LS all ideals any support: dual decomposition, d=1 term <= genuine Lambda, d!=1 terms <= Lambda_0(X/Nd,K) <= 2C1 (X/Nd + K (X/Nd)^{1/2}); C1 = 3(1+2pi).
   (b) bilinear with free coefficients: exact (eta=0).
   (c) polar boxes, T in {all, coprime to q, congruence mod q (q coprime S)}, sharp: diff <= sum_{d!=1} osc_d <= C x^{1/2} tau_{1/2}*(1+..).
   (d) fibre-closed families with error in [c x^gamma, C x^gamma]: diff <= 2C x^gamma (tau_gamma)(1+tau'_gamma).
- Prop E (entire-L collapse): mult graft agreeing off S; L' = L * prod E_p; poles on Re s=0 unless {lam'(pi),lam'(pibar)}={lam(pi),lam(pibar)}. => smooth all-ideal Hecke sums with o(1) error force phi_p = +-theta_p mod pi/(2k). Collapse = Thm D analogue.
- Remark: sharp primitive LS <-> composite max-scale cross-gap (badly approximable type), not collapse; inherits round-10 unclassified status.
- Strength requirement: to see minimal graft need error < Psi_S(X) <= X^{1/A+o(1)}, A ~ c0/2 => gamma(c0) < 2/c0.
Checks to run: LS constants genuine/mult/min graft; smooth Hecke sums genuine vs graft; Psi_S counts; Lambda_opt lower bound.
## turn 3 results
- Rotation issue: reduced tuple zeta=z/g is a DIFFERENT point of X_k (rotation by 1/g). So A_k must use similarity-invariant Z
  (Z in T_k = norm-one torus^k / diag, birational to Mbar^tw_{0,k+2}). Vojta must be on Bl_e Lbar_k^tw (toric heights = block norms exactly).
- k=4 family H_T = s - r_T, r123=i, r12=-6-4i, r13=-6-3i, r23=-6+5i, r1=-7+6i, r2=-3-2i, r3=-4-3i: Im P_i = -420,-180,-240;
  gcd bounded, t bounded, chord/R'^{3/7} -> 12.35 (k4verify.out).
- Rational 2-param subfamily: common real part u, x_j = 1+v_j, x_a+x_b+x_c = 1, u^2 = x_a x_b x_c (conic bundle, rational).
  (X,m,s) -> T_4 rank 3 (dominant) exact (k4subfam.out). Full 4-dim S x s -> P(X_4) rank 7 at rational pts (k4dominance.out).
- So A_4(delta) false for every delta<1/14 (similarity-invariant Z); (L_4^circ) false.
