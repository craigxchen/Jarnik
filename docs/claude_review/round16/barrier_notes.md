# barrier notes (working)
## Read so far
- R12 README, R13 README §0-5, steelman S3/S4, skeptic, fake.md §1-2 + referee_fake §0,1,7.
- Box: A = {k' in {-1,0,1}^d : |k(P)| <= j for patterns P}, base coord free. <k,beta_p> = +-2 k(P).
- u in TC(A,m) iff l_u^m in J(A) (leading-form ideal); nu_A(u) = max{i: l_u^i notin J}.
- Fake Thm 2.1 (corrected): nondeg u* (all lambda.u* != 0), get fakes within theta e^{-W/4}; barrier scale needs (2.2) for relevant lambda with ||lambda||_1 <= r/rho1+1.
- Class-3 bound in Thm2.1: Prob(bad_lambda) <~ (n/(theta' C)) e^{W/4 - tau sigma_lambda/2}. So (2.2) only needed when tau*sigma_lambda/2 <= W/4 + O(log), i.e. sigma_lambda <~ r/2 ("light characters"). OBS 1 (to prove).
- AC(u) (steelman) is BOX-INDEPENDENT: only u, m', d, delta = d^m' eps_0. So AC(u) true => no density->=1/2 subset S of T has near-TC point at u. OBS 2.
- Toy boxes have density < 1/2 (sdc7 L=2: 757/2187 = 0.346): AC not refuted by sdc7 zeros directly. Check others.
- Zonotope certificates cap at zeta <= 2w < m' in regime; lower bounds on nu beyond 2j-ish unavailable.
- Perturbation near rational L in TC: kappa <= 4^m' eta^{o_L}; need o_L N2/P(L) >= d(1+o(1)) for counterexample => P(L) <= m' n/2 (collision-poor).
## Decisions (after reading)
- Densities: sdc7 L2 0.346, pal7 L2 .394, M8 .474, M8' .858, sdc11 L2 .053, pal11 L2 .068 (|A|/3^d from sizes_out.txt).
- New rigorous items planned: (P2.1) AC box-free; (T2.3) refined transfer: only light chars (sigma<~r/2) need Diophantine, |lam.u|>=c e^{-tau sigma/2};
  (P3.1) products: nu additive (self-contained proof), collision-free u has nu>=sum 2min(j,|X_b|);
  (T3.3) cyclic-symmetry circle dichotomy: g of order N>2m acting by rotation on real 2-plane W => kappa_H on circle flat within sqrt(m+1); Galois => all-or-nothing over a.
  Applies to sdc_q (translations, N=q) and Sylvester (Singer cycle, N=d=2^k-1).
- fatpoint already has MIN^char. My T2.3 matches it: barrier test for MIN^char exact.
- J4: M8 TC real pts on collision hyperplanes or 4-term (additive quadruple) hyperplanes, e.g. d0+d7=d1+d6. Quadruples are LIGHT (sigma~0) for Sylvester. Compute: collision-free TC pts on quadruple hyperplanes in M8?
