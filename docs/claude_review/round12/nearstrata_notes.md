# nearstrata working notes (round 12, quantitative half of T=1 programme)

Plan:
- Gamma(u) = (2/N) sum_{x<y} ln^+(1/|u_x-u_y|), u_n = 0, N = binom(n,2). Actual: Gamma(hat u) <= abar = W/(2d) + 2 ln C'.
- Lambda_A(g) = min{kappa_H(u) : |u|_inf=1, Gamma(u) <= g}. MIN: Lambda_A(abar) > eps_0 = C d e^{-W/4}/(m'+1).
- (CL_{theta,gamma}): kappa_H(u) >= gamma G(u)^theta, G = exp(-(1/N) sum ln(1/|u_x-u_y|)) >= e^{-Gamma/2}.
- qualitative: Lambda_A(g) > 0 for all g  <=>  ZE-TC^braid (real zeros of kappa_H on sphere in braid arrangement).
- duality: kappa_V(u) = dist_inf(l_u^m|_A, Poly_<m|_A)/d^m; kappa_H = dist(., V_H^perp)/d^m.
- zonotope lower bound -> uncut cube: kappa >= (m!/(2d)^m) dist(u,TC)^m ; one cut: robust gadget theorem with eta-gap clusters
  (singles, pairs, hexpair), nu >= min(2 sigma, 2L+2t, s+L) + 2a ; distance <= (n-1) eta otherwise.
  Case C odd q=s-L: hexpair + (L-1) singles + perfect matching on q-1, classes <= (q-1)/2 (excess <= L+1 removals).
- barrier: integer forms with zero set in braid but Lambda(g) <= 4^m e^{-(m-2) g N/2} (u_3 small): degree/height data cannot give MIN.
- necessary: rational L in TC_H, P(L) collision pairs, transversal order o >= m'-rho_A: Lambda(g) <= 3^m (2eta)^o with
  (2P/N) ln(n/eta) + Gamma' <= g  => need o N/P < d (approx). Collision-free L => Lambda(g)=0 for g >= g_L.
- Sylvester at m' > 2jk: k-column span collision-free, binary-fraction point has Gamma ~ 3 => MIN false there.
- separation: Markov lemma for general u with Gamma(u) <= g; cheap spans structure (S_w in A, tau_j, window).
- compatibility table: max_{s,e} (m'-2js-2e) N/P_min(s,e) vs d; full inequality at W >= 3 n ln n - 14 n.
- DONE: §1 setup + Lemma 1.1 (rigidity + Step 0 at T=1); §2.1-2.3 (duality, zonotope LB, Thm 2.4 cube, Thm 2.5 one cut).
- NEXT: §2.4 qualitative<->ineffective (Prop 2.6), §2.5 barrier (Prop 2.7), Remark 2.8 literature, Remark 2.10 zonotopes absent in regime.
- then §3 separation, §4 assembly, §5 necessary conds (Sylvester), §6 scripts.
- DONE §2.4 (Lambda_A, ZE-TC^braid, CL, Prop 2.6, barrier Prop 2.7, Remarks 2.8-2.10). Remark 2.9 promises Corollary 3.6 (one-cut CL) and §3.2 P(L) >= (n-m)^2/8 for one-cut strata.
- NEXT §3: Lemma 3.1 pair gap; Lemma 3.2 discriminant (Gamma(hat u) <= abar); Lemma 3.3 Markov for general u; Thm 3.4 cheap spans (S_w, tau_j, window);
  Cor 3.5 growth-regime radii; Cor 3.6 one-cut models satisfy MIN/CL.
- DONE §3 (Lemmas 3.1-3.3, Thm 3.4, Cor 3.5, Cor 3.6). NEXT §4 assembly (Thm 4.1, Cor 4.2 growth, Thm 4.3 hierarchy MIN<=CL<=distance form), minimal open statement.
- DONE §4. NEXT §5 necessary conditions: Prop 5.1 (collision-free strata kill MIN; transversal ceiling => o N/P < d), Prop 5.2 Sylvester at m'>2jk MIN false, Prop 5.3 compatibility at ceil(4w) (computation n_c).
- DONE §5 + compat numbers (n_A=473716/564084/947350; n_B=5185/6398/10280). NEXT: §6 computations list, §7 open points, §0 summary, final consistency read.
- DONE §6, §7, Cor 3.7. NEXT: §0 summary; then full read-through.
- DONE all sections; final read-through done (patched Remark 2.10, cluster diameter, Prop 2.7 reading).
