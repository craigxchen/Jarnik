# Round 16 working notes (reduction of MIN^pair to discrete AC)

## Key derivations parked

D1. Duality chain (H = infinity). g_u = l_u^m |_A. V = V(A,m) = Poly_{<m}|_A^perp.
  d^m kappa_V(u) = dist_inf(g_u, Poly_{<m}|_A) =: Cheb_A(u).
  For any prob. nu on A: dist_{L2(nu)}(g_u, Poly) <= Cheb_A(u)  (sup >= L2).
  Uniform nu: Cheb_A(u) <= sqrt(U) dist_{L2(mu_A)}(g_u,Poly)  [ |c|_1 <= sqrt U |c|_2 ].
  So Psi_A(u) := dist_{L2(mu_A)}^2 satisfies Psi_A <= Cheb^2 <= U Psi_A.
  Psi_A(u) = sum_i Phi_{c_i}(u)^2 / U for c_i ONB of V in l2(A)  (sum of squares form of degree 2m).
  Small-ball forms: sup_A|Q| > delta  <=>  mu_A(|Q|<=delta) < 1. Probability threshold irrelevant
  up to e^{O(d)}: mu_A(|Q|<=delta) <= 1-p => E Q^2 >= p delta^2; and E Q^2 >= delta^2 => sup >= delta.

D2. Product formula. X x Y product measures: gr filtration is tensor; top comp of (l_x + l_y)^m is
  sum_i binom(m,i) T_i(u_X) (x) T_{m-i}(u_Y), so Psi_{XxY}^{(m)}(u) = sum_i binom(m,i)^2 Psi_X^{(i)} Psi_Y^{(m-i)}.
  Ternary single coord: p1 = k (norm^2 2/3), p2 = k^2 - 2/3 (norm^2 2/9).
  Psi_cube^{(m)}(u) = sum_{alpha in {0,1,2}^d, |alpha|=m} (m!/alpha!)^2 (2/3)^{#1} (2/9)^{#2} u^{2 alpha}.

D3. Graded form. theta: R[k]_m / J(A)_m -> H_m(A) := Poly_{<=m}|_A (-) Poly_{<m}|_A is an isomorphism
  (top component of a form phi). Psi_A(u) = q_A(l_u^m mod J(A)_m), q_A pos. def. quadratic form on the
  h_A(m)-dim space. Split: lambda_min(q_A) (conditioning) x algebraic distance of Veronese point to J(A)_m.

D4. H passage. Lambda = V cap Z^A, rank D. Minkowski II + Vaaler (central sections of cube vol >= 2^D):
  prod lambda_i <= det Lambda = det(V^perp cap Z^A) <= U^{r0/2}, r0 = dim Poly_{<m}|_A; lambda_i >= 1.
  => codim(V_H in V) <= r0 ln U / (2 ln H). (theory.md L1.1 had extra U ln2/ln H term.)
  ker(Phi on V) = V(A, m+1). Pointwise passage NOT implied by codim.
  Psi_{A,H}(u) = Psi_A(u) - |P_F g_u|^2, F = V (-) V_H.

D5. Exceptional set. u in G => dist(u,L) >= (1/2) exp(-(1/2) min_{(x,y) in coll(L)} eps_xy).
  So the AC exceptional set is the per-pair braid neighbourhood Br_eps = U {|u_x-u_y| < e^{-eps_xy/2}}.

D6. Monotonicity: B subset A => Cheb_A >= Cheb_B; Psi_A >= (|B|/|A|) Psi_B; TC(A) subset TC(B).
  One cut: slice B = {k(S)=0} x cube_R subset A_{S,L}; J(slice) = (k^3, sigma_S) (SLP, L=0).
  Profile-independent subsets of A are inside S_w = {|k|_1 <= 2j}: rho(S_w) <= 4j < m'. useless.

D7. No-go for cube AC: Q = prod_{|t|<=j}(sigma_S/2... ) vanishes on slab of mass >= 1-1/(2r);
  degree 2j+1 ~ w. Any AC bound in terms of (degree, influences, variance) is vacuous at degree >= ~ sqrt(d ln r).

D8. Hypercontractive thin-cut theorem: E[Q^2 1_A] >= ||Q||_2^2 (1 - 3^{3m/2} mu(A^c)^{1/2}) for deg Q <= m
  (O'Donnell Thm 10.21, lambda = 1/3: ||f||_4 <= (sqrt3 * 3^{1/4})^k ||f||_2 = 3^{3k/4}). So
  Psi_A >= Psi_cube (1 - 3^{3m/2} eta^{1/2}) when mu(A^c) <= eta.  [check constant]

## Scale
 needed: Cheb_A(u) > d^m eps_0 = d^{m+1} C e^{-W/4}/(m+1). Room vs Chebyshev ceiling 2^{1-m} d^m:
 (3/4 - o(1)) n ln n in log.
