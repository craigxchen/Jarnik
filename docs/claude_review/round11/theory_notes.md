# theory notes (round 11) -- working file

## Reformulations
- Bs(A,m) cap real = {u : k -> e^{i<k,u>} on A equals a polynomial of degree < m}.
- Tangent cone: TC = {u : l_u^m in J_m}, J = top-form ideal of I(A); J_m^perp = degree-m part of de Boor-Ron least
  space = {Phi_c : c in V(A,m)}.
- Toric part: Z_tor = union of span_R(v_1..v_s), v_j in Z^{n-1} independent, sum_j width_A(v_j) < m. All
  degenerate (sum of widths of a basis >= m since V(A,m) != 0 => s < n-1).
- Parallelepiped / lines / simplices: every low-dim configuration with "degree-m reach" has width >= m in some slab
  direction unless its directions are relations. Slicing arguments are dead.

## NEW IDEA (valuation certificate with ONE widened slab)
Prop 5.5(ii) kills unique-minimiser certificates only when h_{supp c}(a_p) < m. In Thm 4.1 there is LP slack:
cost <= wW <= mW/4. Widen the slab of ONE prime p (small log p) so that a_p has >= m+1 levels on A: extra cost
m log p, total <= mW/4 + m log p < mW/2 - O(m log n) if log p <= W/4 - O(log n).
Then a c in V(A,m) with a UNIQUE minimiser of v_p(c_k) + <k,a_p> gives F_c(z) != 0 for free (exact pi-adic
valuation), and Thm 4.1's comparison is a contradiction. => unconditional bound, no ZE, no delta!
Reduces to a combinatorial "separation lemma" SL on lattice points:
  exists k* in A, S = {k*} cup {k in A : l(k) > l(k*)}, l = <.,a_p>, integer c in V(S,m), small height,
  with p not dividing c_{k*}  (then unique minimiser automatically since other levels are higher).
char-0 core: no polynomial g, deg g < m, with g = 0 on S \ k*, g(k*) = 1.
p-adic: no g in Q_p[k], deg < m, g(k*)=1, v_p(g(k)) >= l(k)-l(k*) on S\k*.
Consistency: low-rank profiles have a_p in span of few a_q with small coeffs => h_A(a_p) < m automatically, SL
fails, consistent with Prop 5.2.
TO TEST numerically (small models), then think about proof.

## Tropical reading of SL
val_pi(w) = b_p = a_p-direction (pi-adic valuation vector of the actual point w). w in Bs => a_p in trop(Bs)
for EVERY split p. Narrow slab => subgroup t^{a_p} in Bs (Prop 5.2, S={p}) => valuation certificates impossible.
Widened slab => subgroup no longer forced. SL <=> in_{a_p}(linear system) contains monomial w^{k*} (p-adic version).
Semicontinuity: min-initial forms vanish to order m at x_inf = lim_{t->inf} t^{a_p} e (point of orbit of ell-max
face). Necessary: ord_{x_inf}(w^{k*}) >= m  ~ (ell_max - ell* >= m). Not sufficient (dim deficit).

## Height problem
Thm 4.1 needs ln|c|_1 <~ mW/8. SL certificate existence gives no height control (Cramer ~ e^{E m ln T}).
Need short c in Lambda with c_{k*} not = 0 mod p. Siegel counting doesn't give inhomogeneous / congruence.

## Boolean model (box -> {0,1}^M via k = k* + sum t_i v_i, v_i = +-e_x)
X = {t : |<t,u_q> - c_q| <= s}. SL_bool: phi on X, phi(0)=1, phi perp multilinear Poly_{<m}.
Dual: psi = (-1)^{|t|} phi has degree <= M-m and vanishes on B = cube \ X, psi(0) != 0.
Count: |B| <= 2r 2^M e^{-4s^2/M}; dim Poly_{<=M-m} ~ 2^M e^{-2x^2/M} (m = M/2+x)  => m ~ M/2 + sqrt2 s heuristically.
Single slab: deg >= |P_q| - c_q - s (restriction to subcube on P_q).
TC on Boolean slices: u in TC_A => alpha=(<v_i,u>) in TC_X = common zeros of P_phi(alpha)=sum_{|J|=m} mu_J(phi) alpha^J.
I(X) = (slab polys) in C^{cube}; J_0 = (t_i^2, l_{u_q}^{2s+1}); if J_m=(J_0)_m: square-free algebra membership
l_alpha^m in (l_{u_q}^d) in Lambda_M = C[t]/(t_i^2).  Rank-one inverse systems = subcube products (capped).

## Tangent cone negative
Leading-term (Lojasiewicz) approach needs N|v| < 1: fails near toric strata; actual delta can be e^{-W/4}-relatively
close to degenerate subspaces (dense core + outliers). Tangent cone alone insufficient.
