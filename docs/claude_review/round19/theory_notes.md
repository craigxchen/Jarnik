# notes (theory agent, round 19)

KEY IDEA (found before any computation):
 For f a nonconstant rational function on X = Mbar_{0,8}, theta = df is a rational section of
 Omega^1_X(log d) (d = boundary).  div(theta) >= sum_{p=0,1}[sum_{Y int in F_p}(m_Y-1)Y + sum_{D in d, D in F_p} m_D D]
   - sum_{Y int in F_inf}(m_Y+1)Y - sum_{D in d cap F_inf} m_D D.
 Omega^1(log d) has a filtration F_4 c ... c F_8 (pullbacks of Omega^1_{Mbar_{0,k}}(log) along forgetful maps,
 for an ordering of [8]); F_k/F_{k-1} = Lambda_k line bundle, Lambda_k.beta = 32 - 2^{8-k} (k>=5), 16 (k=4).
 Line subsheaf O(div theta) -> Lambda_k nonzero for k minimal => Lambda_k - div theta effective.
 beta >= 0 on effective divisors (Keel-McKernan symmetric effective cone; or beta in cone of psi-monomials).
 => (div theta).beta <= 31 => gain(f) <= 31 - n(f), n(f) = # boundary divisors in F_0+F_1+F_inf.
 Unit pencils: n >= N(u) = #supp div(u) = (1/2)#{x in {0,1}^8: q_e(x) != 0} >= 2^{8-2}/2 = 32 by the
 multilinear Schwartz-Zippel (DeMillo-Lipton) bound for the degree-2 multilinear q_e(x) = sum e_ij x_i x_j.
 => gain <= -1 for every unit pencil, all lambda, all degrees.

TODO: exact checks: Keel consistency of beta, Lambda_k.beta, psi-monomial cone LP, N(u) >= 32 sanity,
 gain + n <= 31 on round-18 enumerations, examples.

beta-RH (exact form): div(df) = sum_{c in C}(F_c - F_c^{int,red}) - F_inf - F_inf^{int,red} (no other zeros: df=0 on V => f const on V).
 => sum_{p in P^1} (L.b - R^int_p) <= 2 L.b + Lambda_k.b,  equality for cross-ratio & 4z(1-z) (k=4).
 => gain <= Lambda_k.b - n'(f) - sum_{p notin {0,1,inf}} (L.b - R^int_p), n' = #D_T with f|D_T constant.
 composites h(g), g unit: gain <= -2^{8-k} (cases j = |h^{-1}{0,1,inf} cap {0,inf}| = 0,1,2).
 psi-monomial cone does NOT contain beta (c2): (D2,D3,D4) pairings all pair-heavy. Use L8.md Thm A(a) (KM, refereed).
 (K+d)^4: (1379,244,150). Not used.
 product unit e = (x1-x2)(sum b_j x_j)-type with 8 vertices gives N = 44 -> method gives only gain <= -13 there.
