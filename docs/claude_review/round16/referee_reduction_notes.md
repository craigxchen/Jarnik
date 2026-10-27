# referee notes (reduction.md)
N1. Thm 2.1 (1)-(5): re-derived, correct. (3) sup in SB not needed; values multiples of 1/U.
N2. Cor 2.2: useful direction AC_2 => MIN_inf lossless. "interchangeable for growth question" / "shift W by O(d)":
    NOT proved; Step 0 level eps_0 fixed by cluster; room = gap to ceiling, not slack of implication. G depends on W too.
N3. Prop 3.1 ok. Cor 3.2 arithmetic ok.
N4. Thm 4.1: BV Thm 2 statement matches my memory. Hadamard ok. counting ok. Phi(V)=J_m^perp ok, ker = V(A,m+1) ok.
    "r_0 = e^{o(d)}" needs ln r = o(d); regime allows ln 4r up to d/316 -> m' ~0.64 d, r_0 = e^{Theta(d)} (check vs 3^d).
N5. Prop 5.7 ok (d_xy = sum log p |a_px - a_py|, a diff = (e_p/2) Delta). Thm 5.1(b),(c) arithmetic ok. Cor 5.6 ok (asymptotic).
N6. MAJOR overclaim: "every probabilistic version same up to e^{O(d)}" — only SB=>Psi=>Cheb direction proved; converse for
    fixed p (most k) and for cube threshold mu(A) NOT proved; toys: cube{-1,0,1}^2 minus 1 pt at m=4 (rho=3): Cheb_A>0 on G,
    but SB_{mu_A}(u;0)>=8/9. Cube form: B=cube minus 3 pts kills degree-3 forms; A 5-6 point set with TC(A,3) small.
N7. (P1) "holds wherever certificates of 5 apply" unproved; certificates give kappa_H lower bounds on G_1, not span+norm.
    V_H as Q-span vs Phi(V) R-space nit. Cube (P1) provable with K = sum_alpha 2^{n_2} = e^{O(m ln d)}.
N8. r_0 = e^{o(d)} needs ln r = o(d) (regime edge ln4r = d/316 gives m' ~ m_max, r_0 ~ 3^d/8).
N9. Remark 4.2 "s <~ 0.15 n" should be ~0.15 m' n (understated, harmless).
N10. n=2901 is exact threshold (README: 2900 for r=n), but not in "d >= 316 ln 4r" (316 ln 11604 = 2957). nit.
N11. Prop 7.1 vacuity: bounds F(D,eps,tau) only; MNV vacuous anyway (rank^{1/(4D+1)}). Q_S = l_{1_S}^{2j+1}+lower, u=1_S not in G.
N12. exact checks: cube formula 63/0, product with irregular Y 21/0; toys for N6 confirmed (A one-cut d=2 m=3 7pts, B 7pts, u=(1/2,1)).
N13. MAJOR: §8 splitting (R_H) <= (R_inf)+(P2) false as written: (P2) loses e^{-o(n ln n)}, (R_inf) has no margin. Need matched margin.
N14. (P1) cube provable: c' = sum a_alpha c_alpha/(m! 2^{n1}), K = sum 2^{n2} = e^{O(m ln d)}.
N15. lit: KS Compositio 161(12) 3089-3139 confirmed; JKSW 2505.23335 confirmed; BV statement shape confirmed (numbering unconfirmed);
     O'Donnell 10.21 constant not confirmed (secondary source: (10q/alpha)^d form).
N16. r <= W/ln5 (sum e_p log p <= W, split primes >= 5) so ln r = O(ln n) when W = O(n ln n): r_0 = e^{o(d)} ok in growth range.
