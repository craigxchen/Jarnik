# notes for t1rigour.md (working; not the write-up)

Plan:
1. closed form m_max >= c* d - kappa (entropy bound binom(N,K) <= N^N/(K^K(N-K)^(N-K))), c*=0.6418, kappa = ln8/g'(c*)
   g(c) = (1+c)ln(1+c) - c ln c ; concave; m-1 <= c d - ln8/g'(c) suffices.
2. m_max >= 4w <= d >= 318 ln 4r.
3. Siegel at T=1: L >= 3 ln 3.
4. W_1 = 24 d ln3/m' + 8 ln(e d C'/m').
5. Step 0: eps0 = C d e^{-W/4}/(m'+1).
6. transversal: kappa(u) <= (2+eta)^m' eta^(m'-rho_F) ; rho-cost.
7. structure thm: h_A(v) >= 2 tau_j(v) (S_w in A), tau_j >= j*ell + e. classes <= prod(ell+1) <= 2^b.
8. separation lemma: gcd bound |delta_x-delta_y| >= e^{-d_xy/2}; Plotkin avg d <= W/2 + W/(2(n-1)).
9. counterexample to cheap-span ZE-TC: unbalanced column S' with (m'-1)/2+16 < s < m'-j; rho(A2(s,j)) <= s+j (level Newton).
10. J1 recheck; Sylvester at T=1, braid fails iff 2jk < m'.
