# survive working notes

## Key idea (found early)
* Cheb_A(u_theta) >= 2^{-m} binom(m,i) Cheb^C_A(M_i)  (translation averaging; m < q needed so 2l-m distinct mod q).
* Cheb^C_A(M_i) = Cheb^C_A(M_i o g_s) for multipliers g_s (box automorphisms) ; M_i o g_s = sigma_s(M_i) (Galois).
* integral c in V: alpha = <c, M_i> in Z[zeta]; conjugates sigma_s(alpha) = <c, sigma_s M_i>; if alpha != 0,
  max_s |sigma_s alpha| >= |N(alpha)|^{1/(q-1)} >= 1. Hence Cheb^C_A(M_i) >= 1/|c|_1.  NO LIOUVILLE LOSS.
* So kappa_V(u_theta) >= 2^{-m} binom(m,i) / (|c|_1 d^m) for every theta, any integral witness c (<c,M_i> != 0).
* kappa_H: same with c' = g_s c, projections tau^a c' all in V_H (same height), max_a |Phi_{tau^a c'}(u_theta)| >= ...
  => kappa_H(u_theta) >= 2^{-m} binom(m,i)/(|c|_1 d^m) whenever |c|_inf <= H. H' = e^{m'W/40} is huge -> passage free.
* Survival <= ln lambda_min(q) <= q tau/2 - m' ln(2d) - ln(C d/(m'+1)), lambda_min = min |c|_1 over integral witnesses.
* No a-priori bound on lambda_min: BV gives e^{r_0 ...}, r_0 = e^{Theta(q)} at regime start (m' ~ 0.63 q).
* Convolution witnesses: Phi_{c1*...*ct} = (m!/prod m_l!) prod Phi_{c_l}; range additivity: need sum_l r_P(S_l) <= j each P;
  cube gadget ratio order/avg-halfrange = 4q/(q+1); needed m'/j ~ 8. Factor 2 short. Gadget = mini-box (self-similar).
* V(A,m) = V(cube,m) cap R^A ; V(cube,m) spanned by tensor mu^alpha, |alpha|>=m (mu^0=delta_0, mu^1=(d1-d-1)/2, mu^2=(d1+d-1)/2-d0).
* Number of patterns: x in q+1 patterns, pair both in (q+1)/2, separated by q+1.

## Regime budget (regime_budget_out.txt)
q=3251: tau_D(c0=.184)=296.1, budget ln|c|_1 < 463237 = 142.5 q ; witness with O(1) entries needs only tau >= 13.3.

## Zonotope obstruction (L2 version)
v integer, S_a = v(P'_a), S^-_a = v(P_a). ||S||_2^2 = (q+1)(K^2 q/(q) ...) : for ternary v: ||S||_2^2 = ||S^-||_2^2 = (q+1)(K^2+N)/4,
N=||v||_1, K = sum v. general integer v: ||S||_2^2 >= (q+1)||v||_2^2/4.  h'(v)=sum_P |sigma_P(v)| >= (q+1)||v||_2^2/(2 max|S|).
Zonotope with m generators in A needs sum_i h'(v_i) <= 4jq.  Rigorous: h' >= (q+1)/2 gives m <= 8jq/(q+1) < 8j <= m'. (barely!)
e_x, e_x-e_y: h' = q+1 exactly.

## toy plan: exact witnesses via fmpz nullspace + LLL ; <c,M_i> in Z[zeta] exactly ; LP Cheb for comparison.
