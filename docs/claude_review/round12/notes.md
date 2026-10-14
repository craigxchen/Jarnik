# working notes (round 12 ZE-TC)

## nu-function
nu_A(u) := max{i : l_u^i not in J_i}, l_u = <k,u>.  TC(A,m) = {u : nu_A(u) < m}  (R2 + monotonicity l^i in J => l^{i+1} in J).
- subadditive: nu(u+u') <= nu(u)+nu(u') (binomial).
- nu(v) <= h_A(v) for integer functional v (takes h+1 values).
- if u in span of integer functionals V (rational subspace L), nu(u) <= rho(pi_L A) (image point set).
- products: A = A1 x A2 (disjoint coords): J = J1 + J2, nu(u1,u2) = nu1(u1) + nu2(u2) (tensor of graded algebras).
- transversal order along L (rational subspace in TC): coefficient of eta^i in Phi_c(u_L+eta v) vanishes for i < m - nu(u_L).
- zonotope lower bound: zonotope k0 + sum [0,L_i]v_i in A with all <v_i,u> != 0 gives nu(u) >= sum L_i.

## one cut A_S(L) = {k in {-1,0,1}^s : |sum k| <= L}
- rho(A_S) <= s+L via Newton along sigma levels; slice Lambda_j has degree <= s-|j| (indicator prod_{Z}(1-k^2) prod_N (k^2-k)/2 for j>=0).
- zonotopes: singles (deg2, range1), pairs e_x-e_y (deg2,range0), hexpair (deg3, range1, 2 coords), triangle (deg3, range0, 3 coords) => zeta = s+L (L<s, (s,L)!=(1,0)).
- conjecture nu_S(u) = min(2|supp u|, 2L+2t(u), s+L), t(u)= s - max value multiplicity.
- => cheap-span ZE-TC FAILS on A_S x cube for m in (s+L, rho(A)] (L <= s-2): u=(generic u_S, q coords) has nu = s+L+2q, cheap cost 2s+2q.
- correct version: ZE-TC^rat: real TC in union of rational subspaces L with nu(L)<m.

## computations so far
- onecut_test: formula nu_S = min(2supp, 2L+2t, s+L) verified mod p on 8157 cases s<=6 (all u in {-1..2}^s for s<=4, random s=5,6). 0 bad.
- twocut_rho + gadget_ip: two-cut blocks (d<=7): rho = gadget zonotope bound for all L>=1 cases; gaps only L=0 (direction (1,-1,1) gadget missing). => two-cut blocks zonotopal (evidence).
- TC equations: span{Phi_c : c in V(A,m)} has dim H_A(m) (degree-m part of inverse system), NOT dim V(A,m).
- Hadamard T=1: sum_p <k,beta_p>^2 = n|k|^2 -> ball |k|^2 < w^2 = 4L^2. H16: |A| = 241 (L=1), 213937 (L=2), 3608377 (L=3).
- Froberg heuristic (froberg_out.txt): needs d>=15 for rho>=8 at L=1 with r=d generic forms.
- zeta upper bound via LP over directions D = A-A with functionals F: sum_v L_v |<v,f>| <= width_A(f).
- KEY: ZE-TC^deg is automatic when zeta(A) >= m (one zonotope c gives Phi_c = m! prod l_{v_i}^alpha_i). One cut: zeta=rho, so automatic for all m. Two cuts: zeta=rho (evidence d<=7) => automatic. So cuts-induction via zonotopes can never reach zonotope-free models; ZE-TC^deg only has content when zeta < m.
- testing ZE-TC^rat with integer u is meaningless (rational u lies on rational line). Need algebraic special points: restrict to 2-dim rational M = span(v1,v2), u = v1 + t v2, gcd_c Phi_c(v1+t v2) over c in V(A, nu(M)); irrational real roots => ZE-TC^rat fails.
- twocut_formula: atom library incomplete (rational lines cheap), not meaningful.
- cyclic_hilbert.py validated on q=7 (M8): H matches direct. mod-p rho >= rho_Q (so mod-p rho<=15 proves rho_Q<=15).
