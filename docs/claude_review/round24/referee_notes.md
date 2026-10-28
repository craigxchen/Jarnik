# referee notes
## First pass observations (theory.md)
- P3.1 table parities consistent; upper bounds from width Thm 2.5 + parity. To recompute ranks m<=5 independently (and try m=6).
- P3.2 (a),(c): ok; Q-rational => support symmetric automatically, so reality vacuous. ok.
- P4.2 lex-min proof ok. P4.3(iii) needs lambda_i>=0, simple profiles, disjoint primes; O_F(1) sign/uniformity: F may vary (Roth degrees grow) -> O_F term = O(log H(F)+deg); need it o(d_i w_i). minor.
- P5.3 prose "weaker than or equal to split ones": false per-section (f=x_a/x_b-1: tw r=1/2|..| < sp r(v)=v_b-v_a). True only vs best over iota-orbit {f, iota^*f}. minor wording.
- P6.2 ISSUE: identity G = prod Pbar^deg g(x) needs per-pair homogeneity (stated "degree deg_j jointly" ambiguous). Non-homogeneous in common factor: "fix P_1=1" changes the object (natural lift P_j ~ zeta_j carries drifting phase x_1=zeta_1/zetabar_1). That's the "drifting phase -> 0" constraint. Repair: g(x)=sum_k x_1^k g_k(x/x_1); some g_k(y)!=0; ord g_k >= t=min_k ord; width of slice <= width of conv(all c); so G-ineq implied by TL for g_k. Conclusion survives.
- iso_ref.py m=1..5 reproduces (2,3),(6,5),(10,9),(12,13),(16,15) via flint mod 2^61-1, parity-split rows. m=6,7 running in bg -> iso_ref67.out
- Thm 7.1 group-words gap: tuple (y1, y1*y2) shares primes; Prop 4.3 needs disjoint. Repair: pull back along hom phi: T^r->T^n; G=phi^*F; RHS_G <= RHS_F (Q_G ⊆ phi(Q_F)); LHS: each nonzero G-Taylor term (a) from some F_J with sum a_i lam_i >= sum_k j_k min_{i in supp w_k} lam_i; choose hypothetical config with lam(w_k)=min (generic directions) -> F-ineq implied by G's (*) at disjoint-prime generators.
- h90_ref.py: 304 tuples/3040 pairs all pass; extra: e_a==e_b iff u_a/u_b in {+-1} (unit class mod sign).
- sec_ref.py: blocks and table m<=4 orders direct; supports/parity m<=7 ok.
- Prop 4.2 has a 3-line proof: W_k filtration, W_{R_i}=all => tensor of all; index > sum w_i R_i forces c=0.
- P2.2 Reading: "distance = |det|/sqrt(Phi0 Phi0)" exact only up to factor 2 (tw) and sqrt(u_b/u_a) (split): minor.
- Inert-place congruence divisibility (twisted pts always in T_5(Z_q), reduce into mu_{q+1}^4): not in L, missing from "exhaustive" list; gain <= sum_q v_q((Cm)!) log q = O(m log m) per section of mB (c_q=0 for q>=2m by character independence on mu_{q+1}^4; weak approx since T_5 rational). Negligible iff log m = o(w).
- Uniformity in F: O_F(1) in (*) includes Taylor const (>= -t log t), log H(F); caps need these o(sum d_i w_i). Prop 4.3(iii)/Thm 7.1(2) state "holds whenever mu<=4(1-eta)" without w_0(F). minor/major.
- P6.2 repair: per-pair multihomog components g_D, choose D* with g_D*(y)!=0 maximizing upper bound; plus common-factor homog decomposition g = sum rho^k g_k.
