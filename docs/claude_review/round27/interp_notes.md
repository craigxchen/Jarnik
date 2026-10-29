# interp lane notes (scratch)
## Definitions recovered (ze.md §1, referee_ze)
- T=1 slab box: A = {k in Z^n_0 : |k_x|<=1 (x<n), |<k,beta_p>|<=w all p}, w = sqrt(8 d ln 4r), d=n-1, D=2d.
  identify k with k' in {-1,0,1}^d (k_n = -sum k'). |A| >= 3^d/2.
- orders m in [ceil(4w), m_max]; m_max: 8 binom(m+d-1,d) <= 3^d. (n=2966: w=471.7, ceil4w=1887, m_max=1906)
- V(A,m) = Ann(Poly_{<m}|_A). V_H = span of integer c in V(A,m') with |c|_inf <= H = U^2 (U = ? find).
- Lambda(S)=sum_p omega_p h_S(beta_p); on A, Lambda(c)<=2w. Rigidity: Z in Bs(V_H) in regime (d>=316 ln 4r, W>W_1).
- Certified for c: P_c(Z^M)=0 if Lambda_M(c) < ord(c) - O(E).
## First observations
- cube form along directions v_1..v_m: P = Z^{a0} prod (Z^{v_i}-1); Lambda(cube) = sum_i Lambda(v_i), Lambda(v) = 2X_v/W,
  pair/character Liouville: X_v >= W/2 - O(log(|v|C)), so Lambda(v) >= 1-o(1): Lambda(cube) >= m(1-o(1)) > 2w for m>=4w.
  => no m'-dim cube (coordinate or general direction) lies in A. Cube forms are exactly Liouville-critical (= product of pair bounds).
- On {-1,0,1}^d: Z^k = prod_x (1 + alpha_x k_x + beta_x k_x^2), alpha=(Z-1/Z)/2, beta=(Z+1/Z)/2-1.
- slab polys g_p(k)=prod_{s in 2Z cap [-w,w]}(<k,beta_p>-s) vanish on A, deg ~ w+1 < m'.
## Idea log
- (S) SHIFT LEMMA: Y in Bs(V(A,m)), p deg<m with p=Y^k on A. For v: q=p(.+v)-Y^v p vanishes on B_v=A cap (A-v); leading form (1-Y^v)p_e.
  If I(B_v)_{<m}=I(A)_{<m} then Y^v=1. Fails for ternary boxes (k_x(k_x+1) of deg 2) and slab boundaries (slab poly with one fewer root).
- binomial-product forms = exactly forms with no nondegenerate zero; cost >= ord(1-o(1)) => no single certified form works (structural).
- grids in A have total degree <= 2w(1+o(1)) (pair Liouville), so "sub-configuration" lower bounds for deg_A(Y^k) reach only ~2w. Need m in (2w,4w].
- I(A) = I(cube)+(g_p) exactly (subscheme of reduced 0-dim scheme is reduced).
- (MS) MULTI-SHIFT CRITERION (proved, primal+dual): L subset Z^d\0, A_l = A cap (A-l). If for all e<m:
   (SHIFT_e)  cap_{l in L} J(A_l)_e  subset  J(A)_e   [equiv: M_e(A) = sum_l M_e(A_l), M_e(S)={Phi^{(e)}_c: c in V(S,e)}]
  then every Y with Y^l != 1 (l in L) is NOT in Bs(V(A,m)); and TC(A,m) subset cup_l l^perp.
  Tangent identity: c in V(A_l,e) => T_l c - c in V(A,e+1), Phi^{(e+1)} = (e+1) <l,u> Phi^{(e)}_c.
  M_e(S) = degree-e part of de Boor-Ron least space of S (= J(S)^perp under differential pairing); D-invariant.
  cube {-1,0,1}^d: cap_l J(A_l) = J(cube) in deg < 2d (works). 1-D interval / grids / lower sets: fails only at e = reg.
- one-cut box: pair directions e_x-e_y free => grids of degree ~d. Regime: all directions cost >=1 => grids <= 2w(1+o(1)).
- TODO: compute SHIFT on small slab boxes (hexagon d=2; Hadamard d=7,8 etc.). Find failure degree vs grid cap.
## Computations (interp_checks/)
- t0: hexagon (d=2, one cut): SHIFT e<=2 holds, fails e=3=reg. V(A,3) = alternating hexagon, P_c = (Y1-1)(Y2-1)(Y1-Y2)/(Y1Y2): ZE^all(A,3) true (zonotope). cube3: SHIFT up to 5, fails 6=2d.
- t1: syl8 L=1 (|A|=57): SHIFT e<=2, fail e=3=reg.
- t2: sdc7 L=2 (|A|=757, H_A matches record 1 7 28 77 161 251 188 43 1): small shifts (98): SHIFT e<=5, FAIL e=6 (188 vs 166), e=7 (43 vs 0), e=8.
       syl8 L=2 (|A|=1037, H 1 7 28 77 161 259 287 175 42): SHIFT e<=5, FAIL e=6 (287 vs 286).
  => criterion certifies ZE^all(A,6) in both; zonotope order zeta: sdc7 6..7, syl8 <=7 (record). Not beyond zonotopes here.
- Product obstructions F = prod_{p in P} l_p^{L_p} prod_{x in X} k_x^2 need covering systems: degree >= 2j log2 n (small L) or ~2d (all L). Failures seen are NOT product type.
- TC-version of criterion: (TC-SHIFT_e): cap_l TC(A_l,e) subset TC(A,e) u (cup l^perp); u notin TC(A,e) => exists l: u notin TC(A_l,e) => u notin TC(A,e+1).
## KEY: exact inductive form (handles H-gap!)
- (SHIFT_e)* : V(A,e) = sum_l V(A_l,e) + V(A,e+1)  [equiv to M_e(A)=sum M_e(A_l)].
- Induction: c_0=delta_k. Given c_e in V(A,e), P_{c_e}(Y)!=0: write c_e = sum c'_l + r. If P_r(Y)!=0 take r; else some P_{c'_l}(Y)!=0,
  take c_{e+1} = T_l c'_l - c'_l in V(A,e+1), P = (Y^l-1)P_{c'_l}(Y) != 0. Heights: |c_{e+1}| <= 2B|c_e| if decomposition has const B.
  => ZE^all_H(A,m) with H = (2B)^m. Rigidity tolerates H'=e^{m'W/40}: need ln(2B) <= W/40 ~ 0.075 n ln n. Siegel-type B plausible.
  => route 5 reduces to (SHIFT-H) for regime boxes: purely algebraic, delta-free, H-gap absorbed.
- Zonotopal algebra (Holtz-Ron 2011; Dahmen-Micchelli): least space of zonotope lattice points = P-space spanned by products.
- t3: sdc7 L=2, e=6 with ALL 15820 shifts (|A_l|>=10): 166 < 188 FAIL (same as 98 small shifts). Genuine failure.
- t4 (2-D random subsets of {0..3}^2, 25 trials): shift reach <= zeta (binomial-product order) always; often <.
- 1-D: criterion reach = max_l |S cap (S-l)|. Example S={0,1,3,10,11,13}: reach 3; zeta probably 2 (checking t5): non-zonotopal certification.
- Sidon-type sets: reach = zeta = 1, but ZE^all likely true to high order (overdetermination) -> invisible to local tools.
- Linkage (round13 (A)): H_A(e) = T_d(e) - H_{A^c}(2d-e) (A^c = cube\A tails). CN: J(A)_e = J(cube)_e for e <= j.
  (SHIFT_e) <=> gr commutes with intersection: J(cup A_l)_e = cap J(A_l)_e.
- HYPOTHESES for reach: (H-zeta) reach <= zeta ~ Liouville cap 2w  vs  (H-peak) reach ~ peak of H_A (toys: failure right at/after peak).
  Toys can't separate (zeta ~ peak there). Regime: zeta<=2w(1+o(1)) << peak(?) .
- t5 (EXACT, sympy): S={0,1,3,10,11,13}: gcd of V(S,m) = (y-1)^m for m<=4; m=5 has non-cyclotomic factor => ZE^all true up to 4.
  shift criterion certifies m<=3; zeta (cyclotomic-product forms in S) = 2.  Hierarchy zonotope 2 < shift 3 < truth 4.
  => the criterion is genuinely NON-ZONOTOPAL (beats every single binomial-product certificate) in this example.
## Plan of write-up (decided)
- Thm A (single-form structural barrier): P with no nondegenerate zero = c Y^a prod (Y^{v_j}-zeta_j)^{m_j} (Baire + UFD);
  edges |v_j|_inf<=2; Lambda(c)=sum m_j Lambda(v_j); char ineq X_v >= W/2 - 2ln(2dC') => every c in V(A,m), m>=4w, W>8ln(2dC'),
  has Zariski-dense nondegenerate zeros. Hypothesis (CG): Lambda(v) > 2w/m for primitive v in {-2..2}^d\0 (full-rank actual; fakes).
- Thm B (multi-shift ZE): exact (primal+dual), tangent, height-graded (SHIFT-H_e): V_{H_e}(A,e) subset sum_l V_{H_{e+1}/2}(A_l,e)+V_{H_{e+1}}(A,e+1).
  Use L_pair = {+-e_x, +-(e_x-e_y)}: Z^l != 1 automatic (distinct points). No full-rank needed.
- Cor C: (SHIFT-H) at cluster's regime box with H_{m'} <= e^{m'W/40} => Z notin Bs(V_{H'}) contradicting rigidity (t1rigour Thm3.1+Rem3.3)
  => n-1 < 316 ln 4r => M = O_C(loglog R).  H-gap -> per-degree growth ln(H_{e+1}/H_e) <= W/40 ~ 0.075 n ln n.
- NOTE: H-gap NOT automatically absorbed: need small-height spanning sets of V(A_l,e), V(A,e+1) (circuits too big). Honest.
- t7 (2-D two-translate S = S' u (S'+(7,3)), S' Sidon 7 pts): EXACT over Q SHIFT e<=2, fail e=3 => ZE^all(S,3);
  zetaQ(S)=2 (complete exact search of Q-forms with only degenerate zeros = +-Y^a prod Phi_k(Y^v)). Non-zonotopal, 2-D.
- pair shifts L_pair (56 shifts) on sdc7/syl8: identical reach (6).
- Sylvester n=2^k>=2^16, m'>2jk: collision-free cheap span K_F (k coordinate columns) subset Bs(V(A,m')) (t1rigour Prop 6.5);
  generic Y in K_F has Y^l != 1 for l in L_pair => contrapositive of Thm B: SHIFT (L_pair) fails at some e<=2jk. matches
  product obstruction prod_{p in Q} k(P_p)^{2j} over separating+covering Q (degree 2j|Q|). At m'=ceil(4w): no known obstruction
  (cheap spans in braid arrangement, t1rigour Prop 6.5(a)).
- Trivial range: (SHIFT_e) e<=j-1, L={+-e_x}, any profile (CN faithfulness J(A)_e=J(cube)_e, J(A_{+-e_x})_e=J(cube_{+-e_x})_e).
- H-gap honest: rational SHIFT for e<m'-1 + last-step (SHIFT-H): V(A,m'-1) subset sum_l V_{H/2}(A_l,m'-1)+V_H(A,m') suffices.
  Only top degree needs heights; "redundancy vs codimension" (round16 Thm4.1 codim bound). Not dissolved.
- Exact route vs fake test: exact ZE^all_H + exact rigidity: fakes (not on Z_P) irrelevant; near-zero cap theorem concerns MIN^pair formulation.
## Final (write-up complete)
- interp.md complete: Thm 2.2 (single-form barrier), Thm 3.1 (multi-shift ZE, 4 parts), Props 3.2-3.5, Cor 4.1, tests, open points.
- Corrected wording: criterion and zonotope certificates are INCOMPARABLE (random 2-D: weaker; Prop 3.3 sets: stronger).
- (SHIFT)_A implies no collision-free TC_H point at ceil(4w) => at least as hard as exact barrier question (honest).
