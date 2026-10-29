# notes (parking derivations)

## N1 (definitions checked)
- slab box (killer klib.build_box): k' in {-1,0,1}^d, k_n=-sum k'; slabs (S,L): |sum_{x in S} k_x| <= L, S subset {0..n-1}
  (canon: min of S, S^c). models.py: patterns omit row 0 (same thing up to complement).
- zeta in k6_zono.py: generators ANY primitive v in {-2..2}^d (multiplicity 1 or 2), NOT only L_pair.
  M8' order-10 zonotope uses e0-e1-e2 (not in L_pair). sdc7/syl8 order-7 zonotopes use only L_pair generators.
  => two zetas: zeta_pair (generators in L_pair: same hypothesis as criterion) and zeta_all.
- reach <= truth_L <= rho (Thm 3.1 + V(A,rho+1)=0).

## N2 (product additivity) A=A1 x A2 (disjoint coordinate blocks; base point shared).
gr(K[A1xA2]) = gr(A1) (x) gr(A2) => M_e(A) = sum_{a+b=e} M_a(A1) (x) M_b(A2).
L_pair(d) = L_pair(d1) u L_pair(d2) u cross {e_x - e_y : x in 1, y in 2} ; cross l = (e_x, -e_y): A_l = (A1)_{e_x} x (A2)_{-e_y}.
sum_l M_e(A_l) = sum_{a+b} (N_a (x) M_b + M_a (x) N'_b)  [cross terms inside N_a (x) N'_b]
=> SHIFT_e(A) fails iff exists a+b=e with SHIFT_a(A1) fails (H_{A1}(a)>0) and SHIFT_b(A2) fails.
=> reach(A) = reach(A1)+reach(A2). zeta additive (projection argument; generators project into L_pair).
=> M8'^t: reach 11t, zeta 10t (rho 11t, L=3). No bound reach <= zeta + c on slab boxes. Not regime: gap fails.

## N3 (ideal of shift-intersection) I(S cap S') = I(S)+I(S') for finite sets (K^S (x)_R K^S' = K^{S cap S'}).
I(A_l) = I(A) + tau_l I(A), J(A_l) = LF(I(A) + Delta_l I(A)) contains J(A) + (d_l J(A)).
Slab boxes: I(A) = I(C) + (omega_p), I(A_l) = I(C_l) + (omega_{p,l}), deg omega_{p,l} = 2L+1-|<l,s_p>|.

## N4 (linkage form) J(A) = J(C):J(B), B = C\A tail (from H_A(e)=T_d(e)-H_B(2d-e) + J(A) subset J(C):J(B)).
Gorenstein duality in R/J(C): M_e(A) <-> J(B)_{2d-e}/J(C). For A_l subset C_l (grid): M_e(A_l) <-> mu_l J(B_l)_{s_l-e},
mu_{e_x}=k_x, mu_{e_x-e_y}=k_x k_y, B_l = C_l \ A_l.
SHIFT_e <=> J(B)_{2d-e} subset sum_l mu_l J(B_l)_{2d-e-deg mu_l} + J(C)_{2d-e}.   (TO VERIFY numerically)

## N5 (symmetry) G = Aut(A) preserving L. For an orbit O of L: cap_{l in O} J(A_l)_e = core_G (largest G-stable subspace).
Obstruction space O_e = cap_l J(A_l)_e / J(A)_e is a G-module. syl8: G = GL(3,2) x {+-1} (base point fixed!),
GL(3,2) simple => 1-dim O_6 is trivial rep => obstruction is a GL(3,2)-invariant sextic mod J(A).
Counting: SHIFT_e fails for rho-isotypic part if mult_rho(M_e(A)) > sum_orbits dim Hom_{Stab(l)}(rho, M_e(A_l)).

## N6 data (t1, t2)
- linkage identity H_A(e)=T_d(e)-H_B(2d-e) verified mod p: syl8 L=1,2.
- syl8 L=2 delta profile (L_pair): [0,0,0,0,0,0,1,34,42] (e=0..8). syl8 L=1: delta=[0,0,0,21], reach 3=rho.
- t2 EXACT: |G|=168 (perms of 7 non-base coords preserving slabs; GL(3,2)), 8 G-orbits of std sextics;
  dim J(A)_6^G = 2, dim J(A_e0)_6^G = dim J(A_{e0-e1})_6^G = 3 and they coincide -> SHIFT_6 fails EXACTLY over Q
  (invariant F in all J(A_l), l in L_pair (G-orbits + sign), F notin J(A)).
  F = orb(k4^2k5^2k6^2)[7] - 2 orb(k3^2 k4 k5 k6^2)[84] + 4 orb(k2k3k4k5k6^2)[21] + 8 orb(k1..k6)[7]  (mod J(A)^G)
- patterns of syl8 = complements of the 7 Fano lines (non-base points = Fano plane), each a 4-subset of {0..6}.

## N7 KEY: syl8 L=2 e=6 obstruction = sigma^6, sigma = k_0+..+k_6 = -k_8 (BASE COORDINATE), exact (t5).
 sigma ranges: A: {-3..3} (7 values) ; A_e0: {-3..2} (6) ; A_{e0-e1}: {-2..2} (5).  => sigma^6 in J(A_l) all l in L_pair.
 sigma^6 notin J(A)_6 (exact). Also: first-order closure J1(A,l)=J(A)+(d_l J(A))+(k_x^2) already contains it (t4, mod p);
 J1 = J(A_l) except 1 dim at (e0,e=6) and (e0-e1,e=5) (t3).
 Why sigma small: Fano: sum_p s_p = 4 sigma => |sigma| <= 7*2/4 -> 3.  General balanced profile (f_x = 1/2):
 sum_p omega_p <k,beta_p> = -2f k_n => |k_n| <= w/(2f) = w = 2L.  Base coordinate = IMPLIED SLAB.
 L=3 (M8'): sigma in {-4..4} on A; A_e0 {-4..3}; A_{e0-e1} {-4..4} (6 pts at +-4: NOT rigid) -> no obstruction. 
 MECHANISM (implied-slab obstruction): integer phi, range R values on A; if range on every A_l (l in L) has <= R-1
 values then phi^{R-1} in cap J(A_l); if moreover phi^{R-1} notin J(A) (e.g. zonotope of order R-1 with all phi(v_i)!=0
 in A) then SHIFT_{R-1} fails: reach <= R-1.  For l with <l,phi> != 0 shortening automatic; for <l,phi>=0 need RIGIDITY
 of the extreme level sets.  L_pair-nondegenerate phi (all coefficients distinct incl. base 0) have range >~ L d: useless.
 Regime Sylvester: h(k_n) = 2L(1-o(1)) (flat-spectrum IP), so if rigid & non-member: obstruction at ~4L = 2w (zonotope scale).

## N8 base-coordinate obstruction (BCO) Prop: sigma=-k_n, h=max_A sigma. (R) no k in A, sigma(k)=h, with k+e_x-e_y in A
 (only top level needed: central symmetry maps -h level of A_l to h level of A_{-l}); (N) sigma^{2h} notin J(A).
 Then SHIFT_{2h} fails (L_pair), reach <= 2h.  (N) <=> pushforward sigma_*: V(A,2h) -> line nonzero <=> exists c in V(A,2h)
 with sum c sigma^{2h} != 0; sufficient: zonotope of order 2h with all sigma(v_i)!=0 (moment = (2h)! prod sigma(v_i)).
 Liouville link: b = sum_p omega_p beta_p, <k,b> = 2 sum_{x!=n} f_x k_x, range <= Lambda(A) <= 2w (same budget as zonotope cap).
 Balanced (f_x = f): <k,b> = -2f k_n => h(k_n) <= w/(2f) = w for f=1/2. Regime heuristic: h ~ 2L - O(sqrt(L ln n)), 2h ~ 2w(1-o(1)).
 Rigidity heuristic in regime: at max level only ~sqrt(L) ln n tight patterns; each blocks 1/4 of pair moves -> rigid w.h.p.

## N9 data t7/t6: BCO explains syl8 L=2 and 3 random boxes (seed105,106 n=6; seed9 n=7: reach==2h, rigid). Not sdc7 L=2
 (h=2 not rigid; delta=22 at e=6, 7 dims spanned by powers phi^6 of ternary phi e.g. (0,1,1,1,1,1,1) range 7 -> 6 on A_l).
 Random boxes: obstruction spaces (partly/fully) spanned by phi^e, phi ternary, often with range-shortening (implied slabs).
 => GENERAL MECHANISM of failure = implied-slab powers phi^{R-1} (R = #values of phi on A), plus a residue.
## N10 SUBSET MONOTONICITY: A' subset A => V(A',m) subset V(A,m) => Bs V(A,m) subset Bs V(A',m).
 So truth_L(A) >= reach*(A) := max_{A' subset A} reach(A',L) >= max(zeta_L (zonotopes with gens in L), reach(A)).
 Removing the rigid top level of sigma may kill the BCO. TEST: syl8 L=2, A' = A cap {|sigma|<=2}.

## N11 zeta_pair exact (t11): syl8 L2: 7, M8': 10, sdc7 L2: 7, sdc7 L3: 10.  sdc7 L3: reach 10 = rho (t1, mod p1).
## N12 IP (t10) Sylvester: rigidity (R) holds iff h = 2L-1 (averaging bound) in all decided cases:
 n=8 L=2 (h3 rigid) L=3 (h4 not); n=16 L=2 (3 rigid), L=3 (4 not), L=4 (7 rigid); n=32 L=3,5 not (h=4,8), L=8 rigid (15);
 n=64 L=8 rigid (15); L=12 h in {21,22}, h_pair=21 (undecided).  h=2L-1 attained by subspace-minus-0 (L power of 2).
 => rigid top level = DESIGN coincidence (Fano line at M8). M8' (L=3): no design -> non-rigid -> no BCO.
 Revised heuristic: generic regime L (not design) -> maximizers have slack -> likely NON-rigid -> BCO absent.
## N13 monomial (staircase) sufficient criterion: for degree-compatible order <, N(S)=standard monomials (greedy
 independent columns). A_l subset A => N(A_l) subset N(A). If N(A)_e = U_l N(A_l)_e then SHIFT_e. (Prop 3.2 is this.)

## N14 t14 linkage form verified (syl8 L2, e=5,6,7): sum_l mu_l J(B_l) dims 259/286/141 vs J(B) 259/287/175; linkage
 identities for A_l in C_l hold.  t15 EXACT: sigma^6 notin J(A)_6, sigma^7 in J(A)_7, sigma^6 in J(A_l) (4 shifts);
 top level = 7 Fano lines.  t10: n=64 L=16 rigid (h=31); L=12 h in{21,22}, h_pair=21.
 Sylvester design: base at origin of F_2^k, V subspace dim m, S=V\0 -> |S cap P_u| = 2^{m-1} or 0: h=2L-1 at L=2^{m-1}.
 Total slack at h=2L-1: n/2-L; a pair move needs slack>=1 on n/4 patterns -> automatic rigidity if L > n/4.
 Subspace configs rigid for m>=2 (affine (k-2)-set {u.x=1,u.y=0} not inside linear V^perp).

## N15 t17: product additivity sanity (mod p): hex x hex reach 6 = 3+3, hex x seg 5 = 3+2, cut3 x hex 7 = 4+3; H convolves.
## N16 t16 random (unbalanced) profiles: h exceeds 2L (n=16 L=2: h=4; n=64 L=6: 14, L=10: 21, L=16: 33); (R) fails at
 n=16 (3 seeds), n=32 L=4 (2 seeds), n=64 L=10; n=64 L=6 undecided (feasibility IP timeout).
## N17 t10: n=128 L=16 rigid (h=31=2L-1, design); n=64 L=26: h in [48,49], h_pair in [47,48] undecided.
