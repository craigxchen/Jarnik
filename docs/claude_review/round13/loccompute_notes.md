# loccompute notes (working; not the write-up)

## Setup decisions
- Work in R = Q[x_1..x_d]/(x_x^3): x^3 in J(A) always (coordinate width 2) and in J_loc (coordinate width powers).
  dim J(A)_i (in R) = N_i - H_A(i), N_i = coeff t^i in (1+t+t^2)^d.
- Models (from zetc §4-5): Syl(2^k,L) = {k in Z^n_0: |k_x|<=1 (x != base), |k(S_a)| <= L}; base unconstrained.
  Paley-12: coords Z_11, base inf, columns i + D, D = {0,1,3,4,5,9}, |k(i+D)| <= L. H(L=2)=[1,11,66,275,880,2265,4785,7525,4345].
  |Syl(16,2)| = 213937, |Syl(16,1)|=241, |Paley12 L=2| = 20153.
- N_i: d=7 max 393; d=8 max 1107; d=11: N_8=14355; d=15: N_8=157950, N_11=827190.
- Certification: J_loc subset J(A) over Q (generators exact over Z, J(P) by fmpz nullspace).
  lower bound on dim_p J_loc (any subset of rows) = sound;  N - rank_p(eval) >= dim_Q J(A). Equality of cumulative => equality.
- Non-equality over Q: exact fmpz rank of generator rows in degree i (< dim J(A)_i exact).
- Isotypic: Z_q acting by cyclic shift on coords, p = 1 mod q. R_1 comps y_chi = sum_l w^{-chi l} x_l; y_chi'' J^{chi'} in comp chi'+chi''.
- Wrong idea avoided: k -> -k acts as (-1)^i scalar on R_i (no splitting).

## Log
- [syl8] colb_syl8_out.txt: Syl(8,1) m=3,4: J=J_w. Syl(8,2) m=8: J_w deficient deg 6 (49 vs 70 in R), 2-col adds nothing
  (2-col images are full 5x5 boxes), 3-col: ONLY the 7 LINE triples {a,b,a+b} matter (image 57 pts, rho 6, new gens deg5:6, deg6:3);
  independent triples have full box image (125). 7 lines x 3 = 21 = the gap. Syl(8,3)=M8' m=11: J = J_w (R=2) in ALL degrees <= 11 !
- m8_lines: pairs (col, half-sum v=(y_a+y_b+y_{a+b})/2) give only 2 per line (63 vs 70); (half,half) nothing. Lines need 3 dirs.
  Mechanism: y_a+y_b+y_{a+b} = 2 k(complement of U) -> parity sublattice of the image.
- IDEA: Sylvester s-dim column subspaces satisfy sum = 2^{s-1} * (0/1 vector): mod 2^{s-1} relations -> maybe 7-col, 15-col families
  needed in higher degree (b0 ~ 2^{m/(2L)}?). Test on Syl16 / projected images P_7.
- [paley L=2 m=8] cyc_paley_L2_out.txt: H_p = zetc's. J(A) in R = [0,0,0,0,0,12,132,1517,10010]; J_w(R=1) = [..,12,132,792,3300];
  columns b=2,3,4 add NOTHING (all 2,3,4-col images are full boxes 5^b). Deficit 725 at deg 7 (= 2L+3), 6710 at deg 8.
- paley_images: b=5: orbits with |P| = 1031/1041 (index-3 sublattice: mod-3 relations, 6 orbits), 3113/3123 (corners cut by ball).
  Column matrix SNF [1^5, 3^5, 6]: 6-dim space of mod-3 relations (min support 5) + the all-columns relation sum y = 6 sigma.
- Congruence lemma (easy): if sum c_q y_q = 0 mod M on A and A contains all e_x, then z = sum c_q beta_q / M is an integer direction.
  Pullback ideal of a family depends only on the Q-span of the family (J(TP) = T^{-*} J(P)). So 'b' = dim of span.
- jfast.py: J(P) via nullspace mod 62-bit prime + ratrec + exact verification (sound).
- paley_b5_out.txt: + mod-3 congruence 5-col families (6 orbits, |P| 1031/1041): deg7 902 (vs 1517), deg8 7546 (vs 10010). Only the
  all-ones-mod-3 family (0,1,2,6,9) gives deg-7 gens (11 per family). Deficit 615 at deg 7 persists.
- Syl(16,1) m=3 (rho=3): J(A)_3 = 560 in R; J_w 46; + lines 155. BALL EFFECT model (<=3 nonzeros): coordinate families matter.
- cyclib patched for non-free orbits (q=15) and per-comp H (global percomp).
- paley_syz: ideal of the 12 quintics in R: dims [..,12,132,792,3300,10120,23353]; free counts 10560 at deg 9 -> 440 NON-KOSZUL
  syzygies in degree 9 (in R) -> candidate deg-7 leading forms (S-polynomial drop by 2: -5 sum c_p y_p^3 + sum d_l x_l).
  Mechanism = 'cut-local' (J(C cap S_T)), invisible to projections when images are full boxes.
- Regime thought: for REAL profile columns, y_p has huge image -> no width power of columns at all; J(A) low degree = (x^3)+integer dirs.
  For Sylvester (integer) profiles, Froberg count (1+t+t^2)^d(1-t^e)^P positive far beyond 4w -> 'expected' no non-Koszul syzygies.
  But Paley shows actual special forms beat the expected count (expected deg-9 coefficient 9295>0 yet 440 syzygies).
- paley_bgen: near-full 5-col images (3113/3123) give only y^5 up to deg 8.
- paley_cert.py: exact invariant-component certificate (fmpz rank of P_0(L_7) rows; exact vanishing invariant polys verified on all A).
- paley_cert_out.txt: EXACT: dim_Q P_0(J_w(R=1))_7 = 72 (fmpz); an exact invariant degree-7 vanishing poly (max|coef| 48, 384 top
  coeffs, verified on all 20153 pts) has top form outside => J(A)_7 not in L_col(4)_7 = J_w over Q. (L_col(4)=J_w exactly:
  all <=4-col images full boxes; J(grid) = (univariate products)).
- paley_bgen: b<=5 all column families: deg7 902 vs 1517 (mod p). b=6,7 running.
- three-cut d<=6: many b0_col=None (even all 3 columns insufficient) -> need coordinates in families; mixed.py: b0_mixed = 2,3,4
  on samples. mixed_syl8_out: M8 b0_mixed=3, M8' b0_mixed=1. m8_cert_out: EXACT dim_Q L_mixed(2)_6 = 49 < 70 <= dim J(A)_6.
- syl16 L=2 m=6 (lines, planes) running; paley_mixed b<=4 running.
- paley_pair_block_out: two-cut (3,3,3) L=(2,2) d=9 (Paley pair cut box): J = J_w all deg <= 8 (exact sandwich) => pair cut
  boxes contribute nothing beyond J_w(A). Paley mixed (cols+coords) b<=3: all images full products -> nothing.
- cuts_three_d6: 953 blocks (= zetc count): b0_col=1: 827, =3: 27, None(all 3 cols insufficient): 99. Never b0=2.
  gap rho - firstfail in 0..4. two-cut d<=7: 161 blocks all b0=1 (J=J_w up to rho+1).
- Paley deficit per component ~66 (= dim R_2): invariant comp J^0_7=137 (mod p kernel count 151-14), L^0_7=72.
- IDEA: regime has m'/rho -> 0; failures in models at i >= rho-4. Track (rho - i_fail).
- syl16 L=2 m=6 run OOM-killed (hilbert_cyc 10GB). New hilbert_lean.py (point-rep column subsets: sound, rank lower bd; rref
  pivots), validated on Paley. Relaunched syl16.py 2 7 (width dirs: column, all-ones, 35 line half-sums, 15 plane quarter-sums).
- paley mixed b<=4: deg7 792 (no change), deg8 5335. paley_cert cong: EXACT dim_Q P_0(L_col(5)_7) = 82 (vs 72 J_w).
- syl16_L1_out: Syl(16,1) (ball model) mixed (cols+coords) b0 = 3 (b<=2: 330 vs 560 at deg 3, mod p).
- line_lattice_out: abstract full line image P_L = [-L,L]^3 cap {sum even} (= Sylvester line image when cosets large, e.g.
  Syl16 L=2): rho(P_L) = 3L, J(P_L) != J_w(P_L) first at degree 3L (L=2..5): non-width line generators at 3L < 8L.
  (regime: 3j < m' ~ 8j). Whether absorbed by J_w(A) in d dims: M8' says yes at L=3 for k=3 (small image). Syl16 L=2 test.
- three-cut ratio table (ratio_table.py): non-local fraction does not trend with rho/(2Lmax); mixed b0 <= 4.
- paley_cert (J_w): rigorous. cong version: ratrec/CRT lifts fail (big denominators); saved poly's top lies IN L_col(5)^0.
  -> paley_cert3.py: Dixon exact solve on pivot rows/cols (E[pr,pc] nonsingular mod p => over Q), verify on all A.
- plane_lattice.py (abstract Sylvester plane image Lambda_3 cap [-L,L]^7, index 32): L=1 (57 pts) J = J_w. L=2 (2605 pts) running.
- paley_cert3_cong_out: RIGOROUS: J(A)_7 not in L_col(5)_7 over Q; dim_Q P_0(L_col(5)_7) = 82 exact; 55 exact vanishing invariant
  polys (Dixon, verified on all A) => dim_Q P_0(J(A)_7) >= 137 (= mod-p upper bound 137 via cumulative 151: exact equality).
- plane_lattice L=2 (R=1): |P|=2605, H_P=[1,7,28,84,210,455,812,756,252], rho=8; J(P)=[..,7,112,960,2751,5005];
  J_w=[..,7,49,232,840,2485] fails deg 6; J_w + 7 lines = [..,7,112,582,1862,4459] fails deg 7 => PLANES needed beyond lines
  (mod p evidence; image ring). Running R=2.
- count heuristic: rho(P_s) ~ (b!/index)^{1/b}(2j+1), b = 2^s-1, index EXACT 2^{(s-2)2^{s-1}+1} (=2,32,131072 for s=2,3,4): s=2: 2.9j, s=3: 4.1j, s=4: 5.9j, s=5: 8.3j, s=6: 11.8j
  -> Sylvester subspace families up to 15 columns may matter below 8j (heuristic only).
- paley_bgen67: L_col(6) = [..,902,7612] (no change at deg 7). b=7 running (slow).
- paley b=7 killed (est. >10h). three_mixed done: b0_mixed 2:6, 3:81, 4:39. three_cert: exact dim_Q L_col(3)_4 = 15 < 18.
- paley_cong_check / paley_nearfull_check: ALL OK -> Prop 4.4(ii) rigorous. plane R=2 same as R=1 (582 < 960 at deg 7).
- two-cut: 318/402 blocks done (d<=5 all, a=0 d<=7, part a>=1), all b0=1; two_a1 still running.
- Syl16 L=2 (H 1h, hilbert_lean): J(A)=[..,15,540,9769]; widths 225 at deg 6; + lines 540 EXACT equality deg<=6;
  widthscan: only columns among 7.17M R=1 dirs have |image|<=5 (exact); mixed pairs all 435 grids -> exact b0>=3.
  deg 7: lines 6226, + planes 9691 (deficit 78), + non-plane 3,4-col (all grids or in planes) 9691. Running 5,6.
  First syl16 run OOM-killed at planes (pullback cache); syl16_planes.py lean (orbit reps).
- syl16_rest56: 5-col non-plane families (180 non-grid orbits) add nothing at deg 7 (9691). 6 running. syl16_planesub.py: 4,5,6-subfamilies of plane (without full plane) to see whether 7-col needed.
- syl16_planesub: in-plane 4,5,6-subfamilies add NOTHING at deg 7 (6226) -> full plane needed. noplane5: non-plane 5-col add
  nothing (6226). rest56 & noplane6 OOM-killed at b=6 (generator accumulation). syl16_inc.py: incremental deg-7 merge, running.
- syl16_inc6: L_col(6)_7 = 6226 (all 328 non-plane 6-col add nothing); + planes 9691.
