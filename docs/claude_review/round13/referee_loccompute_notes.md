# referee notes loccompute
Plan:
- README Loc(b0): families of <= b0 ARBITRARY integer directions (with rho(pi_F A) < ceil(4w)). Author tests only columns(+coords).
  => test whether non-column directions (line half-sums z, sigma, congruence z's) give the "non-local" generators with fewer directions.
- (T1) abstract line image P=[-2,2]^3 cap {sum even}: J(P)_6 vs width powers + all 2-dim integer projections.
- (T2) M8: pairs within line spans in R.
- (T3) Paley: re-derive |A|, widths, Prop 4.4 independently; test families with sigma / congruence directions.
- Check sandwich directions; Syl16 deg-7 numbers: H_A mod p => dim J upper bound => deficits are evidence only.
- T1 done (t1_line.py, t1_line_detail.py): abstract line P=[-2,2]^3 cap even, deg 6: J=18 (exact), widths 9, +all 2-dim integer
  projections (kernel coeffs<=4) = 17 (mod p => >=17 over Q). 2-dim families (y_a, y_b+y_c)=(y_a,z) give 2 each; (y_c, y_b-y_a) 2; many others 1.
  => 8 of 9 line non-width generators are 2-local in README sense. 9th not found.
- Paley: |A|=20153 confirmed; point orbits 1833 (author says 1834); widths: 11 coords w3, 11 cols + sigma w5, 17 orbits w7 (supp 5,7,10). confirmed.
  invariant H by degree mod p gives P0 J_7 <= 137 (matches), J_5<=2, J_6<=12.
- Paley congruences: 6-dim mod-3 space (728 nonzero); e.g. y5+y6-y7+y8-y10 = 3(x6+x9): quotient z is a 2-coordinate SUM.
- paley_run1: reproduced P0 J_w=72, P0 L_col(5)=82 (mod p).
- paley_run2 (pairs from 320 candidates: cols, sigma, coords, x_a+-x_b, width-7 dirs): from J_w 72 -> 77. New: pairs
  (y_0, sigma - x_j), j in {2,6,7,8,10} (j not in D), image 29 pts (non-grid). => 2-direction families carry deg-7 gens.
- pairs from 320 candidates on top of L_col(5): stay 82. triples of {cols, sigma, sigma-x_j}: stay 82.
- M8 (m8_run.py): J(A)=[..7,70,218,315] mod p (matches); J_w(R=2)=7,49,162,294 (matches); lines => equal (matches);
  all 2-dim sublattice families inside line spans (kernel coeff<=3): 63/211/315 at D=6/7/8 -> one per line missing at D=6.
- Paley pair (y_0, sigma-x_2): image 29 pts, J(img) dims [..1,2,5,..] from D=5, rho=8 (=m) -> author's "rho>=m" not contradicted.
- Syl16: |A|=213937, line image 63 (full), plane image 2507: confirmed (syl16_base.py).
- Consequence of T1: in Syl(16,2) 2-dim families inside line spans give >= 225+35*8 = 505 in deg 6 (rigorous via direct-sum arg).
- Paley triples from 45 cands (cols,sigma,sigma-x_j,coords,sigma-y_i) on L_col(5): stays 82 (paley_t3ck_out.txt).
- started: pairs 'ckpq' (incl. sigma-x_a-x_b), triples 'cq', quadruples 'x'.
- Syl16 width-7 scan running (syl16_w7scan.py): author's deg-7 J_w uses only cols, sigma, half-sums, quarter-sums.
- !! Paley quadruples {sigma, sigma-x_a, sigma-x_b, sigma-x_c} = span{k_base, x_a,x_b,x_c} (image 125): P0 L_7 82 -> 87.
  base coordinate excluded from author's "coordinates". (paley_t4x_out.txt)
- three-cut example (S1={0,1,4},S2={0,2,3,4},S3={1,2,3,4},L=1): J(A)=J_w+pairs of {-1,0,1} dirs in ALL degrees<=7 (exact) -> README b0<=2 (author: col None, mixed 3)
- Syl16 plane image ring: J(P)_7=1002 (ub), J_w(P)=421, +lines 771 (plane-only 231, matches author's 3465/15); pairs (y_a,zV),(z_l,zV) +13; triples continuing...
- M8' (m8p_run.py): J(A)=J_w(R=2) all D<=11 exact. confirmed.
- Syl16 width<=7 scan ({-1,0,1}): supp1 w3 (15), supp2 w5 (210), supp3 w7 (1820; vanish in R), supp4 w7 (105 = affine 2-flats not thru 0!),
  supp7 w7 (120 = col - coord in S), supp8 w5 (cols), supp9 w7 (210 = col + coord), supp12 w7 (35 half-sums), supp14 w7 (15 quarter), supp15 (sigma).
  (y+-x)^7 in (y^5) in R -> absorbed. (x_a+x_b+x_c+x_d)^7 for affine 2-flat: NONZERO in R_7, NOT in author's width set. check absorption.
- Paley sigma+5 coords: 87 (nothing more). exact run in progress.
- Syl(16,1): every {-1,0,1} direction has width<=5 on A (|A|=241, ball-bound); pair test infeasible/degenerate, skipped.
- two-cut (1,2,4) L=(1,1): J(A)=J_w(R=2) all D<=10 exact (twocut.py). confirmed.
- Paley triples 89 cands (incl sigma-x_a-x_b): 82 (nothing). started sigma+2coords+col, sigma+3coords+col from 87
- paley_exact.py: EXACT dim_Q P0(L_col(5))_7 = 82 (independent; completeness: 36 near-full orbits J(P)_{<=7}=(y^5) by mod-p ub; 6 congruence exact J).
  + sigma+coords (0,1,5),(0,2,6) -> exact 83, 84 => base-coordinate 4-families give elements of J(A)_7 outside L_col(5) over Q. RIGOROUS.
- three-cut b0_mixed=4 block [0,0,2,1,1,1,0] L=1: pairs of {-1,0,1} dirs suffice all D<=rho+1 (exact).
- plane ring: all <=3-direction families from 22 cands: 771 -> 918 of 1002 (ub). 147 of 231 plane-only gens from <=3 dirs.
- Paley sigma+2coords+col, sigma+3coords+col: stay 87.
- BUG in syl16_inv.py: line orbit reps used coordinate perm on column indices (wrong action). affine 4-sums = (y_a+y_b-y_{a+b})/2 lie in LINE spans -> absorbed by lines. rerunning syl16_inv_fix.py with all 35 lines.
- syl16_inv_fix: J_w+35 lines invariant 416; affine add 0. CONFIRMED absorbed. three_b4: 10 blocks pairs suffice, 0 fail.
