# verify notes
## theory §2 read
- Thm 2.1 (products): checked line by line; OK. zeta_all projection: projected generator may be non-primitive (e.g. 2e_x); fine since slab box = conv cap Z^d so m[0,v'] can be split into m copies (not stated in proof; minor gap).
- Thm 2.3: OK.
- Cor 2.4: "R-1 <= 2h" needs S=-S (or min sigma >= -h); stated for general S subset cube. Minor; slab boxes are symmetric.
- Prop 2.5: OK as stated (f=1/2 exactly); for Sylvester f = n/(2(n-1)) > 1/2.
- TODO verify: h, (R) for syl8 L=1,2,3; syl16 L=2,3,4; M8' 252 legal pair moves; sigma^6 exact in M8.
## plan (own code: vlib.py)
- primes Q1=2147483647 (2^31-1), Q2=1073741827 (both prime; different from primary's P1,P2)
- M_e(S) via one rref of degree-ordered std-monomial eval matrix; pivots scanned incrementally (own code)
- zeta_pair: DFS, state K = feasible base set, K' = K & K[shift by v]; bound sum_x width_x(K)
- zeta_all: own LP (scipy HiGHS) -> rationalised dual, exact integer check over all primitive v in {-2..2}^d up to sign
- theory checks: sigma^6 exact in M8 via J-membership (phi^e|_S in span E_{<e}|_S) with fmpz/fmpq rank
## results so far
- v0_exact_small: 16 small boxes, exact Q (fmpq) vs Q1,Q2, both signs of every L_pair shift: 0 disagreements (incl. non-spanning dims).
- reach_pair mod Q1,Q2: M8' 11 (H 1 7 28 77 161 266 357 386 329 195 63 7); sdc7 L2 6 (166/188); syl8 L2 6 (286/287). Match round 27.
- t1_sigma: M8 (R) holds, top = 7 Fano lines, sigma^7 in J(A)_7, sigma^6 notin J(A)_6, sigma^6 in J(A_l)_6 all 56 signed shifts (EXACT fmpz). M8' h=4, 28 top pts, 252 moves. sdc7 L1,2,3 / syl8 L1 (R) fails. matches theory.
- t2_syl16: n=16 L=2,3,4: h 3,4,7; (R) yes,no,yes; max sigma on A_{e0-e1} 2,4,6 -> matches theory table row.
- zeta: M8' pair 10 (LP 21/2), all 10 (LP 21/2); sdc7 L2 7/7; syl8 L2 7/7; sdc7 L3 10/10. Primary's M8' L_pair witness verified (|c|=2, 540 pts).
- J9a zeta_pair = 10 exhaustive (873100 nodes; LP_pair = 11).
- J9a: reach_pair 11 mod Q2 (e=10 spanned after 33/36 shifts); zeta_all = 10: LP over 155263 primitive candidates = 11 exactly, dual-reduced set = 36 gens, exhaustive DFS 873085 nodes none of order 11.
- t3_product: 4 small products: H convolves, reach adds (both primes), zeta_pair adds. Thm 2.1 confirmed on samples.
- jobs2 launched (further boxes).
- zeta further boxes (all match primary): rand9:8 pair5/all6; rand10:25 4/5; rand9:20 8/8 (LP9, none 9); rand9:2 11/11; rand10:2 7/7; J9b 10/10 (LP 11, 59 gens, none 11); syl16|3,11,.. 10/10 (LP 89/8); had12|J9 L3 11/11; syl16|J10 L2 7/7; had12|J10 L2 7/7.
- EXACT: exact_shift M8' e=6..10 certified over Q (dimV_Q nullities 1337,980,594,265,70,7,0); e=11 all V_Q(A_l,11)=0 exactly.
- t4: random27 (R) on 3/27 (n6 s105,106; n7 s9), reach=2h each; linkage identity holds syl8 L2, L3, sdc7 L2; M8 minus |sigma|=3 reach 7 (141/169).
- had12|J9 L3 reach_pair 10 both Q1,Q2. syl8 L2 cube SHIFT_6 holds (155 classes).
- EXACT: M8' reach_Q = 11 certified (exact_count, counting argument), 247 s.
- EXACT: M8 SHIFT_e over Q for e<=5; H_Q(6)=287; exact lower bound 286 at e=6 + sigma^6 exact => defect exactly 1 => reach_Q(M8)=6 exact.
- syl8 L2 cube (emin 6): SHIFT_6 holds (287, 155 classes), e=7 141/175 (Q1). 
- J9b reach 10 both primes, 215/281 (matches). had12|J9 L3 reach 10 both.
- 27 round-27 random boxes: my reach (Q1,Q2) = csv reach on all 27, 20 < / 7 = vs csv zeta_all.
- J9a' = rsyl:4:1,2,3,4,5,6,7,8,0:3 (base inside hyperplane): |A| 4791, zeta_pair=zeta_all=10 (LP 43/4). reach queued (jobs4).
- process discipline: briefly had 5 heavy procs (my pkill pattern killed my own shell + released jobs3); killed extras by PID; now chain.sh keeps <=2.
- J9a' (base in hyperplane; its face k_8=0 == M8' box exactly, t5_slices): reach_pair 9 (740/741 at e=9) both primes, H 1 8 36 112 266 504 784 1002 1023 741 286 28, rho 11, zeta 10 -> reach - zeta = -1.
  J9a's face k_8=0 is only M8' cap {|k_base|<=1} (1093 pts). No ternary phi with Thm-2.3 range drop (10 -> <=9) on J9a' (t6).
- had12|J10 L2 reach 7 both (1112/1227) matches.
- J9a' cube (150 largest classes) e=9 740/741 both primes: no rescue. sdc7 L3 reach 10 both. syl16|J10 L2 reach 6 (1250/1251) both. DONE.
