# notes (compute lane, round 28)
## Plan / conventions (decided before any new data)
- Box convention = record (round 11/12, killer klib.build_box): k' in {-1,0,1}^d, kfull = k' + (-sum k'),
  slab |kfull(S_p)| <= L. (So |<k,beta_p>| = 2|k(S_p)|; "L" = bound on |k(S_p)|, = j-analogue.)
- rho = top degree of H_A (regularity); reach <= rho always (Thm 3.1 at m = rho+1 would contradict V(A,rho+1)=0).
- Least spaces via ONE rref of the std-monomial evaluation matrix (exponents <= 2; k^3 = k on cube), columns in
  degree order: rows with pivot in degree block e, restricted to block e, span M_e(S). Row-echelon suffices.
- M_e(A_{-l}) = M_e(A_l) (A_{-l} = A_l + l, translation invariance of leading forms) -> shift classes l ~ -l.
- zeta: (a) zeta_all = generators primitive in {-2..2}^d (record/killer class); (b) zeta_pair = generators in L_pair.
  Upper bound: exact rational LP dual (record lpzeta) / MILP; witness verified exactly (support in A, moments).
- truth: only bounds  max(reach, zeta) <= truth <= rho  unless the record decides it.
- primes: P1 = 2147483629, P2 = 2147483587 (both > 2^30).
## Validation log
- v0_validate: slib (mod P1, P2; one sign per shift class) == killer klib exact-Q (both signs) on cube d=3, hexagon,
  10 random boxes d=3,4, every degree, including the non-full dims (10 vs 11, 3 vs 4). 0 disagreements.
- t1_zeta_small: |A| = 757, 1037, 1877, 1765 (record). zeta_all = zeta_pair = 7 (sdc7 L2, syl8 L2), 10 (M8', sdc7 L3);
  LP bounds 7, 7, 21/2, 21/2 -> all exact. (M8' has an order-10 zonotope with ALL generators in L_pair.)
- repro (both primes): sdc7 L2 reach_pair 6 (166/188), cube 6; syl8 L2 pair 6 (286/287), cube 7 (141/175 at e=7).
## Box sizes (survey.py)
- d=8 structured: rsyl9 L2 1837, rsdc9 L2 1367, rhad9 L2 2049; L3 ~4500-5000; L1 = 73 (!).
- d=9 structured: L2 2709-4525; L3 11-13k (infeasible); L1 91.
- random n=9: L1 121-633, L2 1907-3803. random n=10: L1 177-1255, L2 4905-7909 (L2 infeasible: H_A rref too big).
## Observations so far
- 27 killer random boxes (d=5,6): ALL have zeta_all = rho exactly (unfaithful, zonotopal) -> reach > zeta impossible
  there (reach <= rho). So these boxes cannot test the decision rule; "reach - zeta" there is "reach - rho".
  Must separate faithful (zeta < rho) from unfaithful boxes in the trend tables.
- zeta_pair < zeta_all on n=7 seeds 7, 8 (5 vs 6).
- d=8 structured L=1 boxes: all the same 73-point box (H 1 8 36 28), rho = 3: uninformative.
- d=8 cube shifts: 3087 non-pair classes, ~20 min/box/prime -> random d=8 uses the 800 largest classes (lower bound).
- zeta MILP (full candidate set) OOM'd (bad_alloc) on rsdc9 L2 -> replaced by exact DFS with LP-dual reduction
  (only v with sum x_v (c_v - 1) <= B - t can occur; c_v = sum y|lam(v)| >= 1 exact). Fast; exact.
- d=8 results: rsyl9 L2 (both J): rho 8, pair 6, cube 7, zeta 7 (=M8 pattern). rsdc9 L2: rho 7, pair 6, cube 6,
  zeta 7 = rho. rhad9 L2: rho 8, reach 7 (pair=cube), zeta 7.
- **rsyl9 L3 (d=8, |A|=4861): H = 1 8 36 112 266 504 784 1001 1029 756 315 49, rho 11; reach_pair = 11 (both primes);
  LP bound = 11 exactly but NO order-11 zonotope (DFS exhaustive 3.5M nodes; independent MILP on the reduced set
  (= the 36 L_pair classes, c_v = 1) infeasible); zeta = 10 exact.  => second box with reach > zeta (d=8 M8').**
- Truth lemma (for write-up): if H_A(rho) = 1 and zeta_all < rho, then truth <= rho - 1. Proof: V(A,rho) = <c>;
  if all zeros of P_c were degenerate, Thm 2.2 classification gives P_c = cY^a prod (Y^{v_j} - zeta_j)^{m_j}, order
  rho = sum_{zeta_j=1} m_j; Newt(P_c) = a + sum m_j[0,v_j] subset conv(A) and contains the lattice zonotope
  a + sum_{zeta_j=1} m_j[0,v_j] (vertices are lattice pts of conv A, A = conv A cap Z^d) -> zonotope of order rho
  in A, contradiction. So Bs V(A,rho) has nondegenerate points. (v_j primitive in {-2..2}^d as in Thm 2.2 (2).)
- rand:9:1 (d=8, L=1): rho 6, zeta 5 (< rho, faithful), reach 4/4.
- pkill -f chain.sh killed my own shell once (pattern matched own command line): use PIDs.
- d=8 L=3 extra restrictions: J9b and 3 random Sylvester J: reach 10 = zeta 10 < rho 11 (tie); J9b cube (300) no change.
  had12 random J: reach 10 < zeta 11 = rho; sdc11 random J: reach 10 = rho = zeta.
- d=9 structured L2: syl J10 pair 6 cube 7 zeta 7 rho 8; sdc J10 6/6, zeta 7 = rho; had J10 7, zeta 7, rho 8.
- d=9 random L=1 (20): 19 <, 1 = ; all zonotopal.
- Decision reading: UNFAVOURABLE (only M8' and its hyperplane+point lift beat zeta; gap grows with rho/L and d).
