# Round 14 notes (rigcompute)

## Plan
- Universe: labels on exact forced set F_K^ex (free moduli -> collision-free vectors, zero demand WLOG:
  replacing a free-modulus label by a collision-free one never increases any m_(v,s); so existence of a good
  system non-genuine at q0 is decided with zero demand on free moduli).
- Stage 1: exact thresholds F_K^ex (C DFS, exists-free mod n).
- Stage 2: one-modulus deformations (genuine elsewhere all levels): rigid = exhaustive DFS finds no
  non-genuine-class level-1 survivor (level-1 relaxation => certificate valid at all levels);
  flexible = explicit witness with all-level lift, exact integer budget check.
  => q_flex1(K) gives an exact UPPER bound on the joint rigid core.
- Stage 3: joint tests (everything else free): reduction to inclusion-minimal signatures at other moduli.
- Exact budget test: (prod q^a)^2 <= Norm(Z_v) * (2 if s odd).
- Genuine class at q (level 1): c = u L mod n, u unit mod n. System-level: global column sign vector eps
  allowed (zeta_j = pi_j or conj pi_j).
- Symmetry at a single modulus with fixed slack table: ell -> u ell, u in (Z/m)^*, u = 1 mod 4 (preserves unit
  index s exactly, parity since u odd). Column-0 orbit reps.
## Validation (done)
- thresholds.py reproduces round10 forced lists (K=4,5; r=4 first/last) exactly.
- deform14.py count mode reproduces referee table exactly: K4: 47 8/0, 83 36/24, 97 16/0, 101 20/0, 103 48/0,
  107 216/144; K5: 107 18/0, 113 24/0, 127 32/0, 151 36/0, 167 120/96, 199 120/40, 211 884/832.
- r=4 first: q*(6)=1039 (115 forced), q*(7)=2089 (173 forced) [scan to 6000/12000].
## r=4 one-modulus (exact) first flexible q (= exact upper bound for joint core), K=4..7
first {5,13,17,29}: 83,163,269,337 (n=21,41,67,84); q* 199,349,1039,2089
rand {53,61,89,173}: 31,59,83,83 ; bal {97..113}: 41,67,89,131 ; last {353..397}: 31,43,43,59
## joint (everything free) r=4 K=4 first: target none (outside all column-sign classes):
- q0=11,19,31: NO admissible label outside the eps-classes at all (local!).
- q0=47,53: LP relaxation feasible (no Lagrangian certificate). MILP running.
## Key structural observation (to prove in write-up)
- Intrinsic constraint at q: (v,s) can never collide in ANY good system if q^2 > Norm(Z_v)*(2 if s odd).
  These involve only monomials with Norm(Z_v) < q^2, i.e. ||v||_1 < 2 log q / log p_min. => LOCe(q) (every
  admissible label is in an eps-class) is K-INDEPENDENT once K >= k(q) := max{||v|| : Norm Z_v < q^2}.
  first4: LOCe true for q<=31, false at 37,41,43,47,... for every K (monomials with Norm<37^2 have norm<=4).
  => purely local rigidity does not grow with K at fixed P.  One-modulus rigidity beyond it uses 'spent' budgets.
- mech first4 K4: one-mod rigid moduli need norm-3 monomials (q<=43) / norm-4 (q>=47) constraints.
## JOINT results first4 K=4 (target none = outside every column-sign class)
- level-1 CP (cp_joint.py sig_first4_K4.pkl): q<=31 no target label (LOCe); 37..61 INFEASIBLE exact; 73 INFEASIBLE.
  67, 71: level-1 solutions exist but lifts failed (level-2/3 collisions at q=3).
- deep tables at 3(a=3),7(a=2),11(a=2) (sigd_first4_K4.pkl): 67, 71, 83 FOUND + verified (jw_*.pkl).
  => joint strong core = [3,61]; one-mod segment ends 79. joint < one-mod.
- MILP (HiGHS) 300s: no solution found at 47/53 (100k cols) -- superseded by CP.
- LP relaxation feasible (no Lagrangian certificate) at 47, 53.
# notes (pal11cert)
- pal11 = pal_patterns(11) (fplib copy), 23 distinct patterns (11 core size 6, 12 twins size 5), L=2, |A|=12093. All 23 slabs active (k(P) in [-2,2]).
- Aut(A) among signed coordinate permutations (base fixed) = Z_5 x {+-1}: multipliers x->a x, a in squares {1,3,4,5,9}.
  NO Z_11 (shifts are not symmetries: autgroup.py exhaustive). Task's "Z_11" belongs to Paley-12 (paley12_patterns).
- widths: test dirs {-1,0,1}^11 + {-2..2}^11 supp<=4: 129713 dirs; tau in {2,4,6,8,10}; tau<8: 1519 (11 tau=2 coords, 134 tau=4, 1374 tau=6). tau(all-ones)=4 => base incidence <=4.
- pool s<=4,e<=1: 1551 (311 Z5-orbit reps).
- isotan.py p1=1152921504606847081: rank_p K_{<=8} = 12093 = |A| (=> rank_Q = |A|; mod-p nonmembership rigorous). per-degree totals 1,12,78,353,1233,3486,8150,11657,12093 (matches fake H). Z5 blocks: T_8 dims 88,87,87,87,87 = 436.
- brute-force enumeration hopeless (m=3: 1e8 leaves). LP relaxation of width counts: max sum x = 7.3333 for pools supp<=2,<=4 (e<=1), supp<=3/4 e<=2 (8151 gens), supp<=5 e<=1; = 8.259 for supp<=6 e<=1 (10791 gens).
## Joint pipeline final: fj.py (tables via csp mode3 / deepsig, ./cpc C search, extend + verify, check_system.py)
- first4 K4: joint core = 61 exactly (mix INFEASIBLE C<=61 using F<=79 suffices; witness at 67 verified, coherent
  w/ conjugate genuine on <=61, NONE at 67, 79..199).  C and Python searches agree (both 401 nodes).
- one-mod first5 K5: seg 293 (flex 307); first5 K4: seg 127 (flex 131).
- first7 K4 thresholds: 92 forced, q*=599, aborted {499,541,547,557,563,571} (node cap 3e7) -> deform14 handles.
- CAUTION: never pkill -f with a pattern in the same command line (kills own shell). Use ../killpy.sh.
- p2=576460752303423761: rank_p K<=8 = 12093, dim T_8 = 436 (blocks 88,87x4). validate.py (p2): lifted kernel vectors satisfy E_<8 c=0 directly (no symmetry), mu_8(c) members, perturbed/random nonmembers, 436 lifted c independent. OK.
- csearch3 (knapsack bounds) on supp<=6 pool: tight cases (one t_core=2 gen + 7 t=6 gens) => huge tree. LP+28444 informative 2-families: s4e1 still 7.333.
- tau-4 dirs = 110 (e_x +- e_y) + 12 twins + 11 core + all-ones.
- hadflats2 (exact, any real z, sigma = support incl base, t = #core cols touched): min t: s2:6, s3:9, s4:6, s5:7, s6:2 (110 supports = (beta_a+-beta_b)/2), s7:5 ...
  => THEOREM A: any certificate has a gen with sigma>=6 and t<=5 (sum t <= 44 at the 11 core slabs).
  Case analysis: |E|<=2 exceptional (sigma>=6); |E|=1 => t<=2 (beta_c, or span(beta_a,beta_b)); |E|=2 => sigma1+sigma2<=14.
- csearch3 3-way partitioned run on bal:6:1 pool in progress (leaves_b61_{0,1,2}).
- LP + all informative 2-families on b61: still 8.259 (families don't cut fractional).
- csearch4 (stronger filter w_q(i)+(r-1)min_w<=cap) : bal:6:1 pool (10791) -> 0 count-feasible 8-multisets (3 parts, nodes per level ~ 1 22 2836 91815 1.0e6 1.8e6 1.8e5 0 0). EXACT (integer counts).
- excflats: exceptional (t<=5) flats: s6 t2 (110, dim1), s6 t4 (1320 dim1), s7 t5 (9240 dim1), s8 t4 (1980 dim2), s8 t5 (2310 dim1 + 13860 dim2), s12 t1..5 (spans of <=5 core cols). No (6,3),(6,5).
- launched bal:12:1 (36894 all +-1 balanced) search.
## Joint cores so far (exact both sides; check_system OK)
- first4 K4: 61 (witness 67).  first4 K5: 109 (witness 113; <=107 depth0, 109 needs depth level1 tables (3:3,7:2)).
- rand4 K4: 7 (vacuous: n<=2), witness 11.  rand4 K5: 7, witness 11.
## mech_joint first4 K4 (C=[3,61], moduli<=79): needed = outer 73; core moduli 3,7,37,41,47,53 (dropping any -> feasible);
   norms: lifting norm-3 or norm-4 budgets -> feasible; norms 1,2 not needed; no single norm suffices alone.
## LOCe ranges (K-stable): first4 31, first5/6/7 11, rand4/bal4/last4 7 (vacuous).
## first4 K6: >=163 certified with Qsub=163 (depth 1 at Q=163); ladder continuing Qsub=223.
- bal4 K4: joint core 7 (vacuous), witness 11 (check OK).  first4 K6: >=193 (Qsub 223), 197 undecided -> Qsub 293.
- first5 K4: joint core 101 (witness 103) -- n 15 (r=4) -> 25 (r=5) at K=4
- THEOREM 4.3 DONE: all 25 profiles excluded (univrun_out_0_17.txt, univrun_out2.txt, univrun_summary.txt). No 8 real linear forms satisfy pal11 width counts.
- ROUND-13 BUG: prodsearch.py caps base incidence at 2 (tau(1) = 4,6,6,8 actually). Corrected: sdc7 m7 45/24; pal7 m7 2134/300; M8 m8 42/0; M8 m7 23355/1746; M8' m10 1092/1092 (zonotopes, zeta(M8')=10); M8' m11 0. All members zonotopes.
- tcnum.py running (evidence).
- first6 K4 joint: table enumeration infeasible (r=6 label space), killed after 40 min at Qsub=107
- loose K6 (rand4, last4): relaxations feasible at Q=11 at depth levels 0,1 (no certificate beyond 7); 100/100 rand4
  solutions fail lifting because of deep 3-adic collisions (q=3 levels 3-5: 215 cases, q=5 level 3: 53, q=11 level 2: 23).
  OOM when building (3,3),(5,3) deep tables with K=6 (644 monomials, 83k sigs). Witness not built: K6 loose = open [7, q_flex1).
  repair.py from K5 witness diverged (300-500 violated pairs).  first6 K4 joint infeasible to tabulate (r=6).
