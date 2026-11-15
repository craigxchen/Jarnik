# kill_search notes
- Lemma 1.1 re-derived: correct.
- Prop 2.1 re-derived: correct (half-cube type only).
- Thm 3.1 re-derived: correct; typo in third degeneracy case (y_j=-y_i, eps equal gives r_j=-r_i, not -r_j=conj r_i); still degenerate.
- Thm 3.1 scope: cluster through Q only; other clusters on same N not excluded -> check M(N) by full clique.
- 3.2 vertical line 'r<=3' is a dimension count (kill_ff heuristic), not proved; own formula r*ceil(L/2)<=D-1 gives r<=4 at D=2 mod 4.
- Claim "8-pt rank<=4 need 4|D for D<=16 by computation": table has no rank-4 entry at D=11,13,15.
- kcert.py (own Krawczyk, mpmath.iv, radius 1e-25, own gauge: fix best-pivot coords to 8-digit rationals):
  PASS for D9c0 M4 (rank3), D10c0 M5 (rank4), D10c1 M6 rank5 (x2), D12c1 M6 rank5, D12c0 hc rows 0-3 of hc12_lam4 + hc12_lam-4 (all 8 pts rank 4).
  Full max clique over ALL 2^D points of N equals the cluster through Q in every case (4,5,6,6,6,8x5).
  Min nonmember residuals >= 0.015; separations >= 0.0074; all fine.
- rankcodes.py (own DFS over functionals Lambda_r, |Lambda_4|=575, MILP HiGHS): D=8 r3=6, r4=8 (matches lane).
- cubeloc.py: 60 random nondegenerate cube-locus points: cluster through Q = 8 AND max cluster over all 256 points of N = 8.
- balclique.py (networkx, S_D neighbourhood split): 4,7,7,14,6 for D=5..9 c=0; 14 for D=7 c=1. Matches lane.
- circ.py: circle locus clique 6,4,6,4,6 for D=6..10. Matches.
- rankcodes.txt: r<=3: 4 except 6 at D=8,12,16; r<=4: 6@6, 5@7, 8@8, 5@9, 6@10, 5@11, 8@12, 5@13, 6@14, 5@15 (D=16 r4 pending). Fills lane gap at D=11,13,15.
- vert.py + vcert.py: KRAWCZYK-certified rank-4 5-point clusters on the VERTICAL LINE q=i*y at D=10, c=0 (two solutions;
  cluster through Q 5, max on N 5, |y| distinct). REFUTES lane §3.2 'vertical lines give rank <= 3 at c=0'.
- ktypes.py: D=7 M=5: 2 types, one killed by weight-2 span vector; survivor == T_7 (canonical forms equal). Matches.
- rankcodes D=16 r4 = 8; D=10,14 r4 without balance = 6. All match lane.
- MAJOR: six10.py enumerates ALL 6-point rank-4 codes at D=10 (distance>=5; balance-feasible): exactly ONE class up to
  XOR-translation + coordinate/row perms (lane's '141 labelled codes' = this one class, with orientations). Its real span
  contains e5-e4 = ((R1-R2)-(R3-R4))/2 (weight 2 < L+1 = 5) -> impossible by lane's own Lemma 1.1. So NO 6-point rank<=4
  cluster at D=10 c=0: PROVED (exact enumeration + Lemma 1.1). The 'lead' instance #1 is explained by linear algebra.
- D=14: 2 classes of 6-point rank-4 codes; both have span vectors of weight 2 resp 4 < 7 -> killed. No 6-pt rank<=4 at D=14.
- rank<=3 6-point codes impossible at D=10,14 (rank3 max 4), so at D=10,14 every cluster with >=6 points has rank >=5.
- realise.py (own, free orientation + slices, 300 starts): D=8 46 lane types: 35 realised w/ orientation, 40 any. Lane 30.
  Types 9, 12 (lane log: not realised) Krawczyk-certified via real_one.py + kcert.py.
- D=12 cluster column multiset verified half-cube (2,1).
- Plotkin checked by WebSearch (even d, 2d>n).
- ktypes D=9 and D=8 enumerations stopped (process limit / time).
