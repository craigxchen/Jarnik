# notes round30 search
- S1 checked by hand: z_i - z_j = A(B - Bbar) = 2i A Im B; deg Im B <= w - L - 1. OK.
- S2 checked by hand, and generalised: for any nonzero integer vector a in the difference lattice,
  H = prod rho_k^{a_k+} conj^{a_k-}, Im H deg <= |a|_1 - L - 1, real roots >= |sum a| - 1 => a^+ >= L/2, a^- >= L/2.
- Plan: cube-locus exact elimination; even locus D=8; code enumeration rank4 8pt at D=10,12,14 (cube-type multiplicity
  argument: D=10,14 impossible for the half-cube code); balanced clique bounds; numerical types.
- codes4.py (exhaustive S, MILP, balance): rank<=4 code max at c=0: D=8:8, 9:5, 10:6, 12:8, 14:6, 16:8. rank3: 4 unless 4|D (then 6). D=5 r3:4, D=6 r4:6 (true M 3,5: code bound not tight).
  => cube-type 8-point rank-4 needs 4|D; hand proof via multiplicities m_u=(D+lambda chi)/8.
- clique.py balanced code bound: D=5:4, 6:7, 7:7, 8:14 (Plotkin-type unbalanced 4,8,8,16).
- fftypes.py: D=7, 5-point types (balance h=2, dist>=4) up to symmetry: only 2 types; both rank 4. LM 400 starts: only degenerate.
  any 6-point cluster at D=7 contains 5-pt subcluster => M(7,0)=4 iff these 2 types unrealisable (0-dim expected).
- halfcube12.py + cert_hc12.py: D=12, c=0, half-cube code with multiplicities 2 (even u) / 1 (odd u): 8-point rank-4
  cluster found (300 starts: 4 nondegenerate) and Krawczyk-CERTIFIED (60 digits, radius 1e-30, 2 coords fixed to
  rationals). min|Im| 0.026, sep 0.05 resp 0.157/0.031. Not a composition (M(6,0)=5,M(4,0)=6,M(3,0)=3,M(2,0)=4).
- clique: D=9 balanced max 6.
- cube_locus.py: PROVED by hand+sympy: on the cube locus every extra s in {-1,0,1}^4 (s=a-b per root pair) either is
  killed by m=2 (sum t eps != 0) or m1 forces y_i=+-y_j (=> degenerate root coincidence) or y_i=0. So the cluster through
  Q is exactly the 8 sets at every nondegenerate real point of the cube locus: rank exactly 4 on this 2-dim locus.
- halfcube8_sym.py: 8 nondegenerate realisations of the half-cube type at D=8 (400 trf starts): all t->-t symmetric
  after centering (evidence: half-cube locus = cube locus).
- normalised-form solver (ffsolve, exp(u)) is WEAK: D=9 4-pt types 0/21 realised; base-form realise.py: 20/21 realised.
  => run_types results discarded; redo with realise.py. D=7 4pt: 3/5 realised.
- cert: D=9 4-point cluster Krawczyk certified (type 0): M(9,0)>=4.
- D=12 certified cluster: exactly 8 zero-sum sets at 60 digits, min nonzero residual 0.118.
- 29 feasible S at 8 pts (D=8,12,16): each MILP solution antipodal w/ distances D/2 (half-cube form).
- saturated balance kills D=7 type 0 (pairs zero-sum).
- D=7 sole 5-pt type: 3000 base-form starts, no nondegenerate solution (best degen 3.4e-5) -> evidence M(7,0)=4.
- D=10 c=0: 141 labelled 6-pt rank-4 codes (MILP samples): 0 realised (40 starts each).
- hc12 extension: 918 weight-6 extra sets, 0 hits (2 seeds each).
- circle locus clique: D=6..12: 6,4,6,4,6,5,8.
- c=1 D=10 rank 5 (6 pts): 10/15 random types realised; 2 Krawczyk-certified, cluster through Q exactly 6, rank 5.
- CERT: D=10 c=0 5-point rank 4 (rt_D10_c0_M5_1), cluster through Q exactly 5. M(10,0)>=5.
- CERT: D=12 c=1 rank 5 6-point (rt_D12_c1_M6_1), min|Im| .0258 sep .0085.
- c=1 rank5 random types: D=10 10/15, D=12 7/15, D=14 0/15 (20 starts), D=16 pending.
- D=10 c=0 random 5-pt: 25/30; 4-pt 10/10.
- D=8 rank-5 hunt over all 46 5-pt types: 30 realised; 15 hits: 7 rank-3 extensions = cubes; 8 rank-5 near-hits all
  degenerate limits (sep 9e-5; snapped to 8th roots of unity all 6 points coincide). No genuine rank 5.
- D=9 5-pt: 6 types x 600 starts: none.
