# referee notes (rigcompute)
- P3-all* (sharp.md Def 1.1 e-all): odd q<=X not | N, levels a<=A_q (q^a<=X), r_j(q) norm p_j => labels parity coset.
  per-unit budget = |sin|<=1 consequence of (1.2). Joint per-v: m0^2+m2^2<=N, m1^2+m3^2<=2N (same x=<v,phi>) is
  a STRONGER necessary condition; Rig(K) as posed uses per-unit only. -> certificates transfer, witnesses don't.
- check_system.py level convention OK (unit from d mod m, deeper levels d%(n q^b)==0).
- Lemma K: "cannot grow with K" only for q<=p_min^((K+1)/2); beyond, needs separate argument.
- forced sets: my rfree.c DFS agrees with thr_*.json for first4 K4,5,6; rand4/bal4/last4 K4,5; first5 K4,5. free vectors of
  K6/7 r=4 and first6/7/8 verified, coverage complete.
- all 10 jwmix witnesses pass my integer-label all-level verifier (per-unit) and match claimed Q, outside counts.
- per-v angle compat (necessary for (1.2)): FAILS for 8/10 jwmix (only last4 K4,K5 pass). genuine passes always.
- FLEX first witnesses: all verified (genuine on F, my rflex2.py); first8 179 verified with my own lift (no lift stored).
- EXACT angle test (box criterion, Fractions) agrees with float. 8/10 jwmix angle-incompatible (exact).
- onemod_angle (PU = per-unit, ANG = +per-v angle): first4 K4: PU counts match round10 (47 8/0, 83 36/24, 107 216/144);
  ANG: 83 RIGID (12/0), 89..103 rigid, first ANG level-1 flex 107 (180/108). => q_flex1 = 83 is NOT an upper bound in the
  P3-all*-relevant class. loose sets K4: ANG flex persists (rand4 31: 16/8, bal4 41: 16/8, last4 31: 32/24).
- rand4 K4 q=41: 16 outside genuine but 0 outside every eps-class.
- MILP (HiGHS, independent): first4 K4 moduli<=79 level-1, C=F∩[3,61]: INFEASIBLE (283s); Q=67: FEASIBLE, exact recheck ok.
- slab BnB (complete, float w/ margin): last4 K4, last4 K5, rand4 K4 joint witnesses: NO angles phi at norm<=K.
  minimal infeasible subsystem last4 K4: 10 monomials. genuine demands: BnB does not prune (ABORT) - sanity.
- one-mod PU counts match referee_L4K F_4 row: 113 336/240, 127 6592/6080, 167 161088/160320 (90,048 in round13 text = F_5-reference variant).
- flexslab (phi=arg pi, 40 digits): FLEX first witnesses slab-OK: first5 K4 131, first6/7 K4 179, first8 179 (my lift), last4 K4 31, last4 K5 43 (w0,w1);
  slab-violated: first4 K4 83 (all 3), first4 K5 163 (both), rand4 K4 31 (all 3), bal4 K4 41.
- angflex: first4 K4 q0=107: 3 verified witnesses: per-unit ok, angle ok, slab ok at arg pi, non-coherent => Q_core^slab(first4,4) in [61,103].
- Lemma M equality condition verified (F∩[3,Q]=all primes) for first4 K4/5/6, first5 K4/5.
- loc files: first5 K4 LOCe at {3,7,11,23} (23 non-initial); rand4 K4 LOCe at {3,5,7,17}: q=17 (n=4) locally rigid, not vacuous.
  first4 K7 loc rows only q<=113.
- deep MILP types at 3:128, 7:512 match author's deep tables.
