# notes round15 rigSL

## log
- harness: rigSL_checks/sl.py (budgets PU / genuine-phi / phi point / phi box on top of copied round-14 cpc tables),
  slverify.py (independent checker; own Gaussian arithmetic, own dlog, interval arithmetic for general phi).
- rand4 K=4: genuine-phi (budgets log|D_s|) level-1 on F<=47 already non-coherent at Q=11 (24 nodes);
  full witness built (deep 3:3,5:3,7:2, extend, lift) in 208 s; slverify OK: 40 moduli, 1280 slab checks,
  coherent only on q<=7, 34/40 labels outside every eps-class.  => Q_core^SL(rand4,4) = 7 (exact both sides).
- bal4 K4, last4 K4, last4 K5: GB witnesses at Q=11 (verified slverify). rand4 K5: deep cpc aborted 1e6 nodes; bal4 K5 lifts fail (retrying).
- PINNING: first4 pinned for K>=3 (pin.py: global exclusion + Farkas cone + half-space cube). K=2 locally pinned, global not excluded (box cap). first5 K=4 pinned.
  => SL = GB for first4 (K>=3), first5 (K>=4).
- first4 K=4 GB ladder (F<=199): coherent certified Q=61..79 (cpc INFEASIBLE root); witness at 83 (W_first4_K4_gen_Q83.pkl, verified, coherent <=79).
  => Q_core^SL(first4,4) = 79 exactly (modulo independent re-check). Minimal Qsub for Q=79 certificate: 97 (Qsub 79,83,89 feasible).
- MECHANISM of loose witnesses (diag.py): not uncharged cofactors mainly; RE-ROUTING: the fake drops genuine collisions at
  large forced q | D_s (e.g. D_s=107 prime: genuine pays 107, fake pays 11). Budgets |D_s| for loose primes are large and F-smooth.
- indep_cert.py (own enumeration+propagation): first4 K4 Q=79 with Qsub=97 CERTIFIED (root), Qsub=79 non-coherent level-1 found (both impls agree).
- first4 K=5 GB ladder (GB tables, F<=349 = all F_5): cert Q=103..139 (root), witness at 149 (verified; coherent <=139; 33/65 NONE).
  => Q_core^SL(first4,5) = 139.   (PU: 109)
- ext_test: K=4 witnesses violate norm-5 budgets on few monomials with HUGE m (systematic label choices, e.g. (4,-1,0,0) collides at ~15 moduli) -> artefact of DFS extension order.
- launched: first4 K6 ladder Qsub=400 from 193; first5 K4 ladder; ls.py (local search) test rand4 K4.
- first6 K4 PINNED (34 tight pairs). Unit-D pair counts: first4 16/19/21/21 (K=3..6); rand4 2; bal4 2; last4 0 (K=3..6)
  => loose sets are NOT pinned (need >= r+1 pairs positively spanning); SL phi-region positive-dimensional there.
- first5 K=4 GB ladder (F<=293): cert Q=101..113, witness at 127 (coherent <=113). first5 pinned => Q_core^SL(first5,4)=113 (PU 101).
- ls from random start: no witness rand4 K4 in 900s (weak). drv_repair: cpc+extend then ls repair (lift moves added).
- other agents run forced.py/onemod.py on this machine: CPU contention.
