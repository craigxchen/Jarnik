# referee_fake notes

- started; fake.md read in full.
- marcheck: (Mar) small-character margins OK (sdc rho1=1/2, pal >= 1); G_T = Z/4, Z/6 for sdc7/11 all irrelevant; pal G trivial.
- hilb: sdc7 L2 H matches; base 0 vs base 7 same for sdc7. pal7 box DEPENDS ON BASE: base0 861/1877 (L2/L3), base7 807/1781 with rho 7/10.
  => Prop 3.4 (pal7 L3 = M8') only for base = inf. Also actual cluster lies on Z_P(b) for EVERY base b -> joint system.
- tcnum sdc7 m=8 output shows optimizer failures (s=.02: 2e-13, .05: 2.4e-10, .15: 5.2e-7) although zeros with sep .2496 exist;
  the s4.2 table omits exactly these entries. Same AP-shaped minimisers as in overdetermined cases.
- max sep with base at 0 and |u|_inf=1 is 1/ceil(d/2): 1/4 (n=8), 1/6 (n=12), not 2/d.
- sdc7: all 8 base-point boxes give the SAME measure c (support has all |k_x|<=1): base issue void for sdc7.
- sdc11: V(A_0,8) = 5 exact invariant; all 5 measures supported on k_inf = 0 (base coordinate zero) -> forms independent of
  base position; V(core,8) mod p = 5, V(A_0,8) mod p = 5; base-1 box 9553 (job running/killed?).
- sdc7 cert reproduced (sep >= 0.249565, grad, 2nd ff).  sdc11 zeros of author reproduced with own basis.
- NEW: SLSQP max-sep on sdc11 zero set: Krawczyk-certified zeros with sep >= 0.106237, 0.085680, 0.077125, 0.075862
  ((F) holds). Refutes 'none >= 0.02', '<~11% of max', 'cone shrank by order of magnitude'.
- LP zeta check pal11/sdc11 launched (lpzeta_ref.py).
- ZONOTOPES: sdc7 L2, pal7 L2, sdc11 L2 (base inf) contain braid-generator lattice zonotopes of length 7 (verified exactly,
  c in V(A,7), Phi_c = 7! prod <v_i,u>, all v_i pair/coordinate differences) => zeta = 7 = m: m=7 models UNFAITHFUL;
  TC(m=7) inside braid arrangement exactly; kappa >= (5040/(|c|_1 d^7)) sep^7. pal11 LP zeta<=22/3 confirmed independently.
- cheap spans at m=8, n=12: width >= 2(two largest |v_x|) since A contains |k|_1<=2 ternary; all cheap spans in braid arrangement.
