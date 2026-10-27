# referee_cap notes
- regime q=3251 checks: d=q, ln4r=10.166, 316*ln4r=3212.5<=3251 ok; w=514.2, j=257, m'=2057<q.
- Gauss identity derivation re-done by hand: OK.
- record sdc patterns (sdc_sym.py): P_a={a}u{x:chi(a-x)=1}={a}u(a+NQR); Pp_a={a}u(a+QR). matches author.
- MAJOR candidate: "a cap needs kappa_V <= e^-101333" wrong; cap needs kappa_H<=T*; kappa_V-approximant only sufficient.
  cap.md lines 13-14, 94, 200, 442 ("kappa_V level ... exactly what cap needs"). V_H' may be proper (ln H*=(E/2)lnE >> ln H'~2e7..1e8).
- Prop 4.1, 4.2, 3.4, 3.5, 6.1 proofs re-derived by hand: OK. Lemma 3.0 OK.
- toy claim "circles largest for m<q" contradicted by own table sdc7 L=2 m=6 (0.196 vs 0.205-0.264 perm, up to 0.312 unif).
- Lemma 3.1 "no traction" is interpretive (l1 coefficient mass ignores cancellation).
- Remark 6.2(1) sub-box Hoeffding: only a sufficient-condition failure.
- "height not the obstacle" (ln H'~1e8) heuristic.
TODO scripts: regime_ref, gauss_ref, toy_ref (sdc7/sdc11 L=2), margin MILP, (exact N_i rational sdc7 m=6).
- regime_ref: all numbers reproduce (lnT* floor -101332.92 at W=r tau; author uses W=r tau+1/4).
- gauss_ref: identities exact q<=139 OK.
- toy_ref: dims/ratios reproduce; circles NOT largest at sdc7 m=5 (perm 0.405>0.4023) and m=6 (perm .2467, unif .2138 > .1955).
- exact_Ni: sdc7 L=2 m=6: N_i rational, symmetric, Psi constant = sum binom^2 N_i/4^m to 1e-28. Prop 4.1/4.2 confirmed.
- margin_milp: min sum_P|v(P)| = q+1 exactly (q<=43, B=3): Prop 6.1 input confirmed & sharp.

## pass 2 (fresh referee instance, 2026-10-27)
- prior report exists (backup referee_cap_checks/referee_cap_pass1_backup.md). Re-verify independently with new r2_*.py scripts (arb intervals).
- README 7.1 confirms M1: test is on kappa_H; "for the cap direction H is irrelevant since kappa_H <= kappa_V" => kappa_V small sufficient only.
- barrier Lemma 3.4 tau_D third term 2 ln(1000 q n^2/C) c0-free: floor OK. Thm 3.5 also needs tau >= Thm 2.3 thresholds: only lowers T*.
- fake.md: sigma_lambda = sum_S mu_S(|lambda(S)|-1/2) (mu=1 for sdc) -> Prop 6.1 margin form OK.
