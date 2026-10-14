# referee zetc notes
## Findings so far (reading)
- Thm 2.1 slice proof: indicator prod_Z (1-k_x^2) has degree 2|Z|+2|N|, not |Z|+2|N|. Fix: prod_Z (1-k_x) (j>=0): kills k_x=1 on Z; then sum over Z u P = |P| forces P=+1, Z=0. Claim true, proof typo-level wrong.
- Cor 2.2(3) counterexample model (s=2,L=1,d=6,m=8) is ZONOTOPAL (zeta=rho=11>=8): unfaithful in compute.md sense; README ZE-TC was stated for admissible regime profiles. Headline "ZE-TC^cheap from README 8.2 is false" overstated: false for general T=1 boxes only.
- Quantifiers: rat => deg only when m <= rho (else TC = C^d, M = C^d tangent-cheap). "ZE-TC^deg holds at every order" false for m>rho; Lemma 3.2 deg hypothesis "at every order" unsatisfiable literally; should be "every m <= rho".
- §3.4 Rmk(3): slack sections need |Y| <= w/2 (|k(S)|<=|Y|, slab |k(S)|<=L=w/2), not |Y|<=w. Coordinate certificates stop near w ~ m/4, not 2w. Count used worst-case usage 4||c||_1 (upper bound), so "about w/11.4 blocks" is not an upper bound on #blocks. Labelled a count.
- J4 is unrefereed (README §5 J-items), but derivation from refereed Thm 5.2(b) checks out on reading.
- "ZE-TC^deg is what the reduction needs" vs "rat is what near-strata uses": deg lacks height control / transversal order; proved deg results (tower, products) don't feed near-strata step.
## checks done
- t21_check: Thm 2.1 exact s<=5 all L: 8237 cases 0 mismatch. slices rho(Lambda_j)=s-|j| s<=7 exact; corrected indicator verified.
- Cor2.2(3) exact: |A|=567 rho=11 nu=7 for 4 u's. Model zonotopal (zeta=rho=11).
- fib_ref: L=2 ->21, L=3 -> 22 (confirmed).
- small_had_ref: H8 L1 57 rho3; L2 1037; L3 1877; Paley L1 133 rho3; Syl16 L1 241 rho3; Syl16 L3 |A|=3608377 ; Syl16 L2 |A| = 213937.
- paley_ref: Paley12 L2 |A|=20153, deg<=8 spans (22 components full, p=1048433) => rho<=8 confirmed independently.
- z6_ref: |A|=247, HF matches, twin-free, zeta=6 (zono 6, LP 13/2), dim V=3.
- m8usage: sym-M8 block (|B|=869) per-column usage: max 8 (not 11.4), mean 5.64, dist {2:16,4:70,6:112,8:56}.
- J4 derivation re-read: degrees add to m (2+3+3, 4+4), saturation by collision forms => valid deduction from Thm5.2(b) data.
- Running: twocut_ref (all OK so far), threecut_ref, syl16_ref (L=2 D=15 Singer).
- NEW: m8_cheap.py: M8 (rho 8, H matches), u supported on column S (|S|=4): nu = 6 < 8 (exact), upper bound for all u in R^S by inclusion monotonicity + one-cut (s+L=6). Hull R^S dim 4, widths>=2 => cost>=8: ZE-TC^cheap(M8,8) FALSE (zonotope-free, overdetermined, LP-subcritical). README 8.1 "holds in M8 (J4)" only true in deg form.
- M8xM8, m=16: u=(generic R^{S1}, generic R^{S2}) nu<=12<16, hull dim 8, cost>=16: ZE-TC^cheap FALSE in LP-supercritical zonotope-free overdetermined twin-free model. Degenerate points.
- threecut 843/843 OK; twocut all types OK (incl 3,1,3).
- paley sections: 3 classes exact irreducible+indefinite; all-30 running. Syl16 inc running.
