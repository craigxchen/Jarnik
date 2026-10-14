# referee notes
## Reading notes
- Prop 5.2(2) proof: "|u_x-u_y| >= 2/n >= e^{-max eps/2} >= e^{-eps_xy/2}" -- LAST INEQUALITY WRONG DIRECTION (eps_xy<=max => e^{-eps_xy/2} >= e^{-max/2}). G membership needs per-pair gap, worst pair = min eps. Correct (2): min eps_xy >= 2 ln(n/2). Since min eps <= mean <= abar, that needs abar >= 2ln(n/2) (=20.79 at k=16) vs near-extremal abar ~9.64 at k=16 => corrected (2) vacuous at near-extremal W. "window of width 2/(n-1)" claim collapses.
- MIN*=>MIN: e^{-d abar/2}=e^{-W/4}C'^{-d} (writeup says C'^{-2d}); typo, conclusion C<=1 unchanged (C'=max(C,1)).
- Lemma 1.3 proof re-derived: correct.
- Lemma 1.4 identity re-derived: x in 3 or 1 of S_a,S_b,S_{a+b}; correct.
- Lemma 5.0 (nearstrata l.669) needs |eta|<=1/2, |y|=1, |v|<=1: used correctly.
- SLP sl2 proof re-derived: correct. Lemma 3.2 telescoping correct.
- Thm 3.3 boxes are products (decomposable) -> toy scope.
- M8 n=8 not LP-supercritical (smallest n>=32): toy.
## Computations done
- j5_mine.py: all uniform-Q0 numbers match refJ table (o1Q5-d, g*(Q5), R1 exact k<=20). Per-subspace: with sharper mixing bound thresholds even lower (k=12 (1,1): 3.89; k=13,14 (2,0): 2ln15=5.4161). refJ's are conservative. CONFIRMED.
- M8: |A|=1037, H=[1,7,28,77,161,259,287,175,42], dimJ ok; 358 width dirs (7/49/302), tau=h; (J_w)_6=616 and _7=1485 EXACT (fmpz), _8=2940 modp-lower.
- Images: pairs 25 grid rho 8; nonline 125 grid rho 12; Fano 57 pts rho 6, missing 6 pts as stated. J(Y)_e dims Fano [..,6,21,36,45].
- J_w + Fano pullbacks reach J(A) dims in deg 5..8 (modp lower bound, exact Q generators) -> rigorous equality. pairs/nonline add nothing modp (+structural grid argument).
- author's step7 label "None = rho>6" is really rho>=6 (proj_ideal_gens detects rho only if <= maxdeg-1=5); harmless.
- Prop 5.2(2) FALSE as stated: 4 points mutual eps=0 => G empty => MIN^pair vacuous though max eps large. Corrected (2'): min eps >= 2 ln(n/2). Needs abar>=2ln(n/2) i.e. W>=1.9-2.2 Wmin (k=16..20).
- Prop 2.3 verified (rank 45 of 70 conditions, direct TC test). Thm 3.3 verified on 4 new boxes.
