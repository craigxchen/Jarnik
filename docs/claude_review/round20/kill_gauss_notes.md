# notes
## first pass (after reading gauss.md)
- Lemma 1.2, 3.2 algebra re-derived by hand: correct (delta <= W0/(2(W-W0)); Markov reverse; convexity).
- Lemma 3.3: Schwartz-Zippel over C with u-coordinates; uniform in circle. OK.
- Heuristic count for ANY profile at k points: e^{W'(1-(k-1)(1/4+delta/2))}; k=5 -> e^{-2 delta W'} profile-independent.
- phases of q_T don't change class (c_x^2=1 forced by leading coefficient) -> dim count 2-2d stands for single-variable families.
- KILL ANGLES: (A) is A_5 / L5circ at the Liouville exponent -> capped by hardness.md / item 3? (B) check Thm 4.1 numerics; (C) unbalanced profiles; (D) multi-param families k=5 at 7/15; (E) Thm 3.4' Step 1 H_Sigma bound.
## second pass (after overview + gauss.md theory re-derivation)
- Thm 3.4 crude: re-derived fully, correct. Thm 3.4': Markov needs factor 4 (k(k-1)W/(2(M-1)) = 4*E) OK.
- Lemma 1.1 over Q: at split p (embedding via pi), lambda_{E_A+E_Ac,p} = #layers with A_{p,t} in {A,Ac} * log p. OK.
- Thm 4.2: kappa = (k-1)(1+2delta)/4 - 1; k=5 -> kappa = 2 delta. V_E(5) predicts A_5(delta) for EVERY delta>0. k=4: kappa>0 iff delta>1/6 (= CC) consistent w/ Thm 4.1.
- K.gamma-hat list in table starts at k=4 (-3/8 at k=4), consistent with remark (first >0 at k=6). Hand check k=4 = -3/8 OK.
- Thm 4.3(i) overclaim: single genus-0 curve with O of degree>=3 is NOT excluded by (ii) (ii is covering families only). minor.
- Liouville off closed set gives 1/tau; RV gives 1/beta >= 1/tau. So 'never better than Liouville' OK.
- Count dims one-variable polynomial families, general profile: dim <= D(2-(k-1)(1-e)) + k-3; dense needs D(2-(k-1)(1-e))>=1 -> e > 1-2/(k-1) = 1/2 at k=5. So no expected-dim dense polynomial family at k=5 for any subcritical e, any profile.
- TODO: check construct.md Thm B scope (is mean width theorem stated for T_k Gaussian model?), L8.md 6.1 gamma-hat; recompute Thm 4.1 instance + Jacobian rank; kclass; gcdred numbers.
## computations done (kill_gauss_checks/)
- phi_check.py: 4th-root formula makes Im e1=Im e2=0; Im e3 = Phi/(v1+v2+v3). CONFIRMED.
- k4_instance.py: roots reproduce; Im P = -420,-180,-240 const; Jacobian rank 3 (exact, 6 rational pts incl. non-instance (X,m)). CONFIRMED dominance.
- k4_lattice2.py: gcd norm <= 10368000; eta ~ 0.033,0.013,0.0065 at s=1e6,1e15,1e30 (->0 like 1/log s); BUT ratio diam/R''^{3/7} ranges [8.31,142.3] for large s (gauss says [8.3,43.4] -- sample dependent, harmless); max |t| over all pairs incl. anchor 240.
- RESIDUE INTEGRALITY FAILS: 408/1866 pairs have |c_ab|^2 = 2 mod 4 (mixed unit classes), so |t_ab| = |c|/2 not integer (e.g. 1/sqrt2). (L_k^circ) definition's 't_ab is a nonzero integer' false off one unit class. O(1) effect only.
- kclass_own.py: all K.gamma, (K+D).gamma, count exponents, kappa reproduce.
- gcdred_own.py: Lemma 3.2 bound no violations k<=8,M<400; rho_{10,5}=0.00794, E[W0/W]=0.02345; R''>=sqrt R; exponents .606-.702. CONFIRMED.
## key kill idea: invariance lemma
- For bihomogeneous F of bidegree (a,b) in (u,v): F(g zeta) = g^a gbar^b F(zeta); tropical forced divisor d(z) = g^a gbar^b d(zeta) at every prime; archimedean angle Delta invariant. So every forced-divisibility Liouville inequality for z (critical arc) is IDENTICAL to that for zeta (subcritical arc). The 'subcritical' relabelling gains nothing for aux polys / single LW (checked: single-form LW at balanced k=5 gives |Lambda| >= e^{-4(1+eps)w} vs precision e^{-4w}: exactly critical).
- Dirichlet count for approximations to e in T_5 wrt h_B: X^{1-4theta}, theta=1/4+delta/2 -> X^{-2 delta}: A_5 sits strictly between Dirichlet (delta=0) and Liouville (1/30 balanced, 1/10 general). Rational target -> McKinnon-Roth regime; only Liouville known unconditionally.
