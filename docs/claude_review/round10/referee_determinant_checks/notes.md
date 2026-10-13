# Referee notes (working file)

## Re-derivations done (by hand, in head)
- Thm 2.1(a) Hoeffding: ranges 2T|b_x|<=4T, exponent = ln(4r) exactly. OK.
- 2.1(b) closed form m_0 and 0.7 n(n-1) for n>=55. OK (to check numerically).
- 2.1(c) counting Siegel: OK; needs L>=3.3; L>300 holds.
- 2.1(d): 4 sqrt(8/300)=0.653 <0.7. OK.
- Thm 2.3: OK (H H^T = nI, E|S|<=sqrt n, h_A(1)=0 on slice).
- Thm 4.1: G_k formula, valuation, norm bound, Taylor bound, comparison all re-derived. OK.
  SLIP: proof says "Let c in V_H be an integer vector" and uses |c|_1<=UH; only the generators
  (integer kernel vectors with |c|_inf<=H) satisfy this. Conclusion OK by linearity (the proof says so).
- Prop 4.2 OK. Cor 4.3 OK; fakes (round5 Thm2, exact Thm5.1) draw delta uniformly on [0,Delta]^M,
  full rank (twins / S'_a - S_a = e_a - e_0), good set measure > 1/2. Z_P translation-invariant along 1.
- Lemma 3.1 OK. Lemma 3.2 OK but needs (pi/4)Z gap (Lemma 0.1 second form), not (pi/2)Z.
- Prop 5.3(iii): constant wrong: correct is X_v >= W/2 - 2 log(C|v|_1/(2 sqrt2)) (factor 2 on log, and pi/4 gap).
- Prop 5.2 application to actual clusters drops h_A(kappa) (winding term); kappa mod 4 not controlled
  by the slabs, so "k <~ 0.25 sqrt(n/ln 4r)" is unproved as stated.
- Prop 5.5 OK. Lemma 5.4 OK.
- Thm 5.7: logic OK given ZE.

## Points of substance
1. Fake separation is generic: any nontrivial analytic identity in delta is violated by a.e. fake.
   Separation is "a.e. fake", not "every fake"; whether some fake lies on Z_P is equivalent (for
   nondegenerate delta) to a special case of ZE. The real separating input is integrality of
   multi-term forms F_c(z) (outside.md s6 item 1, already listed as outside the barrier).
2. ZE evidence uses V(A,m) (more equations => smaller zero set) not V_H; biased toward ZE.
   Full V(A,m) is NOT reachable by Liouville (successive minima too large), so the relevant object is V_H.
3. Regime is n >= ~2137..4278 (exact) / 300 ln 4r sufficient: astronomically large; Thm 4.1 may be vacuous.
4. Growth-table correction: fair; round8 s2.5 item 5 already says no proved cap covers aux polys.

## To compute
- params: exact m (big-int binomials) for selected n; thresholds; 0.7 n(n-1); Thm 4.1 constant.
- Hoeffding fraction on a skew-double-core profile (evidence).
- small-scale sanity of Liouville+divisibility bounds on actual clusters (exact Z[i]).
- fake test small-scale: g_c(delta) != 0 for fake deltas (trivial) + V_H-proxy ZE test.

## Results of computations
- ref_params: exact m; thresholds 2137/2679/4278 confirmed (fail at 2136/2678/4277). Thm 4.1 const: RHS-16 ln n <= 12.79 (n>=64).
- ref_chain: 240/240 kernel vectors on actual clusters (incl. e_p=2): phase identity, Taylor upper bound, forced-divisibility
  lower bound all hold; all F_c nonzero (no forcing at small n, as expected: LP ratio < 1 there).
- ref_VH_small: at author's evidence sizes V_H = V(A,m) (echelon integer basis sup-norm 2835/87584/1536 << H). So
  author's caveat "V(A,m) not V_H" is moot there.
- ref_fakes: skew double core q=7 (r=14) and round-5 Paley-I b=8 q=7 (r=64): rank 7 = full. Slab-box T=1, m=5:
  200/200 random fake deltas certified off Z_P (leading Taylor term > remainder).
- ref_fullbox: Prop R1 (no-cut box => base locus = sparse strata) consistent numerically.

## Retracted concern
- Prop 5.2 kappa term: for rank k <~ sqrt n, Cramer basis of (colspace)^perp has |v|_1 << 1/Theta, so <v,delta>=0 and
  delta in span{a_p}+R1 with kappa=0. Write-up's sentence is right but the step is unstated.
