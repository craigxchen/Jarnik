# Referee working notes (round 11 theory.md)

## Verified by hand (OK)
- Lemma 1.1: d L >= U(L - ln2) - E_S(L + ln UB). OK. constants 27.73, 41.2 (1.0261 m ln n at n=64), 13.73. OK.
- Lemma 1.2: margin mW(1/2-1/20-3/8)=3mW/40; (80/3)(0.1085+1.3568)=39.09<40. OK.
- Thm 3.1(i): cost <= wW + m log p <= mW/4 + mW/8; ln H_0 <= 2(n-1)ln(2n+1)+2(m-1)ln n <= mW/40 iff W >= 8.68+80 ln n. OK.
- Thm 3.1(ii): R10 Prop 5.5(i) unique minimiser => v_pi(F_c)=eD/2+mu finite. OK (rederived G_k factorisation).
- Prop 2.1-2.4, 2.5 (MVT, m < (n-1)(2n+1)), Prop 3.4 (saturation duality), 3.5, 4.1: OK.
- Cor 3.2 (<= 7 heavy primes) OK. Thm 3.3 proof logic OK *given* SL(n').

## FATAL 1: SL(n) is false (barrier profiles)
- Duplicate column q of p (beta_q = beta_p, e equal): A_p keeps |<k,beta_q>| <= w so h_{A_p}(a_p) <= e_p w < m: SL(P,p) fails.
- Twin a_q = a_p +- e_x +- e_y (x,y != base): l_p = l_q +- k_x +- k_y, width <= w + 4T < 4w <= m. SL fails.
- round5 Paley fake P(H_q,pi,8,a*): 8 identical unflipped copies per label + flipped copy differs in one row.
  => every column except possibly the flipped copy at the base row fails SL. <= 1 column. Admissible:
  rank M-1 (twins), pair separation #differ >= (r+4)/2 (Thm 1 G>=4 at pairs), n-1>=300 ln 4r for M>=3528.
- skew double core: S'_a = S_a xor {0,a}: twins differing in 2 rows; base point != 0 => <= 2 columns escape.
  pair separation exactly M of 2(M-1). Regime M >= 3068 (skew-Hadamard orders, Paley q+1).
- Also abstract: Hadamard profile with each column doubled.
=> SL(n') false at every such n' (infinitely many Paley orders), so hypothesis of Thm 3.3 never holds: vacuous.

## MAJOR 2: parity of m
- |<k,beta_p>| <= m/e_p  <=> |l_p| <= m/2 => at most 2 floor(m/2)+1 levels. m odd => m levels => level-sum
  lemma => no UM_p for ANY profile (author's own Prop 3.5(i)).
- exact: m(n) odd for 4076 of 8137 n in [64,8200]; m odd at n=2137, 3068, 3528, 4096, 8192.
- Bottom-level lemma (m even): UM_p forces l_p(k*) = -m/2 (pushforward = kappa * alternating binomial;
  level j0 has v=0 so -v(binom(m,j*)) > j* impossible unless j*=0).
- So "level obstruction disappears" is false (odd m) / only at bottom level (even m).
- Repair: even m_e = 2 floor(m/2) (all inequalities hold, checked), or widen to lambda*m/2 with lambda log p <= W/8.

## MAJOR 3: evidence boxes are not A_p
- had_sl.py: |<k,beta_p>| <= 2m (|l_p| <= m, 2m+1 levels), cost 2m log p: needs log p <= W/16 variant.
- bool models: l = |t| with M+1 >= 9 > m+1 levels.
- So none of the evidence tests the theorem's A_p (m+1 levels). Need faithful test (ref_had_faithful.py).

## Repair proposal (grouped widening)
- Q set of columns with sum_{q in Q} lambda_q log q <= W/8, widen each q in Q to |l_q| <= lambda_q m/2.
  cost <= wW + mW/8 <= 3mW/8. Thm 3.1 generalises verbatim. Hypothesis SL_Q(n): every admissible profile has
  some (Q, p in Q) with UM_p in V_0(A_Q). Not refuted by twins. Open; must be tested on barrier profiles first.
