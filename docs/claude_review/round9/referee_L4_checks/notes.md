# Referee notes (working, not the report)

## R1 Regime mismatch (main suspicion)
- Thms A-C unconditional need X >= C_L (8 p_j Norm zeta)^L, D <= X: regime (1.4) p_j << X.
- Prop D slack (Cor D.1) needs p_j >= M^(16h+3) > X: construction regime. There Thms A-C vacuous.
- Cor D.2 "or in the setting of Cor D.1": Thm B needs p_j Norm zeta <= D <= X: impossible if p_j > X -> vacuous; Thm C needs X >= C_L(..p_j..)^L: false. So claim unsupported.
- In regime (1.4), aux zeta must be units at all q<=X, q not | N, so Norm zeta <= X => monomial in genuine primes (doc 6.1 bullet 2). Prop D with sigma-permutation: bound 16 log|Z_v| ~ 16 x (genuine height) -> exceeds per-unit height budget; maybe actual m_{e_j} > sqrt(p_j) (case-3 moduli where u^8=1, Re(zeta^4) factors). CHECK numerically.
- So "inverse theorem must use norm>=3 / deep moduli" is shown only in construction regime, not in regime (1.4) which (T') says is all that matters.

## R2 q | N moduli in regime (1.4): Prop D says "any parity-respecting values" at q|N. Does m_v include q|N? check sharp Def 1.1.

## R3 Prop D in construction regime is "certified" in generalized sense (m_v | Im Z_v^16, zeta of height M^h << p_j); coherent with zeta_j up to F1 on Q_3 u S_0 (~1/2 of primes) + half of b_q=3. So "uncertified in every sense" is relative to Q_7-parity notion only.

## R4 minor: Thm B step 2 constant (30 vs 60 log D) -- harmless (sharper accounting gives 30). Q^(B) fraction 2^(1-B) not 2^(-B). Sub-genuine def misses v_p(Norm) even >=2 (zeta = p_j xi): theorems extend.

## R2 resolved: sharp Def 1.1 (e-all) m_(v,s) runs over q not dividing N only. No issue.
## Reduction (T') via stability Cor 3.6(1): checked, eps = 2/c_0, h_0 = 2c_0-3+o(1). OK.

## R5 Twisted coherence remark in Thm C: Def 4.1 needs p'_j > X so h > A -> Thm C bound (A-1)/(L(h_0+h)) < 1/L vacuous. Only C' (GRH) meaningful.
## R6 Bold sentence §1 / summary "no column certified by zeta with p_j not | Norm zeta" false without norm bound:
   p'_j = p_j mod prod_{q<=X, q=7(8)} q gives exact parity match on Q_7 (norm ~ e^{cX}). Theorems themselves have the bound.
## R7 "genuine height" defined as p_j | Norm zeta: zeta = p_j * pi_sigma(j) has p_j | Norm but labels of pi_sigma(j). Should be v_p(Norm) odd (primitive part). Thm B/C extend verbatim (p_j || m still). minor.
## R8 (I1*) at K=2 "false": needs D far from ALL genuine-height coherent systems; never shown. D(sigma=id) in regime (1.4) agrees with G on Q_3, S_0, all b_q=3 (q=7 mod 16 -> half of Q_7) + loose agreement on b_q>=4 => > 1/2. Compute.
## R9 per-column-per-q unit in Def 3.1/(2.1) is a loose notion (associates mod mu_4: gcd(x,m/4)=gcd(y,m/4)), not the F1 'harmless' common-unit op. Fine for A-C (stronger), but makes (I1*)'s conclusion weak.
## R10 §6.3 heuristic content labelled 'precisely'; independence mod squares should be mod squares and 2 (since (2/q)=1 on Q_7). minor.
## R11 additive-label remark: pigeonhole collides only where Legendre symbols agree; count heuristic. minor, label heuristic.

## Computation plan C1..C7 (see report §8)

## RESULTS so far
C1 Lemma 5.1 exhaustive N=16..1200 ALL OK (N=12 also admits a map -> 47 need not be in S_0; cosmetic).
C2 Prop D reproduced: norm<=2 0 fails; norm3 202-446 fails (r=12, X=400..3000); cases 1-3 norm3 0 fails. Case counts match doc.
C3 regime(1.4)-type (X=2e4, p in [29,400]; X=5e4, p in [29,300]): genuine 0 height violations; D(id) norm2 124/144, 93/100 monomials violate; D(sigma) 134/144, 93/100; norm1 too (10/12, 8/10). max log m/budget 8-10.
C4 per-column (Def 3.1) agreement with GENUINE G on Q_7: D(id) ~0.70, D(sigma) median ~0.56, RANDOM parity-respecting ~0.55 (!). => (I1*) threshold 1/2+delta not meaningful for genuine-height certificates; per-column agreement does not control multi-column collisions.
Other: GRH remark in Thm B: need X >> (log D)^4/delta^2 not (log D)^(2+eps). Twisted coherence 'agree on 1/2+o(1) by Thm B' out of scope (p'_j > X > D).
Cor A.1 can be strengthened: all zeta (large ones trivially) once A >= 2 L h_0 + O(1).
