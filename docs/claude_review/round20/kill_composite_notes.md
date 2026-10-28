# notes for kill_composite
## pass 1 (proof re-derivation, by hand)
- Lemma 1.2, 1.3, 1.4, Prop 1.5, Thm M (1),(2), P1, P2, P2', P3, P4, P5, Prop 3.1: re-derived; no error found.
  (P5 twisted: d=1 term identical both sides; diff bound 2 C_E X^g tau^w (1+tau^h) ok.)
  (Prop 3.1: Mellin ok; at most 2(s-1) cancellations ok; L(it,lambda_k)!=0 via FE ok; f = f0 u^{-it_m} ok.)
- F1 (minor): Cor 2.2 says A > max(2,1/g), A >= 4 for lw.md fakes. Round-10 referee Slip 1: lw.md needs A > 4
  (P1 >= (40M/C0)^4). Repeated slip. Correct: lw.md A > max(4,1/g); P3-char A > max(2,1/g).
- F2 (major, classification): §3.3 transfer of referee's max-scale equivalence is valid only for the MINIMAL graft.
  With multiplicativity, primitive support contains mixed d*beta, and the sharp LS needs, for all beta coprime S with
  d beta primitive and all primitive alpha' coprime S:  |sin(psi_d + arg beta - arg alpha')| >= c / max(Nd Nbeta, Nalpha').
  Taking beta ~ sqrt(Ngamma/Nd) for gamma = alpha' conj(beta) gives a geometric-mean gap ~ c/sqrt(Nd Ngamma) on gammas
  with balanced divisors -> LCG-type, heuristically collapsing (convergents must be almost-prime). So with mult. it is
  NOT the max-scale cross-gap. To check numerically.
- F3 (major, scope): position-dependent composite inputs: smooth-radius, sharp-angle sector counts with BA edges have
  genuine error O(log X)?? (Ostrowski). Non-robust, not uniform in position; Thm M needs fixed gamma, Thm P needs
  position-uniform. Not in proposal's list. Check numerically genuine error and minimal-graft violation.
- F4 (major, overclaim): "only non-robust leftover is sharp LS on primitive support" false: smooth sums WITHOUT
  multiplicativity (proposal's open pt 3, called "irrelevant") are unclassified: any finite-support compensated extension
  is invisible to asymptotic O_f statements; effective ones need compensation analysis. Prop 3.1 needs complete additivity.
## computations
- c1: Prop 1.5 lb and Lemma 1.4 ub hold X=1e2..1e5, K=0..100 (Gram eigenvalues). ratio Lambda*/(X+K sqrtX) in [0.70,1.83].
- c2 (F2): p=101, tan psi = sqrt7-2. X*minspacing primitive: genuine ~1.00; minimal 0.23-0.83; multiplicative 1.2e-1 -> 5.6e-5 (X=1e6), 2.3e-4 (4e6).
- c3 (F3): smoothed sector counts, BA edges: |err| <= 0.13 for X up to 1e7; rational edges ~0.0225 sqrt X (71 at 1e7).
- c3b: 30 BA sectors, s=40 fake primes in [1009,1601], X=1.5e6: genuine |err| max .448 mean .141; minimal graft max 8.06 mean 2.30.
