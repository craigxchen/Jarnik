# Referee notes (L5)

## Read so far: L5.md sections 0-3.1, 3.2 up to Lemma 3.4
- Prop 2.1: identity re-derived OK.
- Prop 2.3: OK algebraically (kappa_x = kappa_y => Lambda_y = kappa_x - F_y).
- Cor 2.4: claim ||phi_x-phi_y||_1 >= (||a_x-a_y||_1-4)/2 -- check. Also "with column weights within bounded factor ... >= c r" : HIDDEN ASSUMPTION; margins d_xy >= W/2 means w(a_x-a_y) >= W/2, i.e. weighted; converts to >= c r only with comparable weights. Not valid generally for actual clusters? check stability.md.
- Lemma 3.1: check counting: parity class n_q elements; M-1 used + 4|V_y| forbidden. But pair collision only within same parity class: fine.
- Thm 3.2: bound re-derived OK. (a) needs (T) up to k_max+2c0: stated.
- KEY: IFL (a) is satisfied by GENUINE labels of any distinct split primes p_j <= M^{O(1)} (height bound m_{v,s} <= 2 e^{w(v)/2}); so (a) is not hypothetical; the content is Hall (b). -> author calls IFL "hypothetical object of L4": overstated/misdescribed? check.
- Lemma 3.4: math OK (re-derived). But it's a MODEL (independent uniform candidates); actual candidate sets are translates kappa^lam_y - 2<w_i,f> (structured). Also per-modulus failure 1.5/n_q; summed over q in Q_large: ~ 3 log A expected failures, not o(1). Check how labelled.
- Cor 3.3(2): quantifier: "L5 iff L5_med on such profiles" needs IFL for infinitely many M with uniform constants. loose.

## Sec 3.2 continued
- STATUS para says IFL(a) "is exactly the hypothetical object of L4 (ideal uncertified)". FALSE as description:
  genuine labels zeta_j = pi_j (any distinct split primes p_j <= M^{h0}) are parity-respecting and satisfy (F),(T)
  (even for all sparse v) with beta0 = O(h0) -- the author's own Remark 1.3. So (a) holds unconditionally in regime (1.4);
  the only open input of Thm 3.2 is Hall (b) for a fixed S_y. -> Thm 3.2 is closer to unconditional than stated
  (minor/positive correction). Check genuine Hall numerically? (evidence only)
- "Q_med = where (P) forces collisions": pigeonhole forces only where a parity class has > n_q rows, i.e. n_q < M/2
  (balanced) ... n_q < M. Range M <= n_q < cM (c>1) is in Q_med but no collision forced. minor overstatement.
- L5+ called "strengthened lemma"/"stronger": (T+) => (T), so L5 => L5+; L5+ is WEAKER (stronger hypothesis). wording/logic slip.
  Recommendation (use L5+) still right because L5+ suffices & is implied by L5.
- Cor 3.6 gap: shifted monomial (w - 2 sum c'_y w_i^(y), c') may lie in L when sum c' = 0 (c' twin part of a
  character); ideal bound then not available unless "ideal for all sparse monomials" includes L. For Thm 3.2 no issue
  (c' = +-e_y never in L).
- Remark 3.5 common switching "only at q >= M^{2+eps}": model statement (birthday), labelled? fine as model.

## Sections 4-6 (read)
- Remark 1 of §1 ("unweighted reading false for EVERY system, q=3,5"): WRONG ARGUMENT. q=3,5 give 2 moduli = O(1),
  absorbed by O(||w||+1). Pigeonhole gives SOME difference colliding at each small q (w(q) varies), not one fixed w at
  many q. Correct reason the weighted form is needed: genuine labels can collide at unboundedly many SMALL moduli
  (q | Im pi with Im pi = 3*5*7*...); demo by exact computation. -> unweighted form not valid for actual circles.
- Thm 4.1: re-derived OK (Cauchy-Binet, interlacing, Hadamard). Cor 4.2(1) proof re-derived (E entries <=24/4=6).
  Need: check sharp.md §4 profile matches "columns of H' repeated b times + one flip per row in distinct f(z)".
- Prop 5.3(2) GAP: no collision needs kappa = A nut + 2 kappa' injective (within parity classes), not kappa' injective.
  kappa'=d_M u with u injective can collide when (A nut)_x - (A nut)_y = 2 d_M(u_y-u_x) mod m. Fix: greedy in coset
  e_x + d_M Z/n (size n/g), M-1 forbidden -> still works if n/g >= M. Test by brute force.
- Prop 5.3 / Thm 5.1 assume rank M; rank-deficient actual profiles (relations a1+a2=a3+a4) not covered (Delta_M=0).
- §4.2 "if free labels additive, relation <w,mu> = -eps lambda_y violates (T)": existence of such sparse w not shown.
- Cor 2.4 "column weights within bounded factor" unnecessary: W>=3MlogM-O(M), p_max<=M^h0 give ||a_x-a_y||>=Omega(M).

## More issues
- Lemma 3.1 -> "no counting argument proves L5": Lemma covers only n_q > M-1+4|V| (inside Q_large region), says nothing on
  Q_med where Cor 3.3 localises L5. Single-modulus claim trivial (one modulus contributes <= A_q log q <= log X). OVERSTATED.
- Lemma 3.4: per-modulus 1.5/n; summed over Q_large ~ 3 log A expected failures (not o(1)); Thm 3.2 needs Hall at EVERY q.
  Sharper: s=2 term dominates, P <= O_t(M^{4-2t})/n -> simultaneous Hall whp for t>=3 (model). minor; fixable.
- Hadamard cores (J+H)/2 include the constant all-ones column -> NOT valid profiles (min_x a_xj = 0 fails); the valid
  Hadamard-core profile (nonconstant columns, copies) has rank M-1, Delta_M = 0: Thm 5.1/Prop 5.3 don't apply.
  Also d_M = M/2 is a THEOREM for normalized H (Sherman-Morrison: B^{-1} = 2H^T/M - e1e1^T), not just 5 cases.
- Rank-deficient profiles (actual clusters may have row relations) not covered by §5.
- "cannot help a construction" (non-unimodular): interpretive, overstated.
