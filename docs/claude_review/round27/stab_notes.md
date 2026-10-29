# stab lane notes (parking derivations)

## N1 Cox form (verified by hand, to be checked in stab_checks/cox_form.py)
- cuts indexed by A subset [4] nonempty (side not containing 5). Layer (p,t) with S={j: a_j(p)>=t}:
  if 5 not in S, A=S, gamma_A gets pi_p; if 5 in S, A=S^c, gamma_A gets pibar_p.
- chain condition (Ch): per p, {A: pi_p | gamma_A} chain (decreasing), {A: pibar_p|gamma_A} chain, unions disjoint.
  (Ch) <=> primitivity (proof: pick a in smallest pibar-set; it avoids largest pi-set).
- Pi_a = prod_{A ni a} gamma_A, zeta_a/zeta_5 = Pi_a/conj(Pi_a), t_{a5} = Im Pi_a (integer), |Pi_a|^2 = N_a = e^{d_a5}.
- Pair cofactor: Q_ab = prod_{A ni a, not b} gamma_A prod_{A ni b, not a} conj(gamma_A); zeta_a = G_ab' Q_ab, zeta_b = G_ab' conj(Q_ab).
  Q_ab conjugate-primitive under (Ch). t_ab = Im Q_ab (up to sign), N(Q_ab) = e^{d_ab}.
- Triangle: Q_ab Q_bc = G^{(b)}_{abc} Q_ac, G^{(b)} = prod_{A separates b from {a,c}} n_A (positive integer, 4 cuts).
  With c=5 (Q_{b5}=Pi_b): Pi_a conj(Pi_b) = G_ab Q_ab, G_ab = prod_{A ni a,b} n_A.
  Im: x_b t_a - x_a t_b = G_ab t_ab (x = Re Pi). Re: x_a x_b + t_a t_b = G_ab y_ab.
- Quadratic-form identity: M_{b\a} t_a^2 - 2 y_ab t_a t_b + M_{a\b} t_b^2 = G_ab t_ab^2, disc = -4 t_ab^2.
- Venn form (pairwise coprime blocks): P_a = x_a + i t_a; gamma_A = part of prod dividing exactly {P_a: a in A}.
- Arccot form (t=1): (x_b - x_a) | x_a^2+1 ; cot(theta_a - theta_b) = z_ab integer.
- singleton decoupling: Im(gamma_a R_a) = t_a, R_a = prod_{A ni a,|A|>=2}; x_a = -iota(R_a) t_a mod N(R_a), need |x_a| ~ X^4.

## N2 Equivalences
- Stab(eta,eps): eps-balanced one-class primitive with max|t_a5| <= e^{eta w} => in Z.
- (L5circ)(c,eps) => Stab(eta,eps) if 10 eta + 48 eps < c [t_ab bound |t_ab| <= 2 X^{eta+8 eps}].
- Stab(eta,eps) => (L5circ) with c = eta, eta_0 = eps.
- A_5(1/30 - eta) => Stab(15 eta/2 - C eps, eps).

## N3 MAIN: second-order integrality = tight pencils; Theorem R kills it
- Theorem S: if f1,f2 MW-tight with same support class, c_i(y)=f_i(y)/D(y) in Z[i], |c_i| <= C X^{O(eta)};
  phi = f1/f2 takes values c1/c2 of height X^{O(eta)}; bounded T => finitely many fibres.
- Theorem R (PROVED, all n): for finite S in Z^n, dim A_{MW_n(S)}(S) <= 1.
  Proof: slice x_1; decreasing-g order; peeling dim A_t(S) <= sum_i dim A_{t-i+1}(S_(i));
  critical i: g_(i)+i-1 = K. At most one critical i: improved (2.5)
  K >= (c'-c) + (g_c+g_c')/2 + 2^{1-n}((cmax-cmin)-(c'-c)).
  If i<j critical: C_{g(i)}=[b,b+i-1] strictly inside C_{g(j)}=[a,a+j-1] (extremes have g=g(j));
  pairs (a,b+i-1),(b,a+j-1) give b-a<(j-i)/2 and a+j-1-b<(i+j)/2-1, sum j-1<j-1. Contradiction.
- Consequence: every effective L=pi^*D - tE on Bl_e X with D.beta = 8t has h^0 <= 1 (h^0(mL)<=1): rigid face.
- non-homogeneous Cox combos do not descend; homogeneous => same class => same generalized polytope => Thm R.

## N4 fibrations
- defect pencils: s = D.beta - 8t >= 1, image height (s/2) w; triangle pencil s=4 (fibres = pr_123^{-1}(curve), useless).
- lines through e: exactly tight (s = 2i/conj(Gamma)).
- gluing: fibrewise sections -> H = sum f1^j f2^{n-j} g_j; archimedean bound direction-sensitive.

## N5 counts
- gamma-count X^{15-16}; sequential: Pi_1 X^4, Pi_2 X^0, Pi_3 X^{-2}, Pi_4 X^{-3}.
- singleton elimination: <= X^{11} T^4 trivially.
