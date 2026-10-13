# r8 notes (attempt 2 at L_7, 2026-10-13)

Previous attempt (L7.md, r7/): cub13 exact over quintic field; reduction: need ONE rich addable
family (Theorem 3 there): Q-family V with (B+8K).V >= 0 and ev: U_V -> X having rational section.
For 3-dim families: addable iff ev has degree 1.

New idea: MOBIUS CURVES. x_a(t) = g_a(t), g_a in PGL_2, a = 1..7. Space PGL_2^7/(PGL_2 x PGL_2), dim 15.
pair {a,b} collides at 2 fixed points of g_b^{-1} g_a. Pure class: D_ab = 2 for all 21 pairs.
Triple event {a,b,c}: common fixed point; costs 1 condition. 4-event costs 2.
(B+8K) = 6 tau3 + 11 tau4 + 15 tau5 - 70, dim = 15 - tau3 - 2 tau4 - 3 tau5.
Rich 3-dim types: (tau3,tau4) = (12,0) [B+8K=2], (10,1) [1], (8,2) [0].
Geometry: graphs are (1,1)-curves on P1xP1 = plane sections of quadric Q in P3; event = triple point of
7 planes lying on Q.
Through a general point x=(y_1..y_7) at u=inf: G_a = g_a^{-1}, G_a(z) = (p_a z + q_a)/(z - y_a).
Event e at (z_e,u_e): u_e (z_e - y_a) = p_a z_e + q_a for a in e: points (z_e, u_e(z_e - y_a)) collinear.
Unknowns 24 - 2 gauge = 22, conditions sum(r_a - 2) = 36 - 14 = 22. Finite.
deg(ev) = number of solutions. If 1 (rational construction) -> addable rich family -> L7 false (Thm 3 r7, G1-G3).

## Results so far
- enum_types.py: 23 iso types of 12 distinct triples on [7] with pair mult <= 2 (types12.txt).
- grow.py (continuation adding events): nondegenerate Mobius-12 members found for 16 types
  (0,1,2,3,4,5,6,7,8,9,11,13,17,19,20,22); not found (yet) for 10,12,14,15,16,18,21.
- mono.py monodromy: deg(ev) >= 4..29 for all tested (type1>=4, type11>=4, type2>=5, type5>=10, type4>=29).
  => Mobius-12 families dominant (JX nonsingular) but NOT addable.
- Barrier (free-linear): graph ansatz through a point with r_a <= 2 d_a (linear given event positions)
  gives Q(X)-rational Z_x but max(6n3+11n4+15n5) <= 3 sum|e| <= 6D < 10D needed for rich. Never rich.
- Class beta in graph ansatz: forced d_a = 3 all a, 35 triple events, 21 pure pairs; dim 43-35 = 8.
  NEXT: try continuation (grow) for class beta itself (round 4 failed with other solvers).

## KEY IDEA (G(2,4) model) -- to verify
For 7 planes H_a in P^3 (arrangement A over Q), Phi_A: G(2,4) -> X, R |-> (R cap H_a)_a. Inverse: given x,
L in Hom(C^2,C^4) with l_a(L v_a) = 0 (7 linear eqs on 8 unknowns) -> unique line R_x (if rank 7; fails if
a point lies on >= 6 planes). So X ~ G(2,4) birationally over Q.
Pair D_ab <-> sigma_1(l_ab) (lines meeting axis l_ab = H_a cap H_b); triple D_abc <-> lines through P_abc (codim 2)
BUT coaxial triple {a,b,c} (common axis l): sigma_1(l) contacts are TRIPLE events: free richness.
Mobius curves = rulings of quadrics Q; event at P_abc iff P_abc in Q: LINEAR in Q.
STAR arrangement: H_O plane, triangle of lines l1,l2,l3 in H_O, two more planes through each l_i:
labels O,a1,a2 (l1), b1,b2 (l2), c1,c2 (l3). 3 coaxial triples (each 2 events), choose 6 of the 8 points
P_{ai bj ck}. Quadrics through 6 points: P^3 -> 3-dim family; class pair 6, triple 12 => (B+8K)=+2 RICH.
Through general x: R_x unique, Q through R_x + 6 pts unique (9 linear conds) -> deg ev = 1 -> ADDABLE over Q.
=> with Thm 3 of r7 (G1-G3): L7 FALSE?! Must verify exactly.

## VERIFIED (exact, star_run.py seeds 2,5): star arrangement, chosen 6 = (1,3,5),(1,3,6),(1,4,5),(1,4,6),(2,3,5),(2,4,5)
(labels 0=O,1=a1,2=a2,3=b1,4=b2,5=c1,6=c2). Member through random rational x: Q-rational Mobius curve,
contacts: O a1 a2 x2, O b1 b2 x2, O c1 c2 x2, six triples, pure pairs a2b1,a2b2,a2c2(x2),b1c2,b2c2.
pair 6 triple 12 (B+8K)=2, 18 distinct contact pts, no contact at infinity.
CAVEAT: 6-subsets that are union of two faces (e.g. H_a1 u H_b1 contains them) give rank-2 quadric. Avoid.
The 8 cube points P_{ai bj ck} = base locus of net <H_a1 H_a2, H_b1 H_b2, H_c1 H_c2>.
NEXT: re-read round4 Thm 2(a) + proof; re-derive Thm 3 (iterated addition) rigorously incl. genericity.

## MAJOR CORRECTION to r7 Theorem 3 (iterated addition)
Lemma 4.2 of r7 assumes ruled surface Y = C x P^1 with section Sigma of class h + delta F (Sigma^2 >= 0).
(G1) "iota defined along Sigma" FAILS at every contact of C with a divisor not covered by the family
(e.g. forgetful f_i: contacts D_S with i in S, |S|>=3; V_A: the 26 non-coaxial non-event triples and
pairs with D.gamma_A = 0). Honest surface = pullback family over C with reducible fibres; Sigma^2 = -K_{X/B}.C.
Forgetful: Sigma^2 = -psi_i.C + sum_j D_ij.C = -9 lambda on C ~ lambda beta. Irreducible sections in
|Sigma + mF| need m >= -Sigma^2. Switching components at a special fibre: +1 to Sigma^2 but class + e_abc
(e: D_abc -1, D_ab,D_bc,D_ac +1), i.e. loses triple contacts. => Thm 3 invalid; addable+rich NOT sufficient.
Last-step necessary condition: -K_{Xt/B}.gamma_final >= m > 0.
V_A is rich + addable: so the r7 "reduction to one rich addable family" was wrong (V_A would refute otherwise).

## lambda_A(beta) (elamp.py, mod 2^61-1, 2 curves per class agree):
E_Lambda.mu = 18 (pure Mobius), sum_T E.mu_T = 605, beta = sum_T mu_T - 32 mu  => E.beta = 29,
lambda_A(beta) = 7 - 2*29 = -51 (K_B = -2 f^*H for the double cover of the web). forgetful: -9.
Kapranov projection from general pt: 7 - 4*15 = -53; from p_i: 7 - 4*14 = -49.
Product chart (birational): lambda_k(beta) >= 2 c_k > 0.
E.f_j (O,a1,a2,b1,b2,c1,c2) = 4,3,2,2,2,3,2 ; sum = 18 = E.mu (consistency: mu == sum_j f_j numerically). lambda_A(f_j) = 2-2E.f_j all < 0.

## ISSUE: lambda_A(beta) = -51 was (-K_X + pi^*K_B).beta; true -K_F differs by sum_E (r_E - 2) E for divisors
contracted to codim-r centres in U (vector field blowup: F_X = sigma^*F_U(-E) for transverse centres).
Intrinsic: lambda_X(C) = Sigma_U^2 - k(C), k = # contacts with the 42 uncovered divisors.
Sigma_U^2 via event sections: Sigma^2 = 2 Sigma.s1 - s1.s2 - s1.s3 + s2.s3 (P^1-bundle sections).
NEXT: lamA_direct.py computing this exactly mod p on Mobius test curves.

## lamA_direct.py FAILED consistency: along general Mobius test curves Q_t passes through rank-2 quadrics
(R_t in H_a1 or H_c1, planes with 4 event points) -> leaf surface has degenerate fibres; event-section formula
invalid. Local corrections: blowup of radial singular centre (sigma_2(P_e)): F_X = sigma^*F(+E);
transverse regular centre: F_X = sigma^*F(-E). Too many divisor types (41) to settle now.
DECISION: report lambda_A(beta) sign as UNDETERMINED; naive (-K_X + pi^*K_B).beta = -51.
r7a Thm 3 with V_A: forgetful steps must precede (lambda_{f_i}(tau gamma_A) in {-2,-3,-4} exact), then
V-steps need lambda_A(mu) >= 0 and lambda_A(beta) >= 0 (necessary); undecided.

## kfol.py: c_1(F_G) = -2 sigma_1 on G(2,4) (two random pencils, exact). dR.py: d_R(beta)=21, d_R(mu)=8,
d_R(f_j) = 2,1,1,1,1,1,1, d_R(gamma_A)=2. lambda_A = Phi^*(-2 sigma_1) + sum_E c_E E over 41 exceptional
boundary divisors (if no other exc. divisors). lambda_A(gamma_A)=2 => c=+1 at the 6 event divisors (radial).
lambda_A(beta) = -42 + sum_41 c_E ; c_E <= m_E (vanishing order of generator along centre).
>= 0 needs some c_E >= 2. lambda_A(mu) = -16 + 2 sum_{9 coaxial pairs} c_E.
