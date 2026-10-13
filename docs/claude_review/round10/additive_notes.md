# Notes for additive.md (working file)

## Plan
1. Collapse theorems: relation lattice <=> lattice coset; rigidity; AG(H) decides; rational conic => order; mod-N small lifts.
2. Tangent-range lattice geometry: Prop 3.1 (monomial gaps => mu <= W/4 + log(sqrt2 C) off one matching) [class A_0];
   Thm 3.2 (effective height h*: M <= C^2 h* + 2, mu <= 2 log(2 C h*) off chords parallel to i v*), Cor (R1-arch on h* <= M^b).
3. Relation p-adics: Lemma 5.3 generalised (Newton polygon + leading residue cancellation), isolation lemma.
4. Barrier map.

## Key derivations
- w(b_xy) <= W/2 - mu_xy where b_xy = a_x + a_y - e (since |a+a'-e| <= e - |a-a'|).
- arg A_{b_xy} = (theta_x+theta_y)/2 - (pi/4)(k_x+k_y).
- |u_xy| <= |g_xy||gamma_s| <= 2|g_xy| since B | (z_x - z_y) as rational integer.
- Prop C.2 with v radial: bound |v| C (N^{1/4} eta + C/2) + 1 <= C^2 h(v) + 2, h(v) = |v|(1 + N^{1/4} eta/C), C := D N^{1/4}.
- mixed gap: |u| >= 1/(|v| sin(eta + D)) >= N^{1/4}/(C h(v)).
- Minkowski: h* <= (2 pi/C)^{1/2} N^{1/8} + 1.

## AG ladder (for §2.5)
- rank bound: AG(2(M-1)H) => (H+1)^(rho-2) <= 2^M M^2 L^2 (from Thm 2.4(2) proof with kappa = M-1-rho).
  => rho <= 2 + (W/2 + M log 2 + 2 log(ML))/log(H+1). With H = M^A, W = c0 M log M: rho <= 2 + (c0/(2A)) M + O(M/log M).
- generic fakes: P(|sum c w| < 1) ~ 1/(|c|^2 C^3 N^(1/4)) per c (Jacobian heuristic) => AG(M^A) holds w.h.p. for A < c0/4.
- relation-imposed: k relations = 2k equations on M-1 positions, generic solvable iff k <= (M-1)/2 i.e. rho >= (M-1)/2;
  then AG fails at height ~ (2ML)^(2/(rho-2)) ~ M^(c0 (1+o(1))). So A > c0 needs overdetermined (non-generic) positions.
- CM-field circles: rank 2d modules, but AG fails at small height (other embeddings make small sigma_1).
