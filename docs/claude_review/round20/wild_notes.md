# wild notes (round 20 wildcard)
## Screening (turn 1)
Record grep: Zilber, Chabauty, Skolem, Frey, Dobrowolski, Hilbert symbol absent. Masser/Habegger only in passing.
Screened candidates:
- Pila-Wilkie/o-minimal: definable set in exponent space is semilinear (slab) -> algebraic part = everything; algebraicity of phi invisible. dead (item 4 shape).
- Chabauty/Skolem: no rank deficiency (T 1-dim, Gamma rank r~M); archimedean inequality not an equation; equalities are S-unit eqs (item 5).
- Erdos-Turan/Bilu equidistribution of roots: (X-1)^D shows integer polys may have all roots at one point; F=prod(conj(P)X-P) ratio 2^M -> trivial.
- Elekes-Szabo/sum-product: circle group is the exceptional group case. dead.
- Hilbert/quartic reciprocity: residue symbols of pi_j mod varpi = [varpi/pi_j]_4 -> subsumed by genuine residues (P3-G), item 1.
- Frey/modular: conductor ~ rad, squarefree => abc-vacuity row 15.
- Strip geometry of numbers: M <= 2(2/lambda_2+1), lambda_2 >= (C^3 sqrt R)^{-1/2} -> R^{1/4}; near-axis only (Chan-type).
## Main candidate: unlikely intersections (BMZ / Habegger bounded height), S-uniform
Natural varieties: X_M cap torus = coset (V^oa empty); chord linear space L_M: fixed-shape triangles anomalous -> V^oa empty;
M_{0,k} Plucker: fixed 5-point shapes anomalous -> V^oa empty.
Better: the Q-line ell_z = span(X,Y) (z = X+iY) contains u=(z_1:...:z_M). Curve! Not in coset iff no two z_a,z_b R-proportional.
Coordinates in Gamma rank r; if r < M-1 point in subgroup codim >= 1 -> BMZ bounded height for curve ell, constant c(ell) depends on h(ell).
Heights: h_sum(u) = sum_x h(z_x/z_1) ~ sum d_{1x}/2; h(ell) = log max |Im(conj z_a z_b)| - log gcd <= (3/4)W + log C.
Full circle (all lattice pts, r primes): ratio h_sum/h(ell) ~ (M-1)/4 ; short arc ~ (M-1)/3 (to verify).
Needed: sharp-constant line BMZ:  h_sum(u) <= (n/4 + o(n)) h(ell)  for u in ell cap (subgroup codim>=1)?  -- check plausibility / counterexamples.
Rank condition not automatic: near-extremal omega(N') in [M/4, 2M/eps]; fakes have r >= M-1 (not algebraic anyway).
## Turn 2 decisions
Main mechanism = S-uniform unlikely intersections (BMZ 1999 bounded height on curves; Habegger BHC; Remond cones).
Results to prove:
 Lemma P (projection anomaly): V subset torus T, phi: T->T' homomorphism, f=phi|V with dim f(V) < dim T' and dim f(V) < dim V
   => every fibre component Y of f is anomalous (Y in coset of ker phi; dim Y = dimV - dim f(V) > max(0, dimV - dimT')). Cover => V^oa empty.
 Thm A: V^oa empty for (a) X_M cap torus (coset), (b) chord space L_M (M>=4; triangle shapes), (c) M_{0,k} in intrinsic torus (k>=5; forgetful fibres to M_{0,4}/M_{0,5}).
 Thm B: carrier curve = Q-line ell_z=span(Re z, Im z); ell_z in coset iff two coords R-proportional; cluster: h(z)=W/2, h(ell)<= (3/4)W+log C => ratio >= 2/3 - o(1).
 Prop C: unconstrained (point on its own line, in subgroup of codim n-1): u=(1,g,g^2), g=A+i: h(u)/h(ell_u) -> 1. So any UI-type inequality h(u) <= kappa h(ell)+c(n) needs kappa>=1; cluster needs kappa<2/3. Capped.
   With circle constraint: full circle ratio 1/2; constraint = equal norms = Ptolemy/(L_k) (items 3,7).
 Also rank deficiency absent for full-rank clusters.
Checks: wild_checks/anomaly_ranks.py ; wild_checks/line_examples.py (rank-1 examples, full circles, 10-pt cluster n=1176852625).
## Turn 3 results (exact, cluster_heights.out)
F1=(1,g,..,g^{n-1},g^K): proj ratio >= K/(K+n-1) (|minor|<=|u_a||u_b|, G>=1). (1,g,g^2): ratio exactly 1.
F2=(1,g^{K-n+1..K}): anchor-0 sum ratio >= n(K-n+1)/(2K) -> n/2. BUT min-anchor sum ratio ~1/2 (anchor in block).
Cluster: h(z)=log R exactly; h(ell)<=(3/2)log R+log C; every-anchor sum >= n((log R)/2 - log C). Data: G=1, ratio = bound.
Full circles ratio exactly 1/2 (G=1). 10-pt cluster rank 5 < 9 (rank-deficient).
Mixtures alpha*h + beta*h^(0)_Sigma: max(F1,F2) >= cluster lower bound (2/3 alpha + n/3 beta). Proven by cases alpha>=(n-3)beta / alpha<=n beta.
S_min (min-anchor sum) NOT refuted; it encodes all-pairs margins (pair data) + covers unequal norms; unconstrained cluster-like tuples need size-equality as extra condition -> rarer than circle clusters. Record as open non-standard, not UI-type.
Uniform-from-one-n0: subcluster W' >= (n0/(n0+1))W - O(1) by Thm 3.11 gcd_B bound => (UI-l) at one n0 gives uniform bound.
Full rank: if rank B = n, z in no proper algebraic subgroup (gamma_j indep mod torsion).
