# The auxiliary-polynomial threshold tau(M) (lead author, 2026-10-08)

## Reduction (claimed; to be verified)

Let X_M = {x_1^2+y_1^2 = ... = x_M^2+y_M^2} in P^{2M-1}, Y = diagonal line {z_1 = ... = z_M}.
Write z_j = x_j + i y_j. Suppose F is a homogeneous polynomial of degree D with coefficients in
Z[i] in the 2M variables, not vanishing identically on X_M, whose restriction to X_M vanishes to
order m along Y, with **m > 2D**. Then the uniform theorem holds (for every C, a bound M(C)
independent of R and arc position):

1. Liouville. For lattice points z_1..z_M on |z|=R inside an arc of length C sqrt(R),
   F(z) is in Z[i] and, by Taylor expansion along Y on the compact X_M(C),
   |F(z)| <= K R^D (C R^{-1/2})^m = K C^m R^{D - m/2} -> 0. So F(z) = 0 once R >= R_0(C,F).
2. Grid lemma. F is not identically zero on the product of circles C_N^M (by homogeneity: the
   slice t = N of the cone {u_j v_j = t} determines all slices, and (S^1)^M is Zariski dense in the
   torus). So F' = F * prod_{i<j} |z_i - z_j|^2 is nonzero on C_N^M and vanishes on S^M for every
   cluster S with R >= R_0. Bezout on the conic, coordinate by coordinate, gives |S| <= 2 deg F'.
3. Small R: finitely many circles, finitely many points.

## Combinatorial reformulation

On the circle, z_j = R e^{i theta_j}; a homogeneous F of degree D restricted to X_M is a Laurent
polynomial sum_k c_k N^{(D - ||k||_1)/2} z^k with ||k||_1 <= D, ||k||_1 = D mod 2. Writing
z_j = z w_j along the diagonal, the order of vanishing along Y decouples by s = sum k and by the
shell t = ||k||_1 (the N-weights separate shells). So

  tau(M) := sup over (s, t) and nonzero c on S_{s,t} = {k in Z^M : sum k = s, ||k||_1 = t}
            of  ord(c)/t,
  ord(c)  := max m such that sum_k c_k P(k) = 0 for every polynomial P of degree < m
             (equivalently the Laurent polynomial sum c_k w^k vanishes to order m at w = 1).

The reduction says: **tau(M) > 2 for some M implies the uniform theorem.** tau(M) is the
pseudo-effective threshold of pi^*H - t E on Bl_Y X_M (sup of ord_Y F / deg F), which also
dominates the Ru-Vojta beta constant (beta <= threshold); the referees noted that beta_M > 2
would give the theorem and is open, so tau(M) > 2 is the weaker, sufficient condition.

## Known data

* Vandermonde / Ramana: c = alternating sum over the S_M-orbit of v = (-s', ..., s'), M = 2s'+1,
  t = s'(s'+1), ord = C(M,2): ratio (2s'+1)/(s'+1) -> 2 from below.
* Exact computation (mod p rank, lead author): general measures M=2..5 (t <= 9), and alternating
  measures M <= 7 (t up to minimal + 12): the maximum ratio is exactly the Vandermonde value
  (1, 1.5, 1.5, 1.667, 1.667, 1.75 for M = 2..7). Scripts: scratchpad/lead/shell_reg*.py,
  alt_reg.py, sym_reg.py.
* Dimension count (generic) gives ratio < 2 M^{-1/(M-1)} < 2; special structure beats it.
* Movable-curve (BDPP) upper bounds: lines/conics through a point of Y give only O(M); nothing
  proves tau <= 2 yet. Counting heuristics suggest no movable rational curves with contact ratio
  >= 1/2 exist for M >= 5, so tau(M) > 2 is NOT excluded.

## Questions

(A) Prove tau(M) <= 2 for all M (a barrier: no single-place auxiliary polynomial proves the
    theorem), or determine tau(M) exactly for small M.
(B) Construct c with ord(c) > 2t for some M (this proves the uniform theorem).
(C) Generalisations that keep the Liouville argument valid: several shells with N-dependent
    coefficients; polynomials whose values at cluster points carry forced Gaussian divisibility
    (e.g. products over pairs inherit prod g_ij) so that a ratio below 2 already suffices once the
    forced divisor is counted.
