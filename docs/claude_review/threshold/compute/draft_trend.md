## 4. The trend in t at fixed M: proofs

Notation: `S = S^M(p,n)`, `t = p + n`, `Delta_a(q) = {x in Z_{>=0}^a : sum x = q}`.

**Lemma 4.1 (peeling; upper.md Lemma 2.1).** If `A ⊂ Z(L_1) ∪ ... ∪ Z(L_r)` for affine
functions `L_i`, and `A_i = A ∩ Z(L_i) \ (Z(L_1) ∪ ... ∪ Z(L_{i-1}))`, then
`reg(A) <= max_{A_i nonempty} (reg(A_i) + i - 1)`.

*Proof.* Induction on `r`: interpolate `g` on `A_1` by `P_1` (degree `reg A_1`), then
`(g - P_1)/L_1` on `A \ Z(L_1)` by `P'` (induction), and take `P_1 + L_1 P'`. `[]`

**Lemma 4.2.** (a) `reg(A x B) <= reg(A) + reg(B)`. (b) `reg Delta_a(q) = q` for `a >= 2`,
`= 0` for `a = 1`. (c) `reg` is monotone under inclusion and invariant under affine bijections.

*Proof.* (a) `C^{A x B} = C^A ⊗ C^B` is spanned by products. (b) `<=`: peel
`Delta_a(q)` by the hyperplanes `x_1 = 0, 1, ..., q` (in this order): the `c`-th piece is
`Delta_{a-1}(q - c)`, of regularity `<= q - c` by induction, so the bound is
`max_c (q - c + c) = q`. `>=`: `Delta_a(q)` contains `q + 1` collinear points, and a polynomial of
degree `< q` cannot vanish at `q` of them without vanishing at the last. (c) Clear. `[]`

**Proposition 4.3 (shell regularity is t + O_M(1)).** Let `M >= 3`, `p >= n`, `t = p + n`.

1. If `n = 0`: `reg S = t`.
2. If `n >= 1`: `reg S <= max( t + F - 1, p + F + M - 1, n + 2^M - 3 )`, where
   `F = 2^M - 2 - 2M` is the number of sign vectors with at least two `+` and two `-`.
3. (line bound) `reg S <= reg Q^M(p,n) <= p + (M-1) n = t + (M-2) n`.
4. If `M >= 4` and `n >= 1`: `reg S >= t`. If `M = 3`: `reg S >= p`.

Consequently `max_s reg(S_{s,t}) / t -> 1` as `t -> infinity`, for every fixed `M >= 3`.

*Proof.* 1. `S^M(p,0) = Delta_M(p)`.
2. For a sign vector `eps in {±1}^M` put `H_eps = {sum_j eps_j k_j = t}`. Since
`sum eps_j k_j <= sum |k_j| = t` with equality iff every nonzero `k_j` has sign `eps_j`,
`S ∩ H_eps ≅ Delta_a(p) x Delta_b(n)` (`a`, `b` = numbers of `+`, `-` in `eps`), which is
nonempty iff `a, b >= 1`, and `S` is covered by these `2^M - 2` hyperplanes. By Lemma 4.2 the
piece has regularity `<= p[a >= 2] + n[b >= 2]`: `<= t` for the `F` vectors with `a, b >= 2`,
`<= p` for the `M` vectors with `b = 1`, `<= n` for the `M` vectors with `a = 1`. Peel in this
order (Lemma 4.1; the pieces `A_i` are subsets of `S ∩ H_eps`).
3. Restrict a Laurent polynomial `f = sum c_k w^k`, supported in `Q^M(p,n)`, with
`ord_{w=1} f = m`, to the line `w = 1 + lambda v`:
`f(1 + lambda v) prod_i (1 + lambda v_i)^n = sum_k c_k prod_i (1 + lambda v_i)^{k_i + n}` is a
polynomial in `lambda` of degree `<= sum_i (k_i + n) = p + (M-1) n` that vanishes to order `>= m`
at `lambda = 0`, and it is nonzero at `lambda = 1` if `v = w_0 - 1` with `f(w_0) != 0`. So
`m <= p + (M-1) n`; apply Lemma 2.1 (this is upper.md Remark 3.3). The bound for the shell follows
by monotonicity.
4. `S` contains `{(x, p - x, -y, -(n - y), 0, ..., 0)} ≅ Delta_2(p) x Delta_2(n)`, a
`(p+1) x (n+1)` grid, whose regularity is `p + n` (the measure
`(-1)^{x+y} binom(p,x) binom(n,y)` annihilates every monomial `x^i y^j` with `i < p` or
`j < n`). For `M = 3`, `S` contains `p + 1` collinear points. The limit statement follows from
2, 4 and 1 (for `M = 3` use the shells with `n = 1`, `reg >= p = t - 1`). `[]`

**Corollary 4.4 (finite verification for M <= 5).** Let `V(M) = 2 - 1/ceil(M/2)`.

* `M = 3`: by 3, `reg S <= t + n <= 3t/2` for every shell.
* `M = 4`: a shell with `reg S > 3t/2` has `n > t/4` (by 3) and then, by 2
  (`F = 6`: bound `max(t + 5, p + 9, n + 13)`), `t <= 12`.
* `M = 5`: a shell with `reg S > 5t/3` has `n > 2t/9` (by 3) and then, by 2
  (`F = 20`: bound `max(t + 19, p + 24, n + 29)` with `p < 7t/9`, `n <= t/2`), `t <= 28`; and
  `s = t - 2n < 5t/9`, i.e. `s <= 15` for `t <= 28`.

Hence the computed range (M = 3: t <= 30; M = 4: t <= 26; M = 5: t <= 24 all s, and
t = 25..28 with s <= 15) covers every shell that could exceed V(M), and the shell threshold
`tau_sh(M) = sup reg(S_{s,t})/t` is

```text
tau_sh(3) = tau_sh(4) = 3/2,      tau_sh(5) = 5/3          (computer-assisted theorem).
```

The upper bounds come from modular ranks (rigorous direction), the lower bounds from the
Vandermonde (alternating-orbit) measures. For `M = 6` the same argument would need
`t <= 73` (since `F = 50`), beyond the computed range.
