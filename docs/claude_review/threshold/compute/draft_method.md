## 2. Method: the S_M-isotypic reduction

Throughout, `A` is a finite S_M-stable subset of `{k in Z^M : sum k = s}` (a shell or a slice).
`S_M` acts on points by permuting coordinates and on functions/measures and polynomials by
the induced actions. The pairing `<c, P> = sum_{k in A} c(k) P(k)` is S_M-invariant.

**Lemma 2.1 (duality).** `max{ord(c) : 0 != c in C^A} = reg(A)`, where
`reg(A) = min{d : Poly_{<=d}|_A = C^A}` and `ord(c) = max{m : c ⊥ Poly_{<m}}`.

*Proof.* `W_m := Poly_{<m}^⊥` (annihilator in `C^A`). `W_m != 0` iff
`Poly_{<m}|_A != C^A` iff `m - 1 < reg(A)`, i.e. iff `m <= reg(A)`. `[]`

For a partition `lam` of `M` fix the row-reading tableau `t` (positions `0..M-1` filled row by
row), its row group `R_t`, column group `C_t`, and

```text
a_t = sum_{r in R_t} r,      b_t = sum_{c in C_t} sgn(c) c,      e_lam = a_t b_t.
```

Standard facts (Young symmetrisers; e.g. Fulton–Harris §4.1, James §4):

* (Y1) `e_lam^2 = h_lam e_lam` with `h_lam = M!/f^lam != 0`; for every finite-dimensional
  `C[S_M]`-module `U`, `dim e_lam U = mult(S^lam, U)`; in particular `e_lam U != 0` iff `U`
  contains `S^lam`.
* (Y2) (Young's rule) for an orbit `O` with value multiplicities `mu`, `C^O = Ind_{S_mu}^{S_M} 1`
  and `mult(S^lam, C^O) = K_{lam,mu}` (Kostka number = number of semistandard tableaux of shape
  `lam` and content `mu`).
* (Y3) the adjoint of `a_t` (resp. `b_t`) for the invariant pairing is `a_t` (resp. `b_t`),
  because `g -> g^{-1}` preserves `R_t`, `C_t` and signs. Hence `<a_t b_t z, P> = <z, b_t a_t P>`.

**Lemma 2.2 (isotypic reduction).** Put `E_lam = e_lam C^A` (dimension
`N_lam = sum_O K_{lam,mu(O)}`). For `d >= -1` let
`H_lam(d) = dim{ (<x,P>)_x : P in Poly_{<=d} }` = rank of the pairing `E_lam x Poly_{<=d}`.
Define `reg_lam(A) = min{d : H_lam(d) = N_lam}`. Then

1. `reg_lam(A) = max{ord(c) : 0 != c in the lam-isotypic component of C^A}`;
2. `reg(A) = max_lam reg_lam(A)`, and `h_A(d) = sum_lam f^lam H_lam(d)` (Hilbert function).

*Proof.* `W_m` is S_M-stable. By (Y1), the lam-isotypic part of `W_m` is nonzero iff
`e_lam W_m != 0`. Since `e_lam^2 = h e_lam`, `e_lam W_m = W_m ∩ E_lam`. And
`W_m ∩ E_lam != 0` iff the pairing `E_lam x Poly_{<m}` is degenerate on the `E_lam` side iff
`H_lam(m-1) < N_lam`. So the lam-part of `W_m` is nonzero iff `m <= reg_lam(A)`; this is 1.
Item 2: `W_m != 0` iff some isotypic part is nonzero; the Hilbert-function identity is the
isotypic decomposition of `Poly_{<=d}|_A` (each copy of `S^lam` contributes `f^lam` dimensions
and one dimension of `e_lam`-image). `[]`

**Columns.** For each orbit `O` and each semistandard tableau `T` of shape `lam` whose content is
the value multiset of `O` (values ordered as integers), let `k_T` be the point carrying the
entry of `T` in box `b` at position `t(b)`, and `c_T = a_t b_t 1_{k_T}`. There are exactly
`N_lam` of them. In every run the final rank of the generator matrix below equals `N_lam`,
which proves (rank mod p <= rank over Q) that the `c_T` are linearly independent, hence a
basis of `E_lam` (this is the transposed form of James's semistandard basis theorem; we do not
need the theorem, the rank check proves it in each case).

**Rows (higher Specht polynomials).** For a standard tableau `S` of shape `lam` let `i(S)` be
the Ariki–Terasoma–Yamada index tableau (read the columns of `S` bottom-to-top, left to right;
`1` has index 0; `k+1` has the index of `k` if it stands to the right of `k` in this word,
otherwise one more) and `charge(S) = sum of indices`. Put
`F^S = b_t a_t x_t^{i(S)}` (`x_t^{i(S)} = prod_b x_{t(b)}^{i(S)(b)}`), a homogeneous polynomial
of degree `charge(S)`. For a monomial `g` in the power sums `p_2..p_M` of weighted degree `w(g)`,
`b_t a_t (g x^{i(S)}) = g F^S` and therefore

```text
<c_T, g x^{i(S)}> = g(k_T) F^S(k_T) = g(v_O) F^S(k_T),        degree w(g) + charge(S).        (2.1)
```

The generator matrix `G_lam(d)` has the rows (2.1) of degree `<= d`. Only the restriction of
`g` to the orbit representatives matters, so `g` runs over a greedy basis (in weighted-degree
order) of `C[p_2..p_M]|_{orbits}`; `p_1 = s` is constant on `A`.

**Rigour.** Every row of `G_lam(d)` is a pairing with a polynomial of degree `<= d`, and
`rank_p <= rank_Q`. Hence `rank_p G_lam(d) <= H_lam(d)` and

```text
reg_lam(A) <= reg_lam^{(p)}(A) := min{d : rank_p G_lam(d) = N_lam}     (rigorous upper bound).
```

Equality needs (a) completeness: the rows (2.1) span the full pairing space. This holds by the
ATY theorem (the `F_T^S` form a basis of the coinvariant algebra and, for fixed `S`, span a copy
of `S^lam`; with `T = t` fixed the `F^S` therefore form a basis of `b_t a_t` of the coinvariant
algebra, and graded Nakayama gives `b_t a_t Poly = Sym . span{F^S}`), and it is confirmed
independently by the cross-check `sum_lam f^lam h_lam^{(p)}(d) = h_A^{(p)}(d)` against a direct
(non-equivariant) Hilbert function (§2.4); and (b) `rank_p = rank_Q`, supported by two primes
`P1 = 2^31-1`, `P2 = 2^31-19` and proved for the certified cases of §5.

**Parity split (s = 0).** For `s = 0` the involution `nu : k -> -k` preserves `A` and commutes
with `S_M`, so `E_lam = E^+ (+) E^-`. Columns: one orbit of each pair `{O, -O}` and all
self-negating orbits, with `c_T^eps = c_T + eps nu c_T`. For a homogeneous `P` of degree `d`,
`<c_T^eps, P> = (1 + eps (-1)^d) <c_T, P>`, so `E^+` only sees generators of even degree and
`E^-` only those of odd degree. `dim E^eps` is obtained exactly: `K_O` for each pair, and for a
self-negating `O`, `rank_p(F-rows of parity eps on the columns of O)`; these two ranks are lower
bounds for `dim E(O)^+` and `dim E(O)^-`, which add up to `K_O`, and the program checks that the
two ranks add up to `K_O`, which forces equality. Then
`reg_lam = max(reg_lam^+, reg_lam^-)`. This halves the column count; results were checked to
coincide with the unsplit computation (`test_parity.py`, 20 shell/slice cases).
