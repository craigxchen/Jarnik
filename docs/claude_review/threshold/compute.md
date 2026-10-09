# Exact isotypic computation of shell vanishing orders (tau-type data), M = 3..9

Task: for `M = 3..9` and the shells `S_{s,t} = {k in Z^M : sum k = s, |k|_1 = t}`, compute the
maximal order `ord(c)` of a nonzero measure `c` on the shell (largest `m` with
`sum_k c_k P(k) = 0` for all polynomials `P` of degree `< m`), over ALL measures, decomposed into
`S_M`-isotypic components. Report the best ratio `ord/t`, the components attaining it, whether
anything beats the Vandermonde ratio or 2, and the trend in `t`. Problem statement:
`../tau_problem.md`.

* Scripts and raw data: `scratchpad/tau/compute/` (§8).
* Full per-shell tables: `../compute_tables_shell.md` (shells) and `../compute_tables_slice.md`
  (full slices).

Conventions.

* *Theorem / Lemma / Proposition / Corollary* means proved here, or cited with the place where
  the proof is written. *Data* means a finite exact-modular computation.
* Modular ranks give **rigorous upper bounds** on `ord` (rank mod p <= rank over Q). Lower
  bounds are rigorous only where an exact certificate or an explicit construction is cited (§5).
* `V(M) := 2 - 1/ceil(M/2)` is the Vandermonde ratio: `3/2, 3/2, 5/3, 5/3, 7/4, 7/4, 9/5` for
  `M = 3..9`.
* Partitions are written in exponent notation: `21^3 = (2,1,1,1)`, `1^M` = sign, `M` = trivial.

## 0. Verdict

**No computed shell beats the Vandermonde ratio, and none comes near 2.**

* **Coverage.** Every shell `S_{s,t}` with `0 <= s <= t` in these ranges was computed:
  `M = 3: t <= 30`; `M = 4: t <= 26`; `M = 5: t <= 24` (plus `t = 25..28` for `s <= 15`);
  `M = 6: t <= T6`; `M = 7: t <= T7`; `M = 8: t <= T8`; `M = 9: t <= T9`, plus the Vandermonde
  shell `(p,n) = (10,10)`, `t = 20`, `s = 0`.
* **Size.** That is 1002 (1002 with both primes) shells. Every shell was computed with all `S_M`
  components and two primes `P1 = 2^31-1`, `P2 = 2^31-19`, and the two primes give identical
  per-component rank profiles in every case.
* **The maximum ratio is V(M).** For every `M`, `max reg(S_{s,t})/t = V(M)` over the computed
  range, with M9_CAVEAT. It is attained at the Vandermonde shell by the sign component alone,
  with `reg = C(M,2)`; every other component there is at most `C(M,2) - 1`. There are also a
  few ties (§3.1). No case exceeds `V(M)`, so no case exceeds 2.
* **Theorem 4.4 (computer-assisted).** `tau_sh(3) = tau_sh(4) = 3/2` and `tau_sh(5) = 5/3`,
  where `tau_sh(M)` is the supremum of `ord/t` over ALL shells and all `t`. The proof combines a
  new elementary bound, Proposition 4.3, which reduces each of these `M` to finitely many shells,
  with the data.
* **Trend (Proposition 4.3).** At fixed `M`, `ord/t` decreases beyond the Vandermonde shell and
  tends to 1, because `t <= reg(S^M(p,n)) <= t + 2^M - 3` for `M >= 4` and `p, n >= 1`. In the
  data the excess `max_s reg - t` is eventually periodic for `M = 3, 4, 5`, with maximum
  `1, 3, 5` respectively, and still grows slowly for `M >= 6`. So the shell quantity is
  maximised at small `t`, near the Vandermonde shell.
* **Full slices.** These carry the invariant that the reduction actually needs. Data on
  616 (616 with both primes) slices (`M <= 9`, `t <= 24/20/18/15/13/13/12`, plus the `M = 8`
  Vandermonde slice): `reg Q^M(p,n) <= floor(phi_M(p,n))` holds in every case. This is an
  independent check of `../upper.md` Theorem 3.1 far beyond its `M <= 6, t <= 7` range.
  * `reg = C(M,2)` at the Vandermonde slices for `M = 3..8` (`upper.md` Theorem 4.4(a)).
  * Unlike shells, slices keep the ratio `V(M)` for arbitrarily large `t`. They reach it exactly
    at multiples of the Vandermonde point, as products allow.

## 1. What is computed

For a finite set `A` in a hyperplane `{sum k = s}` let
`reg(A) = min{d : every function A -> Q is the restriction of a polynomial of degree <= d}`.

* By duality (Lemma 2.1), `reg(A) = max{ord(c) : c != 0 on A}`.
* So the requested quantity for a shell is `reg(S_{s,t})`, and the ratio is `reg(S_{s,t})/t`.
* Since `k -> -k` maps `S_{s,t}` onto `S_{-s,t}`, only `s >= 0` is computed.
* Shells are written `S^M(p,n)`, with `p = (t+s)/2` (positive mass) and `n = (t-s)/2`
  (negative mass).

Two refereed notes (`../verify.md` §3.2 and `../upper.md` Remark 1.3) showed that the invariant
relevant for the reduction is the full slice `Q^M(p,n) = union_i S^M(p-i, n-i)`, not the shell,
and that `tau_sh <= tau`. The task asks for shells. §6 adds slices, computed with the same code.

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

### 2.4 Validation

* **Direct cross-check.** `sum_lam f^lam h_lam(d)` was compared, degree by degree, with the
  Hilbert function of the whole shell, computed directly from all monomials in `k_2..k_M`
  (`direct.py`, mod `P1`). The result: 132 of 132 shells agree. This covers `M = 2..5`
  with `t <= 8`, `M = 6` with `t <= 8`, `M = 7` with `t <= 7` and `M = 8` with `t <= 6`, all `s`.
  * Each `h_lam^{gen}` is at most `h_lam`, so equality of the sums forces
    `h_lam^{gen} = h_lam` (mod `P1`) for every `lam`.
  * So in these cases the generators are complete, independently of the ATY theorem.
* **Previous independent computations.** These are the direct, non-equivariant shell
  regularities of `../verify.md` (E1: `M = 3, 4, 5`) and `../upper.md` (`upper_checks/log_M*.txt`:
  `M = 3..6`). All PREV_N overlapping shells agree (`compare_prev.py`).
* **Two primes.** Every shell and slice was run with `P1` and `P2`. The full per-component rank
  profiles `h_lam(d)` agree in every case.
* **Parity split.** It was checked against the unsplit computation on 20 shells and slices
  (`test_parity.py`).
* **Sizes.** `sum_lam f^lam N_lam` equals the shell cardinality (closed formula) in every run.
* **Columns.** The final rank equals `N_lam` in every run, so the `c_T` are independent.

## 3. Results for shells

### 3.1 Maxima per M

| M | V(M) | shells computed (t range) | two primes agree | max reg/t | attained at (s,t): components | any > V(M) | any > 2 |
|---|---|---|---|---|---|---|---|
| 3 | 3/2 | 255 (t = 1..30) | 255/255 | 3/2 | (0,2): 1^3 | no | no |
| 4 | 3/2 | 195 (t = 1..26) | 195/195 | 3/2 | (0,2): 21^2; (2,4): 1^4; (0,6): 1^4 | no | no |
| 5 | 5/3 | 200 (t = 1..28) | 200/200 | 5/3 | (0,6): 1^5 | no | no |
| 6 | 5/3 | 120 (t = 1..20) | 120/120 | 5/3 | (0,6): 21^4; (3,9): 1^6; (0,12): 1^6 | no | no |
| 7 | 7/4 | 89 (t = 1..17) | 89/89 | 7/4 | (0,12): 1^7 | no | no |
| 8 | 7/4 | 80 (t = 1..16) | 80/80 | 7/4 | (0,12): 21^6; (4,16): 1^8 | no | no |
| 9 | 9/5 | 63 (t = 1..14) | 63/63 | 7/4 | (0,12): 31^6 | no | no |

(MAXIMA_NOTES pending)

### 3.2 The Vandermonde shells, component by component

* M = 3, (p,n) = (1,1), t = 2, |S| = 6, reg = 3, C(3,2) = 3 (two primes agree); components: 3: 0, 21: 2, 1^3: 3
* M = 4, (p,n) = (3,1), t = 4, |S| = 40, reg = 6, C(4,2) = 6 (two primes agree); components: 4: 3, 31: 4, 2^2: 4, 21^2: 5, 1^4: 6
* M = 5, (p,n) = (3,3), t = 6, |S| = 340, reg = 10, C(5,2) = 10 (two primes agree); components: 5: 6, 41: 7, 32: 8, 31^2: 8, 2^21: 8, 21^3: 9, 1^5: 10
* M = 6, (p,n) = (6,3), t = 9, |S| = 4340, reg = 15, C(6,2) = 15 (two primes agree); components: 6: 10, 51: 10, 42: 11, 41^2: 12, 3^2: 11, 321: 12, 31^3: 13, 2^3: 12, 2^21^2: 13, 21^4: 14, 1^6: 15
* M = 7, (p,n) = (6,6), t = 12, |S| = 65226, reg = 21, C(7,2) = 21 (two primes agree); components: 7: 15, 61: 15, 52: 16, 51^2: 16, 43: 16, 421: 17, 41^3: 18, 3^21: 17, 32^2: 17, 321^2: 18, 31^4: 19, 2^31: 18, 2^21^3: 19, 21^5: 20, 1^7: 21
* M = 8, (p,n) = (10,6), t = 16, |S| = 1264032, reg = 28, C(8,2) = 28 (two primes agree); components: 8: 19, 71: 20, 62: 21, 61^2: 21, 53: 21, 521: 22, 51^3: 23, 4^2: 21, 431: 22, 42^2: 23, 421^2: 24, 41^4: 25, 3^22: 22, 3^21^2: 23, 32^21: 24, 321^3: 25, 31^5: 26, 2^4: 24, 2^31^2: 25, 2^21^4: 26, 21^6: 27, 1^8: 28

(VAND_NOTES pending)

### 3.3 Trend in t: max over s of reg/t

| M | t : max_s reg (ratio; `*` = equals V(M)) |
|---|---|
| 3 | 1:1 (1.000), 2:3 (1.500*), 3:3 (1.000), 4:5 (1.250), 5:5 (1.000), 6:7 (1.167), 7:7 (1.000), 8:9 (1.125), 9:9 (1.000), 10:11 (1.100), 11:11 (1.000), 12:13 (1.083), 13:13 (1.000), 14:15 (1.071), 15:15 (1.000), 16:17 (1.062), 17:17 (1.000), 18:19 (1.056), 19:19 (1.000), 20:21 (1.050), 21:21 (1.000), 22:23 (1.045), 23:23 (1.000), 24:25 (1.042), 25:25 (1.000), 26:27 (1.038), 27:27 (1.000), 28:29 (1.036), 29:29 (1.000), 30:31 (1.033) |
| 4 | 1:1 (1.000), 2:3 (1.500*), 3:4 (1.333), 4:6 (1.500*), 5:6 (1.200), 6:9 (1.500*), 7:9 (1.286), 8:10 (1.250), 9:11 (1.222), 10:13 (1.300), 11:13 (1.182), 12:14 (1.167), 13:15 (1.154), 14:17 (1.214), 15:17 (1.133), 16:18 (1.125), 17:19 (1.118), 18:21 (1.167), 19:21 (1.105), 20:22 (1.100), 21:23 (1.095), 22:25 (1.136), 23:25 (1.087), 24:26 (1.083), 25:27 (1.080), 26:29 (1.115) |
| 5 | 1:1 (1.000), 2:3 (1.500), 3:4 (1.333), 4:6 (1.500), 5:7 (1.400), 6:10 (1.667*), 7:10 (1.429), 8:12 (1.500), 9:14 (1.556), 10:14 (1.400), 11:15 (1.364), 12:17 (1.417), 13:17 (1.308), 14:19 (1.357), 15:19 (1.267), 16:20 (1.250), 17:21 (1.235), 18:23 (1.278), 19:23 (1.211), 20:24 (1.200), 21:25 (1.190), 22:27 (1.227), 23:27 (1.174), 24:28 (1.167), 25:29 (1.160), 26:31 (1.192), 27:31 (1.148), 28:32 (1.143) |
| 6 | 1:1 (1.000), 2:3 (1.500), 3:4 (1.333), 4:6 (1.500), 5:7 (1.400), 6:10 (1.667*), 7:11 (1.571), 8:12 (1.500), 9:15 (1.667*), 10:15 (1.500), 11:17 (1.545), 12:20 (1.667*), 13:20 (1.538), 14:21 (1.500), 15:22 (1.467), 16:24 (1.500), 17:24 (1.412), 18:26 (1.444), 19:27 (1.421), 20:28 (1.400) |
| 7 | 1:1 (1.000), 2:3 (1.500), 3:4 (1.333), 4:6 (1.500), 5:7 (1.400), 6:10 (1.667), 7:11 (1.571), 8:13 (1.625), 9:15 (1.667), 10:16 (1.600), 11:17 (1.545), 12:21 (1.750*), 13:21 (1.615), 14:23 (1.643), 15:24 (1.600), 16:27 (1.688), 17:27 (1.588) |
| 8 | 1:1 (1.000), 2:3 (1.500), 3:4 (1.333), 4:6 (1.500), 5:7 (1.400), 6:10 (1.667), 7:11 (1.571), 8:13 (1.625), 9:15 (1.667), 10:16 (1.600), 11:18 (1.636), 12:21 (1.750*), 13:22 (1.692), 14:23 (1.643), 15:24 (1.600), 16:28 (1.750*) |
| 9 | 1:1 (1.000), 2:3 (1.500), 3:4 (1.333), 4:6 (1.500), 5:7 (1.400), 6:10 (1.667), 7:11 (1.571), 8:13 (1.625), 9:15 (1.667), 10:16 (1.600), 11:18 (1.636), 12:21 (1.750), 13:22 (1.692), 14:24 (1.714) |

(For `M = 5` and `t = 25..28` the maximum is over the computed `s <= 15` only. By Corollary 4.4,
the omitted shells cannot exceed `5/3`.)

The excess `max_s reg(S_{s,t}) - t`, for `t = 1, 2, ...`:

| M | max_s reg(S_{s,t}) - t for t = 1, 2, ... |
|---|---|
| 3 | 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 |
| 4 | 0 1 1 2 1 3 2 2 2 3 2 2 2 3 2 2 2 3 2 2 2 3 2 2 2 3 |
| 5 | 0 1 1 2 2 4 3 4 5 4 4 5 4 5 4 4 4 5 4 4 4 5 4 4 4 5 4 4 |
| 6 | 0 1 1 2 2 4 4 4 6 5 6 8 7 7 7 8 7 8 8 8 |
| 7 | 0 1 1 2 2 4 4 5 6 6 6 9 8 9 9 11 10 |
| 8 | 0 1 1 2 2 4 4 5 6 6 7 9 9 9 9 12 |
| 9 | 0 1 1 2 2 4 4 5 6 6 7 9 9 10 |

(TREND_NOTES pending)

### 3.4 Which components attain the maximum

(COMPONENT_NOTES pending)

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

## 5. Exact certificates (lower bounds)

(CERT_SECTION pending)

## 6. Full slices Q^M(p,n) (the invariant the reduction needs)

(SLICE_SECTION pending)

## 7. Consequences for the problem

(CONSEQ_SECTION pending)

## 8. Files

(FILES_SECTION pending)
