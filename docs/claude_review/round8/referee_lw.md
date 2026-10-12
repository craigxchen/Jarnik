# Referee report on `round8/lw.md` (Lang–Waldschmidt and the growth rate)

Referee: adversarial, round 8. Scripts: `round8/referee_lw_checks/`. Each runs from that directory with
`python3 <script>`, and its output is stored next to it as `*_out.txt`. The repository and `lw_checks/` were not
modified.

## 0. Verdict

**The result survives as a barrier, with two major corrections, both about the scope of the statement.**
The write-up claims no growth improvement, so there is no growth claim to test against the round-5 fakes. Its
exclusion claims (Theorem L, Theorem 3.1) use the Lang–Waldschmidt hypothesis, which is exactly the input the
construct/lower/sharp fakes violate. That is consistent with the barrier map.

I re-derived every proof the conclusions rest on:
* the valuation lemma (Lemma 2.2);
* the skew double core uncertainty principle (Theorem 5.1);
* the union bound of Theorem 6.1;
* the cited but unrefereed Liouville margin `Phi(c) >= 2(M-1) + n` (`exact.md` Theorem 4.1(iii)), which
  Theorem 6.1 needs;
* `flipped.md` Theorem L.

None of these proofs contains a mathematical error. Independent exact computations confirm:
* Lemma 2.2 on a wider pool of representations, including units, `1+i`, `3` and a split prime outside `J`;
* the exact minimum `s = M/2 + 1` at `M = 8, 12`, now proved rigorously by a rank certificate over every
  candidate zero set;
* the census of light characters;
* the random two-basis examples;
* the arithmetic of the union bound.

What does not survive as stated:

1. **(major, statement of the conjecture)** The displayed (LW) is the *max form*
   `C^n B^(-(n-1+eps)) prod A_k^(-(1+eps))`, and §1.2 calls it "the Lang–Waldschmidt conjecture" while citing
   Waldschmidt 2004, Conjecture 2.5. That conjecture, and Lang 1978 pp. 212–217, are stated in the *product form*
   `C(eps)^n B / (|b_1 ... b_n| A_1 ... A_n)^(1+eps)` (for rational integers in Conj. 2.5). The two forms are
   incomparable at the same `eps` when `n` is unbounded.
   * Theorem 6.1, and with it outcome (B) "LW(eps) at fixed `eps` gives no growth improvement", is proved only
     for the max form, and only for lists that contain `i`.
   * The standard (product-form) conjecture is not covered, and neither are homogeneous instances (no `Log i`).
     The union bound of the proof diverges in both cases (§3.4).
   * The write-up lists these as caveats, but its §0 bold verdict and the task summary state (B) for "LW(eps)".
2. **(major, scope)** "No growth improvement follows from Lang–Waldschmidt through single linear forms" holds
   *relative to the angle-model input class*. That class is the exponent profile, all monomial (Lemma 0.1)
   gaps and Haar-generic angle properties. Point-level residues come in only through the unrefereed
   `exact.md` Theorem 5.1. The skew double core fakes are not known to satisfy the prime-level parity law of
   `sharp.md` (P3-char / P3-all). They do not model additive data either: for an actual cluster,
   `Log((z_x + z_y)/(2 z_u)) = O(Delta)` is a cluster-controlled form that is not a character form. So the
   refutation of outcome (A) is: no argument that adds one LW form to these inputs works. It is not a statement
   about LW plus every valid input for actual circles.

Minor issues (§4) include:
* an arithmetic slip in the proof of Corollary 3.3 (the claimed error bound `W/(2 log M)` is false; the
  corollary survives with a larger `M_1(eps)`);
* census counts at `M = 8` that miss positive-dimensional families at `s = 8`;
* an imprecise "the product form is the stronger statement";
* a constant `4/eps` that can be lowered to `2/eps + o(1)`, though the `1/eps` itself is intrinsic to the
  union-bound method;
* a diagonal argument missing from the "`eps`-uniform barrier" reading of Proposition 7.2;
* the restriction of Lemma 2.2 to `Q(i)`;
* the scope of the Theorem L consequence: Paley orders only, since the fakes at other `M` are row-deleted
  sub-profiles.

The corrected statement is in §5.

## 1. What was checked, and how

| item | method | result |
|---|---|---|
| Lemma 2.1 (heights in `Q(i)`) | re-derived; exact check on 26,244 numbers (8 places incl. `1+i`, `3`, `4+-i`; 4 units) — `r4_valuation.py` | correct |
| Lemma 2.2 (valuation lemma) | re-derived; exhaustive representations `K <= 2` (8 places, `\|b\| <= 3`), `K = 3` (6 places, `\|b\| <= 2`), 8 targets — `r4_valuation.py` | correct; tight (`v = (0,1)`, `(1,-1)` attain `w(supp v)` exactly) |
| Cor. 2.3, Def. 2.4, Thm 3.1, Thm 3.2 | re-derived | correct (constant `9` needs `eps <= 1.34`) |
| Cor. 3.3 | re-derived; numbers in `r5_union_bound.py` | statement OK for `M >= M_1(eps, C_0)`; proof's inequality false (§4.1) |
| §4: `flipped.md` Theorem L, Props 2.1–2.2, (L.1)–(L.2), `LW => (H_(1+eps))` | re-derived; `check_two_log_identity.py` re-run (`referee_lw_checks/rerun_flipped/`) | correct; 2,600 identities, output byte-identical |
| Theorem L kills construct/lower/sharp | profiles re-read in `construct.md` §2, `lower.md` §3.1, `sharp.md` §4: `F/W ~ 1/(b+1)` or `1/b`, every label has unflipped columns | correct at Paley orders (§4.8) |
| Theorem 5.1 (`s(c) >= M/2`) | re-derived both cases; independent code (`r1_dcore.py`): identities (5.1), `Sigma_a = R(a) mod q`, random `c`, adversarial Case-2 vectors (`p_k = 0 mod q`, `k <= d`) for `q = 7..59` | correct; Case-2 vectors give `s >= q+3` as claimed |
| exact minimum `M/2+1` (`M = 8, 12`) | **rigorous**: every set of `#types - M/2` types plus `1` has rank `M` mod `2^31-1` (1,001 and 74,613 sets) — `r2_census.py` | confirmed; no `c` with `s <= M/2` |
| census (`14, 28, 140, 525`; `22, 0, 220, 1771`), lightest = the `2q` types | independent method (signed maximal minors, exact verification) — `r2_census.py` | confirmed, but see §4.2 on positive-dimensional families |
| `exact.md` Thm 4.1(iii), `Phi(c) >= 2(M-1)+n` | re-derived line by line (amplitude `>= 2`; ternary: `c_0 = 0`, `c_0 != 0` with `n = 4` and `n >= 6`) | correct |
| `exact.md` Lemma 5.0, `outside.md` Thm A.1 | re-derived | correct as used |
| Theorem 6.1 | re-derived (a), (b), (c); arithmetic in `r5_union_bound.py` | correct **for the max form with `i` listed** (§3) |
| §7: pair supports, Lemma 7.1, random 2-basis data | re-derived; `r6_random_bases.py` recomputes the reported light `c` and tests the Lemma 7.1 input for `q = 19..47` | correct (`s = 2, 4, 4, 4` reproduced) |
| conjecture statement and references | web search (Waldschmidt 2004 Conj. 2.5 as quoted in arXiv:2205.07841, Conj. 5.1; Stewart's lecture notes) | §1.2 misattributes the max form (§2.1) |

## 2. Major issues

### 2.1 (major) The conjecture is not stated in its standard form; (B) is proved only for a variant

`lw.md` §1.2 says: "The Lang–Waldschmidt conjecture asserts: for every `eps > 0` there is `C(eps) > 0` with
`LW(eps, C(eps))`". Here (LW) is the max form

```text
|Lambda| >= C^n B^(-(n-1+eps)) prod_k A(alpha_k)^(-(1+eps)).                                   (max form)
```

Its sources are Lang 1978, Waldschmidt 2000 Ch. 1, and Waldschmidt, *Open Diophantine problems*, Mosc. Math. J.
4 (2004), Conjecture 2.5. Conjecture 2.5 is the product form, stated for positive rational integers
`a_1, ..., a_n`. As quoted in arXiv:2205.07841 (Conjecture 5.1, which attributes it to Lang 1978,
pp. 212–217, and to Conj. 2.5):

```text
|a_1^(b_1) ... a_n^(b_n) - 1|  >=  C(eps) max_j |b_j| / |b_1 ... b_n a_1 ... a_n|^(1+eps).   (product form)
```

The algebraic-number version uses `A_j = max(H(alpha_j), e^(|log alpha_j|))` and the constant `C(eps)^(nD)`. This
is also the form that `flipped.md` §4 quotes. I could not open the 2000 book or the 2004 paper; both are behind a
network block. The form is confirmed by two independent secondary sources. The write-up itself admits that the
2000 book has the product form.

The two forms are **incomparable** at the same `eps` when `n` is unbounded:
* For balanced coefficients the max form is stronger: `B^(-(n-1+eps))` against `B^(-(n-1)-n eps)`.
* For unbalanced coefficients the product form is stronger. For `b = (B, 1, ..., 1)` it gives `B^(-eps)`, against
  `B^(-(n-1+eps))`.

So the sentence "the product form is the stronger statement" is not right as written. What is true is: product
at `eps` implies max at `n eps`, which is useful for bounded `n` only.

Consequences:
* **Exclusion results.** Theorems 3.1 and L and the corollaries use a bounded number `K + 1` of logarithms
  with coefficients `<= poly(M)`. They hold under either form, with the error term
  `(K + (K+1) eps) log beta_K` in place of `(K + eps) log beta_K`. No change in substance.
* **Single-form immunity.** Theorem 3.2 and Corollary 5.2 hold for both forms. The product RHS is
  `<= B^(-eps) prod A^(-(1+eps))`, because `prod |b_k| >= B`.
* **Simultaneous fakes (Theorem 6.1, outcome (B)).** These are proved only for the max form, and only for
  lists containing `i` (Definition 1.1). For the product form, the sum over coefficient vectors in part (c)
  becomes `sum_c B_c prod |c_k|^(-(1+eps))`. When the `Log i` coefficient is small, this has a divergent
  one-dimensional tail (`c = (B, 1, ..., 1)` contributes `B^(-eps)`). The same happens for homogeneous
  instances (`sum_B B^(K-1) B^(-(K-1+eps))`). Closing either gap needs the density gain
  `1/(Delta ||lambda||_inf)` together with a lower bound `||lambda||_inf >~ B_c / poly`. The proof does not
  supply this bound: cancellations among overlapping `e_k` can make `v` small while `B_c` is large.

The write-up records all of this in its caveats. But its §0 bold sentence ("For each fixed `eps`, `LW(eps)` ...
gives no growth improvement at all") and summary items (4) and (B) attribute the result to "LW(eps)". They must
say "`LW` in the max form, for linear forms whose list of numbers includes `i`". Whether the standard product-form
conjecture at a fixed `eps` admits simultaneous fakes at `W = O(M log M)` is **not proved**. I expect it to be
true, since the divergences are union-bound artefacts, but it is open as far as this write-up goes.

### 2.2 (major, scope) The barrier is relative to the angle-model class

Outcome (A) is "a conditional growth theorem obtained by finding, in every profile with `W = O(M log M)`, one
linear form that LW forbids". It is refuted only for profiles that survive all the *other* inputs the argument
may use. The witness is the skew double core, and what is proved about it is:
* parts 1–2 of Theorem 6.1: the class `A_0` (exponent profile, all Lemma 0.1 gaps, Haar-generic angles),
  proved self-containedly. I re-derived the Liouville margin it borrows from `exact.md` Thm 4.1(iii);
* the class `A_1(Q_0)` of point-level residues, but only through `exact.md` Theorem 5.1, which is unrefereed
  (the write-up says so).

Not covered:
* **Prime-level parity.** The parity law of `sharp.md` is satisfied by actual circles; `referee_construct.md`
  showed that the construct fake violates it. Nobody has checked P3-char or P3-all data for the double core.
  If P3-all forced `W >> M log M` for the double core, a "single LW form + P3" argument would need another
  immune witness. P3-all is listed as open in `round5/README.md`.
* **Additive data.** For actual clusters, `beta = (z_x + z_y)/(2 z_u)` satisfies `|Log beta| = O(Delta)`
  (the argument is `(delta_x + delta_y)/2 - delta_u`, and `|beta| = cos((theta_x - theta_y)/2)`). It is a
  cluster-controlled linear form (`n = 1`) that is not a character form. The angle models have no `z_x`, so they
  say nothing about such forms. With `n = 1`, LW is Liouville, a theorem. Its height is about `e^W` unless the
  gcd is large, and controlling gcds is the "single-place auxiliary polynomial" route already barred in
  `threshold/README.md`. So this is not a hole in the barrier map, but the verdict's "through single linear
  forms" must mean "single linear forms in the logarithms of the primes of `N` and of `i` (character forms)".
  That is the task's scope, and §1.1 of `lw.md` says as much.

The write-up phrases (B) carefully in Remark 6.2 ("no argument whose inputs are ..."). The §0 verdict, §8 (A) and
the summary should carry the same qualifier.

## 3. Re-derivations of the load-bearing steps

### 3.1 Lemma 2.1 and the valuation lemma (Lemma 2.2): correct

* **Lemma 2.1(1).** `2h = sum_v log^+ ||alpha||_v + 2 log^+|alpha|`. The finite part is `P-`, and
  `2 log|alpha| = P+ - P-` by the product formula, so `2h = max(P+, P-)`.
* **Lemma 2.1(3).** For `alpha = (A+Bi)/D`, the primitive minimal polynomial is
  `(D^2, -2AD, A^2+B^2)/g`. Its Mahler measure is `max(D^2, A^2+B^2)/g = e^(2h)`, and `2|A|D <= A^2 + D^2`
  gives `H <= 2 e^(2h)`.
* **Lemma 2.1(4).** `ord_(pi_j) U = e_j = -ord_(conj pi_j) U`, so `P+ = P- = ||e||_w`.
* **Lemma 2.2.**
  * Subtracting the two valuation equations at `pi_j`, `conj pi_j` forces some irrational `alpha_(k_1)` with
    `x != y`.
  * If `|x| + |y| = 1`, adding the equations forces a second index with a nonzero contribution. Its
    coefficient `b_(k_1)` is nonzero, because zero-coefficient terms were dropped.
  * Each such index contributes `>= 1/2` to `tau_j`. Rational numbers contribute `|z|/2`, from
    `max(P+,P-) >= (P+ + P-)/2`.
  * Hence `sum_k log H(alpha_k) >= sum_k 2h(alpha_k) >= w(supp v)`.
* Units are irrelevant: only `H >= e^(2h)` is used, and `h` is unit-invariant.

*Exact check* (`r4_valuation.py`):
* Lemma 2.1 holds on all 26,244 pool elements (4 units, exponents in `[-1,1]` at 8 places).
* For the targets `v in {(1,0),(0,1),(1,1),(1,-1),(2,1),(2,-1),(2,0),(3,-2)}` over `J = {5, 13}`, no
  representation with `K <= 2` numbers (`|b| <= 3`) or `K = 3` (6 places, `|b| <= 2`) beats `w(supp v)`.
* The bound is attained exactly by single numbers: `u_13`, with `H = 13`, and `u_5/u_13`, with `H = 65`.

Scope remark: the lemma is stated for `alpha_k in Q(i)^*`. The conjecture is stated for all algebraic numbers
(constant `C(eps)^(nD)`). Representations by higher-degree numbers, for example Kummer roots `u_j^(1/m)`, are not
covered. With naive heights, `log H(alpha) >= D h(alpha) - (1/2) log(D+1)`, and the valuation count at the places
above `pi_j` should go through. But this is not written (minor, §4.7).

### 3.2 Theorem 5.1 (skew double core, `s(c) >= M/2`): correct, and the minimum `M/2 + 1` is now rigorous

Notation: `inf` is the border. Then

```text
<c, K(.,a)> = (c_a - c_inf) + Sigma_a,     <c, s2_a> = -(c_a - c_inf) + Sigma_a,     Sigma_a = sum_y chi(a-y) c_y.
```

I checked these against the matrix `K = I + S` built independently: rows `inf` and `a` of `s2_a` are negated
relative to `K(.,a)`.

* **Case 1.** `R != 0`, where `R(X) = sum_y c_y (X - y)^d` has degree `<= d`, so `|Z_u cap Z_v| <= d` and
  `s >= 2q - q - d = (q+1)/2`.
* **Case 2.** `R = 0`. All `binom(d,k)` are units mod `q`, so `p_k = 0` for `k <= d`. With
  `sum_y y^k = 0 mod q` for `k <= q-2` this gives `sum_y f_y y^k = 0`, `k <= d`, for `f = c - c_inf`. If
  `f != 0 mod q`, the Vandermonde argument gives `|supp f| >= d+2`, hence `s >= q+3`. If `f = 0 mod q`,
  zero-sum gives `(q+1) c_inf = 0 mod q`, so `c = 0 mod q`, contradicting primitivity.

*Checks:*
* `r1_dcore.py`, for `q = 7, 11, 19, 23, 31, 43, 47, 59`:
  * the identities and `Sigma_a = R(a) mod q` on random `c`;
  * `s >= M/2` on random sparse `c`;
  * adversarial Case-2 vectors (`c_y = g(y) + q * noise` with `deg g <= (q-3)/2`): `s >= q+3`, as the proof
    says;
  * column characters: `s = M/2 + 1`, values `{M/2, (M-4)/2, 2}`, so `K <= 3` and `beta_K = M/2`;
  * label pairs: `M/2 + 3`; row pairs: `M`.
* `r2_census.py` (R2): every set of `2q - M/2` types together with `1` has rank `M` mod `2^31 - 1`. Since rank
  mod `p` is at most rank over `Q`, no nonzero zero-sum `c` has `s(c) <= M/2`, for `M = 8` (1,001 sets) and
  `M = 12` (74,613 sets). With the columns this proves **`min s = M/2 + 1` exactly**. The write-up's own
  minimum used float SVD to *find* hyperplanes; the bound is now certified independently of floating point.

### 3.3 `exact.md` Theorem 4.1(iii) (needed by Theorem 6.1(a)): correct

`Theorem 6.1` uses `Phi(c) = sum_j |<c, s_j>| >= 2(M-1) + n` from the unrefereed `exact.md`. I re-derived it:
* `u = K^T c`, `v = -K c~` (with `c~` equal to `c` with `c_inf` negated), and `u_inf = v_inf = 0`.
* `||u||_2^2 = ||v||_2^2 = M ||c||_2^2`.
* Defect identity `n ||u||_1 = ||u||_2^2 + sum |u_a|(n - |u_a|)`.
* Amplitude `>= 2`: `c = K u/M`, so `||u||_1 >= MA`, and `2MA >= 2M - 2 + n`.
* Ternary `c`: all `u_a, v_a` are even.
  * `c_inf = 0`: `|u_a - v_a| = 2` on the support forces a middle value unless `n = 2`.
  * `c_inf != 0`, `n >= 6`: `|u_a - v_a| = 4` on the `n/2` rows with `c_a = -c_inf`.
  * `n = 4`: every empty row gives exactly one `|.| = 2`, because `Sigma_y` is odd. So `D >= 4(M-4)`.

All steps check. Part B of `lw_checks/check_dcore_families.py` spot-checks it for `q <= 47`.

### 3.4 Theorem 6.1: correct for its stated scope; where the `1/eps` comes from

**Parts (a) and (b).** These are `exact.md` Theorem 5.1 for the class `A_0` with the primes in a window of
log-width `1/(8q)`. I re-derived the following:
* the window-count probability `4 eps(n + pi/Delta)` (Lemma 5.0);
* `2M P'^(-1/4) <= C_0/20`;
* the sum `4.3 sum_(n>=2) (n + pi/C_0)(C_0/20)^n`. It equals `0.0588` at `C_0 = 1` and is smaller for smaller
  `C_0` (`r5_union_bound.py`). It is `< 1/16`;
* the deterministic handling of rational non-integral characters through the injective winding `k = e_1`;
* `1.33(Pi - 1) <= 1.33/(200M)`.

**Part (c).** The non-character part is `(4/pi) L` per tuple, from Haar-uniformity of `<v, phi>` for `v` not in
`span L` (`outside.md` Thm A.1(ii) proof). The character part is `4 M L/Delta` per tuple: density
`<= 1/(Delta ||lambda||_inf)`, `||lambda||_inf >= 2/M`, at most `2 ||lambda||_1 Delta/pi + 2` windows. Further:
* `h' >= w(supp v) >= W/4 + x_0`, because every `j in supp v` lies in some `supp e_k`.
* `sum_(c : B_c = B) B^(-(K+eps)) <= K 2^K B^(-1-eps)`.
* `Z(s) <= 4r/(P'^s - 1)` once `2r <= P'^s - 1`.

All correct. Torsion entries, zero coefficients and repeated `alpha`'s are all inside the enumeration of `K`-tuples
of nonzero `e_k`, as the proof says. Numerically, the character part is `< 1/8` at `P' = P_1` already for
`M = 8`, for every `eps` tested (`r5_union_bound.py`).

**On the constant.** The `1/eps` is real for *any* union bound over representations. The number of ways to
group a support of size `m = M/2 + 1` into numbers `alpha_k` is at least the Bell number
`B_m = e^((1-o(1)) m log m)` (`log B_51 = 111.7` against `(M/2) log M = 230` at `M = 100`; `r5_union_bound.py`).
Only `e^(-eps W/4) = e^(-eps (M/2) log P)` is available against it, which needs `eps log P >~ log M`. The proof's
split `(eps/2, eps/2)` gives `4/eps`. Using `theta eps h' + (1-theta) eps W/4` with `theta -> 1` gives
`2/eps + o(1)`. Numerically the threshold falls to `9.5, 23.7, 47.4` times `M log M` at `M = 10^7`, for
`eps = 1/4, 1/10, 1/20`. This changes constants only. The write-up's count `exp((W/4) log r/log P)` counts
ordered tuples, but it has the right order `e^((M/2) log M)`.

**Scope.** The proof covers exactly Definition 1.1: max form, `alpha_n = i` listed. It does not cover the product
form or homogeneous instances (§2.1).

### 3.5 Theorem L and its use against the refereed fakes: correct (at Paley orders)

I re-derived Lemma 1.1, (1.2)–(1.4), Props 2.1–2.2 (including the nonvanishing via the axis gap, and `|b_3|`),
the averaging in (L.2) (`P_x N_x <= (M-1)^2/4`), and the passage from LW to `(H_(1+eps))`: with three
logarithms, `A(alpha) <= e^pi H(alpha)` for `|alpha| = 1` and `A(i) = e^(pi/2)`.

The profiles are as follows:
* `construct.md` §2: `b` unflipped copies plus one flipped copy per row, so `F/W -> 1/(b+1)`.
* `lower.md` §3.1: capacity two.
* `sharp.md` §4: the flip column `(x, 0)` is one of the `b` copies of label `x`, so `F/W ~ 1/b`.

All are full Paley cores at `M = q+1`, and every label keeps unflipped columns (`|F_j| = 0 < M/4`). So
Theorem L applies with `F <= (1/(2(1+eps)) - delta) W` for small `eps`. The single-label character `H(., a)` also
works through Theorem 3.1: its support is block `a` plus all `M` flipped columns, `sigma ~ 1/(b+1)`, values
`{M, M-2, 2}`, `K = 3`. The re-run of `check_two_log_identity.py` (`referee_lw_checks/rerun_flipped/`) passes
2,600 identities with byte-identical output.

## 4. Minor issues

### 4.1 Corollary 3.3: the displayed error bound is false; the corollary survives with a larger `M_1`

The proof says the error terms are `<= W/(2 log M)`. But under its own hypotheses the first term is already
`(W/(log M)^3)(log M)^2 = W/log M`. The full error is

```text
E = W/log M (1 + 3/log M + (9 + log(1/C(eps)))/(log M)^2) + (log M)^2 + log C_0,
```

which is `1.43, 1.26, 1.17` times `W/log M` at `M = 10^4, 10^6, 10^9` (`r5_union_bound.py`, with
`log(1/C) = 0`). The conclusion `w(supp v) >= W/(4(1+eps)) - W/log M` needs `E <= (1+eps) W/log M`. That holds
once `3/log M + (9 + log(1/C(eps)))/(log M)^2 + (log M)^3/W <= eps`, i.e. for `M >= M_1(eps, C(eps), C_0)` with
`log M >~ 3/eps`. Replace "`<= W/(2 log M)`" by "`<= (1+eps) W/log M` for `M >= M_1`".

### 4.2 The census counts only zero-set-determined characters; at `M = 8` it misses infinite families

The census of §5.2 is restricted to "primitive, up to sign, determined by their zero sets". That restriction
matters:
* At `M = 8`, 154 of the `(M-2)`-subsets are degenerate (`rank(Z u {1}) = 6`). Each spans a 2-dimensional
  kernel whose generic member has `s = 8`.
* So there are **infinitely many** primitive characters with `s = 8 = M`. The figure "525" counts only the
  isolated ones.
* At `M = 12` the degenerate spans have generic `s = 11` (330 spans) or `12` (11,011 spans). So the counts
  `22, 0, 220, 1771` for `s = 7..10` are complete (`r3b_families.py`; minimiser rank confirmed over `Q`).

This does not affect Theorem 5.1 or the minimum. It does affect the §8 heuristic "the next level has `O(q^2)`
members": the count of characters at support excess `x` is infinite from `s ~ M` on. Any version of (N) must
weight by coefficient size, which (N) already does through `||lambda||_inf <= T`.

### 4.3 "The product form is the stronger statement" (§1.2)

This is false as a statement about a fixed `eps`; see §2.1. Product at `eps` implies max at `n eps`. Neither
implies the other at the same `eps` for unbounded `n`.

### 4.4 Proposition 7.2 and the "`eps`-uniform barrier"

Proposition 7.2 gives fakes for every `eps in (0,1]` once `M >= M_3(eps, mu, C_0)`. So a fake of order `M` is
certified only for `eps >= eps_M`, with `eps_M -> 0`. To conclude that no argument from the *full* conjecture,
with ineffective `C(eps)`, proves growth, one more step is needed. Each fake must satisfy `LW(eps, C_M(eps))` for
**every** `eps > 0` with some `C_M(eps) > 0`. Then `C(eps) := min(1, min_(M : eps_M > eps) C_M(eps)) > 0` works
for the whole family. That step is true: in the union bound of Theorem 6.1(c) the factor `C^(K+1)` makes the
failure probability tend to 0 as `C -> 0`, for fixed `M` and `eps`, and one sums over `eps = 1/m` with
`C_m -> 0` fast. But it is not stated. The same remark applies to reading Theorem 6.1 as a statement about
`C = C(eps)`.

### 4.5 The lead's multi-basis claim: false for `t >= 3` by pairs; for `t = 2` only small-`M` evidence

The body of §7 says "false for `t >= 3`". For `t = 2`, the claimed `(1 - 1/2) W = W/2` is *consistent* with
pairs. It is contradicted only by the `M = 8, 12` random examples:
* `s = 2, 4, 4` of `14, 22, 22` types, reproduced by `r6_random_bases.py`;
* "small" in the lead's statement could absorb these.

Summary item (5) ("false: pair characters always have support about `W/2`") should say "false for `t >= 3`; for
`t = 2` refuted at `M = 8, 12` by coincidences, asymptotically undecided".

### 4.6 Constants and normalisations

* The constant `9` in Definition 2.4 and (3.1) needs `(1+eps)(pi + log 2) <= 9`, i.e. `eps <= 1.34`. State
  `eps <= 1`.
* "Any true instance of `LW(eps, C)` has `C <= 1`" uses dummy terms `1 * Log 1`. Formulations that require the
  `alpha_k` to be multiplicatively independent, or `Lambda` to be a minimal relation, exclude these, and then
  `C(eps) > 1` is allowed. Theorem 6.1 then needs `LW_J(eps, C)` for `C > 1`. The union bound absorbs a factor
  `C^(K+1)` if `Z(eps/2)` is made smaller, and torsion entries are then excluded by the same formulation. This is
  harmless but should be said.

### 4.7 Lemma 2.2 is proved for `alpha_k in Q(i)^*` only

See §3.1. The conjecture covers all algebraic numbers. Single-form immunity (Theorem 3.2, Cor. 5.2) is proved
only against `Q(i)`-representations. For a barrier against LW as stated for `Qbar`, one needs the analogue
`sum_k log H(alpha_k) >= w(supp v) - O(sum_k log D_k)` for representations in any number field. It is plausible
by the same valuation count at the places above `pi_j`, but not written.

### 4.8 Theorem L's consequence holds at Paley orders

Theorem L (`flipped.md` §0) needs the full set of rows of a Hadamard matrix. The construct/lower/sharp fakes are
full cores only at `M = q + 1`. At other `M` they come from monotonicity (`lower.md` Lemma 5.2: keep `M` rows of
a larger core), and these row-deleted sub-profiles are not configurations "as in Section 0". So "under LW none of
them is a genuine circle for large `M`" is proved along the Paley orders, which suffices for the growth question.

For a bounded number `m` of deleted rows, Theorem 3.1 still applies. Take a combination of `m + 1` columns
vanishing on the deleted rows: it has `K <= m + 1 + 2^(m+2)` distinct values and
`log beta_K <= (m/2) log m + log M` (Hadamard). For gaps `m -> infinity` this is not checked.

### 4.9 Smaller points

* §3 table, "fully flipped core ... excludes iff `F < W/2`": only "if" is proved. "Only if" would need every
  other character of that profile to be spread.
* "Plotkin-optimal profile" and "best known fake" rest on `exact.md`, which is unrefereed. The optimality is
  for the point-level classes only (§2.2).
* Lemma numbering: `check_valuation.py` calls the valuation lemma "Lemma 2.3"; the text calls it 2.2.
* The section numbers of Lang 1978 and Waldschmidt 2000 are unverified (the write-up says so). Secondary
  sources place Lang's statement at pp. 212–217.

## 5. Corrected statement

**Setting.**
* *Max-form LW over `Q(i)`*, `LW_max(eps, C)`: for all `n` and `alpha_k in Q(i)^*`,
  `|Lambda| >= C^n B^(-(n-1+eps)) prod_k A(alpha_k)^(-(1+eps))`, with
  `A(alpha) = max(H(alpha), e^(|Log alpha|))` and `H` the naive height.
* *The Lang–Waldschmidt conjecture*: the product form (Lang 1978, pp. 212–217; Waldschmidt 2004, Conj. 2.5 for
  rational integers; Waldschmidt 2000, Ch. 1), with `C(eps)^n B/(|b_1 ... b_n| A_1 ... A_n)^(1+eps)`.

**Theorem A (single forms; proved).**
1. For every representation of a character form `sum_j v_j Log u_j + b Log i` by numbers of `Q(i)^*`, the
   logarithmic naive heights sum to at least `w(supp v)`.
2. Hence, under `LW_max(eps, C)` or the product form, a single linear form can exclude an `M`-point cluster on an
   arc `C_0 sqrt R` only through a character with `w(supp v) < W/4 + log(1/C_0)`.
3. It does exclude the cluster when a character with `K = O(W/(log M)^3)` distinct values, coefficients
   `e^(O((log M)^2))` and `w(supp v) < W/(4(1+eps)) - (1+eps) W/log M` exists, for `M >= M_1(eps, C(eps), C_0)`.
4. `flipped.md` Theorem L is correct. Under LW, the genuine versions of the construct, lower and sharp
   (P3-char) fake profiles do not exist at Paley orders `M = q + 1 >= M_0(eps, C)` (ineffective).

**Theorem B (skew double core; proved).**
* For the Paley skew double core of order `M = q + 1` (`q = 3 mod 4` prime), every nonzero zero-sum `c` has
  `s(c) >= M/2`.
* The exact minimum is `M/2 + 1` at `M = 8, 12`, certified by rank. The minimisers are exactly the `2q` type
  vectors.
* With primes in `[P, P e^(1/(8q))]`, every character has `w(supp v) >= W/4 + (1/2) log P - 1/16`.
* So no single LW inequality, in either form, for any `eps >= 0` and any representation by numbers of `Q(i)^*`,
  beats the arc's a priori bound for this profile.

**Theorem C (simultaneous fakes for fixed `eps`; proved for the max form with `i` listed).**
* For `0 < eps <= 1`, `0 < C_0 <= 1` and `M = q + 1 >= M_2(eps, C_0)`, there are angle models on the skew double
  core with `M` points on an arc `C_0 sqrt R`, every Lemma 0.1 gap, and `LW_max(eps, C)` (`C <= 1`) for every
  linear form in the logarithms of numbers of `Gamma_J` *whose list contains `i`*.
* Their size is `W <= (max(8, 4/eps) + o(1)) M log M + 8 M log(1/C_0)`. The constant `4/eps` can be lowered to
  `2/eps + o(1)` by the same method.
* Hence, within the input class {exponent profile, all monomial gaps, Haar-generic angle properties, those LW
  instances}, no argument proves `M <= (min(1/4, eps/2) - eta) log R/loglog R`. With the improved constant this
  becomes `min(1/4, eps)`.
* Point-level residue data can be added if `exact.md` Theorem 5.1 (unrefereed) holds.

**Not proved.**
1. Simultaneous fakes for the product-form (standard) conjecture at fixed `eps`.
2. Fakes for homogeneous instances.
3. Fakes for representations by numbers of higher degree.
4. Anything `eps`-uniform: the full conjecture with `eps -> 0` remains open, as the write-up says.
5. Anything involving the prime-level parity law (P3-char/P3-all) or additive data.

The question "does LW give `M = o(log R/loglog R)` for actual circles?" is therefore answered "not through single
forms on top of the angle-model inputs, and not through `LW_max(eps)` at a fixed `eps`". That is a barrier, not a
growth theorem. The decisive open cases are the ones the write-up lists in §8, (S_mu) and (N), together with
items 1 and 2 above.

## 6. Files (`round8/referee_lw_checks/`)

| script | checks | output |
|---|---|---|
| `r1_dcore.py` | independent skew double core; (5.1); `Sigma_a = R(a) mod q`; `s >= M/2` on random and adversarial Case-2 `c` (`q = 7..59`); column / label-pair / row-pair supports | `r1_dcore_out.txt` |
| `r2_census.py` | R2: rigorous `min s = M/2 + 1` (`M = 8, 12`) by rank mod `2^31 - 1` over all `(#types - M/2)`-sets; R3: independent census by signed maximal minors with exact verification | `r2_census_out.txt` |
| `r3b_families.py` | positive-dimensional families of characters (degenerate spans), generic support | `r3b_families_out.txt` |
| `r4_valuation.py` | Lemma 2.1 on 26,244 numbers; Lemma 2.2 by exhaustive representations (`K <= 3`) with units, `1+i`, `3`, `17` | `r4_valuation_out.txt` |
| `r5_union_bound.py` | Theorem 6.1(a) sum; character part at `P_1`; improved split (`2/eps`); Corollary 3.3 slip; Bell-number count | `r5_union_bound_out.txt` |
| `r6_random_bases.py` | the random two-basis light characters of §7 (`s = 2, 4, 4, 4`); Lemma 7.1 input for `q = 19..47` | `r6_random_bases_out.txt` |
| `rerun_flipped/` | copy of `round4/flipped_checks/check_two_log_identity.py` (+ `fcommon.py`), re-run | `my_two_log_out.txt` (identical to the stored output) |

Exact arithmetic throughout:
* `r2_census.py` evaluates minors in floating point only to propose integer kernel vectors (`|minor| < 6e5`).
  Every vector is then verified by exact integer products, and every reported support is recomputed in Python
  integers.
* The rigorous minimum (R2) uses modular ranks only.
* `r5_union_bound.py` is arithmetic on the proof's bounds (floats in log space), not a claim about data.
