# A fourth-moment gap for Paley characters through square-root support

The finite [order-twelve all-character obstruction](paley_twelve_full_character_height_obstruction.md) does not extend by enumeration. A uniform character-sum argument does, however, exclude cheap saturated characters throughout an **unbounded support range** in growing normalized Paley matrices. For each fixed `0<epsilon<1/3`, every zero-sum ternary character of support

```text
4 <= h <= (1/3-epsilon) sqrt(q)
```

has normalized Hadamard transform norm bounded below by a constant strictly greater than one, depending only on `epsilon`, once `q` is sufficiently large. With `b>=5` repeated columns and one distinct physical one-entry flip per row, every such character has physical coefficient sum above half the number of physical primes. The same is true for all integer characters of amplitude at least two and for row-pair characters. Thus any cheap saturated character in this model must eventually be ternary and have support larger than `(1/3-epsilon)sqrt(q)`. This is a support-range obstruction, not an all-character theorem or an endpoint arc result.

## Fourth moment of a ternary Paley transform

Let `q` be a prime with `q=3 mod 4`, `M=q+1`, and `H` the normalized Paley Hadamard matrix from the [quartic-slack note](general_four_row_quartic_slack_obstruction.md). Its first column is all ones, every other column has zero row sum, and `HH^T=MI`. Let `c` be a zero-sum vector in `{0,1,-1}^M` with support `h`. Write

```text
T_a = sum_x c_x H_(x,a),
L(c) = (1/M) sum_a |T_a|,
E_k = (1/M) sum_a |T_a|^k.
```

The constant column has `T_infinity=0`. Orthogonality gives `E_2=h`. For distinct rows `x_1,...,x_4`, let

```text
P(x_1,...,x_4) = sum_(all columns a) product_(i=1)^4 H_(x_i,a).
```

The [uniform Paley quartic estimate](general_four_row_quartic_slack_obstruction.md) gives

```text
|P(x_1,...,x_4)| <= 3 sqrt(q)+5.                       (1)
```

Indeed the nonconstant-column sum has absolute value at most `3sqrt(q)+4` by the squarefree cubic or quartic quadratic-character bound and at most four diagonal replacements; the constant column adds one. The needed complete-character estimate is stated in Kim--Yip--Yoo, [Lemma 2.1](https://link.springer.com/article/10.1007/s11139-024-00888-5).

Expand `T_a^4` and sum over columns. Terms involving two distinct rows with odd multiplicities vanish by Hadamard orthogonality. The single-row and two-row terms contribute `M[h+6 binom(h,2)]`. The remaining terms use four distinct rows, so the **exact** identity is

```text
E_4 = 3h^2-2h
      +(24/M) sum_(four-row subsets Q of supp(c))
               (product_(x in Q) c_x) P(Q).            (2)
```

Put `delta_q=(3sqrt(q)+5)/(q+1)`. Equation (1) and `24 binom(h,4)<=h^4` give

```text
E_4 <= 3h^2-2h+24 binom(h,4) delta_q
    <= 3h^2-2h+h^4 delta_q.                            (3)
```

Interpolation between the first, second, and fourth absolute moments gives `E_2^3<=L(c)^2 E_4`. Hence

```text
L(c) >= h^(3/2) / sqrt(3h^2-2h+h^4 delta_q)
     = sqrt(h) / sqrt(3-2/h+h^2 delta_q).             (4)
```

This uses only the uniform four-row Weil bound, rather than counting individual certificates. The bound is deliberately worst-case over the signs and locations of the support.

For completeness, fix `epsilon in (0,1/3)` and put `alpha=1/3-epsilon`. Since `delta_q sqrt(q)->3`, sufficiently large `q` satisfies `delta_q h<=1-2epsilon` whenever `h<=alpha sqrt(q)`. Set `H_epsilon=ceil(6/epsilon)`. For `h>=H_epsilon`, dividing (3) by `h^3` gives

```text
E_4/h^3 <= 3/h-2/h^2+delta_q h <= 1-3epsilon/2.
```

For the finitely many `4<=h<H_epsilon`, increase `q` until `delta_q h<=1/8`. The function `3/h-2/h^2` decreases for `h>=4` and equals `5/8` at `h=4`, so `E_4/h^3<=3/4`. Thus the following **uniform** gap holds across the entire stated support range:

```text
L(c) >= gamma_epsilon
     = min(sqrt(4/3), (1-3epsilon/2)^(-1/2)) > 1.    (5)
```

An explicit conservative version needs no implicit threshold. For `q>=4099`, one has `delta_q<=4/sqrt(q)`. If `4<=h<=sqrt(q)/16`, then `delta_q h<=1/4`, so (3) gives `E_4/h^3<=7/8` and

```text
L(c) >= sqrt(8/7) > 17/16.                           (6)
```

The range is nonempty at the first admissible prime `q=4099` and grows without bound.

## Consequence for distinct row-assigned flips

Take `b>=5` physical copies of each of the `q` nonconstant Paley columns, giving `r=bq` physical columns. Flip one entry in one physical column at each row, and assume the `M` chosen physical columns are distinct. For saturation, also assume each column label receives at most two assignments; then at least three unflipped copies of every label remain. For any zero-sum integer `c`, let

```text
v_j = (1/2) sum_x c_x s_(x,j),
B(c) = sum_(physical j) |v_j|,
A = max_x |c_x|.
```

The zero row sum makes every `v_j` integral. A selected flip changes only its own physical coefficient by absolute value at most `|c_x|`, yielding

```text
B(c) >= (bM/2)L(c)-sum_x |c_x|.                      (7)
```

Hadamard inversion gives `L(c)>=A`. Consequently, if `A>=2`, then

```text
B(c) >= MA(b/2-1) >= M(b-2) > b(M-1)/2=r/2.       (8)
```

If `A=1` and `h=2`, row orthogonality gives `L(c)=1`, and only its two supported rows can change coefficients. Therefore

```text
B(c) >= bM/2-2,
B(c)-r/2 >= b/2-2 > 0.                             (9)
```

For `A=1` and `4<=h<=alpha sqrt(q)`, equations (5) and (7) give

```text
B(c)-r/2 >= (bM/2)(gamma_epsilon-1)+b/2-h > 0     (10)
```

for all sufficiently large `q`, uniformly in `b>=5` and every admissible flip assignment. In the explicit range of (6), `sqrt(8/7)-1>1/16` and `M>=2sqrt(q)` make the right side of (10) greater than `b/2` for every admissible `q>=4099`.

There is a sharper restriction on the possible **norm-one** characters at any order. If a ternary `c` has `L(c)=1`, equality in Hadamard inversion at every support row forces each active transform label to agree with the signed character on all `h` support rows. Hence every active `|T_a|=h`, their number is exactly `s=M/h`, and `h` divides `M`. The `h` by `s` submatrix on these rows and labels is rank one with area `M`; in Hadamard trade terminology it is precisely a minimal rectangular trade by Ó Catháin--Wanless, [Theorem 8](https://arxiv.org/html/1502.02353v1). This is a dictionary for the equality case, not a classification of Paley trades. The [general equality-case identity](paley_twelve_full_character_independent_audit.md) then reads `B(c)=bM/2+h-2k`, where `k` counts support rows assigned to active labels. Label capacity two gives `k<=2s=2M/h`, so

```text
B(c)-r/2 >= b/2+h-4M/h.                              (11)
```

In particular, every norm-one character with `h>=2sqrt(M)` is expensive. A cheap norm-one Paley character, if one exists, must have `h|M` and lie between the square-root lower cutoff proved above and `2sqrt(M)`. This refinement does not control high-support characters whose norm is only slightly greater than one.

The unflipped copy of each label also makes this the **saturated** rational-row character lattice. If `lambda` is rational, has zero row sum, and `lambda^T S` is integral, comparing a flipped copy with an unflipped copy of the same label forces `2lambda_x` integral for each row `x`. Conversely, every zero-sum integer `c=2lambda` gives integral `c^T S/2` because every sign column is congruent to the all-ones column modulo two. Thus (8)--(10) apply to all saturated characters in their stated amplitude and support classes.

For actual distinct split rational primes with `ell<=log p_j<ell+1/r` and `ell>1`, the conjugate-primitive composite attached to any covered nonzero character has exact height `V=sum_j |v_j|log p_j`. Since `B(c)` is integral and strictly above `r/2`, it has `B(c)-r/2>=1/2`; meanwhile `W_0=sum_j log p_j<r ell+1`. Hence

```text
V-W_0/2 > ell(B(c)-r/2)-1/2 >= (ell-1)/2 > 0.       (12)
```

For each **fixed** `b`, the [nearby-prime construction](strict_obtuse_prime_box_countermodels.md) supplies such clusters. As in the prior Paley notes, this is a primitive-height obstruction within a sign-profile family. This argument alone neither proves endpoint localization nor rules out ternary characters beyond its square-root support range.

The later [high-support bound](hadamard_ternary_support_spectral_gap.md)
and [prime-field polynomial theorem](paley_primefield_polynomial_character_gap.md)
close the remaining supports for every capacity-two assignment at
sufficiently large prime Paley orders. The latter also classifies all
norm-one characters. A separate
[canonical skew-Hadamard proof](paley_canonical_all_character_half_height.md)
gives the exact minimum for the diagonal assignment at every Paley order.
None of these theorems controls simultaneous angular localization.

Root and a fresh Astra instance independently audited the fourth-moment
expansion, support-uniform epsilon argument, saturation, norm-one capacity
bound, and actual-prime height comparison. The
[companion checker](check_paley_sqrt_support_fourth_moment_gap.py) verifies
12 moment fixtures, a norm-one order-twelve rectangular trade and its
physical flip identity, and the explicit square-root-support cutoff.

The [bounded checker](check_paley_sqrt_support_fourth_moment_gap.py) verifies the exact second/fourth-moment normalization at four Paley orders, the order-twelve norm-one rank-one rectangle and its switch, the physical flip identity, and the explicit `q=4099` constants. The unbounded theorem rests on (1), rather than on these finite fixtures.
