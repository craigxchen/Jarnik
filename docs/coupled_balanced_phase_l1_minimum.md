# The complete balanced-cut phase lattice has no quarter-cost vector

Let `m=2k` and let `C` be the `k`-subsets of `{1,...,m}`.  We may use all
cuts or one representative of each complementary pair: if

```text
sigma_i(S) = +1  (i in S),   -1  (i not in S),
```

then `v(S^c)=-v(S)`, so the normalized absolute average is unchanged.
Consider

```text
v(S) = sum_i lambda_i sigma_i(S),
sum_i lambda_i = 0,
v(S) in Z for every S in C.
```

The exact minimum over nonzero `v` is

```text
min  (1/|C|) sum_S |v(S)| = (m-2)/(2(m-1)).              (1)
```

For `m=6,8,10`, this is respectively `2/5, 3/7, 4/9`; in particular it is
strictly larger than `1/4`.  Thus the proposed coupled-phase fractional
packing target, which would require normalized average cost below `1/4`, is
unavailable in the complete balanced-cut row-difference space.

## Integral lattice normal form

Put `n_i=2 lambda_i`.  Since `sum lambda_i=0`,

```text
v(S) = sum_{i in S} n_i.                                (2)
```

Comparing two `k`-subsets that differ by one label shows that all
`n_i-n_j` are integers.  Hence `n_i=a_i+c` for integers `a_i` and a common
rational `c`.  Writing `A=sum_i a_i`, the zero-sum condition gives
`c=-A/m`; integrality in (2) is equivalent to `A` even.  Therefore every
allowed vector has the form

```text
v(S) = sum_{i in S} a_i - A/2,                           (3)
lambda_i = (m*a_i-A)/(2*m),
```

and adding a common integer to all `a_i` changes nothing.  Normalize
`min(a_i)=0` and write `Delta=max(a_i)`.  The denominator of the unique
sum-zero coefficient vector `lambda` is the torsion denominator
`D_lambda`; the corresponding phase targets lie in the grid
`pi/(2*D_lambda)`.

## Exact minimum

If `Delta=0`, (3) is zero.  If `Delta=1`, let `2j` entries of `a` equal one
(`A` must be even).  For a uniform balanced cut, the intersection `X` with
those entries is hypergeometric, and

```text
E|v| = F_j
     = [j(k-j)/k] * binom(k,j)^2 / binom(2k,2j).         (4)
```

The sequence is symmetric under `j -> k-j` and increases from the edge:
`F_1=(m-2)/(2(m-1))`, while `F_j>=F_2>1/2` for every interior value when
it exists.  Hence the smallest `Delta=1` vectors have exactly two ones or
exactly `m-2` ones.

For `Delta>=2`, the conditional two-label swap bound for balanced subsets
gives

```text
E|v| >= m*Delta/(4*(m-1)) >= m/(2*(m-1))
       > (m-2)/(2*(m-1)).                               (5)
```

This is the same conditional-swap estimate used in the complete-cut height
classification in [balanced_pair_products.md](balanced_pair_products.md),
with `W=1`.  Equations (4)--(5) prove (1) for arbitrary rational `lambda`,
without a finite-range assumption on its coefficients.

The minimizers are precisely

```text
v = +/- (sigma_a + sigma_b)/2,   a != b.                (6)
```

For the plus sign, the unique sum-zero coefficients are

```text
lambda_a=lambda_b=(m-2)/(2m),
lambda_i=-1/m  (i not in {a,b}),                        (7)
```

so `D_lambda=m`.  The alternative row-difference vectors
`(sigma_a-sigma_b)/2` have the smaller denominator `D_lambda=2`, but their
average cost is `m/(2(m-1))`, which is larger than the minimum and still
near `1/2`.

## Exact finite audit

The checker [check_coupled_balanced_phase_l1.py](check_coupled_balanced_phase_l1.py)
enumerates normalized integer profiles with levels `0,1,2` at `m=6,8,10`.
It finds the predicted minima and all minimizers:

```text
m=6:  min=2/5, minimizers=30, D_lambda=6
m=8:  min=3/7, minimizers=56, D_lambda=8
m=10: min=4/9, minimizers=90, D_lambda=10
```

The enumeration is only a sanity check.  The lattice normal form and the
`Delta>=2` swap bound establish the all-coefficient theorem.  Therefore a
coupled-phase argument based on this complete balanced sign matrix cannot
obtain the required normalized average below `1/4`; changing from all cuts
to complementary representatives does not alter that conclusion.
