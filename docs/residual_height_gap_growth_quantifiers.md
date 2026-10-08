# Which residual-height gaps would improve the endpoint growth bound?

This is a conditional quantifier audit. None of the arithmetic lower bounds
posed below is proved here. The current general prose bound remains
`M <= (1+epsilon) log R/loglog R`, and the radius-independent bound remains
unproved.

This note concerns the **prime-disjoint extraction with corrections**.
The new [exact nested-layer extraction](exact_nested_profile_residual_extraction.md)
has no corrections and aggregate residual height `O(k^2 W/M)`, while
allowing prime overlap on two nested orientation chains. For that different
class, a lower bound exceeding `log w` by an unbounded factor would improve
growth, and `w/log w` would give `O(loglog R)`. No such arithmetic lower
bound is currently proved. The comparisons below remain valid for the
prime-disjoint class and must not be transferred by dropping its hypotheses.

For one fixed profile dimension, the extraction would turn a logarithmic
residual-height lower bound larger than `sqrt(w log w)` by an unbounded
factor into an improvement in the order of growth. A bound `w/log w`
would give `O((loglog R)^2)`. A fixed positive multiple of `w` would
give uniformity. Sublinear lower bounds, even an entire family indexed by
all fixed dimensions, do not by themselves force uniformity through these
extraction comparisons.

## 1. The height convention and the exact extraction loss

Use the normalization of
[endpoint_full_profile_quantifiers.md](endpoint_full_profile_quantifiers.md):
subdivide to a fixed `C_0<sqrt(2)`, remove the whole cluster's Gaussian
gcd, and retain its cardinality `M`. Write

```
W = log(R^2),
lambda = log(sqrt(2)/C_0) > 0,
k = m+1,       b = 2^m,       w = W/b.
```

The selected tuple has `m` nonanchor rows. The scale `w` refers to the
original normalized conductor, before optional empty-block removal. The
exact primitive rows and disjoint-support blocks satisfy the conditions
in the linked extraction note, including the allowance that correction
factors share primes with incident blocks. Define the individual
**logarithmic** residual height

```
H = max( max_i log max(1,|K_i|,|Y_i|),
         max_(i<j) log |t_ij| ).
```

Thus a proposed lower bound `H>=sqrt(w)` concerns the logarithm of a
residual or correcting factor, not its ordinary magnitude.

For `m>=2` and `M>=16k^2`, the audited extraction gives

```
W >= 4lambda(M-1),
epsilon_profile <= 8kb/sqrt(M),
H <= 4k W/sqrt(M) = (4kb/sqrt(M))w.                  (1)
```

In particular, the factor `b` occurs when comparing height to one block;
it cancels from the absolute height upper bound in `W`. It must reappear
when an arithmetic lower bound is evaluated at `w=W/b`.

An aggregate version is also available. With the products from the
extraction note, put

```
H_sum = log K_all + log Y_all + log B_all.
```

Its exact upper bound from that note is at most

```
[(2+sqrt(2))/2] m k^2 W/sqrt(M)
 < 2k^3 W/sqrt(M) <= 4k^3 W/sqrt(M).                 (2)
```

The extra powers of `k` matter in growing dimension. The following
formulas use `H_r` and `r`, where `(H_1,r)=(H,1)` or
`(H_3,r)=(H_sum,3)`, and the convenient common bound is
`H_r<=4k^r W/sqrt(M)`.

## 2. A precise hypothetical theorem and its inversion

Suppose that, for some fixed `m>=4`, there are fixed positive constants
`delta_m,w_0(m)` and a positive function `f_m(w)` such that every genuine
profile in the extracted class with

```
w >= w_0(m),
epsilon_profile <= delta_m,
H_r <= delta_m w
```

obeys the hypothetical arithmetic lower bound

```
H_r >= f_m(w).                                      (3)
```

The small-height hypothesis may be omitted if the intended arithmetic
theorem does not need it. It is included here to permit a theorem proved
only near the extracted limiting regime. Robustness on a fixed relative
profile neighborhood is essential: a theorem only at exactly equal block
weights cannot be substituted for (3).

Set

```
Q_m = max(16k^2,
          [max(8k,4k^r)b/delta_m]^2,
          1 + b w_0(m)/(4lambda)).                   (4)
```

Whenever `M>=Q_m`, (1) guarantees every hypothesis of (3), including
its height threshold. Combining the two height inequalities gives the
explicit necessary comparison

```
f_m(W/b) <= 4k^r W/sqrt(M),
M <= [4k^r W/f_m(W/b)]^2.                            (5)
```

For `W>=b w_0(m)`, a bound valid without assuming `M>=Q_m` is

```
M <= max(Q_m, [4k^r W/f_m(W/b)]^2).                  (6)
```

Fixed arc subdivision and the passage back from the normalized radius
only affect fixed constants in the growth conclusions below. In
particular, the comparison with the existing rate is with
`W/log W`, since `W=2log R`.
For an arbitrary nonmonotone `f_m`, returning to the original radius
uses the eventual monotone envelope of the resulting bound; one must
not simply substitute the original `W` into a nonmonotone expression.

The restriction `m>=4` is substantive: the known full two- and three-row
integer families have bounded residual heights as `w` tends to infinity.
They preclude an eventually unbounded universal lower bound of the stated
kind at those row counts. No such lower bound at four or more rows is
being asserted.

## 3. Fixed-dimension growth consequences

In this section `m`, the theorem's constants, and any positive multiplier
in `f_m` are fixed. All conclusions are conditional on (3).

| Proposed eventual lower bound, up to a fixed positive multiplier | Bound supplied by (6) | Comparison with `W/log W` |
|:---|:---|:---|
| `loglog w` | `O(W^2/(loglog W)^2)` | No improvement |
| `sqrt(w)` | `O(W)` | No improvement |
| `w^alpha`, `0<alpha<1` | `O(W^(2-2alpha))` | Improves the order exactly when `alpha>1/2` |
| `sqrt(w)(log w)^beta` | `O(W/(log W)^(2beta))` | Improves the order when `beta>1/2` |
| `sqrt(w log w loglog w)` | `O(W/(log W loglog W))` | A strict improvement |
| `w/log w` | `O((log W)^2)` | `O((loglog R)^2)` |
| `w/loglog w` | `O((loglog W)^2)` | `O((logloglog R)^2)` |
| `c(m)w`, `c(m)>0` | A constant independent of `W` | Uniformity |

For `alpha=1/2`, a pure power gives only `O(W)`, not the current
`O(W/log W)` bound. At the finer borderline
`f_m(w)=c sqrt(w log w)`, (6) has the same order as the existing result;
constants might improve a leading coefficient but do not improve growth.

More generally, the sufficient regular threshold for a strict order
improvement at one fixed dimension is

```
f_m(w)/sqrt(w log w) -> infinity.                     (7)
```

Indeed, write `f_m(w)=sqrt(w log w)L(w)` with `L(w)->infinity`.
Then (5) reads exactly

```
M <= 16 k^(2r) b W / [log(W/b) L(W/b)^2],             (8)
```

which is `o(W/log W)`. A lower bound only on a subsequence of `w`
would not yield a bound for all large radii. There is no smallest
unbounded factor `L`, so no literally weakest lower-bound function in
this class. A simple concrete target is the fifth line of the table.

Likewise a hypothetical bound `f_m(w)=c_m w/L(w)` gives

```
M <= 16 k^(2r) b^2 L(W/b)^2/c_m^2.                   (9)
```

This explains why an almost-linear but genuinely sublinear gap can give
an extremely slowly growing bound without giving a constant bound.

## 4. Optimizing the dimension does not erase its losses

For a family of hypothetical theorems, keep every `delta_m`, `w_0(m)`,
and `f_m` explicit. Formula (6) gives the valid envelope

```
M <= inf_(m>=4 with W/2^m>=w_0(m))
       max(Q_m, [4(m+1)^r W/f_m(W/2^m)]^2).           (10)
```

If the index set is empty, (10) supplies no bound. The same comparison
can be used directly with a chosen `m=m(M,W)`, but the admissibility
requirements must then be checked:

```
M >= 16(m+1)^2,
max(8(m+1),4(m+1)^r)2^m/sqrt(M) <= delta_m,
W/2^m >= w_0(m).                                    (11)
```

For a fixed neighborhood width these give
`2^m=O(sqrt(M)/(m+1)^r)`. If errors are required to tend to zero,
replace this `O` by `o`. Thus the near-boundary choice is
`m<=0.5 log_2 M-r log_2 log_2 M-h(M)` with `h(M)->infinity`,
subject also to the height thresholds. The cases `r=1,3` agree with
the individual and aggregate ranges in the extraction note.

For example, if

```
f_m(w) = c_m w^alpha L(w),
```

then the second term of (10) is exactly

```
16 (m+1)^(2r) 2^(2alpha m) W^(2-2alpha)
 / [c_m^2 L(W/2^m)^2].                              (12)
```

When `f_m=f` is the same nondecreasing function for every allowed `m`,
the expression in (5) worsens as `m` increases: its numerator grows and
its denominator decreases. If neighborhood and height thresholds are
also uniform, the smallest dimension admitting the theorem is optimal.
There is no free gain from setting `m` close to half of `log_2 M`.

A useful gain from growing `m` must instead be supplied by the arithmetic
theorem's actual dimension dependence. In (12), the coefficient `c_m`
must compensate for `(m+1)^r 2^(alpha m)` and the changed argument of
`L`; its profile width and threshold costs in (11) must also permit that
choice. No uniformity in those constants follows from a theorem stated
separately for each fixed `m`.

For a linear bound, however, one fixed dimension is enough regardless
of how badly the coefficients at other dimensions behave. From
`f_m(w)=c_m w` and (6),

```
M <= max(Q_m, [4(m+1)^r 2^m/c_m]^2).                 (13)
```

Thus even a very small `c_m>0` gives uniformity when it belongs to one
fixed applicable theorem with a fixed positive neighborhood width and
finite height threshold.

## 5. Why sublinear gaps alone leave arbitrarily slow counts open

The following are compatibility examples for the numerical inequalities,
not constructions of lattice-point clusters or Gaussian profiles.

Take

```
M_n=n,       W_n=exp(n^2),       R_n=exp(exp(n^2)/2).
```

For every fixed `m`, the extracted relative error tends to zero and
`w_n=W_n/2^m` tends to infinity. The strict-gap condition
`W_n>=4lambda(M_n-1)` also holds. The permitted height upper bound is
`4(m+1)^r exp(n^2)/sqrt(n)`. For the largest of the three listed
sublinear examples `f(w)=w/log w`,

```
f(W_n/2^m) / [4(m+1)^r W_n/sqrt(M_n)]
 = sqrt(n)/[4(m+1)^r 2^m (n^2-m log 2)] -> 0.        (14)
```

The functions `sqrt(w)` and `loglog w` are still smaller for large `w`.
Thus each of those hypothetical lower bounds is compatible with this
unbounded, very slowly growing count scale. The same calculation works
uniformly for the dimensions admitted by (11) when widths are fixed and
the displayed functions have no additional dimension multipliers.

There is a general diagonal version which allows arbitrary dimension
dependence. Suppose every fixed `m` has an eventual lower-bound function
with

```
f_m(w)/w -> 0.                                      (15)
```

For each integer `n`, choose `W_n` so large that `W_n>=exp(n^2)`,
all the finite height thresholds for `4<=m<=n` are passed, and

```
f_m(W_n/2^m)/(W_n/2^m) <= 1/n
                       for every 4<=m<=n.            (16)
```

This is possible because it imposes only finitely many eventual
conditions for each `n`. Put `M_n=n`. For every dimension that could be
extracted from an `n`-point cluster, (16) gives

```
f_m(W_n/2^m) <= W_n/(n2^m)
             <= 4(m+1)^r W_n/sqrt(n).               (17)
```

Dimensions whose profile-width requirements fail are unavailable and
cannot produce a contradiction; all the available dimensions satisfy
(17). Consequently no collection of comparisons (1), (3), and (15),
even for every fixed `m`, rules out all unbounded counts. An additional
arithmetic relation, a non-sublinear gap at a fixed dimension, or a
different extraction mechanism would be needed to deduce uniformity.
The diagonal argument is about these comparisons alone and does not
assert existence of any underlying arithmetic configuration.

## 6. A concrete disjoint-core target and the distinction from earlier audits

A concrete growth-oriented target is: for one fixed `m>=4`, prove (3)
on a fixed positive profile neighborhood with

```
H >= c sqrt(w log w loglog w),       c>0,
```

and with a finite height threshold, allowing all correction and prime-power
features supplied by the extraction. This would give the general bound

```
M = O_C(log R/[loglog R logloglog R]),
```

with constants also depending on that one arithmetic theorem. This target
is much weaker than a positive linear exponent gap, although it remains
unproved. An arbitrarily slower divergent multiplier in (7) would give
an arbitrarily smaller strict improvement.

[weaker_intrinsic_finite_tests.md](weaker_intrinsic_finite_tests.md)
instead weakens the coefficient of a conductor-linear finite-test bonus
while retaining a positive fixed binomial gap; its conclusions are still
conditional uniform bounds. It does not invert sublinear residual-height
gaps as above.

[reciprocal_divisor_conditional_growth_audit.md](reciprocal_divisor_conditional_growth_audit.md)
obtains a similar-looking growth rate under the additional existence of
an admissible nonconformal endpoint-to-endpoint map. Its integer divisor
argument is different and that extra map is not available for arbitrary
arcs. The displayed target here would apply to general endpoint clusters
only if its new arithmetic assertion were proved for all extracted
profiles; no existing conditional-map result proves it.
