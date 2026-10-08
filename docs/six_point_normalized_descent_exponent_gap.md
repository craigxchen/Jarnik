# What polynomial-height square descent changes at the endpoint

Polynomial-height squareclass representatives remove the uncontrolled
coefficient constants from the finite branch descent. They do not change
its power of the target height. In degree twelve the elementary norm
comparison remains compatible with an isolated endpoint solution, even
if the target height greatly exceeds the coefficient height.

## 1. The normalized norm comparison

Let `F=Delta_u` be the controlled squarefree degree-twelve binary form,
of coefficient height `U^O(1)`. Let `(s,t)` be primitive, put
`H=max(|s|,|t|)`, and assume the necessary endpoint window

```text
F(s,t)=d^2!=0,       |d|<=C_end U^a H.                (1)
```

Above a polynomial threshold in `U`, the branch approximation theorem
places `r=s/t` (or the reciprocal chart, with denominator `H`) near a
unique simple real root `alpha`, with

```text
0<|r-alpha|<=C_end^2 U^c H^-10.                      (2)
```

Let `K` be the field of this root, of degree `m<=12`. Use the integral
root `theta=a0 alpha`, where `a0` is the nonzero leading coefficient,
and the normalized square descent

```text
a0 s-theta t = eta beta^2,
beta in b^-1,
max_sigma max(|sigma eta|,|sigma eta|^-1)<=U^C,
Norm(b)<=U^C.                                      (3)
```

The theorem in
`six_point_polynomial_squareclass_representatives.md` establishes (3)
with a finite list of polynomially many choices, independently of `H`.
One can absorb
fixed endpoint constants by writing `C(C_end) U^C` below. Root bounds,
separation, and (2) give

```text
|sigma_0 beta| <= C(C_end) U^C H^-9/2,
|sigma beta| <= C(C_end) U^C H^1/2  (sigma!=sigma_0).
```

The first embedding is real. Complex embeddings in the remaining
product are counted with both conjugates. Taking the field norm yields

```text
|Norm(beta)| <= C(C_end) U^C H^((m-10)/2).            (4)
```

Since `(beta)b` is a nonzero integral ideal, its norm is an integer at
least one. Thus (3) also gives

```text
|Norm(beta)| >= U^-C.                               (5)
```

For `m<=9`, combining (4)--(5) gives a polynomial bound on `H`:

```text
H^((10-m)/2) <= C(C_end) U^C.                        (6)
```

This recovers precisely the low-degree branch consequence already
proved by integer polynomial evaluation in
`six_point_endpoint_parameter_gap.md`. For `m=10`, the power of `H`
is zero; for `m=11` it is one half; and for `m=12` it is one. None of
these three cases gives an upper bound on `H` from this comparison.

Equivalently, in degree twelve the eleven ordinary embeddings and the
fixed ideal denominator force only

```text
|sigma_0 beta| >= U^-C H^-11/2,
```

whereas the endpoint condition supplies an upper bound of order
`U^C H^-9/2`. These overlap by a factor `H`, apart from polynomial
coefficient costs. Even an additional assumption that `H` exceeds every
fixed power of `U` would not reverse these powers of `H`.

## 2. Finite twists do not bound the first solution

The normalized twists have controlled equations and form a finite list
of size `U^O(1)`. This is useful input for any subsequent arithmetic
argument, but the list contains possible equations, not a bound on
their solutions. An isolated point belongs to one member of the list;
there is no pigeonhole contradiction. As the central vector varies,
the list and its number fields vary as well.

The residual equations are substantial: `eta beta^2` must lie in the
rational two-plane `span_Q(a0,theta)`, with the original primitivity and
square-norm conditions. An estimate exploiting these equations could
improve the result. Discarding them and retaining only (4)--(5) cannot.
Fixed-curve finiteness, or a gap between successive approximations to a
fixed branch root, does not by itself bound the first isolated solution
uniformly in these moving coefficients.

## 3. The fair scale does not make coefficient powers negligible

For a full fair extracted six-core,

```text
log N=31w+o(w),       log U=7w+o(w),
log X=8w+o(w),       X=least endpoint cotangent.
```

The controlled-frame radius comparison is only

```text
U^-C2 H^20 <= N <= U^C2 H^20.
```

Consequently it gives the interval

```text
(31-7C2)/140 + o(1)
    <= log H/log U <= (31+7C2)/140 + o(1),            (7)
```

with the trivial lower bound zero when the displayed lower endpoint
is negative. In particular `log H=31w/20+o(w)` does not follow: the
coefficient loss has the same exponential scale as `N` and `H`.
Nor does the extraction give `H` larger than every fixed power of `U`.

Even a future theorem `H<=U^C` would require its exponent and the
frame conversion costs to be compared quantitatively with the fair
scales. A concrete sufficient target in the intrinsic cotangent variable
is the already identified estimate `X<=U^(c+o(1))` with `c<8/7`.
An unspecified polynomial exponent, a finite twist count, or a bound
exponential in a power of `U` does not establish that target.

The new normalization removes a real coefficient-height obstruction.
The remaining obstruction exhibited here is the unchanged norm exponent,
together with the need to control isolated solutions of the two-plane
equations at the actual moving-weight scale. This is not a theorem
excluding stronger square-descent or approximation methods.
