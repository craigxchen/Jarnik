# A bounded test of the effective polynomial model

The proposed polynomial model has real polynomials `F_i` of a common
degree `D`, with

```text
deg(F_i-F_j)=D,
F_i-F_j divides F_i F_j+1       in R[t].                     (1)
```

It arises when `h_i=f_i+i c_i`, with nonzero constant `c_i`, and
all primitive pair imaginary parts are nonzero constants: put
`F_i=f_i/c_i`. The difference in (1) has no real zero, because a
real common value at such a zero would satisfy `F_i^2+1=0`.

Condition (1) by itself admits arbitrarily large families, even in
degree two. In that degree there is a complete elementary
classification: the families are pencils with a common complex
linear factor. They therefore fail the full uniform independent-block
profile. This is a useful separation of the algebraic division
condition from the stronger faithful endpoint model, not a
counterexample to endpoint uniformity.

## 1. Arbitrarily many exact degree-two polynomials

For distinct nonzero real constants `a`, set

```text
F_a(t)=a(t^2+1)+t.
```

Every pair satisfies

```text
F_a-F_b=(a-b)(t^2+1),
F_a F_b+1
 =(t^2+1)[ab(t^2+1)+(a+b)t+1].                              (2)
```

The differences have degree two and no real zeros. The exact complex
factorization is

```text
F_a+i=(t+i)(a t+1-i a).                                    (3)
```

The factors for distinct `a` have no further common roots. Dividing
the primitive transition by `N(t+i)=t^2+1` gives

```text
Im((F_b+i)(F_a-i)/(t^2+1))=a-b.                            (4)
```

Thus all imaginary parts in this model are constants, both at the
anchor and after primitive pair reduction. For integer `a`, all
displayed polynomials are Gaussian integral polynomials.

The common factor in (3) is the obstruction to using this as a full
uniform profile. For `m` nonanchor rows the complex polynomial lcm
has degree `m+1`: one common linear factor and `m` distinct private
linear factors. Meanwhile each anchor numerator has degree two and
its angle is of order `t^(-2)` along the positive real axis. The
normalized arc length on the polynomial least-radius construction
therefore has order

```text
t^((m+1)/2-2)=t^((m-3)/2).                                  (5)
```

It grows when `m>=4`. In particular, the example does not supply the
small radius required by the faithful full-pattern construction.

## 2. Classification when D=2

Let at least three distinct real quadratic polynomials satisfy (1)
for every pair. Then there is one quadratic polynomial `P`, positive
on the real axis, such that

```text
F_i=F_0-a_i P,
P divides F_0^2+1,                                          (6)
```

with distinct real constants `a_i` (and `a_0=0`). Conversely every
such family with all degrees equal to two satisfies (1).

To prove this, fix `F=F_0`. The positive quartic `F^2+1` has four
distinct nonreal roots. A repeated root would also be a root of
the linear polynomial `F'`, whose root is real, which is impossible.
Thus over the real numbers

```text
F^2+1=c P Q,
```

where `P,Q` are distinct monic irreducible quadratics and `c>0`.
Every difference `F-F_i` has degree two and divides this product,
so it is a nonzero scalar multiple of `P` or `Q`.

Suppose both types occur, say

```text
F_i=F-aP,       F_j=F-bQ,       a b !=0.
```

The pair condition says `aP-bQ` divides `F_i^2+1`. Since it is
coprime to `P`, it divides

```text
cQ-2aF+a^2P.
```

Both have degree at most two, and the divisor has degree exactly two
by hypothesis. It follows that `F` lies in the real linear span of
`P,Q`, say `F=uP+vQ`. Consequently

```text
(uP+vQ)^2-cPQ=-1.                                          (7)
```

Factor the binary quadratic form on the left over `C`. If its two
linear factors are independent, their evaluations at `(P,Q)` are
polynomials with nonzero constant product. Both evaluations must be
constant, forcing `P,Q` to be constant, a contradiction. If the two
linear factors coincide, its discriminant vanishes. Since `c>0`,
this gives `c=4uv`, and the left side of (7) is the real square
`(uP-vQ)^2`, which cannot equal minus one.

Thus the two types cannot both occur. All differences use one fixed
factor `P`, proving (6). Its roots are a conjugate pair; exactly one
of them is a zero of `F_0+i`. Hence all the `F_i+i` share the same
complex linear factor. Their other linear factors are distinct.

For arbitrary real coefficients that linear factor is asserted only
over `C`. If the `F_i` have rational coefficients, choose `P` monic
over `Q` from a nonzero difference. The degree-one gcd of `P` and
`F_0+i` is then computed over `Q(i)`, so the common linear factor
is Gaussian-rational. No rationality of that factor is assumed in
the real-coefficient classification.

This proves that the degree-two division model produces only a common
pencil, even when no common factor was assumed beforehand. Such a
pencil has only a full cut and private factors, rather than every
positive-weight cut of a uniform independent-block profile.

## 3. The precise remaining higher-degree issue

In larger degree a difference can choose many real factors of
`F_0^2+1`. If two such degree-`D` divisors have common factor `C`,
write them as `CP` and `CQ`, with `gcd(P,Q)=1`. Their pair condition
reduces to a divisibility whose quotient can have degree `deg C`.
The degree-two proof used `deg C=0`; it therefore does not extend
by simply repeating the same argument. Uniform profiles have
substantial intersections, precisely where that quotient is allowed
to be nonconstant.

No proof has been obtained that the faithful full-profile polynomial
model is impossible for large `m`, and no such polynomial construction
has been produced here. The concrete conclusions are the exact
arbitrary-size pencil example and the complete quadratic
classification, both with all primitive divisions retained.
