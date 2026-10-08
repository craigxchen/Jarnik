# Fixed rational maps with quasi-integral rational specializations

This note records a necessary condition for a **fixed** rational function to
have very small reduced denominators at rational parameters of growing
height. It is a fixed-map statement only. It does not apply to moving maps,
and it does not say that either allowed pole type realizes a prescribed
Gaussian norm profile.

The degree-seven frame in
[scaled_frame_specialization_denominator.md](scaled_frame_specialization_denominator.md)
is an elementary application of the same principle.

## Statement

Let `f=P/Q in Q(t)` be nonconstant. Choose homogeneous
forms `P(X,Y),Q(X,Y) in Z[X,Y]` of the same degree `d>=1`,
with no common geometric zero, after clearing coefficient denominators. Let

```text
t_n=[a_n:b_n] in P^1(Q),   gcd(a_n,b_n)=1,
H_n=max(|a_n|,|b_n|) -> infinity.
```

Discard the finitely many terms at which `Q(a_n,b_n)=0`, and write the
finite value `f(t_n)` in lowest terms with positive denominator `q_n`.
If

```text
log q_n = o(log H_n),                                  (1)
```

then the pole divisor of `f` over the algebraic closure of `Q` is one of the
following:

1. one rational projective point, with multiplicity `d`; or
2. two real quadratic-conjugate projective points, each with multiplicity
   `d/2`.

This is a necessary classification. It is not a converse, and it makes no
claim about the numerator, the target value, or any norm/arc condition.

## Resultant reduction

Put `g_n=gcd(P(a_n,b_n),Q(a_n,b_n))`. The two homogeneous resultant
Bezout identities, evaluated at the primitive pair `(a_n,b_n)`, imply

```text
g_n | Res(P,Q).                                        (2)
```

Indeed `g_n` divides both `Res(P,Q) a_n^(2d-1)` and
`Res(P,Q) b_n^(2d-1)`, and `gcd(a_n,b_n)=1`. Therefore

```text
log |Q(a_n,b_n)| = o(log H_n).                          (3)
```

after changing the error by a fixed constant. An exact rational pole is one
fixed primitive projective point, so it cannot occur infinitely often while
`H_n -> infinity`. Values at zeros of `P` likewise occur only finitely often
unless the projective point is fixed; they cause no issue in the eventual
subsequence.

## Pole-factor estimate

Factor `Q` over the algebraic closure of `Q` as

```text
Q(a,b)=c product_alpha L_alpha(a,b)^(m_alpha),
```

where `L_alpha` is a linear form for the projective pole `alpha`, and
`sum_alpha m_alpha=d`. Passing to a subsequence, compactness of
`P^1(R)` gives a limit of `t_n`. If that limit is not a
real pole, every factor has modulus bounded below by a positive constant
times `H_n`, contradicting (3). A real
point can be close to at most one fixed pole at a time, since distinct poles
have a fixed positive separation.

Nonreal poles have the same lower bound on every real primitive point:
the complex linear map `(a,b) -> L_alpha(a,b)` is invertible over
`R` up to a fixed constant. Thus the approached pole is real.

### Rational pole

For a rational pole `alpha=[u:v]`, take `u,v` primitive integers and
`L_alpha(a,b)=va-ub`. Away from the exact pole,

```text
|L_alpha(a_n,b_n)| >= 1.                               (4)
```

Every other pole factor has modulus at least a positive constant times
`H_n` along a subsequence approaching `alpha`. If alpha has multiplicity
`m`, then

```text
|Q(a_n,b_n)| >> H_n^(d-m).                              (5)
```

Condition (3) forces `m=d`. This includes `alpha=infty=[1:0]`, whose
factor is `b`: along a sequence tending to infinity, `|b|>=1` and all
other factors have modulus at least a positive constant times `H_n`.

### Irrational real pole

Let `alpha` be an irrational real pole of multiplicity `m`. On a
subsequence approaching it, `|b_n|` is comparable with `H_n`. Roth's
theorem gives, for every `epsilon>0`,

```text
|alpha-a_n/b_n| >>_(alpha,epsilon) |b_n|^(-2-epsilon)  (6)
```

for all sufficiently large reduced denominators. Hence

```text
|a_n-alpha b_n| >> H_n^(-1-epsilon).                   (7)
```

All other pole factors have modulus at least a positive constant times `H_n`, so

```text
|Q(a_n,b_n)| >> H_n^(d-2m-epsilon*m).                   (8)
```

Since `Q` has rational coefficients, every Galois conjugate of `alpha`
is a pole with the same multiplicity `m`. Equations (3) and (8) imply
`d<=2m`. A conjugacy orbit of degree at least three would give
`d>=3m`, impossible. Thus `alpha` is quadratic. Because it is real,
its conjugate is also real; the two poles contribute degree `2m`, so
`d=2m`, and there can be no additional pole factor. This proves the
second alternative.

The approximation input is K. F. Roth, “Rational approximations to
algebraic numbers,” *Mathematika* 2 (1955), 1–20,
[DOI: 10.1112/S0025579300000644](https://doi.org/10.1112/S0025579300000644).
The number-field and one-place forms are also stated in Section 1,
printed pages 31--32, of [Lang's Integral points on curves](https://www.numdam.org/article/PMIHES_1960__6__27_0.pdf).

## Fixed source reparametrization

If `phi in PGL_2(Q)` is fixed and `g=f composed with phi`,
then `phi` bijects `P^1(Q)`, and

```text
h(phi(t))=h(t)+O_phi(1).                               (9)
```

The pole divisor is pulled back by the fixed automorphism `phi`: pole
multiplicities and residue degrees are unchanged, and rational points and
real quadratic-conjugate pairs remain of those respective types. Applying
the statement to `g` therefore gives the same necessary classification in
any fixed rational source coordinate. This invariance is only for a fixed
reparametrization; it does not permit a parameter-dependent `phi_n`.

## Scope

The conclusion uses only fixed `P,Q`, primitive rational parameters, and
subheight reduced denominators. It does not handle moving coefficient forms,
uncontrolled coefficient heights, or a claim that a permitted pole divisor
produces a desired endpoint cluster. In particular, the result supplies a
fixed-map obstruction, rather than a global lattice-point growth estimate.
