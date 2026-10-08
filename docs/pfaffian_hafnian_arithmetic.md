# Pfaffian--Hafnian factorization and actual denominator cancellation

The reciprocal-sum Hafnian has a common quadratic factor that must be
removed before calling it a rational invariant. After this normalization,
its matching terms have exact common denominator `L_X`. The
Pfaffian--Hafnian identity gives an exact integer factorization, but the
new positive factor has height `27W/4+o(W)` in the critical profile.

There is also a concrete cancellation phenomenon beyond the earlier
termwise-denominator audit: on actual Gaussian circle tuples, an
arbitrarily large power of 19 can disappear from the reduced Hafnian
denominator even though 19 divides none of the primitive chord residues.
These tuples have a fixed positive angular span, not endpoint span.
Neither result proves the desired eight-point inequality.

## 1. The common quadratic factor

Use the primitive common-unit model and notation of
[eight_point_pfaffian_content.md](eight_point_pfaffian_content.md):

```text
h_ij=x_ij+i t_ij,       x_ij>0,
q_ij=x_ij^2+t_ij^2=product_p p^|a_i(p)-a_j(p)|,
X=product_(i<j)x_ij,       T=product_(i<j)|t_ij|,
L_X=lcm_M product_(ij in M)x_ij,       B_X=X/L_X | T.
```

The signs of `t_ij` follow the chosen ordering. Write `r_(p,l)` for
the number of marked rows in a threshold layer and define

```text
O=product_p p^kappa_p,       kappa_p=sum_l (r_(p,l) mod2).
```

This positive integer need not be squarefree. For a perfect matching
put

```text
a_M=sqrt(product_(ij in M)q_ij/O).
```

Every `a_M` is a positive integer. At a layer the number of crossing
matching edges has parity `r mod2` and is at least that parity. The
minimum can be achieved simultaneously at every layer of one prime
by sorting the allocations and matching adjacent ranks. Consequently
`O` is exactly the gcd of the 105 integers `product_M q_ij`.
The parity assertion then shows that each quotient is a square.

Let

```text
S_ij=2sqrt(z_i z_j)/(z_i+z_j)
    =sec((theta_j-theta_i)/2)=sqrt(q_ij)/x_ij,
tau_ij=t_ij/x_ij.
```

Consistent square roots along the proper arc make every `S_ij`
positive. The normalized Hafnian is rational:

```text
Hf(S)/sqrt(O)=sum_M a_M/product_M x_ij.                  (1)
```

In general `Hf(S)` itself is not rational. Squaring it does give
a rational invariant.

## 2. Exact matching denominator and positive numerator

The least common denominator of the 105 separate terms in (1) is
exactly `L_X`, even after reducing every term. To prove the lower
bound at a conductor prime, retain only the positive-weight edges
of a matching maximizing `sum_M v_p(x_ij)`. Such edges join rows
with equal allocation at `p`, since unequal allocations make both
primitive pair coordinates `p`-units.

These retained edges can be completed to a matching minimizing
`sum_M |a_i-a_j|`: remove their equal-allocation pairs, then match
the remaining allocations in adjacent sorted order. Removing two
equal entries preserves the parity of every threshold count and
does not change the minimum total distance. The resulting matching
therefore has `v_p(a_M)=0` and the maximum real-part denominator
exponent. No added edge can raise its real weight, by maximality.
At a nonconductor prime all `q_ij` and `O` are units, so any maximum
matching works. This includes the prime two, since the primitive
circle has odd norm.

It follows that

```text
K=sum_M L_X a_M/product_M x_ij in Z_(>0),
Hf(S)=sqrt(O) K/L_X.                                    (2)
```

The integer terms defining `K` have gcd one. For each prime the
constructed matching has maximum denominator valuation and
`v_p(a_M)=0`, so its cleared term is a unit. This includes
`v_p(L_X)=0`: complete the empty positive-edge matching by
adjacent sorted allocations to obtain `v_p(a_M)=0` again.
The least-common-denominator assertion alone would not establish
this numerator primitivity. Nor does primitivity of the separate
terms assert that `gcd(K,L_X)=1`.

The normalized reciprocal-sum kernel in
[sum_kernel_denominator_budget.md](sum_kernel_denominator_budget.md)
uses permutation products and denominator `L_X^2`. Formula (2)
is its matching analogue with the quadratic factor retained.

There is a useful case where the final sum also has no cancellation.
Suppose `p` does not divide `T`. By the exact local opposition
description in the Pfaffian-content note, the edges on which
`p|x_ij` form a matching. If it contains at least three edges,
there is a unique perfect matching containing all of them. Its
norm numerator `a_M` is a unit by the same minimum-norm completion
argument; every other cleared term is divisible by `p`. Thus

```text
p does not divide K   if p does not divide T
and at least three of the x_ij are divisible by p.       (2a)
```

This is an actual statement about the reduced rational sum, not
only its termwise denominator. In this case the denominator of
`Hf(S)^2` has exact exponent
`max(0,2v_p(L_X)-v_p(O))`. With fewer pole pairs, a complementary
Hafnian remains in the leading coefficient and can cancel, as in
Section 5.

## 3. The exact nonlinear integer factorization

The primary Pfaffian--Hafnian identity is

```text
Pf((z_j-z_i)/(z_i+z_j)^2)
 =product_(i<j)(z_j-z_i)/(z_i+z_j) * Hf(1/(z_i+z_j)).
```

See [Ishikawa--Kawamuko--Okada, Theorem 1.1](https://arxiv.org/pdf/math/0408364).
For eight rows, diagonal rescaling by the square roots of the `z_i`
gives exactly

```text
Pf((tau_ij S_ij))=product_(i<j)tau_ij * Hf(S).           (3)
```

There is no missing sign: the relevant powers are `i^4=i^28=1`.
Each Pfaffian matching term divided by `sqrt(O)` is the rational
number `a_M product_M t_ij/product_M x_ij^2`. Therefore

```text
J=L_X^2 Pf((tau_ij S_ij))/sqrt(O) in Z,
J=(T_signed/B_X) K,       |J|=(T/B_X)K.                 (4)
```

Both identities are exact, including conductor prime powers and all
real-part cancellations. If the rows follow their angular order,
`T_signed>0`, so `J` is positive. The divisibility `B_X|T`
ensures integrality of the factor on the right.

This retains the strong nonlinear vanishing: the Pfaffian equals
a product of all twenty-eight tangents times a positive Hafnian.
It does not turn its cleared numerator into a product only of old
primitive residues; `K` is a new factor.

## 4. The actual endpoint scale

On a shrinking arc, `Hf(S)=105+O(Theta^2)>105`. In the critical
fair eight-row profile with inherited weight `W` and
`log T=o(W)`,

```text
log O=W/2+o(W),
log L_X=7W+o(W),
log K=(27/4)W+o(W).                                    (5)
```

The first identity follows because half of the 256 threshold
patterns have odd size. Constant layers have even size, so removing
a common Gaussian factor does not change `O`. The other identities
follow from `log x_ij=W/4+o(W)`, `B_X|T`, and (2).
Equation (4) consequently gives `log|J|=27W/4+o(W)`.

Thus the exact Pfaffian vanishing does not create a small nonzero
integer in this normalization. Its positive Hafnian factor uses the
remaining height. For the rational square, the exact rational expression is

```text
Hf(S)^2=O K^2/L_X^2.                                   (6)
```

The reduction can depend on both the conductor factor `O` and
the numerical cancellation in `K`. The latter is not controlled
prime by prime by `T`, as the next construction shows.

## 5. An actual prime-power cancellation away from every old residue

Let `p=19` and fix five integers

```text
(a_1,...,a_5)=(38,2,22,4,10).
```

The six-row rational function

```text
F_6(x)=Hf(1/(1+u_i u_j)),       u=(38,2,22,4,10,x),
```

is exactly

```text
F_6(x)=P(x)/(721355647620465 product_(j=1)^5(1+a_j x)),

P(x)=2297907909463+113054407839984x
     +1602230240264968x^2+7550950497889920x^3
     +10293635455067632x^4.                             (7)
```

Every displayed denominator is a 19-unit at `x=12`. Modulo 19,

```text
P(x)=13+2x+14x^2+4x^3+14x^4,
P(12)=0,       P'(12)=15 !=0.                           (8)
```

For every `e>=1`, ordinary one-variable Hensel lifting gives an
integer `x_e=12 modulo19` with `P(x_e)=0 modulo19^e`.
It may be chosen even and positive by adding a multiple of `19^e`.
The lifting step is elementary: if `P(x)=0 modulo19^e`, replace
`x` by `x+c19^e` and solve
`P(x)/19^e+c P'(x)=0 modulo19`.

Choose an even positive `b_e=18 modulo19` such that
`v_19(20b_e+1)=e`; among the even representatives modulo `19^e`,
at most one next digit is forbidden. Use the eight half-angle rows

```text
H_i=X_i+i,       X=(20,b_e,38,2,22,4,10,x_e).
```

All `X_i` are even, so every `H_i` is conjugate-primitive and has
odd norm. They yield actual Gaussian integer circle points

```text
G=product_i H_i,       z_i=bar(G) H_i/bar(H_i),
N(z_i)=product_i(X_i^2+1).
```

With primary Gaussian prime representatives, each `H_i` has unit
factor `+i` or `-i`. Its quotient by its conjugate therefore has
unit factor `-1`, independent of the row. The `z_i` belong to one
common literal unit class. A common Gaussian divisor may be removed
without changing any of the quantities used below.

The eight residues of `X_i` modulo 19 are distinct. Since 19 is
inert, all their norms are 19-units. Moreover the only pair with
`1+X_i X_j=0 modulo19` is the first pair `(20,b_e)`.
Consequently all twenty-eight primitive residues are 19-units,
and

```text
v_19(T)=0,       v_19(L_X)=e.                            (9)
```

To see the cancellation, write
`F_8(X)=Hf(1/(1+X_i X_j))` and expand on whether a matching uses
the first pair. Exactly,

```text
F_8(X)=F_6(x_e)/(1+20b_e)+(19-adically integral terms).
```

The first quotient is 19-integral by (7)--(8); all other terms
have unit denominators. Hence `F_8(X)` is 19-integral. On these
actual directions,

```text
Hf(S)=sqrt(product_i(X_i^2+1)) F_8(X),
Hf(S)^2=product_i(X_i^2+1) F_8(X)^2.
```

Thus the reduced denominator of `Hf(S)^2`, and also of
`Hf(S)^2-105^2`, has zero 19-valuation despite (9). In (2),
`O` is a 19-unit, so this also proves

```text
v_19(K)>=e,       v_19(gcd(K,L_X))>=e,       19∤T.       (10)
```

This rules out any universal primewise assertion that this final
denominator cancellation is supplied only by powers of the old
primitive residue product and a fixed coefficient denominator.
It is not a counterexample to a numerical height bound involving
`T`, since the other primes and the archimedean size of `T` can grow.

All these tuples lie in a proper arc: the `X_i` can be chosen at
least two, and all half-angle arguments are between zero and
`arctan(1/2)`. But the rows with `X=2` and `X=4` retain a fixed
positive angular separation while the radius grows. Thus they do
not satisfy a shrinking endpoint-arc hypothesis. The construction
isolates a real special-value cancellation that a successful
endpoint argument would have to constrain globally.

The least radius also grows after common-factor removal: the moving
row against a fixed `X_i+i` has primitive norm tending to infinity,
since its common divisor divides that fixed row. This avoids using
artificial common Gaussian inflation to produce the examples.

## Exact verification and remaining scope

[check_pfaffian_hafnian_arithmetic.py](check_pfaffian_hafnian_arithmetic.py)
checks the common quadratic factor, exact term denominator, integer
factorization (4), the polynomial certificate (7), and successive
Hensel cancellation examples with their actual Gaussian primitive
pair coordinates.

The Pfaffian--Hafnian identity is a valid additional nonlinear
compatibility equation. Its current exact normalization does not
prove a residue-growth exponent or a uniform point count. The
unproved input would have to constrain the new positive numerator
`K` or its reduced denominator using the actual shrinking arc,
without excluding the valid cancellation mechanism (7)--(10).
