# An integral-closure criterion for fixed factor families

This note isolates the reusable part of the primitive-content argument
from [six_factor_primitive_content_bound.md](six_factor_primitive_content_bound.md).
It gives a fixed-family criterion with a constant independent of the
integer parameters and the resulting circle radius.  It does not assert
that an arbitrary endpoint cluster belongs to a fixed factor family.

## 1. Fixed product rows and their content ideal

Put `O=Z[i]`, `K=Q(i)`, and let `X=(X_0,...,X_d)`.  Fix nonzero
homogeneous linear forms

```text
F_s=A_s+iB_s in O[X],        bar(F_s)=A_s-iB_s,
```

and fixed orientation patterns.  More generally than the binary case,
allow fixed nonnegative multiplicities and write

```text
W_j=eta_j product_s F_s^e_js bar(F_s)^(n_s-e_js),     (1)
```

where `eta_j` is a Gaussian unit. Assume that every row polynomial is
nonzero. All rows have the same degree `D` and, at every real parameter
vector,

```text
|W_j(X)|=rho(X)=product_s |F_s(X)|^n_s.               (2)
```

Let `I=(W_1,...,W_m)` in `O[X]`. Fix a nonzero homogeneous polynomial
`P in O[X]` of degree `D`. On the projective chart `X_k!=0`, put

```text
A_k=K[X_0/X_k,...,X_d/X_k],
I_k=(W_1/X_k^D,...,W_m/X_k^D),
P_k=P/X_k^D,
```

where the redundant ratio `X_k/X_k` is omitted.

The content condition is

```text
P_k belongs to the integral closure of I_k in A_k
for every k=0,...,d.                                  (IC)
```

This condition involves only the fixed coefficient matrix, the row
patterns, and `P`. It does not involve a specialized parameter or its height.
The stronger identity `P_k in I_k` on every chart is an especially
simple certificate; `(IC)` permits the strictly larger integral closure.

## 2. Uniform primitive-content lemma

**Lemma.** If `(IC)` holds, there is a fixed nonzero `C in O`, depending
only on `I` and `P`, such that for every primitive parameter vector
`x in O^(d+1)` and every specialization with `P(x)!=0`, a Gaussian gcd
`G(x)` of the values `W_j(x)` satisfies

```text
G(x) divides C P(x),
|G(x)| <= |C| |P(x)|.                                (3)
```

Here primitive means that the coordinates of `x` generate the unit
ideal.  The same conclusion applies to primitive rational-integer
vectors, since for every Gaussian prime above a rational prime some
coordinate is a unit.

To prove the lemma, integral dependence on chart `k` supplies a monic
equation

```text
P_k^n+b_1 P_k^(n-1)+...+b_n=0,       b_r in I_k^r.    (4)
```

There are finitely many charts and finitely many coefficients in chosen
representations of the `b_r`.  Choose one nonzero `C in O` clearing all
their `K`-coefficient denominators.

Fix a Gaussian prime `pi` and choose `k` for which `x_k` is a `pi`-adic
unit.  Put

```text
h=min_j v_pi(W_j(x)),       q=v_pi(P(x)).
```

Evaluate (4) at the integral ratios `x_l/x_k`.  Division by `x_k^D`
does not change any valuation.  Every product of `r` generators of
`I_k` has valuation at least `rh`, while the fixed coefficient
denominators cost at most `v_pi(C)`.  Hence

```text
v_pi(b_r(x)) >= rh-v_pi(C).                           (5)
```

If `h>q+v_pi(C)`, every term after `P^n` in (4) has valuation strictly
larger than `nq`, which is impossible in a nonarchimedean field.  Thus
`h<=q+v_pi(C)` for every `pi`, proving (3), including all prime powers.

The degree equality makes (3) invariant under common parameter
scaling.  A rational parameter point may therefore be cleared by any
common denominator and then made primitive; different clearings give
the same primitive Gaussian row tuple up to a unit.

## 3. Necessity over ordinary integer parameters

The content condition is also necessary; Gaussian parameters are not
needed for the converse.

**Theorem.** The following statements are equivalent.

1. Condition `(IC)` holds.
2. There is a nonzero `C in O` such that `G(x)` divides `C P(x)` for
   every primitive `x in Z^(d+1)` with `P(x) product_j W_j(x)!=0`.

The lemma proves `1 => 2`.  We prove the converse by retaining the
field of definition and one fixed split rational prime.

Suppose `(IC)` fails on chart `k`.  The finite valuative criterion gives
a Rees valuation `nu` with

```text
a=nu(I_k) > b=nu(P_k).                               (6)
```

The existence and finiteness of these valuations, and their detection
of integral closure, is Rees's valuation theorem; see
[Rees, *Valuations Associated with Ideals (II)*](https://doi.org/10.1112/jlms/s1-31.2.221).
Let `Y` be the normalized blow-up of `I_k`, and let `E` be the prime
divisor giving `nu`.  Shrink a nonempty open subset of `E` so that `Y`
and `E` are smooth, every nonzero pulled-back row has its exact
`E`-order, and `P_k` has exact order `b`.  This also ensures that none of
these functions vanishes identically on the transverse curves below.

There is a closed point `y` in this open set whose residue field `L` is
a finite separable extension of `K`.  This uses the density of such
closed points on a smooth finite-type scheme; see
[Stacks Project, Lemmas 33.25.5--6](https://stacks.math.columbia.edu/tag/04QM).
After base change to `L`, take `y` to be rational. Since both `Y` and
the divisor `E` are smooth there, a local equation for `E` has nonzero
differential and can be completed to smooth local coordinates

```text
(t_1,...,t_d),       E=(t_1),                         (7)
```

with the coordinate map etale at `y`.  In these coordinates the
pullbacks have the form

```text
W_j/X_k^D=t_1^a_j u_j,    min_j a_j=a,
P_k=t_1^b u,                                         (8)
```

where all displayed `u_j,u` are units at `y`.  The standard etale-local
description and the corresponding isomorphism on completed local rings
are recorded in
[Stacks Project, Section 29.37](https://stacks.math.columbia.edu/tag/02GH).

Spread this finite collection of data over the integers of `L` after
inverting finitely many primes, including the denominators of `y`, the
coordinate maps, and the unit values in (8). The coordinate map remains
etale and those units remain invertible in the residue disk at every
remaining prime. By Chebotarev there is a rational prime
`p` outside that finite set which splits completely in the Galois
closure of `L/Q`; see
[Lenstra--Stevenhagen, *Chebotarev and his density theorem*](https://pub.math.leidenuniv.nl/~lenstrahw/PUBLICATIONS/1996d/art.pdf).
Because `L` contains `Q(i)`, this `p` is `1 mod 4`.  Fix the resulting
embedding `L -> Q_p` and the corresponding Gaussian prime `pi` above
`p`.

For every positive integer `n`, prescribe the etale coordinates

```text
(t_1,...,t_d)=(p^n,0,...,0).
```

Form the fiber product of this etale coordinate neighborhood with
`Spec Z_p` via the prescribed coordinate vector. Its special fiber has
the selected point of `E`. Hensel lifting gives a point `y_n in Y(Z_p)`
with exactly those coordinates. One may use the smooth lifting statement
[Stacks Project, Lemma 15.13.3](https://stacks.math.columbia.edu/tag/0D49).
Equations (8), with the units still `p`-adic units, give exactly

```text
v_pi(W_j(y_n))=n a_j,          v_pi(P(y_n))=nb.       (9)
```

The blow-down of `y_n` is a point `z_n` of the parameter chart over
`Z_p`.  Approximate its `d` affine coordinates by ordinary rational
integers modulo `p^M`, with `M>n max(b,a_1,...,a_m)`.
Polynomial evaluation is continuous, so these approximations preserve
all the finite valuations in (9).  Homogenize by taking `x_k=1`.  The
result is a primitive vector `x_n in Z^(d+1)`, all `P(x_n),W_j(x_n)` are
nonzero, and

```text
v_pi(G(x_n))-v_pi(P(x_n))=n(a-b) -> infinity.         (10)
```

Thus statement 2 cannot hold.  This proves the equivalence using
ordinary integer parameters at one fixed rational prime.  The conclusion
here is the exact divisibility failure.  A single local defect need not,
without additional global height control, make the Archimedean ratio
`|G(x_n)|/|P(x_n)|` unbounded.

This distinction occurs even in fixed homogeneous circle-product rows.
Put `S=X^2+Y^2` and take

```text
W_1=S^2,       W_2=-S^2,       P=X^4+Y^4.
```

For every nonzero ordinary integer pair, the Gaussian gcd is `S^2`
up to a unit and

```text
|G|/|P|=(X^2+Y^2)^2/(X^4+Y^4)<=2.
```

Nevertheless `(IC)` fails on `Y=1`: at `X=i`, the row ideal has
order two and `P(i,1)=2` has order zero. Exact divisibility fails on
primitive ordinary integer pairs as well. Starting with `x_1=2`, lift
inductively to `x_n^2+1=0 mod 5^n`, with `x_n=2 mod 5`. The next
digit is available because `2x_n` is a unit modulo five. Then
`v_5(G(x_n,1))>=2n` while `P(x_n,1)=2 mod 5`. No fixed Gaussian
constant can absorb this local defect, despite the uniform norm ratio.
Thus the equivalence concerns exact divisibility, not every possible
method of bounding primitive content in absolute value.

## 4. The finite valuative test

On every chart, let `nu` run through the finitely many Rees valuations
of `I_k`.  The valuative criterion used above says

```text
(IC)  iff  nu(P_k) >= nu(I_k):=min_j nu(W_j/X_k^D)   (11)
            for every such nu and every chart k.
```

Equivalently, normalize the blow-up of the row ideal.  It is enough to
compare the vanishing order of `P_k` with the row ideal along the
finitely many prime divisors that govern its pullback.  For the product
rows (1), each test is the explicit valuation-arrangement inequality

```text
nu(P_k) >= min_j sum_s[
 e_js nu(F_s/X_k)+(n_s-e_js)nu(bar(F_s)/X_k)].       (12)
```

Thus normalization of one fixed blow-up, or a standard integral-closure
computation on its affine charts, certifies the uniform arithmetic
content bound.  The exceptional rational primes introduced by reducing
the coefficient arrangement are absorbed in the single constant `C`;
they cannot acquire parameter-dependent extra exponents.

Set-theoretic control of the base locus is insufficient.  It gives only
`P in radical(I)`, and hence at best `P^N in I` for an uncontrolled
`N`.  Even with equal homogeneous degrees, take
`I=(X^2,Y^2)` and `P=XZ` in three variables.  The projective base locus
is `[0:0:1]`, where `P` vanishes, but on the chart `Z=1` the exceptional
valuation has

```text
nu(P)=1 < 2=nu(I).
```

The exponent one in (11) is exactly what is needed to compare primitive
content with a chord-scale polynomial.

Failure is therefore equally precise: a Rees divisor with positive
defect is not merely a formal obstruction.  The construction above
turns it into ordinary primitive integer vectors with unbounded content
defect.

## 5. Combining content with a chord certificate

The content lemma becomes an endpoint obstruction when it is paired
with an archimedean identity.  Fix `q` pairs of rows, repetitions allowed,
and suppose that on a parameter domain `Omega`, with `P(x)!=0`, one has

```text
product_(l=1)^q |W_r_l(x)-W_s_l(x)|
    >= c (|P(x)| rho(x))^(q/2)                       (13)
```

for a fixed `c>0`.  This is homogeneous because `P` and the rows have
the same degree.  Exact factorizations of chord differences are one way
to establish (13); the lower bound on any remaining modulus ratios is a
separate archimedean requirement.

Divide the evaluated rows by their actual Gaussian gcd `G` and put
`R=rho/|G|`.  The product of the `q` primitive chord lengths normalized
by `sqrt(R)` is at least

```text
c (|P|/|G|)^(q/2).
```

Consequently one selected chord, and hence every circle arc containing
all the rows, satisfies

```text
arc length/sqrt(R)
    >= c^(1/q) sqrt(|P|/|G|)
    >= c^(1/q)/sqrt(|C|).                            (14)
```

This proves a positive normalized-span lower bound, uniform in the
radius, for every fixed family satisfying both `(IC)` and (13).

The two requirements are logically independent.  The known four-point
Pell family in
[four_point_bonus_counterexample.md](four_point_bonus_counterexample.md)
has primitive content `G=1` while its normalized span tends to zero.
Since every nonzero Gaussian-integer value of `P` has modulus at least
one, no such `P` can satisfy (13) with fixed `c>0` along that sequence,
regardless of its integral-closure behavior.  A proposed chord
polynomial must vanish on the Pell locus or fail the balanced chord
inequality there.  The failure is archimedean, rather than hidden in
primitive content, so the criterion retains this counterexample.

## 6. The signed six-factor family is an instance

For the five signed rows of
[six_factor_primitive_content_bound.md](six_factor_primitive_content_bound.md),
take `X=(u,v,T)` and

```text
P=2u^2 v(u+v)(2u+v)(3u+v).
```

The exact local proof in that note gives the stronger divisibility

```text
G divides (1+i)^6 5 P.                              (15)
```

for every primitive integer triple with `P!=0`.  The necessity theorem
therefore proves, without a genericity assumption, that `P_k` belongs to
the integral closure of the five-row ideal on all three projective
charts.  The signed-cover checker is an exact finite certificate for the
local valuation arrangement behind this membership.

The two-chord identity has `q=2` and `c=4` in (13).  Since
`|(1+i)^6 5|=40`, formula (14) recovers exactly

```text
arc length/sqrt(R) >=1/sqrt(10).
```

## Scope

The lemma is uniform over all rational parameters of one fixed
homogeneous factor family, and the constant is effectively determined
by finitely many integral-dependence equations.  It supplies neither a
finite list of factor families containing arbitrary lattice-circle
clusters nor an extraction theorem from a varying cluster.  Those are
separate requirements for the global endpoint problem.
