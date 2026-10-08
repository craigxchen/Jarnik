# Exact Newton support for the eight-point cut orders

This gives the exact Newton-support formulation of the weighted-certificate
question in the continuation status. The multiplicity bound in Section 4
now rules out a strict generic cut-order gain in every degree. It does not
prove a uniform arc bound.

## 1. Polynomial coordinates on the phase torus

Treat the eight unit-circle phases as nonzero complex variables `z_i`.
The complex moment matrix has rows with exponents `-2,-1,0,1,2`.
For a five-subset `J`, its minor is

```text
V_J / product_(j in J) z_j^2,
V_J=product_(j<k in J)(z_k-z_j).
```

After multiplying all minors by the same nonzero factor, its polynomial
coordinates, indexed by the complementary triple `I`, are

```text
R_I(z)=(product_(i in I) z_i^2) V_(I^c).             (1)
```

The fixed complementary-minor signs can be incorporated into the
coordinate names. They affect no support or valuation argument below.
Every `R_I` has total degree 16 and individual variable degree at most
four. Since the Vandermonde contains ten factors,

```text
R_I(1/z)=R_I(z)/product_i z_i^4.                     (2)
```

This is the common inversion symmetry of the first two harmonics.

Let `F` be a homogeneous polynomial of degree `d` in these coordinates,
and assume its pullback `f(z)=F(R(z))` is not identically zero. Then
`f` has total degree `16d`, individual degrees at most `4d`, and

```text
f(1/z)=f(z)/product_i z_i^(4d).                      (3)
```

Thus its exponent support is centrally symmetric under
`alpha -> 4d*1-alpha`. Put `beta=alpha-2d*1`; all these centered
exponents have coordinate sum zero, and their convex hull is centrally
symmetric about zero.

## 2. Every generic cut order is a support minimum

For a subset `S`, substitute `z_i=epsilon x_i` on `S` and `z_j=y_j`
off `S`, with generic nonzero parameters. If `s=|S|` and
`a=|I intersect S|`, the raw coordinate has order

```text
w_S(I)=2a+binomial(s-a,2).                          (4)
```

The second term counts the pairs of small phases in the complementary
five-subset. For minority sizes one through four, the common minimum is

```text
m_1=0, m_2=1, m_3=3, m_4=5.                       (5)
```

Divide all coordinates by `epsilon^(m_s)` to make their minimum
valuation zero. This is the primitive local normalization. After all
exact cancellations in `f` have already been performed, the generic
order of `F` on this normalized cut is exactly

```text
nu_S(F)=min_(alpha in supp f) sum_(i in S)alpha_i-d*m_s.  (6)
```

Indeed different surviving monomials in the coefficient of the least
power of `epsilon` remain different monomials in the independent `x,y`
parameters. That coefficient is a nonzero polynomial. Formula (6) also
gives a universal lower order for arbitrary unit specializations; special
cancellations can increase the order. Fixed coefficient denominators and
exceptional primes must be retained when applying this to integers.

Inversion identifies the complementary cut. Use all eight singleton
cuts, 28 pairs, 56 triples and one representative of each of the 35
balanced cuts. Equivalently, sum over all nonempty proper subsets and
divide by two, with their corresponding baselines.

## 3. The complete order budget

For each chosen cut write

```text
h_S=max_(beta in supp f-2d*1) sum_(i in S) beta_i.
```

Central symmetry and (6) give

```text
nu_S(F)=2d*|S|-d*m_(|S|)-h_S.
```

Therefore

```text
sum_S nu_S(F)=373d-sum_S h_S,                       (7)
373=8*2+28*3+56*3+35*3.
```

The fair eight-point lattice height has logarithm `53w+o(w)`. Hence a
degree-`d` pullback could have cut-order sum exceeding the archimedean
budget precisely if

```text
sum_S h_S < 320d.                                  (8)
```

For a single coordinate `R_I`, its Vandermonde support gives
`sum h_S=320` and `sum nu_S=53`. For a nonzero product of `d`
coordinates, support functions add under multiplication, so equality
remains exact with `320d` and `53d`. The theorem below shows that arbitrary
cancellations in sums cannot produce (8) either. The earlier proposed
strict weighted-order certificate therefore cannot exclude a fair
configuration in this way.

## 4. The all-polynomial multiplicity obstruction

At the diagonal point `z_1=...=z_8=1`, each coordinate (1) vanishes
to order exactly ten. Thus every degree-`d` pullback `f` has
vanishing order at least `10d`, unless it is identically zero.

The general support-multiplicity theorem proves

```text
mult_1(f) <= (sum_S h_S)/32,                        (9)
```

and hence rules out (8) in every degree. Coordinate products attain
equality. The proof applies to every nonzero homogeneous polynomial in
the stated individual-degree box with reciprocal support, including the
complete pullback after cancellations between different row weights.
It uses separate polarization, a signed coefficient-indicator transform,
and an elementary distance-to-support bound on a rectangular grid.
See [the proof](polarized_coefficient_support_multiplicity.md) and
[its independent audit](polarized_support_multiplicity_independent_audit.md).

Combining (7) and (9) gives the sharp ceiling `sum_S nu_S(F)<=53d`.
This concerns universally prescribed generic cut orders. Extra divisibility
at actual arithmetic values and numerical height comparisons are not
bounded by this conclusion.
