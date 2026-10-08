# Direct value-sharing and polynomial `abc` estimates retain degree dependence

Consider distinct real polynomials `F_1,...,F_m` of one even degree `D`, with

```text
deg(F_i-F_j)=D,
F_i-F_j divides F_i F_j+1,
deg lcm_C[t](F_i+i)<=2D.                               (1)
```

This note audits the natural attempts to apply polynomial value-sharing
theorems and Mason--Stothers.  The hypotheses encode pair-dependent sharing
of half of one fiber, rather than sharing of fixed complete fibers.  The
direct `abc` identities have enough distinct zeros that their degree bounds
are weaker than the identities already known.  No degree-independent bound
on `m` follows from these standard reductions.

This is an applicability result, not a construction with arbitrarily many
rows and not a theorem that no stronger polynomial argument exists.

## 1. The exact sharing pattern

Put

```text
h_i=F_i+i,       G_ij=gcd(h_i,h_j).
```

The divisibility in (1), including multiplicities, gives

```text
deg G_ij=D/2,
F_i-F_j=c_ij G_ij bar(G_ij),       c_ij in R*.        (2)
```

Thus `F_i` and `F_j` share a degree-`D/2` sub-divisor of the fiber over `-i`.
They do not generally share the complete degree-`D` fiber.  Moreover the
shared divisor `G_ij` changes with the pair.

All polynomial maps `F_i:P^1->P^1` share the point over infinity, with total
ramification `D`, but this is only one shared target.  At the finite targets
`-i` and `+i`, equation (2) supplies pairwise overlaps rather than equality of
preimage divisors.  Consequently uniqueness theorems whose hypotheses say
that two meromorphic functions share several fixed values, ignoring or
counting multiplicity, do not apply: even for a single pair, `-i` is only a
half-shared fiber, and across the family there is no common fixed list of
complete shared fibers.

The degree-four eight-block example makes the distinction concrete.  Its four
rows have pairwise intersections of size two in their `-i` fibers, while the
eight-root union has four private roots and four roots incident to triples.
No two `-i` fibers are equal, yet all conditions in (1) hold with
`deg lcm(h_i)=8=2D`.

## 2. The direct Mason--Stothers identity

For each pair define

```text
Q_ij=(F_i^2+1)/(F_i-F_j).
```

Then

```text
F_i^2+1=(F_i-F_j)Q_ij,       deg Q_ij=D.              (3)
```

The additive Mason triple is

```text
F_i^2+1-F_i^2=1.                                      (4)
```

Its three terms have gcd one.  Mason--Stothers therefore gives

```text
2D<=deg rad(F_i(F_i^2+1))-1<=3D-1.                    (5)
```

The factorization in (3) does not add the roots of `F_i-F_j` and `Q_ij` a
second time: together they are exactly the roots of `F_i^2+1`, with
multiplicity.  Counting just those factors as the radical for Mason would
incorrectly omit the roots of `F_i`.  The correct bound (5) is compatible with
every `D>=1` and contains no occurrence of `m`. A contradiction from this
identity would require an upper bound below `2D+1` for this radical; no
such bound has been deduced from (1) here.

Rewriting (3) as a rational-function unit equation does not change this
count.  Its set of allowed places includes the moving roots of `G_ij` and of
the complementary half-fibers, so the size of the `S`-set grows with `D` and
depends on `(i,j)`.  Function-field unit-equation bounds consequently remain
degree- or support-dependent.

## 3. Summing pair identities does not create a fixed small radical

For three rows, telescoping gives

```text
c_ij G_ij bar(G_ij)+c_jk G_jk bar(G_jk)
 =c_ik G_ik bar(G_ik).                                (5)
```

This is an `abc` identity among three degree-`D` real polynomials.  Before
applying Mason, divide all three terms by their common gcd `W`; the resulting
terms are pairwise coprime, since a common divisor of two also divides the
third.  If `w=deg W`, all three degrees are exactly `D-w`.  When `D>w`,
their combined radical degree is at most `3(D-w)`, so Mason gives only
`D-w<=3(D-w)-1`.  When `D=w`, all three terms are constants and the
nonconstant hypothesis of Mason fails.  Neither case gives a contradiction.  The
lcm hypothesis controls the union of the oriented factors
`G_ij` inside the `h_i`, but (5) also contains their conjugates.  Since
`h_i` and `bar(h_i)` are coprime for each row, no uncharged cancellation turns
that union into a bounded set of places.

A proposed many-term `abc` argument would have to track both the linear
dependencies among its terms and their common factors. The available support
bound has degree proportional to `D`. No degree-independent many-term
conclusion is proved or ruled out by the two direct calculations above.

## 4. The known linear moment bound from real coefficients

After division by its real leading coefficient, `h_i` is monic and every
nonconstant coefficient is real.  Equivalently, if its roots are counted with
multiplicity, their power sums satisfy

```text
Im sum_(alpha root of h_i) alpha^r=0,
1<=r<=D-1.                                            (6)
```

These are genuine constraints absent from an abstract block design.  On the
layer-incidence matrix of the lcm, (6) supplies a real nullspace of dimension
at least `D/2`.  Combined with `deg lcm(h_i)<=2D`, this is exactly the mechanism
behind the existing bound

```text
m<=3D/2.                                              (7)
```

Standard value-sharing theorems do not exploit (6), because they work with
preimage divisors of target values rather than the reality of all intermediate
coefficients.  Conversely, using only the linear moment equations (6) gives a
nullity proportional to `D`, while the available root-coordinate space also
has dimension proportional to `D`.  The quartic four-row family shows that
these equations do not collapse every nontrivial incidence pattern.

## 5. Why this is not a polynomial `D(1)` tuple

A polynomial `D(1)` tuple asks that products of two members, plus one, be
squares.  Condition (1) instead asks that the pair-dependent difference
`F_i-F_j` divide `F_iF_j+1`.  The quotient need not be a square, and the
divisor changes with the pair.  Passing to the identity

```text
F_iF_j+1=(F_i-F_j)H_ij                                (8)
```

does not supply a common squareclass or a fixed auxiliary polynomial shared by
all pairs.  Results whose finiteness mechanism uses those features cannot be
imported without proving an additional reduction, and (1) alone provides no
such reduction.

## 6. Precise remaining possibility

The standard routes leave one genuinely stronger polynomial question open:
can the pair-dependent half-fiber overlaps (2), the global support bound
`deg L<=2D`, and all real moment equations (6) be combined nonlinearly to bound
`m` independently of `D`? The direct pairwise and triple Mason estimates
audited here do not supply the needed compatibility among changing overlap
divisors.

Even a positive answer would apply only to polynomial families satisfying one
common finite-degree model.  An arbitrary sequence of integer endpoint
clusters does not automatically specialize from polynomials of uniformly
bounded degree, and the degree `D` here is itself allowed to grow.  A transfer
to the original lattice problem would therefore require a separate bounded-
degree or controlled-specialization theorem.
