# Sharp separation of multiplicative rectangles

This is a finite theorem for arbitrary Gaussian lattice points, without
prime-pattern, unit-class, or Pell-template assumptions. It does not prove
the uniform endpoint count.

The later [higher product rigidity theorem](short_arc_product_rigidity.md)
excludes product collisions through length eight at `C=1/2`, using
equal-radius polynomial identities. The sharp two-factor constant
proved here remains stronger for two-factor relations.

## Theorem

Let nonzero Gaussian integers `z0,z1,z2,z3` have common modulus `R` and
lie in a circular arc of angular width `Delta < pi`. Suppose

```text
z0 z3 = z1 z2,       z0 != z1,       z0 != z2.
```

The four entries need not all be distinct. Then

```text
Delta sqrt(R) >= 2 sqrt(2).
```

Consequently any such arc of length less than `2 sqrt(2) sqrt(R)` is
multiplicatively Sidon: equality between products of two points, allowing
repetitions, implies equality of the two unordered pairs.

## Proof

The Euclidean domain `Z[i]` gives a rank-one factorization

```text
z0 = a c,   z1 = a d,   z2 = b c,   z3 = b d
```

with all four factors nonzero Gaussian integers. Explicitly, take
`c = gcd(z0,z2)`, write `z0 = a c`, `z2 = b c`, and cancel `c` in the
product identity. Since `a,b` are coprime, `a` divides `z1`; writing
`z1 = a d` and cancelling `a` gives `z3 = b d`. Units in the chosen gcd
do not affect this argument.

Equal moduli imply `|a|=|b|=A`, `|c|=|d|=B`, and `AB=R`.
The two inequalities in the hypothesis give `a != b` and `c != d`.
For distinct equal-modulus Gaussian integers `a,b`,

```text
|a-b|^2 = 2 (|a|^2 - Re(a conjugate(b)))
```

is a positive even integer. Thus `|a-b| >= sqrt(2)`, and likewise
`|c-d| >= sqrt(2)`.

Choose real arguments `phi_j` of all four points in the given interval.
The product identity says `phi0+phi3-phi1-phi2` is an integer multiple of
`2 pi`. Its absolute value is at most `2 Delta < 2 pi`, so it is zero.
Put `u=phi1-phi0`, `v=phi2-phi0`. The four arguments are therefore

```text
phi0, phi0+u, phi0+v, phi0+u+v.
```

Their span is exactly `|u|+|v|`, regardless of the signs of `u,v`.
Moreover the chord bound `|exp(i t)-1| <= |t|` gives

```text
|u| >= |d-c|/B >= sqrt(2)/B,
|v| >= |b-a|/A >= sqrt(2)/A.
```

It follows that

```text
Delta >= sqrt(2) (1/A + 1/B) >= 2 sqrt(2)/sqrt(AB).
```

For the Sidon conclusion, write an equality of pair products as
`z0 z3 = z1 z2`. If either excluded equality `z0=z1` or `z0=z2`
holds, cancellation identifies the unordered pairs. Otherwise the
proved lower bound applies. This also handles repeated factors.

## Sharpness

For positive integers `n`, put

```text
a = n + (n+1)i,     b = (n+1) + n i,
R = |a|^2 = |b|^2 = 2n^2+2n+1.
```

The three distinct Gaussian points `a^2, ab, b^2` have radius `R`
and satisfy `a^2 b^2 = (ab)^2`. With `z0=a^2`, `z1=z2=ab`,
`z3=b^2`, their angular span is

```text
Delta = 4 arcsin(1/sqrt(2R)).
```

Indeed `|a-b|=sqrt(2)`, so the angular difference between `a,b` is
`2 arcsin(1/sqrt(2R))`, and squaring doubles the span. Consequently
`Delta sqrt(R)` tends to `2 sqrt(2)` as `n` tends to infinity.
No larger universal constant is possible in the stated theorem.

## A stronger count for fixed templates on sufficiently short arcs

Under the hypotheses of the
[fixed quadratic-unit template theorem](quadratic_unit_template_bound.md),
if `C < 2 sqrt(2)`, there are at most **six** persistent rows.
This sharpens that theorem's count of eight for this range of `C`.

Indeed, its Section 7.4 proves that seven rows force active width at
most four. Width less than four allows at most four rows by its cube
bound. At width four, any nonbinary width partition allows at most six
mutually nonadjacent vectors, as proved there. Thus seven rows would
give seven distinct vertices of the binary four-cube.

Partition its sixteen vertices into eight complementary pairs. This
partition alone does not force a collision among seven vertices; the
nonadjacency constraint is essential. Here is a direct way to retain it.
For a subset containing vertices of both parities, suppose the smaller
parity part has `t` vertices. In a putative seven-vertex independent set,
`1 <= t <= 3`. Any one vertex has four neighbors. Any two distinct
same-parity vertices have at most two common neighbors, so their union
of neighbors has at least six vertices. Any three same-parity vertices
have at least seven neighbors: if one pair is antipodal their union
already has eight; otherwise translate one vertex to zero. The other
two have weight two, intersect in one coordinate, and their three
neighborhoods have pairwise intersections of size two and common
intersection of size one, giving union size seven. Therefore the two
parity parts together have at most `5,4,4` vertices for `t=1,2,3`,
respectively. Seven independent vertices must all have the same parity.

Each parity class of the four-cube consists of four complementary
pairs. Seven of its vertices contain at least three full complementary
pairs. Choose two, with exponent vectors `v0,v3` and `v1,v2`.
Their coordinatewise sums agree. In the original template, this makes

```text
z0(n) z3(n) / (z1(n) z2(n)) = c0 c3 / (c1 c2),
```

a fixed complex number of modulus one. Endpoint angular widths tend
to zero, so the left side tends to one; the fixed quotient is therefore
exactly one. Division by the common Gaussian gcd preserves this
identity. All four selected rows are distinct. For sufficiently large
indices their primitive angular width is less than `pi`, and the
rectangle theorem forces `C >= 2 sqrt(2)`, a contradiction.

This argument also explains why the prefactors cannot simply be
ignored: their cancellation follows from an exact exponent identity
and persistence along shrinking angular widths. It is not valid for
arbitrary isolated configurations with varying prefactors.

## Counting consequence and remaining gap

At `C=1/2`, nonempty arcs have `R>=1` and angular width at most
`1/(2 sqrt(R)) < pi`. If an arc contains `M` lattice points, its
`M(M+1)/2` unordered pair products are all distinct. Its ordered
multiplicative energy is exactly `2M^2-M`.

These products lie on a circle of radius `R^2` in an arc of length
at most `2C R^(3/2)`. Relative to the new radius, this is exponent
`3/4`, not exponent `1/2`; applying an endpoint estimate to that arc
would be invalid. Nor does the Sidon property alone bound cardinality:
abstract torsion-free groups have arbitrarily large Sidon subsets.

In particular, this removes exact two-versus-two multiplicative
collisions from a possible counterexample, but supplies neither the
intrinsic eight-point conductor bonus nor a uniform bound on the
remaining collision-free configurations. The proof above is in prose;
no Lean formalization of this theorem is claimed.
