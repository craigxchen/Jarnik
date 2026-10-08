# Even cut occupancy does not divide the direct Runge remainder

This note gives an exact local formula for the direct-source monic
model, valid in every degree. At a mixed binary source cut, even
occupancy does not force even one factor of its odd prime into the
cleared Runge remainder. The formula identifies an additional residue
congruence that would be needed to obtain such divisibility.

This is an obstruction to deriving divisibility from the indicated
local cut and squareclass data. It does not exclude extra congruences
forced by a large actual endpoint tuple. The finite global fixture below
is a literal primitive binary circle tuple, but is not asserted to meet
the endpoint constant two.

## 1. The exact remainder congruence

Use an even selection of `2k` distinct roots from the
[direct source model](direct_projection_monic_height_bridge.md):

```text
Z=2N,       d_i=2N/f_i,
c_i=d_i b_i^2,       Z-c_i=d_i a_i^2,
product d_i is a square.
```

Let `P(X)=product_i(X-c_i)`, and define the polynomial part of its square
root at infinity and its cleared remainder by

```text
s_a=[t^a] product_i(1-c_i t)^(1/2),
S(X)=sum_(a=0)^k s_a X^(k-a),
q=2^(2k-1),       R(X)=q^2(P(X)-S(X)^2) in Z[X].    (1)
```

The coefficient denominators and notation are those of the
[even-product lemma](runge_even_product_square_bound.md).

Suppose `p` is odd, `p^e|Z`, and exactly `j>=1` of the roots are
divisible by `p^e`, while the other roots are `p`-units. Let `U` be
the list of these unit roots and set

```text
F_k(U)=[t^k] product_(u in U)(1-u t)^(1/2).
```

Then the following congruence is exact:

```text
R(Z) = -[q F_k(U)]^2  (mod p^e).                   (2)
```

Indeed `P(Z)` is divisible by `p^e`. In `S(Z)`, all positive powers
of `Z` disappear modulo `p^e`, leaving `s_k`. Setting the divisible
roots equal to zero in this coefficient leaves precisely `F_k(U)`.
All coefficient denominators are powers of two, so this reduction is
valid at every odd prime, independently of its size or of `e`.

In particular, if `F_k(U)` is a `p`-unit, then

```text
v_p(R(Z))=0.                                      (3)
```

No content division of the integer remainder polynomial can recover a
factor `p` in this case. The integer square root of `P(Z)` can be
divisible by a large power of `p` while `qS(Z)` remains a `p`-unit;
there is consequently no common enlarged integer lattice for both
quantities in the usual rounding step.

## 2. Binary allocations and both orientations

At a binary source prime `p^e||N`, fix either Gaussian orientation and
compare each selected allocation with the actual anchor. A row matching
the anchor has `v_p(d_i)=e`, and hence `p^e|c_i`. An opposite row has
`v_p(f_i)=e`, `v_p(d_i)=0`, and ordinary primitivity makes both `a_i`
and `b_i` units at `p`, so `c_i` is a unit. Therefore `j` in (2) is
exactly the number of selected rows matching the anchor.

The squareclass relation gives

```text
v_p(product d_i)=e j is even.                      (4)
```

Thus odd `e` forces even `j`; even `e` imposes no occupancy parity.
Changing Gaussian orientation complements every allocation, including
the anchor, and leaves matching, opposite, and (2)--(4) unchanged.
The half-angle parity factor is already retained in `Z=2N` and `d_i`.
The prime two is not a source-conductor prime for a Gaussian-primitive
tuple; its coefficient clearing is the explicit factor `q` in (1).

For example, if `k=2,j=2`, write the two unit roots as `u,v`. Then

```text
F_2(u,v)=-(u-v)^2/8,
R(Z)=-(u-v)^4  (mod p^e).                          (5)
```

The even occupancy is compatible with a unit remainder whenever the
two opposite roots are distinct modulo `p`.

## 3. A variable-degree obstruction, including local square conditions

For every `k>=2` and every mixed even occupancy

```text
2<=j<=2k-2,       j even,
```

put `l=2k-j`. The polynomial `F_k(u_1,...,u_l)` is nonzero: its
coefficient of `u_1^k` is `(-1)^k binom(1/2,k)`. This coefficient is
a unit modulo every prime `p>2k`.

More strongly, for any split prime `p>8k`, choose `l-1` distinct
nonzero quadratic residues for `u_2,...,u_l`. As a polynomial in
`u_1`, `F_k` has degree `k` and nonzero leading coefficient. There
are at most `k` roots and at most `l-1` already used values to avoid.
Since `(p-1)/2>k+l-1`, another nonzero quadratic residue can be chosen
so that all the unit roots are distinct and `F_k(U)` is a unit.
Thus the obstruction persists for arbitrary `k`, even with the stronger
restriction that every unit root is a quadratic residue.

These local data can satisfy all of the individual norm equations and
the square-product relation; the unit roots are not independent of
those equations by assumption. To check compatibility, work over `Q_p`
with `p=1 mod 4` and any integer `Z` of valuation `e`. For an opposite
row choose a unit square `u=b^2`. Since `-u` is a nonzero square modulo
`p`, `a^2=Z-u` has a `p`-adic solution. Take `d=1`, so

```text
a^2+b^2=Z,       c=u,       Z-c=a^2.
```

For each matching row choose a distinct rationally parametrized
`p`-adic point `a^2+b^2=1` and take `d=Z`, so `c=Z b^2`.
The product of all the labels is exactly `Z^j`, a square because `j`
is even. The matching roots can be made distinct and nonzero.

At an opposite row, one of `a+ib,a-ib` has valuation `e` and the
other valuation zero. Choosing the sign of `a` selects either
orientation, so all opposite rows may use the same orientation.
Matching rows have valuation zero at both. Multiplying the phases by a
local anchor with complementary valuations `0,e` gives integral local
circle points with binary allocations. These are local compatibility
data, not a constructed large integral endpoint tuple.

If instead all `2k` roots are divisible by `p^e`, homogeneity does give
`p^(2ke)|R(Z)`: simultaneous scaling of `Z` and every root scales the
remainder by its total degree `2k`. This full-occupancy case does not
extend to the mixed even occupancies above.

## 4. A literal primitive binary fixture

Take the source half-angle rows

```text
1,       3+2i,       4+7i,       10+11i,       9+32i.
```

They come from the products `1,D,AD,BD,ABD` with
`A=2+i,D=3+2i,B=4+i`, of pairwise coprime norms `5,13,17`.
One exact primitive physical realization is

```text
N=1105,
(32,9), (4,33), (-24,23), (-12,31), (-32,9).
```

Every point has ordinary content one. The actual least circle norm is
`1105`, and every source allocation is binary. At the first point as
anchor, the four nonanchor rows have

```text
c=(680,1666,1210,2048),
d=(170,34,10,2),       product d=340^2,
Z=2210,
S(X)=X^2-2802X+1701512,
P(Z)=367200^2,
R(Z)=-1264902967296=4  (mod 5).
```

Exactly two selected roots are divisible by five. This verifies that
the local failure can occur on a global primitive binary circle tuple
with an exact even squareclass relation. This five-point tuple is not
claimed to satisfy the endpoint constant two.

## 5. The precise extra condition still needed

To obtain a factor of `p` from a mixed occupancy using this remainder,
one must prove the extra congruence `F_k(U)=0 mod p`. Higher divisibility
requires corresponding higher precision and cancellation information
in (2). Occupancy parity, the two Gaussian orientations, and the local
norm-square relations do not imply even the first congruence.

This is more specific than the generic fixed-polynomial ceiling in
[the Hermitian polynomial barrier](hermitian_polynomial_barrier.md):
it computes the particular moving-degree remainder and exhibits its
nonvanishing at mixed even cuts. The
[short squareclass separation theorem](norm_squareclass_separation.md)
uses additional global cut and angle information, and is not contradicted.
Whether many actual endpoint rows force the required extra residue
congruences remains open; no improvement of the general bound follows.

The [checker](check_direct_monic_remainder_local_obstruction.py) verifies
(2) at prime powers, variable-degree unit examples, full-occupancy
homogeneity, and the exact global fixture.
