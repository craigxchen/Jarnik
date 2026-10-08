# Accidental conductor in the ordered affine index

The forced endpoint-versus-interior factor and the conductor-supported
affine index are genuinely different.  This note classifies that difference
at a primitive two-level split prime.  Equal-level residue collisions can
make `gcd(Q,N)` contain the entire prime power even when the ordered gap is
zero.  The phenomenon is locally unbounded, and it occurs in actual integer
circle data.  No example here lies in the critical `C=1/2` arc class.

The forced factor itself has the index-free description proved in
[the pair--triple denominator note](ordered_pair_triple_gap_denominator.md).
In that notation

```text
K=product_p p^H_p
 =Norm(gcd_G(z_0,z_4)) Norm(gcd_G(z_1,z_2,z_3))
   /gcd(N,Norm(gcd_G(z_0,z_4)) Norm(gcd_G(z_1,z_2,z_3))).       (1)
```

Thus the ordered sampling argument measures `log K` exactly.  It does not
measure all of `log gcd(Q,N)`.

## Complete two-level formula

Fix an odd split prime `p=pi conjugate(pi)`, let `p^e || N`, and assume the
primitive allocation has exactly the levels zero and `e`.  Write its bits as
`b_0,...,b_4`.  For a same-bit pair put

```text
r_ab=v_pi(z_a-z_b)-t_a >= 0.
```

The depths within either bit class are ultrametric: the minimum of the three
depths on any triple is attained at least twice.

For any row set `W` containing at least one same-bit pair define

```text
mu(W)=min {r_ab : a,b in W and b_a=b_b}.
```

For a constant four-set define

```text
nu(W)=min_(T subset W, |T|=3) sum_({a,b} subset T) r_ab.
```

The triangle formula gives the following complete rules.

* The full five-row triangle gcd has valuation `mu(01234)`.
* A mixed four-row window `W` has gcd valuation `mu(W)`.  A constant
  four-row window has gcd valuation `e+nu(W)`.
* If rows `1,2,3` are mixed, let `a,b` be their unique repeated-bit pair.
  Then `v_p(D_123)=r_ab`.  If they are constant, then
  `v_p(D_123)=e+r_12+r_13+r_23`.

These statements include pure triangles: in a mixed set a pure triangle has
valuation `e` plus three nonnegative depths, and hence cannot beat the
available mixed triangle realizing the relevant minimum.

Substitution in

```text
v_p(Q)=v_p(g)+v_p(D_123)-v_p(g_L)-v_p(g_R)               (2)
```

now classifies every two-level allocation.  In particular, all cases with
zero ordered gap are as follows.

If the three interior rows are mixed, both windows are mixed and

```text
v_p(Q)=mu(01234)+r_ab-mu(0123)-mu(1234).                 (3)
```

If the three interior rows are constant, primitivity and `H_p=0` force
exactly one endpoint to be the opposite bit.  When row zero is that singleton,

```text
v_p(Q)=mu(1234)+(r_12+r_13+r_23)-mu(123)-nu(1234).       (4)
```

When row four is the singleton, the reversed formula is

```text
v_p(Q)=mu(0123)+(r_12+r_13+r_23)-nu(0123)-mu(123).       (5)
```

The only remaining primitive case with constant interiors has both endpoints
in the opposite class.  It is the forced pattern and has `H_p=e`; it is not
an `H_p=0` case.

## The accidental part is locally unbounded

Take the bit word `01100`.  Give its two-element class the depth
`r_12=k`, and give the three pairs in the other class depths

```text
r_03=r_04=r_34=0.
```

These are valid ultrametrics.  Equation (3) gives

```text
H_p=0,                    v_p(Q)=k.                     (6)
```

This is literal local norm-torus data, not merely a formal depth table.  Under
the split identification

```text
Q(i) tensor_Q Q_p = Q_p x Q_p,
```

choose units `u_1=1,u_2=1+p^k`, and choose three units in distinct residue
classes for rows `0,3,4`.  The two local coordinates

```text
(p^t u, p^(e-t) u^(-1))                                 (7)
```

have product `p^e`.  Their first-coordinate differences have exactly the
prescribed equal-level excesses, while unequal-level differences have the
minimum valuation.  Consequently (6) is realized on the local norm conic.
Taking `k>=e` makes `gcd(Q,N)` retain all of `p^e`, although `K` retains none.

This local family is not asserted to lift to five Gaussian integers in a
uniform short arc.  It proves only that the triangle-gcd and local norm-conic
conditions permit an arbitrarily large accidental index valuation.  Any
exclusion at `C=1/2` must use additional global or quantitative Archimedean
information.

## Actual integer-circle fixtures

There are also literal global examples with `K<gcd(Q,N)`.  Consider the
ordered primitive tuple

```text
(446,23), (439,82), (343,286), (329,302), (302,329).
```

All five norms are

```text
N=199445.
```

Its triangle data give

```text
g=6,       g_L=g_R=6,       Q=220,       gcd(Q,N)=5.
```

At `pi=1+2i` above five, the allocation and nonzero collision depth are

```text
(t_0,...,t_4)=(0,0,1,0,1),       r_13=1,
```

while `r_01=r_03=r_24=0`.  Thus (3) gives `H_5=0` and `v_5(Q)=1`.
Moreover

```text
Norm(gcd_G(z_0,z_4))=353,
Norm(gcd_G(z_1,z_2,z_3))=1,
```

so (1) gives `K=1`.  This is not a critical short-arc example: its endpoint
chord alone is longer than `(1/2)sqrt(R)`, as follows exactly from
`|z_4-z_0|^2=114372` and `16(114372)^2>N`.

The rational-quartic family supplies a second kind of actual accidental
overlap.  At parameter `n=12`, in geometric order `(3,0,4,1,2)`, its primitive
tuple has

```text
Q=2937,       gcd(Q,N)=89,       K=1.
```

For `pi=5+8i`, `v_89(N)=2` and

```text
(t_0,...,t_4)=(0,1,1,1,2),       H_89=0,
v_89(Q)=1.
```

This is a three-level accidental collision, so it lies outside the
two-level classification above while confirming that the distinction is not
only local.  The family lies on an endpoint arc of some fixed normalized
constant; it gives no `C=1/2` counterexample and its fixed polynomial
resultants do not establish unbounded accidental exponent.

The [exact checker](check_ordered_index_accidental_conductor.py) exhausts
two-level bit patterns on many finite ultrametric residue tables, verifies
(2)--(5), checks the local norm-torus construction at finite precision, and
checks both global fixtures from their actual Gaussian coordinates.
