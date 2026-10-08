# Triangle squareclasses at good unbalanced cuts

This note computes the ten triangle-field squareclasses at an exact
singleton or pair allocation cut.  A minority-pair prime splits four of the
triangle fields locally, but the other six fields and all fields at a
singleton cut can have variable local class.  These facts retain the full
prime-power height and do not produce a new size obstruction.

Fix six distinct local nodes satisfying

```text
sum_i u_i r_i^j=0,             -2<=j<=2,             (1)
```

where `u` is the ordinary primitive central vector.  Let `p` be an odd
rational prime split in `Q(i)`, choose one orientation above `p`, and assume
`p` does not divide `2L`.  At such a good circle prime, unequal allocations
have the exact minimum-difference valuation and equal-allocation nodes have
distinct reduced unit parts.

For a balanced triple `I`, write

```text
S_I=sum_(i in I)u_i,
D_I=-S_I product_(i in I)u_i.                         (2)
```

The exact triangle formula in
[`six_point_triangle_splitting_fields.md`](six_point_triangle_splitting_fields.md)
says that `D_I` is the squareclass of the discriminant of the binary
triangle form `M_I`.  Notice that a cut prime divides some of the `u_i`, so
it is an exception to the simple finite-field reduction statement for the
raw cleared frame.  We compute its characteristic-zero class in `Q_p`
directly; no simple-root assertion modulo `p` is being imported.

## 1. The local central vector

For six nodes the unique central vector has the barycentric expression

```text
u_i=c r_i^2/product_(j!=i)(r_i-r_j),                 (3)
```

with one common scalar `c`.  This is immediate because `1/P'(r_i)`
annihilates powers zero through four, and multiplication by `r_i^2`
shifts those equations to (1).  At the exact good cuts below, primitive
normalization makes `c` a `p`-adic unit.

Suppose first that the minority pair is `T={1,2}` and

```text
r_1=p^a x_1, r_2=p^a x_2,       a>0,
r_j=y_j for j in M,              M={3,4,5,6},        (4)
```

where the `x` residues are distinct nonzero units and the four `y` residues
are distinct nonzero units.  Put

```text
Y=product_(j in M)y_j,       P_M(X)=product_(j in M)(X-y_j).
```

Reduction of (3) gives

```text
u_1/p^a = c x_1^2/[(x_1-x_2)Y],
u_2/p^a =-c x_2^2/[(x_1-x_2)Y],
u_j       = c/P_M'(y_j)                    modulo p. (5)
```

In particular, the minority coefficients have valuation `a`, the majority
coefficients are units, and

```text
-(u_1/p^a)(u_2/p^a) is a square modulo p.            (6)
```

For a singleton `T={1}`, write `r_1=p^a x` and let `M` contain the five
majority nodes `y_j`.  Formula (3) instead gives

```text
u_1/p^(2a)=-c x^2/Y,
u_j=c y_j/P_M'(y_j)                         modulo p, (7)
```

where now `Y=product_M y_j`.  Thus the singleton coefficient has valuation
`2a` and all five majority coefficients are units.  Equations (5) and (7)
are the residue refinements of the known valuation formula; the rectangular
Vandermonde kernels remain one-dimensional and have only unit entries.

## 2. A minority pair guarantees four local splittings

There are four balanced partitions for which the two minority labels lie
on the same side.  If `I=T union {j}`, then `S_I` reduces to the majority
unit `u_j`.  Equations (2) and (6) show that

```text
v_p(D_I)=2a,          p^(-2a)D_I is a square mod p.  (8)
```

If instead `I` consists of three majority labels and omits the fourth one
`k`, then `S_I=-S_(I^c)` reduces to `-u_k`.  Moreover

```text
product_(j in M) 1/P_M'(y_j)
 =1/product_(j<k)(y_j-y_k)^2,                        (9)
```

because there are four majority nodes.  Hence

```text
v_p(D_I)=0,           D_I is a square mod p.         (10)
```

For odd `p`, (8)--(10) say that `D_I` is a full square in `Q_p`, including
the even valuation in (8).  Thus the four noncrossing triangle quadratics
split over `Q_p`.

No analogous assertion holds for a crossing partition.  Take
`I={1,j,k}`, let `{l,m}` be the other two majority labels, and put

```text
b_I=v_p(S_I).
```

Then exactly one minority coefficient occurs in (2), so

```text
v_p(D_I)=a+b_I.                                      (11)
```

The majority cofactor kernel does not force `b_I=0`.  If
`lambda_j=1/P_M'(y_j)`, direct subtraction gives

```text
lambda_j+lambda_k=0 mod p
 iff y_j+y_k=y_l+y_m mod p.                          (12)
```

This additive symmetry is compatible with four distinct majority nodes.
When `b_I=0`, the remaining unit squareclass is, up to a square,

```text
p^(-a)D_{1jk} =
 -(y_j+y_k-y_l-y_m)/[(x_1-x_2)Y]
                                      in F_p^*/(F_p^*)^2. (13)
```

It is not fixed.  Thus a generic crossing class has valuation parity `a`,
but cancellation in `S_I` changes it to `a+b_I`; the value of `b_I` and the
unit squareclass both depend on higher residue data.

The complementary triple always has the same squareclass.  Globally this
already follows from hyperbolicity:

```text
D_I D_(I^c)=-S_I^2 product_i u_i,
-product_i u_i is a rational square.                 (14)
```

## 3. Singleton cuts leave every individual class variable

At a singleton cut, a triple containing the minority label has

```text
v_p(D_I)=2a+v_p(S_I),                                (15)
```

whereas its complementary majority triple has valuation `v_p(S_I)`.
These parities agree, as required by (14).  But the sum of the two majority
unit coefficients occurring in `S_I mod p` can vanish, and when it does not
vanish the unit in (2) can be either a square or a nonsquare.  Formula (7)
also verifies the local version of (14): using

```text
product_(j in M) P_M'(y_j)
 =product_(j<k)(y_j-y_k)^2
```

for five majority nodes shows that `-product_i u_i` is a square modulo `p`
after removing the even valuation `2a`.  This relates complementary fields
but does not force either one to split.

## 4. Examples, bad-prime cost, and scope

The local central equations themselves realize all the variable cases.
At the Gaussian-split prime `p=13`, the pair-cut nodes

```text
(r_1,...,r_6)=(13,26,1,2,3,4)
```

have primitive central vector

```text
(-1495,234,-759,10350,-26730,18400).
```

The four noncrossing classes are `Q_13`-squares.  Among crossing triples,
there are examples with odd discriminant valuation, with even valuation and
nonsquare unit, and with a square local class.  The singleton-cut nodes

```text
(13,1,2,3,4,7)
```

have primitive central vector

```text
(169,-165,2592,-8019,7040,-1617)
```

and likewise give ramified, unramified nonsplit, and split classes
among their triangle fields.  The accompanying checker verifies these
claims from exact rational barycentric weights.

Prime powers are fully retained in (8), (11), and (15): the cut height is
the actual exponent `a`, and the extra subsum exponent is the actual `b_I`.
The clean formulas require an exact two-level profile and the good
difference rule.  At primes supported on `2L`, or when private factors
perturb the profile, the existing central-vector estimates give additive
valuation errors whose total size is controlled.  They do **not** preserve
valuation parity or residue squareclass.  A single private factor at a prime
supporting a cut of arbitrarily large height can destroy the conclusion,
and discarding that prime would discard the entire cut height.  Therefore
no `o(w)` transfer of the squareclass statements is asserted here.

A pair cut therefore makes its primes split simultaneously in four moving
quadratic fields.  This is a genuine local correlation, but it yields no
new deterministic size or divisibility gain.  The even factor `p^(2a)` in
(8) is already the product of the two known minority coefficient factors;
(11) and (15) merely restate the coefficient and subsum valuations.  Common
splitting in finitely many quadratic fields permits arbitrarily large prime
powers and a positive-density set of prime supports, while the fields move
with `u`.  Any growth conclusion would require an additional global
distribution or height input, not the rectangular Vandermonde reduction or
quadratic reciprocity alone.
