# A rational four-row quartic construction and five lattice points

The eight supports

```text
{1},{2},{3},{4},{123},{124},{134},{234}
```

admit four simultaneous quartic constant-imaginary numerators over
`Q(i)`. Their primitive pair numerators also have constant imaginary
parts, and their polynomial lcm has degree eight. Clearing the eight
linear factors gives five integer lattice points on circles of
radius comparable to `n^8`, with an arc containing them of length
`O(n^4)`.

This is a fixed five-point family with a fixed implicit constant;
it does not contradict, or prove, a uniform bound on point counts.
It does show that the arithmetic obstruction for the symmetric
real family does not extend to all realizations of this support
design.

## 1. Exact rational root data

Set

```text
r_1=2,       r_2=-51/22,       r_3=11/3,       r_4=-1/17,
z_j=(1-r_j^2+2i r_j)/(1+r_j^2).                            (1)
```

The four `z_j` are the missing-row triple roots:

```text
z_1=-3/5+4i/5,
z_2=-2117/3085-2244i/3085,
z_3=-56/65+33i/65,
z_4=144/145-17i/145.                                      (2)
```

Define two rational numbers and the four private roots by

```text
B=-sum_j Im z_j=-107712/232609,
C=-sum_j Re z_j Im z_j=28929719808/54106946881,

w_j=(Re z_j Im z_j+C)/(Im z_j+B)+i(Im z_j+B).               (3)
```

All denominators in (3) are nonzero. The monic rows are

```text
h_j(t)=(t-w_j) product_(l!=j)(t-z_l).                      (4)
```

The integer coefficient table in Section 4 also specifies every
root as `(A+iB)/d`, avoiding long separate fractions for the `w_j`.

## 2. Complete exact verification

Direct rational multiplication of (4) gives real coefficients of
`t^4,t^3,t^2,t`, and the following constant imaginary parts:

| Row | `tau_j=Im h_j` |
|---|---|
| 1 | `16334372208/45504135625` |
| 2 | `206448/235625` |
| 3 | `-5077968/11183125` |
| 4 | `5592048/185485625` |

These four constants are nonzero and distinct. All eight roots
are nonreal and distinct, and the entire set is disjoint from its
complex conjugate. Hence, up to nonzero constants,

```text
G_ij=gcd(h_i,h_j)=product_(l notin {i,j})(t-z_l),
deg G_ij=2,
L=lcm(h_1,h_2,h_3,h_4)=product_l(t-z_l)(t-w_l),
deg L=8,             gcd(L,bar L)=1.                       (5)
```

Every pair polynomial

```text
d_ij=h_i bar(h_j)/(G_ij bar(G_ij))                         (6)
```

has degree four and imaginary part exactly `tau_i-tau_j`, a
nonzero rational constant. Indeed its numerator's imaginary part
has degree at most four and leading coefficient `tau_i-tau_j`,
while its real monic denominator has degree four. The division
in (6) is exact by (5).

Thus `F_j=Re h_j/tau_j` are rational quartics with all pair
differences of degree four and

```text
F_i-F_j divides F_i F_j+1 in Q[t].                         (7)
```

The accompanying
[exact certificate](check_four_row_rational_quartic.py) uses only
standard-library rational arithmetic. It checks (1)--(6), all
root and conjugate separations, and the degree-four circle
differences below. No floating-point computation or bounded search
is used to certify this example.
The root orchestrator independently checked the two-conic identity,
the osculating-quadratic calculation, the integer factor and parity
data, and ran the complete exact certificate successfully.

## 3. How the rational parameters were obtained

The common-circle and two-conic necessary conditions from
[the eight-block real model](four_row_quartic_eight_block_family.md)
suggest placing the triple roots on the unit circle. Under its
rational parameter `z=(1-r^2+2ir)/(1+r^2)`, the second conic

```text
x^2-y^2-2Uxy-Vy=0
```

becomes

```text
r^4+(4U-2V)r^3-6r^2+(-4U-2V)r+1=0.                      (8)
```

It therefore suffices at this stage to seek four rational numbers
whose product is one and whose second elementary symmetric sum
is minus six, followed by the actual coefficient verification.

Fix `r_1=2`, write `p=2r_2`, and eliminate `r_4` using the
product condition. The quadratic for `r_3` has a rational root
precisely when the following discriminant is a square:

```text
Y^2=p^4+11p^3+30p^2-4p+1.                                (9)
```

The degenerate antipodal configuration gives `(p,Y)=(-1,5)`.
Its osculating quadratic

```text
q(p)=(-37p^2-214p+23)/40
```

satisfies the exact identity

```text
q(p)^2-[p^4+11p^3+30p^2-4p+1]
 = -21 (p+1)^3 (11p+51)/1600.                            (10)
```

The other intersection is `p=-51/11`, with `Y=665/121`.
Solving the remaining quadratic gives `r_3=11/3` and
`r_4=-1/17`, exactly (1). The resulting example is nondegenerate,
as the independent checks in Section 2 establish. This derivation
does not assert that every rational solution of (8) yields a
nondegenerate four-row system.

## 4. Integer factors, parity, and five circle points

For each root use the integer linear factor `d n-A-iB` in the
following table. Each `d` is positive and odd.

| Root | `d` | `A` | `B` |
|---|---:|---:|---:|
| `z_1` | 5 | -3 | 4 |
| `z_2` | 3085 | -2117 | -2244 |
| `z_3` | 65 | -56 | 33 |
| `z_4` | 145 | 144 | -17 |
| `w_1` | 364033085 | 59073189 | 122657188 |
| `w_2` | 1885 | -1637 | -2244 |
| `w_3` | 984115 | 2144984 | 43923 |
| `w_4` | 1483885 | -1069488 | -861101 |

Call these eight Gaussian integer factors `q_zj(n),q_wj(n)`.
For every even integer `n`, all eight have odd norms: in each
row of the table, `A` and `B` have opposite parity. Therefore
the ramified Gaussian prime `1+i` occurs in none of the factors,
rows, or pair products. This progression removes the possible
45-degree rotation in canonical primitive reduction.

Let `P(n)` be the product of the eight integer factors, and let
`H_j(n)=q_wj(n) product_(l!=j) q_zl(n)`. The five points

```text
Z_0(n)=bar(P(n)),
Z_j(n)=bar(P(n)) H_j(n)/bar(H_j(n)),       1<=j<=4,         (11)
```

are Gaussian integers: each denominator in (11) is a product
of four factors of `bar P`. Their moduli are all `|P(n)|`.
The positive integer leading constants in the `q` factors
cancel in `H_j/bar H_j`, so their angles relative to `Z_0`
are exactly the angles of `h_j/bar h_j`.

Equivalently, in the rational polynomial model,

```text
bar L * h_j/bar h_j - bar L
 = 2i tau_j * bar L/bar h_j,                              (12)
```

whose right side has degree four. Clearing the fixed denominators
preserves that degree. Thus all five points in (11) are distinct
for sufficiently large even `n`, their radius is `asymp n^8`,
and they lie on an arc of angular span `O(n^(-4))`, hence of
length `O(n^4)=O(sqrt(radius))`. All constants depend only on
the displayed fixed rational data.

Distinct block polynomials, and blocks versus their conjugates,
have nonzero fixed resultants. Their common prime factors at
integer specializations are consequently bounded by fixed
integers. Each of the eight blocks has norm of order `n^2`,
so the specified cut incidences persist at leading height after
removing these bounded overlaps. No claim of exact coprimality
at every even integer is required for (11).
Every one of the eight factors has different orientations at some
of the five points, and no block is conjugate to another. Hence
the five point polynomials have constant common gcd over `Q(i)[n]`.
A polynomial Bezout identity with fixed denominators bounds the
common Gaussian divisor after specialization. Primitive reduction
therefore also leaves radius comparable to `n^8` and the same
endpoint order of arc length.

## 5. Scope

Each three-row subtuple has all seven nonempty cut supports once.
Thus some rational full-seven-cut quartic triples do extend to
a fourth row. The maximality of the different specified triple
in [quartic_extension_maximality.md](quartic_extension_maximality.md)
remains valid with its original scope.

The rank-zero obstruction in the symmetric family also remains
valid: this construction lies outside that slice. Neither result
classifies all quartic systems, and this fixed five-point example
does not improve the global upper bound or settle uniformity.

On four nonanchor rows, this design has only the eight odd-cardinality
supports. All six two-row cuts and the full four-row cut are absent.
It therefore does not realize the full fifteen-block profile extracted
from a hypothetical unbounded endpoint cluster. The degree computation
in Section 4 gives a fixed multiple of `sqrt(radius)`, not an
`o(sqrt(radius))` arc.
