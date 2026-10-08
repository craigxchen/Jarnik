# Pair-support packing and two opposite-phase extraction limits

The support-packing estimate below repackages the existing exact pair
chord bound in [the affine allocation note](linear_allocation_affine_rigidity.md).
It is not a new arithmetic separation theorem: arbitrary row units and
phase winding are already handled by that bound. Its conditional use here
is to test a proposed extraction of two-block exchanges from a general
endpoint tuple.

Two actual endpoint families delimit that proposal. A four-point Pell
family can be written with exactly one negative anchored block phase,
while its maximum two-block packing mass is two. A separate, already
primitive five-point family has exactly two negative anchored phases
after a common-unit rewrite, yet has no two-block row exchange at all.
Both keep the literal Gaussian block heights and the endpoint scale.

## 1. Pair supports and a conditional packing bound

Use the literal factorization

```text
z_i=d u_i product_(j=1)^k Gamma_j^((1+s_ij)/2)
                      bar(Gamma_j)^((1-s_ij)/2),
R=|d| product_j |Gamma_j|,
```

where the `z_i` are distinct Gaussian integers of common modulus `R`,
`d` is a nonzero Gaussian integer, and `u_i` are Gaussian units. No
phase sector, disjoint support, or individual conjugate-primitivity is
needed. For two rows `i,h`, let `T_ih` be their set of differing block
columns. Orient those factors from row `i` and let `K_ih` be their
Gaussian integer product. Factors outside `T_ih` are common, so

```text
z_i=g u_i K_ih,       z_h=g u_h bar(K_ih),
|K_ih|=product_(j in T_ih)|Gamma_j|.
```

The chord numerator `|K_ih-(u_h/u_i)bar(K_ih)|` is one of
`2|Im K_ih|`, `2|Re K_ih|`, and
`sqrt(2)|Re K_ih +/- Im K_ih|`. The relevant integer is nonzero
because the source points are distinct. Thus the numerator is at
least `sqrt(2)` (at least `2` if the units are common). Since the
chord is at most `R Delta` for a containing angular arc of width
`Delta`, the **raw support height** obeys

```text
product_(j in T_ih)|Gamma_j| >= sqrt(2)/Delta,       (1)
```

with `2/Delta` in the common-unit case. The factorization can have
shared Gaussian primes: any rational cancellation in a primitive pair
cofactor only makes that cofactor smaller than this raw product.

Give row-pair supports nonnegative weights `w_ih` with each block load
`c_j=sum_(ih:j in T_ih) w_ih<=1`, and let `P=sum_ih w_ih`. Multiplying
(1) and using `|Gamma_j|>=1` yields

```text
(sqrt(2)/Delta)^P <= product_j |Gamma_j|^c_j
                    <= R/|d|.                          (2)
```

At `Delta<=C/sqrt(R)`, this is

```text
R^(P/2-1) <= (C/sqrt(2))^P/|d|.                      (3)
```

For a packing mass `P>=5/2` and `C>sqrt(2)`, it gives the uniform
radius bound `R<=(C/sqrt(2))^10`; for `C<=sqrt(2)` it excludes
`R>1`. With common row units, replace `sqrt(2)` by `2`, giving
`R<=(C/2)^10` for `C>2`. The exponent is maximal at `P=5/2` and
decreases for larger `P`. A bounded radius gives a bounded point count:
a radius-`R` integer circle has at most `4 floor(R)+2` points.

A two-block exchange graph simply keeps row pairs with
`|T_ih|=2`; a complete graph on `q` exposed block directions has
packing mass `q/2`. This application adds no extraction theorem:
general endpoint tuples need not contain such supports or a packing
of mass above two. The general Roth certificates in
[the coupled-product note](coupled_block_roth_packing.md) concern
different, potentially smaller integral monomials in the row-difference
span and have their own `P>4` height threshold.

## 2. A synchronously flipped opposite group

Let `H_1,...,H_r,G_1,...,G_b` be nonzero Gaussian integers and set
`L=product_(ell=1)^b bar(G_ell)`. Consider the `r+1` rows

```text
z_0=d u bar(L) product_(j=1)^r H_j,
z_j=d u L bar(H_j) product_(k!=j) H_k,    1<=j<=r.
```

Every non-anchor row flips the whole `G` group and exactly one `H_j`.
The synchronized `G` columns are merged literally into the single block
`bar(L)`, with no reduction of shared prime factors, before applying
the two-column exchange criterion.
The circle radius is `R=|d||L|product_j|H_j|`. Every pair of rows
supplies a two-block support among these `r+1` merged directions,
regardless of their phase branches or row units. The graph is complete
and has fractional edge packing mass

```text
P=(r+1)/2
```

by giving every edge weight `1/r`. With one common unit, `r>=4` gives
the uniform bound `R<=(C/2)^10` for `C>2`, even for several
opposite-phase blocks inside `L`,
shared prime support, and arbitrary phases in blocks held fixed across
the displayed rows. If `r<=3`, the displayed profile has at most four
points. The anchor-to-leaf phase can be `2(arg L-arg H_j)`, a signed
cancellation rather than a one-sided sum.

## 3. Primitive four-point Pell boundary example

Use the Pell progression and Gaussian blocks from
[the original four-point construction](four_point_bonus_counterexample.md):

```text
U+V sqrt(5)=(9+4 sqrt(5))^n,  n=1 mod 10,
U^2-5V^2=1,
P'=1+2i,
A=(U+V)+2Vi,
B=2V+(U-V)i,
C=(2U+8V)+(3U+V)i,
D=P' bar(A)=(U+5V)+2Ui.
```

Take the four literal blocks `A,B,C,G=bar(D)`. All of `A,B,C,D` lie
in the open first quadrant; `G` has a negative anchored phase and the
other three blocks have positive anchored phases. The common-unit star
rows are

```text
z_0=bar(P') A^2 B C,
z_1=P' bar(A)^2 B C,
z_2=P' a bar(B) C,
z_3=P' a B bar(C),                 a=N(A).
```

They are exactly `z_0=ABC G`, with each `z_j` for `j=1,2,3`
flipping `G` and its own block. Their norm is

```text
R^2=5 a^2 b c,
b=N(B), c=N(C),
a=1+10V^2+2UV,
b=1+10V^2-2UV,
c=13+130V^2+38UV.
```

The original construction proves `5,a,b,c` pairwise coprime on this
progression and `A,B,C` coprime to their conjugates. Since
`N(D)=5a`, `D` and `G` are also conjugate-coprime. The complete row
gcd is one: at primes of `A`, row 0 uses `A^2` while row 1 uses
`bar(A)^2`; at primes of `B` and `C`, rows 2 and 3 respectively use
the opposite orientation; at the prime over five, row 0 uses
`bar(P')` and the other rows use `P'`. The overlap between `D` and
`A` is therefore harmless and is retained exactly.

For pairs in the order `01,02,03,12,13,23`, primitive pair cofactors
and their *exact* imaginary coordinates are

```text
q_ij:    bar(P')A^2, bar(P')AB, bar(P')AC,
         bar(A)B, bar(A)C, bar(B)C,
Im q_ij:     -2,       1,        -1,
              1,       3,        -2.
```

The identities follow directly from `U^2-5V^2=1` and show all four
points are distinct. They also give a fixed endpoint constant. Since
`U<=3V`, `V>=4`,

```text
10V^2<=a<=17V^2,
 4V^2<=b<=11V^2,
130V^2<=c<=257V^2.
```

Each cofactor satisfies `|Im q_ij|/|q_ij|<=1/(sqrt(40)V^2)`.
The pair geodesic angle is `2 asin(|Im q_ij|/|q_ij|)`, hence is at
most `2/(sqrt(10)V^2)` using `asin x<=2x` here. All four points lie
in a common short neighborhood; their containing angular span obeys
the same maximum pair bound. Moreover

```text
sqrt(R)<= [5*17^2*11*257]^(1/4) V^2 < 45V^2.
```

Thus `Delta sqrt(R)<(2/sqrt(10))*45<29`, so the points lie on an arc
of length `29 sqrt(R)` for every index in the progression. Their
radius tends to infinity.

After conjugating `G` for the exchange chart, the four selected
directions are `D,A,B,C`, with arguments in a common first-quadrant
sector. Their exchange graph is `K_4`, whose maximum packing mass is
two. In the original anchored chart, row 0 is the only row retaining
the negative block, while the other three rows flip singleton positive
blocks. Those singleton sets form an antichain, so there is no
same-negative-state monotone edge. The endpoint assumption therefore
does not force the one-sided witnesses that would bound all three
positive block phases individually.

The Pell example tests the packing threshold and the failure of a
simple one-opposite-to-one-sided reduction. It does not assert that
large endpoint tuples can be reduced to this star profile.

## 4. An actual primitive five-point tuple with no exchange

The [six-factor affine-shape family](five_point_affine_shape_cubic_family.md)
provides a stronger incidence obstruction. Let `T` be a positive
multiple of `600`, put `n=T/60`, and write

```text
b_a=60/a in {60,30,20,15,12,10},
h_a=1+i b_a n,             kappa_a=-i h_a=b_a n-i,
S_0={4,6},       S_1={1,3,6},       S_2={1,4,5},
S_3={2,3,5},     S_4={1,2,3,4}.
```

The original exact points are

```text
z_j=(-1)^|S_j| product_(a in S_j) h_a
                  product_(a notin S_j) bar(h_a).
```

There are six blocks, so for every row

```text
product_(a in S_j) kappa_a product_(a notin S_j) bar(kappa_a)
 = -(-1)^|S_j| product_(a in S_j) h_a
                    product_(a notin S_j) bar(h_a).
```

Hence the same actual tuple has the **literal common row unit** `-1`:

```text
z_j=- product_(a in S_j) kappa_a
        product_(a notin S_j) bar(kappa_a),
R=product_a |kappa_a|.
```

Orient each block from the actual anchor `S_0`. The anchored blocks
for `a=4,6` are `kappa_a` and have negative arguments; the other four
are `bar(kappa_a)` and have positive arguments. Thus exactly **two**
anchored nonconstant phases are negative, and every phase has
absolute value `atan(1/(b_a n))<pi/2`. No individual phase winds.

The linked note proves that all five points are distinct, their
complete Gaussian gcd is one, the six block norms are pairwise
coprime, and their least containing arc has length less than
`5 sqrt(R)` for every `T>=600`. In particular `R->infinity` along
this endpoint family. Its normalized arc constant tends to
`2 sqrt(5)`, so it is **not** an obstruction to an extraction theorem
restricted to sufficiently small endpoint constants such as `C<=2`.
Sections 5--6 of the linked note already prove a positive lower bound
for the primitive normalized arc constant throughout this six-factor
rational shape family; no vanishing-constant variant is inferred here.
Its exact row Hamming distances are

```text
01=3, 02=3, 03=5, 04=4,
12=4, 13=4, 14=3,
23=4, 24=3, 34=3.
```

There is **no** two-block row support. The six sign columns are
pairwise distinct even modulo complement, so merging synchronized
columns creates none. Column conjugation or permutation preserves
these distances. This is an actual Gaussian-primitive, common-unit
endpoint counterexample to the naive rule that one or two opposite
anchored phases must expose a local exchange after block merging.

The full pair-support packing is critical as well. Every row-pair
support has at least three of the six blocks, so summing the six
vertex-load inequalities gives `3P<=6`, or `P<=2`. Equality holds:
put weight `1/2` on each of the four pairs `01,02,14,24`, whose
supports are respectively

```text
{1,3,4}, {1,5,6}, {2,4,6}, {2,3,5}.
```

Each block then has load exactly one. Thus the exact pair-support LP
optimum is `P=2` throughout this unbounded primitive endpoint family.
The failure of the strict exponent gap in (3) is real, not an artifact
of considering only two-block exchanges.

It also resists a natural exact integral block-coordinate change.
For `U in GL_6(Z)`, form *unreduced integral Gaussian monomial blocks*

```text
Lambda_v=product_a kappa_a^((U_av)_+)
                  bar(kappa_a)^((-U_av)_+).
```

Their total literal radius product is

```text
product_v |Lambda_v|
 =product_a |kappa_a|^(sum_v |U_av|).                (4)
```

Each row of a unimodular matrix has absolute load at least one. If
`U` is not a signed permutation, some row has load at least two;
because every `|kappa_a|>1`, (4) is strictly larger than the actual
radius `R`. No exact same-radius sign factorization using these
unreduced blocks and a Gaussian common multiplier `d'` can exist:
it would require `R=|d'|product_v|Lambda_v|` with `|d'|>=1`.
The remaining signed permutations preserve the empty exchange graph.

This last statement concerns literal integral Gaussian monomial
blocks. Rational transforms with denominators, exact gcd reductions,
or moving correction factors have different height accounting and
are outside this counterexample. The five-point incidence example
does not prove that every larger endpoint tuple lacks a different,
more global extraction mechanism.
