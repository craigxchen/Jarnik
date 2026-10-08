# One and two relations among four signed products

A single four-term relation with bounded coefficients is possible even
when all three signature groups are conjugate-primitive, have disjoint
rational-prime norm supports, and grow at the same rate. Two independent
relations obey a coefficient-height product bound. The bound has the
correct order for Gaussian coefficients on the explicit family below.

This settles the single-circuit question left in
[the signed-witness note](signed_block_flip_witness_height.md). It does
not supply two sufficiently small relations from an actual endpoint
profile, and does not improve the general lattice-arc count.
For four actual thin corrected witnesses, the stronger
[matching-norm argument](four_signed_witness_matching_norm_gap.md)
now supplies a positive height exponent. Its repeated opposite-edge
divisors still need a common-support construction from general rows.

## 1. A full three-group family

For a positive integer `t=1 mod15`, put

```text
A=t+i(t+5),        B=(t+1)+i(t+4),        C=(t+2)+i(t+3),
T0=ABC,           T1=A bar(B) bar(C),
T2=bar(A) B bar(C), T3=bar(A) bar(B) C.
```

The exact polynomial identity is

```text
-T0-6T1+8T2-3T3=0.                               (1)
```

One certificate is the following coefficient table, in ascending powers
of `t`:

| Product | Constant | Linear | Quadratic | Cubic |
|---|---|---|---|---|
| `T0` | `-55-50i` | `-71-19i` | `-24+6i` | `-2+2i` |
| `T1` | `55-50i` | `51-41i` | `16-14i` | `2-2i` |
| `T2` | `25-70i` | `29-59i` | `12-18i` | `2-2i` |
| `T3` | `-25-70i` | `-1-69i` | `8-22i` | `2-2i` |

The indicated linear combination vanishes in each column. All three
factors have coprime coordinates of opposite parity: their coordinate
gcds are `gcd(t,5)`, `gcd(t+1,3)`, and one. Thus they are
conjugate-primitive. Their odd norms are

```text
nA=2t^2+10t+25, nB=2t^2+10t+17, nC=2t^2+10t+13.   (2)
```

The differences are eight, four and twelve. The first two give gcd one
by oddness; the third could only contribute three, but `nA=1 mod3` on
the progression. Hence all three norms are pairwise coprime. Every
factor is a nonunit and `log nA,log nB,log nC=2log t+log2+o(1)`.
Distinct sign products cannot be associates because a changed block
has disjoint conjugate support.

This is not an endpoint circle family. Let `rho=sqrt(nA nB nC)` and
align the four points as `(-T0,T1,T2,T3)`. With

```text
s=2t+5, deltaA=atan(5/s), deltaB=atan(3/s), deltaC=atan(1/s),
```

their arguments, relative to `-pi/4`, are

```text
deltaA+deltaB+deltaC,   deltaA-deltaB-deltaC,
-deltaA+deltaB-deltaC, -deltaA-deltaB+deltaC.
```

They are in the displayed decreasing order. The exact angular span is
`2(deltaA+deltaB)`. Therefore the arc length `L` satisfies

```text
rho ~ 2sqrt(2)t^3,
L ~ 16sqrt(2)t^2 ~ 8sqrt(2)rho^(2/3),
L/sqrt(rho) -> infinity.                            (3)
```

The example disproves a coefficient gap based only on one signed
circuit and its support data; it is not a counterexample to the goal.

## 2. Two independent circuits cost a full block

Now let `A,B,C` be any nonunit conjugate-primitive Gaussian integers
with pairwise coprime norms, and define `T0,...,T3` as above. Suppose
`c,d in Z[i]^4` are independent over `Q(i)` and

```text
sum_j c_j Tj=sum_j d_j Tj=0.
H(c)=max_j |c_j|,       H(d)=max_j |d_j|,
Mij=c_i d_j-c_j d_i.
```

Reduction at each oriented block gives the six exact divisibilities

```text
bar(A)|M01,  A|M23,
bar(B)|M02,  B|M13,
bar(C)|M03,  C|M12.                                 (4)
```

For example, modulo `bar(A)` the last two products vanish, whereas
`T0,T1` are units. The two coefficient rows both kill the column
`(T0,T1)`, so their determinant `M01` vanishes modulo `bar(A)`.
This argument works in the full quotient ring, including prime powers:
a column with a unit entry can be completed to an invertible matrix.
The other five assertions follow identically. This is the existing
[coefficient-minor principle](auxiliary_form_route.md) on these particular
two-versus-two cuts, not a new general lattice argument.

Rank two ensures some `Mij` is nonzero, while
`|Mij|<=2H(c)H(d)`. Hence

```text
2H(c)H(d) >= min(|A|,|B|,|C|).                       (5)
```

If both coefficient vectors lie in `Z^4`, their minors are ordinary
integers. Divisibility by a conjugate-primitive block then forces its
entire norm to divide that minor, so the stronger bound is

```text
2H(c)H(d) >= min(nA,nB,nC).                          (6)
```

Thus, if all three log norms equal `w+o(w)`, a common upper bound `H`
for the two heights is at least `exp(w/4-o(w))` for Gaussian
coefficients, or `exp(w/2-o(w))` for integer coefficients. More usefully,
if the first relation is subpower, the second must have height at least
`exp(w/2-o(w))`, or `exp(w-o(w))` in the integer case. The coefficient
vectors may have zero entries; independence is essential.

The Gaussian order is attained by the family (1). Its first vector is
`c=(-1,-6,8,-3)`, of height eight. Since `A-2B+C=0`, a second relation is

```text
d=(0,-bar(A),2bar(B),-bar(C)),
sum_j d_j Tj=bar(A)bar(B)bar(C)(-A+2B-C)=0.           (7)
```

It is independent of `c` because `d0=0` and `c0` is nonzero, and its
height is `2|B|~2sqrt(2)t`. Bound (5) gives the matching lower order
`H(d)>=constant*t`. The three nonzero coefficients in (7) contain
exactly the oriented blocks demanded by singleton divisibility.
No sharpness assertion for the integer-coefficient bound (6) is made.

## 3. The available row identities do not supply the two circuits

Actual endpoint rows do give a cheap four-term identity:

```text
Y_j(P_i-bar(P_i))-Y_i(P_j-bar(P_j))=0.                (8)
```

Conjugation returns its negative, not a second independent relation.
The terms also have differing supports. Write the core of `P_i` as
`DE`, that of `P_j` as `DF`, with `D` shared and `E,F` exclusive.
At a block of `E`, the four oriented valuation pairs are

```text
(e,0), (0,e), (0,0), (0,0),                         (9)
```

at each Gaussian prime of core depth `e`, apart from subpower
corrections. These are not the two-versus-two signed states in (4).

Consider the restricted conversion that multiplies all four terms of
(8) by one common Gaussian monomial and writes each resulting term as
a Gaussian integer coefficient times a full signed product on the union of
the core supports. Let the common valuation shift at this prime and
its conjugate be `(a,b)`. Each absent term must contain one full
orientation afterward, so `a+b>=e` and `a,b>=0` in the clean model.
After removing the four chosen target orientations, the total exponent
left in the four coefficients is

```text
4(a+b)-2e >= 2e.
```

Summing over the exclusive blocks shows that the sum of coefficient
log moduli is at least `log N(EF)`. To retain correction overlaps,
let `S` be their total depth at this prime and its conjugate across
the four terms, including the coefficients `Y_i,Y_j`. An originally
absent term forces `a+b>=e-S`; the total remaining coefficient depth
is therefore `4(a+b)-2e+S>=2e-3S`. The sum of `S log p` is `o(w)`,
so allowing the corrections changes the bound by only `o(w)`.
For a fixed extracted
`m`-row profile with `r=2^(m-1)` incident blocks per row,
`log N(E)=log N(F)=rw/2+o(w)`. At least one converted coefficient
therefore has modulus at least

```text
exp(rw/4-o(w)).                                      (10)
```

This common-monomial conversion cannot provide the subpower signed
circuits needed for (5). It does not cover arbitrary nonlinear
recombinations. Likewise the
[five-star lift](four_row_integer_five_star_lift.md) supplies only its
established first-minimum bound, and a pair bracket identity changes
the products and their supports rather than supplying a second relation
on this fixed tuple.

The missing implication remains the full mathematical task: produce
additional compatible arithmetic information from actual many-point
endpoint clusters. Equations (5)--(6) are conditional coefficient bounds,
not a proof that every endpoint cluster supplies their small input data.

## Verification

The [exact checker](check_four_term_signed_circuit_height.py) verifies
the full polynomial identity (1), primitive and norm-support conditions
on a parameter progression, both independent relations, all six
divisibilities in (4), their integer-coefficient strengthening, and the
valuation cost in (9). The family and bounds are proved above for all
parameters; finite tests are supplemental.
