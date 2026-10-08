# Changing the row degrees in the invariant-relation test

For equal positive row degree, degree two gives the smallest proved
successive-minima coefficient exponent once there are at least six rows.
Degree one has an injective balanced evaluation map and a smaller space;
degrees above two keep the same balanced rank while increasing the
primitive numerical height. Four rows have exact threshold equality in
degrees one, two, and three. None of these comparisons supplies a strict
height margin for the uniform bound.

The statements concern the same actual integer full-cut system as
[invariant_relation_lattice_threshold.md](invariant_relation_lattice_threshold.md).
The weighted formulas at the end give a finite test for unequal degrees;
they are not a classification of all possible degree choices.

The later [arbitrary-degree scalar obstruction](arbitrary_degree_scalar_height_obstruction.md)
settles that specific comparison for every admissible positive integer
multidegree: `A(d)>=2 min(s(d),r(d)-1)`. It requires no full-rank claim
for the weighted scalar map. It leaves individual irreducible numerical
zeros in the scalar kernel outside its conclusion.

The finite table and unequal-family counts are checked exactly by
[check_invariant_degree_variation.py](check_invariant_degree_variation.py).

## 1. Equal row degree and the complementary coefficient threshold

Let `m=2q>=4`, and let `V_d` consist of invariants of degree `d>=1`
in each row. Its dimension is

```text
r_d=[x^(md/2)](1+x+...+x^d)^m
       -[x^(md/2-1)](1+x+...+x^d)^m.                    (1)
```

At a balanced cut the scalar still multiplies
`product_inside Y_i^d product_outside P_j^d`.
Swapping the two directions gives

```text
c_(S^c)=(-1)^(md/2)c_S.
```

The sign has no effect on divisibility. For a numerical zero relation,
primitivity of the actual rows gives

```text
M_(S,d) | c_S,
M_(S,d)=N(H_S)/N(gcd_G(H_S, product_outside K_j^d)).
```

Complementary norm supports are disjoint. Thus a nonzero balanced scalar
forces polynomial coefficient sum

```text
log C_Q >=2(1-eta)w-2d sum_i log|K_i|
        >=2(1-eta)w-2dm sigma.                          (2)
```

The leading coefficient threshold remains `2w`; it does not scale by
`d`. The correction loss does increase with degree.

## 2. Exact balanced rank

For `d=1`, evaluation on all balanced cuts is injective. Indeed every
multilinear weight-zero polynomial has an expansion in the monomials
`product_(i in S)x_i product_(j outside S)y_j`, with `|S|=q`, and the
balanced evaluations recover exactly these coefficients. Consequently

```text
rank C_1=r_1=binomial(m,q)/(q+1).                        (3)
```

The dimension is the Catalan number, also given by the noncrossing
matching basis. This basis description is recalled in
[Patrias--Pechenik--Striker, Section 2](https://www.mat.univie.ac.at/~slc/wpapers/FPSAC2022/55.pdf).

For every `d>=2`, the balanced rank is instead

```text
rank C_d=s=binomial(m,q)/2.                              (4)
```

For degree two this was proved by the positive-definite squared-matching
incidence Gram matrix. Choose one degree-one invariant `L` whose value
is nonzero at every balanced cut. Such an integer invariant exists: in
an integral matching basis take coefficients `1,2,4,...`; every cut
has a nonzero basis value in `{0,1,-1}`, and its largest nonzero weighted
term cannot cancel. Multiplication by `L^(d-2)` sends degree-two
evaluations to degree-`d` evaluations by an invertible diagonal map on
the finite set of unordered cuts. This proves (4), including odd `d`.
All constants in this construction depend only on `m,d`.

## 3. Primitive evaluation heights and the critical minimum

Each graph now has `md/2` brackets. Its least internal-edge count at a
cut of size `a` is `d max(0,a-q)`, attained by the `d`th power of a
perfect matching chosen for that cut. Therefore, with

```text
A_m=m 2^(m-2)-(m/2)binomial(m,q),
```

the primitive evaluation height on `V_d` is

```text
h_d=(d/2)A_m w+o(w).                                    (5)
```

This is the same actual gcd argument as in the degree-two note: graph
values attaining each local minimum bound any extra common factor by
`product_edges b_ij^d`. It does not estimate cancellation in a sum.

For `m>=6,d>=2`, evaluation on `ker C_d` is nonzero. Multiply a
two-triangle degree-two invariant by any nonzero matching product to
the power `d-2`. Its numerical value is nonzero at distinct directions,
and all balanced scalars vanish. Thus the relevant outside-relation
successive-minima denominator is `s`, yielding coefficient exponent

```text
E_(m,d)=d A_m/(2s),       m>=6,d>=2.                     (6)
```

For degree one the balanced kernel is zero. The first relation in
the rank-`r_1-1` numerical kernel instead gives

```text
E_(m,1)=(A_m/2)/(r_1-1).                                (7)
```

These are the guaranteed upper exponents from the lattice method,
to be compared with the lower threshold `2` in (2). A bound at an
exponent larger than two cannot force a contradiction below that
threshold. For `m>=6`, (7) is larger than (6) at `d=2`, because
`s>2(r_1-1)`. For `d>=2`, (6) increases linearly with `d`.

Four rows need separate treatment. Here `r_d=d+1`, `A_4=4`, and the
balanced rank is two for `d=1` and three for `d>=2`. The kernel is
zero for `d=1,2`. For `d>=3`, the product of all three perfect-matching
invariants, times any matching to power `d-3`, is a nonzero numerical
element of the balanced kernel. Hence the correct exponents are

```text
E_(4,1)=2,       E_(4,2)=2,       E_(4,d)=2d/3 for d>=3.
```

In particular degree three also has exact equality with the arithmetic
threshold. Applying a denominator equal to the balanced rank in the
degree-two four-row case would incorrectly give `4/3`.

| `m` | `d` | dimension | balanced rank | primitive height / `w` | guaranteed exponent |
|---|---:|---:|---:|---:|---:|
| 4 | 1 | 2 | 2 | 2 | 2 |
| 4 | 2 | 3 | 3 | 4 | 2 |
| 4 | 3 | 4 | 3 | 6 | 2 |
| 4 | 4 | 5 | 3 | 8 | 8/3 |
| 6 | 1 | 5 | 5 | 18 | 9/2 |
| 6 | 2 | 15 | 10 | 36 | 18/5 |
| 6 | 3 | 34 | 10 | 54 | 27/5 |
| 6 | 4 | 65 | 10 | 72 | 36/5 |
| 8 | 1 | 14 | 14 | 116 | 116/13 |
| 8 | 2 | 91 | 35 | 232 | 232/35 |
| 8 | 3 | 364 | 35 | 348 | 348/35 |

## 4. Exact formulas for unequal positive degrees

Let the row degrees be positive integers `d_i`, with even total `D`,
and assume `max d_i<=D/2` so that the invariant space is nonzero.
Write `d(S)=sum_(i in S)d_i`. The dimension and the number of unordered
weighted-balanced cuts are

```text
r(d)=[x^(D/2)]product_i(1+x+...+x^d_i)
        -[x^(D/2-1)]product_i(1+x+...+x^d_i),
s(d)=(1/2)#{S:d(S)=D/2}.                                (8)
```

Only those weighted-balanced cuts have a scalar specialization.
Complementing again identifies the scalars up to sign, and the exact
coefficient threshold is

```text
log C_Q>=2(1-eta)w-2 sum_i d_i log|K_i|.                 (9)
```

For a cut `S`, the smallest internal-edge count is
`g(S)=max(0,d(S)-D/2)`. To see attainability when `d(S)>D/2`, send all
outside stubs across the cut. If `g=d(S)-D/2`, give an inside vertex
at least `max(0,d_i-g)` crossing edges. The sum of these requirements
is at most the total outside degree: for one positive requirement use
`d_i<=D/2`; for two or more use their total degree at most `d(S)`.
The remaining inside degrees have total `2g` and maximum at most `g`,
so admit a loopless multigraph. The minority-side case follows by
complementation. Consequently the primitive evaluation height is

```text
h_d=A(d)w+o(w),
A(d)=D 2^(m-3)-sum_S max(0,d(S)-D/2).                   (10)
```

If all `d_i` are even, the balanced rank is exactly `s(d)`. Replace
row `i` by `d_i/2` identical clones, square perfect-matching monomials
on the clones, and use the Gram argument on the resulting admissible
balanced clone cuts. Its matrix is a principal submatrix of a positive
definite balanced-cut Gram matrix. Matchings using two clones of the
same row vanish on every admissible cut, so deleting those zero columns
does not affect the rank. If there are no weighted-balanced cuts,
the rank and the scalar restriction are both zero.

For other unequal degree vectors, (8) only gives the upper bound
`rank C_d<=min(r(d),s(d))`; full rank is not asserted here. Even when
the rank is known, the outside-relation index also depends on whether
evaluation on its kernel is nonzero, as explained in the main note.
Thus formulas (8)--(10) are a finite comparison procedure, not a blanket
claim that every unequal choice is worse.

For a tractable exact family, the degrees `(2,2,2a,2a)`, `a>=2`, factor
out `Delta_34^(2a-2)` from every invariant. The remaining space is the
three-dimensional equal-degree-two four-row space. There are only two
weighted-balanced partitions, so its scalar kernel is one dimensional,
with a nonzero matching-product value at distinct directions. Formula
(10) gives `A=4` independently of `a`, and the critical exponent is
`4/2=2`, again exact threshold equality. Similarly `(1,1,a,a)` with
`a>=2` factors out `Delta_34^(a-1)`; its two-dimensional residual space
has zero balanced kernel and exponent `2/(2-1)=2`.

These computations leave no positive exponent margin for the displayed
equal-degree choices or these unequal four-row families. They do not
exclude more detailed arithmetic input, higher jets with new numerical
control, or other nonuniform degree constructions.
