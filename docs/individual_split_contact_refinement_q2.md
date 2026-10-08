# Individual degree-two split contacts: two bounded reductions

The aggregate example in
[the `q=2` construction](scaled_split_contact_profile_q2.md) supplies
contact degrees `14,14,14,8` with common pole degree seven. It does
not supply seven separate `Q`-defined degree-two contacts in each
of the first three directions: its first three degree-fourteen
numerators are irreducible over `Q`.

This note gives two exact statements about that stronger contact
question. First, degree eight cannot hide a nonsplit residual point
after seven rational points have been selected in an oriented fiber.
Second, a natural restricted polynomial ansatz has at most two such
points even in one required fiber. Neither statement resolves the
stronger problem for arbitrary rational functions of degree at most
eight, nor transfers to arbitrary integer circle arcs.

## 1. Seven points force a completely split degree-eight fiber

Let `f in Q(t)` have degree `n<=8`, and put `K=Q(i)`. In the frame

```text
U_f=((f,f-1),(1,1)),
```

seven distinct reduced degree-two split contacts for `f^2+1` are
equivalent to seven distinct `K`-rational points in the fiber `f=i`.
Indeed every degree-two split residue field over `Q` is `Q(i)`, and
each conjugate pair contributes exactly one point to this oriented
fiber. Conversely each such point gives its conjugate contact pair;
it cannot be a `Q`-rational point because its image is `i`.

The same equivalence holds for the other three directions, with
oriented target values `1+i`, `1/2+i`, and `-2+i` and contact counts
`7,7,4` respectively.

Every fiber divisor of the degree-`n` morphism `P^1_K -> P^1_K`
has degree `n`, with multiplicity. If seven distinct `K`-rational
points occur in one fiber, subtract their reduced sum. The remainder
is an effective divisor defined over `K`, of degree `n-7`. Therefore:

* `n<7` is impossible;
* if `n=7`, the fiber consists of exactly seven simple rational points;
* if `n=8`, the remainder has degree one and is itself a rational
  point. There are either eight distinct simple rational points, or
  seven distinct rational points with exactly one double point.

Thus the three first oriented fibers must split completely over
`Q(i)` at both allowed degrees seven and eight. The degree-eight
case permits one extra rational point or one double point, but no
nonsplit leftover factor. Conjugation gives the same conclusion for
the opposite oriented fibers.

## 2. A sharp obstruction for monic integral polynomial frames

**Proposition.** If `F in Z[t]` is monic and nonconstant, the equation
`F(t)=i` has at most two distinct solutions in `Q(i)`.

Every rational-over-`Q(i)` root of the monic polynomial `F(t)-i`
is integral over `Z[i]`. As `Z[i]` is integrally closed in `Q(i)`,
such a root belongs to `Z[i]`; write it as `r=x+iy`, with integers
`x,y`. The coefficients of `F` are real, so `F(bar r)=-i` and `y!=0`.
Polynomial divided differences over `Z[i]` give

```text
(F(r)-F(bar r))/(r-bar r)=2i/(2iy)=1/y in Z[i].       (1)
```

The rational elements of `Z[i]` are integers. Hence `y=1` or `y=-1`.

Suppose two roots `r=x+epsilon i` and `s=u+epsilon i` have the same
imaginary sign `epsilon in {1,-1}`. A second divided difference gives

```text
r-bar s=(x-u)+2epsilon i  divides  F(r)-F(bar s)=2i
                                                in Z[i]. (2)
```

Taking Gaussian norms shows that the positive integer
`(x-u)^2+4` divides four. It must equal four, so `x=u` and `r=s`.
There is at most one root of each imaginary sign, proving the bound.

The bound two is attained:

```text
F=t^2+t+1,
F(i)=i,       F(-1-i)=i.                              (3)
```

In particular no monic integral polynomial `f` of degree seven or
eight can realize the seven individual split contacts even for the
first norm direction.

This also excludes the explicitly larger rational-function ansatz

```text
f(t)=F(h(t)),
F in Z[t] monic,  1<=deg F<=8,
h(t)=(at+b)/(ct+d) in PGL_2(Q).                       (4)
```

The map `h` is a bijection of `P^1(Q(i))`; its pole cannot lie over
the finite target `i`. Consequently `f=i` still has at most two
`Q(i)`-rational points. The frame's common pole degree is `deg F`.

Monicity and integrality in the proposition are substantive. General
rational-coefficient polynomials need not have Gaussian-integer
roots in this fiber, so (1) and the divisor argument (2) need not
apply. In particular this argument does not exclude every rational
function of degree at most eight, and does not assert that the
original lattice profiles must belong to ansatz (4).

For example `F=(t^3+7t)/6` has three distinct points `i,2i,-3i`
over `i`, since `F-i=(t-i)(t-2i)(t+3i)/6`. This explicitly prevents
extending the two-point conclusion to all rational polynomials.

Independently of these individual contacts, the aggregate example
also has a proved [specialization denominator obstruction](scaled_frame_specialization_denominator.md).

## Reproduction of the bounded checks

`python3 docs/check_individual_split_contact_refinement_q2.py` checks
the sharp example and the complete finite list of Gaussian integers
that can divide `2i`. The symbolic proof above, rather than a search
over rational functions, is the obstruction.
