# An incidence bound for the faithful polynomial division model

This note proves a count bound for each fixed polynomial degree, with
the least-radius degree condition retained. The strongest general
bound here is `m<=3D/2`, obtained by combining root incidences with
the real Newton sums of the actual polynomials. It does not give a
count independent of the degree, and it does not transfer a
polynomial statement to the integer endpoint problem.

For the additional equal-degree full Boolean block hypotheses, the
later [Kummer-cover theorem](boolean_section_kummer_uniform_bound.md)
does give a degree-independent bound of 35 nonanchor rows. That
stronger structural assumption is not imposed in this note.

Let distinct polynomials `F_1,...,F_m` in `R[t]` have common positive
degree `D`. Assume, for every distinct pair,

```text
deg(F_i-F_j)=D,
F_i-F_j divides F_i F_j+1.                                  (1)
```

Put `h_i=F_i+i`, and let `L` be their monic lcm in `C[t]`. The
faithful least-radius condition in this model is

```text
deg L <= 2D.                                                (2)
```

The normalization from `f_i+i c_i`, with nonzero real constants
`c_i`, is `F_i=f_i/c_i`; constant scaling does not change any gcd
or lcm degree.

## 1. The bound

For `m>=2`, the degree `D` is even. Write

```text
r_D = binomial(D,D/2).
```

Under (1) and (2), the following bounds all hold:

```text
m <= 3D/2,
m <= max(3, 2+(r_D-1)^2).                                   (3)
```

More precisely, if

```text
g = deg gcd(h_1,...,h_m),
```

then `0<=g<=D/2`, and one may replace `r_D` in the second bound
of (3) by

```text
r_(D,g) = binomial(D-g,D/2-g).                              (4)
```

There is also the bound `m<=2D/g-1` when `g>0`.

Thus degree four gives `m<=6`, with no common-factor hypothesis.
If the common factor has degree two, it gives `m<=3`. The earlier
quadratic classification gives `m<=3` in degree two, consistently
with (3). The ordinary-line argument alone would have given 27 for
quartics, or 6 if a common linear factor were assumed; the moment
argument removes that latter assumption.

## 2. Exact primitive gcd degrees

Since `F_i-F_j` divides `F_i^2+1`, it has no real zero. This already
forces its degree `D` to be even.

For later use, (1) gives the exact identity up to a nonzero real
constant

```text
F_i-F_j = constant * G_ij bar(G_ij),
G_ij = gcd(h_i,h_j),       deg G_ij=D/2.                    (5)
```

Here conjugation acts on coefficients. To verify multiplicities,
fix a nonreal root `alpha`. The two polynomials `F_i+i` and
`F_i-i` are coprime. If `F_i-F_j` has positive order `e` at
`alpha`, its divisibility into `(F_i+i)(F_i-i)` says that exactly
one of the latter factors vanishes there to order at least `e`.
Suppose it is `F_i+i`, of order `a>=e`. Then

```text
ord_alpha gcd(F_i+i,F_j+i)=min(a,e)=e.
```

If instead `F_i-i` vanishes, the corresponding contribution is to
the conjugate gcd. Because the difference is real, its roots pair
with their conjugates with equal multiplicity. This proves (5).

In particular the global gcd `H=gcd_i h_i` divides every `G_ij`,
so `g<=D/2`. Also `H` and `bar H` are coprime, and their real
product divides every difference.

## 3. Only finitely many lines through any coefficient vector

Regard the polynomials as coefficient vectors in `R^(D+1)`. Fix
`F_i`. Any affine line through `F_i` and another member has
direction `F_i-F_j`, up to a nonzero real scalar. This direction
is a degree-`D` real divisor of `F_i^2+1`.

The polynomial `F_i^2+1` has degree `2D`, no real zeros, and hence
exactly `D` irreducible real quadratic factors when multiplicities
are counted. A degree-`D` real divisor chooses `D/2` of those
quadratic factors. The number of distinct choices is at most
`binomial(D,D/2)`: repeated factors can only identify choices.
Consequently at most `r_D` lines through `F_i` contain other
members of the configuration.

For the refinement, each direction contains `H bar H`. Removing
these `g` real quadratic factors leaves `D-g` quadratic factors
in `(F_i^2+1)/(H bar H)` and requires the choice of `D/2-g` more
factors. This proves the line-count bound (4), including repeated
roots. It does not assert that all these divisors occur.

## 4. The real ordinary-line argument

We use the following elementary finite incidence fact: a finite
noncollinear set of real points has a line containing exactly two
of its points.

One can reduce to the plane by a generic real linear projection,
which preserves the existing collinearities and creates no new
ones among this finite set. In the plane, minimize the positive
distance `h` between a point `P` and a line `ell` determined by
two other points. If `ell` contains at least three points, let `Q`
be the perpendicular foot of `P` on `ell`. Among the three points,
choose two `A,B` on the same closed side of `Q`, with `A` farther
from `Q`. They satisfy

```text
0<|AB|<=|AQ|<|AP|.
```

The distance of `B` from the determined line `PA` is
`h |AB|/|AP|`, strictly between zero and `h`, a contradiction.
Thus the minimizing line is an ordinary line.

Apply this fact to noncollinear coefficient vectors, obtaining an
ordinary line through `F_a,F_b`. Each other vector lies on one
of at most `r-1` other determined lines through `F_a` and one
of at most `r-1` other determined lines through `F_b`, where
`r` is either line-count bound above. Each ordered pair of such
lines contains at most one common point: if they coincided, they
would be the excluded line through `F_a,F_b`. Therefore

```text
m-2 <= (r-1)^2.                                            (6)
```

This part uses (1) but does not require the lcm condition.

## 5. The lcm condition excludes long collinear pencils

If the coefficient vectors are collinear, write

```text
F_j = F_0+a_j P,
```

where the real constants `a_j` are distinct and `deg P=D`.
For `m>=2`, (1) implies `P | F_0^2+1`. Set

```text
H = gcd(F_0+i,P).
```

Because `P` is real, has no real roots, and divides the product
of the coprime conjugate polynomials `F_0+i,F_0-i`, exactly half
its degree lies in each conjugate factor. Hence `deg H=D/2`.
For every pair,

```text
gcd(F_i+i,F_j+i)
 = gcd(F_i+i,(a_i-a_j)P)
 = gcd(F_0+i,P) = H.                                      (7)
```

Thus the quotients `(F_j+i)/H` are pairwise coprime. Their lcm
is their product, even if an individual quotient has further
factors in common with `H`. It follows exactly that

```text
L = H * product_j ((F_j+i)/H)       up to a constant,
deg L = D/2 + m D/2 = (m+1)D/2.                            (8)
```

Combining (8) with (2) gives `m<=3`. Combining this with (6)
proves the ordinary-line bound in (3) and its refinement (4).

## 6. Incidence Gram matrices and an exact common-factor bound

For every distinct complex root `alpha` of `L`, create one
coordinate for each integer layer `1<=r<=ord_alpha L`. Let `A`
be the zero-one matrix whose row `i` has a one at `(alpha,r)`
exactly when `r<=ord_alpha h_i`. It has `ell=deg L` columns.
Identity (5), including multiplicities, gives

```text
A A^T = (D/2)(I_m+J_m).                                    (9)
```

The right side is positive definite, so `rank A=m`. This already
gives `m<=ell<=2D`.

For completeness, removing the `g` columns belonging to the
global gcd gives row weight `D-g` and pair overlap `D/2-g`.
Cauchy--Schwarz applied to its column sums yields

```text
ell >= g + m(D-g)^2 / [D-g+(m-1)(D/2-g)].                  (10)
```

Combining (10) with `ell<=2D` and simplifying gives
`m<=2D/g-1` when `g>0`. These incidence observations were supplied
and independently checked by the parallel symmetry route.

## 7. Newton sums give the stronger general bound

Put `D=2n`. A polynomial `h_i=F_i+i`, divided by its real leading
coefficient, is monic with all coefficients except its constant
coefficient real. The first `D-1` Newton sums of its roots, counted
with multiplicities, are therefore real. For `1<=j<=D-1`, define
a vector on the layer coordinates of `L` by

```text
v_j(alpha,r) = Im(alpha^j).
```

The value repeats for every layer at the same root. Newton's
identities say exactly

```text
A v_j = 0       (1<=j<=D-1).                               (11)
```

We show that these vectors span a real space of dimension at least
`n`. The elementary Vandermonde fact needed is the following:
if `q<=n` distinct nonreal numbers `alpha_1,...,alpha_q` have
no conjugate pair, the `q` columns

```text
(Im alpha_l^j)_(1<=j<=D-1)
```

are linearly independent over `R`. Indeed a relation with real
coefficients `c_l` gives

```text
sum_l c_l alpha_l^j - sum_l c_l bar(alpha_l)^j = 0
                                      (0<=j<=2q-1).
```

The case `j=0` is automatic, and the remaining exponents are
within `1,...,D-1`. This is the square Vandermonde system on the
`2q` distinct nodes `alpha_l,bar(alpha_l)`. Its coefficients
`c_1,...,c_q,-c_1,...,-c_q` must all be zero.

Each individual `h_i` has no real root or conjugate pair of roots,
because `F_i` is real and `h_i,bar(h_i)` differ by `2i`. Moreover
it has at least `n+1` distinct roots. If it had `q<=n`, its
nonzero positive root multiplicities would give a dependence among
the `q` columns just shown to be independent, by (11) for that
one polynomial. Choose any `n` distinct roots from one `h_i`.
The preceding Vandermonde fact shows that the full layer moment
matrix has rank at least `n`.

Consequently (11) gives an `n`-dimensional subspace of `ker A`.
Together with (9), this proves

```text
m = rank A <= ell-n <= 2D-D/2 = 3D/2.                     (12)
```

The argument does not require `gcd(L,bar L)=1`: the selected
nodes come from a single row, where conjugate-coprimality is
automatic. Nor does it require simple roots: the repeated layer
columns and the Newton multiplicities agree exactly.

The moment argument was independently audited by both the root
orchestrator and the parallel symmetry route. It is additional
actual-polynomial information: an abstract incidence design need
not have root coordinates satisfying (11).

## 8. Scope and the unresolved strengthening

The division condition without (2) admits arbitrarily many
collinear quadratic examples; see
[the exact pencil and its classification](polynomial_heron_pencil_test.md).
Equation (8) identifies precisely why those examples fail the
least-radius condition as the number of rows grows.

The bounds here still grow with `D`. A balanced full-cut polynomial
model can itself have degree growing exponentially with its number
of rows, so these bounds cannot establish a uniform count in that model,
let alone in the integer problem.

The potential further input is the conjunction of (5), the common
root union, and its complete moment system (11). The pair incidence
data alone allow large balanced designs. The moment constraints
rule out some such designs, as (12) proves, but the lower bound
`rank{v_j}>=D/2` still spends one root coordinate per independent
constraint. No degree-independent consequence is claimed here.

For a separate positive example and a simultaneous obstruction, see
[the three-row quartic construction](three_row_quartic_construction.md)
and [its polynomial maximality certificate](quartic_extension_maximality.md).
The latter proves that one specified full-cut quartic triple has no
fourth compatible polynomial of any degree. An independent full read
and execution of `check_quartic_extension.py` confirmed its exhaustive
six-candidate list and the survival of only the three original anchors;
this is distinct from the fixed-degree count proved above.
