# Five-row elliptic heights: an exact compatibility theorem

Canonical-height comparison on a smooth five-row invariant zero curve
does not, by itself, contradict the full-profile proportions. Its
anticanonical degree is five, and each of the ten boundary divisors
has degree one on the curve. Thus the predicted heights `10w` and
`2w` agree exactly with the degree ratio.

There is a concrete fixed smooth elliptic section with ten distinct
rational boundary points and infinitely many rational interior points
realizing these proportions. More strongly, along an infinite
subsequence the ten boundary intersection ideals have disjoint
good-prime supports, each with logarithmic norm `H/5+O(sqrt(H))`;
the bad-prime and archimedean terms remain bounded. These are actual
arithmetic points on the five-row moduli surface, not independently
assigned formal heights.

They are **not** asserted to satisfy the oriented Gaussian block
profile or the small primitive-residue hypotheses. In particular,
the good-prime ideals may contain inert rational primes. That
distinction is central to the scope of this result.

The finite exact certificate is
[check_five_row_elliptic_height_compatibility.py](check_five_row_elliptic_height_compatibility.py).
The proof of the infinite subsequence uses only compact groups and
fixed-curve height comparison, not an elliptic logarithm estimate.

## 1. Geometric and height normalizations

Let `S` be the split del Pezzo surface of degree five, realized as
the blow-up of `P^2` at four rational points in general position.
Write `H_0` for the plane line class and `E_1,...,E_4` for the
exceptional classes. The anticanonical class is

```text
L=-K_S=3H_0-E_1-E_2-E_3-E_4,       L^2=5.
```

Its six-dimensional section space is the five-row invariant space
of degree two in each row. The ten boundary curves are the four
exceptional curves and the six proper transforms of joining lines.
They satisfy

```text
L.D_e=1,             sum_(e=1)^10 D_e=2L.               (1)
```

The blow-up, anticanonical embedding, moduli interpretation, and
boundary description are given in
[Bauer--Catanese, Sections 1--2](https://link.springer.com/article/10.1007/s12215-020-00483-9).
The second identity in (1) also follows immediately by summing the
classes `E_i` and `H_0-E_i-E_j`.

Let `C` be a smooth anticanonical section containing none of the
boundary curves. Adjunction gives genus one. Each `D_e` cuts out
one rational point `B_e` on `C`, because the boundary line and the
defining hyperplane are rational and have intersection degree one.
Different `B_e` can coincide if `C` passes through a boundary node;
our explicit example below avoids every such coincidence.

All heights here are **absolute logarithmic heights over Q**. Set

```text
H(Q)=h_L(Q),
q(Q)=hat h_[O](Q)=lim_(r->infinity) 4^(-r) h_x([2^r]Q)/2.
```

Thus `q` is the canonical height for the degree-one divisor `[O]`,
and the pairing used below is

```text
<P,Q>=(q(P+Q)-q(P)-q(Q))/2.
```

In particular `<P,P>=q(P)`. This normalization is half the
canonical height normalized directly against `h_x` in some sources.
The fixed-curve height machine and the quadratic canonical-height
properties are recalled in
[Silverman, slides 8--10 and 16](https://legacy.slmath.org/attachments/workshops/301/HtSurveyMSRIJan06.pdf).
For explicit Weierstrass/local height bounds and their factor-of-two
normalization, see
[Cremona--Prickett--Siksek, Theorem 1 and Section 4](https://johncremona.github.io/papers/CPS.pdf).

## 2. An exact conditional comparison, with the varying data retained

Choose a rational origin `O` on `C`. There is a point `S_0 in C(Q)`
such that

```text
L|_C ~ 4[O]+[S_0].
```

The following canonical translates are valid Weil heights for the
indicated divisor classes:

```text
h_e^can(Q)=q(Q-B_e),
H^can(Q)=4q(Q)+q(Q-S_0).                                (2)
```

For the chosen actual height functions write finite comparison
constants `C_e,C_L` satisfying

```text
|h_(D_e)(Q)-h_e^can(Q)|<=C_e,
|H(Q)-H^can(Q)|<=C_L.                                   (3)
```

For a fixed curve and fixed metrics these constants are independent
of `Q`. They are not silently uniform when the curve varies.
The exact canonical difference is

```text
h_e^can(Q)-H^can(Q)/5
 =-2<Q,B_e-S_0/5>+q(B_e)-q(S_0)/5.                     (4)
```

The notation `S_0/5` belongs to `C(Q) tensor Q`; no rational
five-division point is asserted. Let

```text
B=max(1,C_L,max_e C_e,q(S_0),max_e q(B_e)).
```

Since `H^can>=4q(Q)` and
`sqrt(q(B_e-S_0/5))<=sqrt(q(B_e))+sqrt(q(S_0))/5`,
Cauchy--Schwarz in (4) gives the fully specified bound

```text
|h_(D_e)(Q)-H(Q)/5|
 <=(6/5)sqrt(B(H(Q)+B))+(12/5)B.                       (5)
```

This is a rigorous conditional uniform proposition: in any family
where `B<=c(J+1)`, with `J` the logarithmic height of the defining
relation, it yields

```text
h_(D_e)(Q)=H(Q)/5+O_c(sqrt((J+1)H(Q))+J+1).             (6)
```

Even granting this coefficient-uniform data bound, `J=o(w)` and
`H=10w+o(w)` give `h_(D_e)=2w+o(w)`, exactly the expected profile,
with no saving in the leading exponent.

There is a simple uniform bound on the *ambient* boundary-point
height: `h_L(B_e)<=J+O(1)`. Indeed, if a fixed rational boundary
line is spanned by integer vectors `u,v` and the hyperplane is
given by an integer linear form `g`, its intersection is
`g(v)u-g(u)v`, whose coordinates have the stated height. This
does not alone justify replacing every canonical height and every
comparison constant in `B` by `O(J+1)`. Such a replacement requires
control of a Weierstrass model, its discriminant, the origin and
boundary sections, and the height/metric changes. We make (5)
unconditional in these explicit quantities and (6) conditional
on their displayed bound, rather than hiding these dependencies.

## 3. A fixed smooth section with ten distinct boundary points

Start with the plane cubic

```text
E: y^2=x^3-2,             P=(3,5).                     (7)
```

The point `P` has infinite order. Direct counting gives

```text
#E(F_5)=6,                #E(F_7)=7.
```

Both primes have good reduction. Prime-to-p torsion injection rules
out every prime divisor of a rational torsion order: primes other
than five and seven would have to divide both counts, five would
have to divide seven, and seven would have to divide six. Hence
the rational torsion group is trivial. The needed reduction
statement appears in
[Silverman, slide 64](https://www.math.brown.edu/~jhs/Presentations/WyomingEllipticCurve.pdf).

Blow up the four points

```text
O, P, 3P, 7P.                                          (8)
```

No three are collinear: three points of a Weierstrass cubic are
collinear precisely when their group sum is zero, and each sum
of three distinct indices in `{0,1,3,7}` is positive. These four
rational points form a projective frame, so a rational projective
change sends the resulting surface to the fixed split del Pezzo
model. This change is fixed once and for all.

The proper transform of `E` is smooth, isomorphic to `E`, and is
anticanonical. Its intersections with the four exceptional curves
are the points in (8). Its intersection with the proper transform
of the line through `aP,bP` is `-(a+b)P`. Thus the ten boundary
indices are

```text
K={0,1,3,7,-1,-3,-7,-4,-8,-10}.                        (9)
```

They are all distinct. In particular this section avoids all
boundary nodes; it has no coincident-boundary obstruction of the
special two-pentagon pencil.

On the cubic, `H_0|_E~3[O]`. Hence

```text
L|_E~9[O]-([O]+[P]+[3P]+[7P])~4[O]+[-11P].            (10)
```

Put `a=q(P)>0` and evaluate at `Q_N=NP`, omitting the finitely many
boundary indices. Equations (2) and (10) give

```text
H(Q_N)=5N^2 a+22N a+121a+O(1),
h_(D_k)(Q_N)=(N-k)^2 a+O(1),          k in K.            (11)
```

Consequently, with `w_N=H(Q_N)/10`,

```text
h_(D_k)(Q_N)=2w_N+O(sqrt(w_N))       for all ten k.      (12)
```

The defining invariant relation has fixed rational coefficients,
so after clearing denominators its logarithmic coefficient height
is constant, hence `o(w_N)`. Its zero curve is smooth and
irreducible, and contains no boundary component. On the open
moduli chart these are actual distinct rational five-point
configurations. The original invariant is irreducible as well:
the open preimage is a bundle with irreducible projective-change
and row-scaling fibers, and the smooth section has multiplicity one.
Any additional
divisor factor invisible there would be a bracket boundary
component. No such component is present. Thus a fixed-curve or tiny-coefficient global
height-gap assertion cannot contradict (12).

The square-root error scale is genuine in this example. The exact
canonical difference is

```text
-2(k+11/5)N a+(k^2-121/5)a,
```

whose linear coefficient is nonzero for every integer `k` in (9).
It cannot be replaced by `O(1)` for these fixed divisor heights.

## 4. Disjoint finite-prime contacts realize the same heights

Write the x-coordinate of a nonzero multiple in lowest terms as

```text
x(rP)=u_r/v_r^2,       u_r in Z, v_r>0, gcd(u_r,v_r)=1.
```

The denominator is a square: a negative p-adic valuation of `x`
in the integral equation (7) satisfies `2v_p(y)=3v_p(x)`.
Let `Sigma` be a fixed finite set containing two and three and
every prime dividing some `v_(k-l)`, for distinct `k,l in K`.
Enlarge it, if necessary, to account for the fixed projective
coordinate change and the boundary models, including primes where
the four blow-up centers fail to remain in general position. Define

```text
D_(k,N)=product_(p notin Sigma) p^v_p(v_(N-k)).          (13)
```

At such a good prime, `p|v_r` precisely when `rP` reduces to `O`.
If `p` divides two different ideals in (13), then `(k-l)P` also
reduces to `O`, contradicting the definition of `Sigma`. Thus

```text
gcd(D_(k,N),D_(l,N))=1        for k!=l.                 (14)
```

These are exactly the good-prime intersection multiplicities
with the boundary points, including every prime-power depth.
Indeed the parameter `-x/y` at the identity has valuation
`v_p(v_r)` when `rP` reduces to the identity. They are not merely
radicals. Transport to the del Pezzo model
preserves this interpretation outside the fixed enlarged set.

We now choose an infinite subsequence for which all omitted local
terms are bounded. Fix `N_0=11`, which is not a boundary index.
For each `p in Sigma`, the point `N_0P` differs from all ten
boundary points in `E(Q_p)`. Choose an open subgroup `U_p` small
enough that the compact coset `N_0P+U_p` avoids all of them.
The compact group `E(Q_p)` has finite quotient by `U_p`; therefore
there is an integer `M>0` with `MP in U_p` for every such prime.
For

```text
N=N_0+Mj,
```

the p-adic x-coordinates of all `(N-k)P` are bounded at the
finitely many bad primes. Hence the `Sigma`-part of each `v_(N-k)`
is bounded independently of `j`.

The real cubic (7) has one connected component, a circle group.
The point `MP` is nontorsion, so its integer multiples are dense.
Choose a compact real interval with nonempty interior avoiding
the ten boundary points. Infinitely many positive `j` put `Q_N`
in this interval. Along these `N`, all ten real x-coordinates
`x((N-k)P)` are bounded as well.

For these multiples,

```text
h_x((N-k)P)=2 log v_(N-k)+O(1),
h_x((N-k)P)/2=q((N-k)P)+O(1).
```

The first equality follows directly from bounded real x-coordinate
and the reduced fraction `u/v^2`; the second is the fixed-curve
canonical-height comparison. Removing the bounded bad-prime parts
therefore gives the stronger exact asymptotic

```text
log D_(k,N)=(N-k)^2 a+O(1)
           =H(Q_N)/5+O(sqrt(H(Q_N))).                  (15)
```

All ten ideals are pairwise coprime by (14). Thus not only the
global Weil heights but also their disjoint good-prime intersection
parts realize the proposed leading proportions. The archimedean
and bad-prime contributions stay bounded. No quantitative rate
of real equidistribution or elliptic logarithm bound was used.
Extending these rational heights to `Q(i)` does not introduce
a factor two: the ideal `(D)` has norm `D^2`, and the absolute
height normalization divides its logarithm by `[Q(i):Q]=2`.

## 5. The exact limits of this dictionary

For the actual Gaussian full profile on five retained rows, the
primitive invariant height is `10w+o(w)`. A pair boundary groups
the two complementary oriented core blocks and carries `2w+o(w)`.
These are proved, with finite correction and archimedean bounds, in
[five_row_del_pezzo_arithmetic.md](five_row_del_pezzo_arithmetic.md),
Sections 2--3. Equations (12) and (15) show that this aggregate
height picture is fully compatible with actual rational points
on a fixed tiny-coefficient irreducible relation curve.

Two arithmetic conditions remain outside the comparison:

* The ideals `D_(k,N)` can contain inert primes `p=3 mod4`. The
  oriented Gaussian conductor blocks use odd split primes. Small
  primitive residues would require the inert-prime contribution
  here to be `o(w)`, a claim not supplied by canonical heights or
  by the subsequence construction.
* Each boundary point records one **unoriented** contact. It does
  not separate `H_S` from `H_(S^c)`, determine their Gaussian
  orientations, or produce a common rational chart in which all
  row corrections and primitive imaginary parts are small.

In particular, a rational moduli point can always be cleared to
a Gaussian lattice circle, but that fact supplies none of these
height or arc controls. No endpoint family is claimed here.

A concrete next arithmetic question is the split-versus-inert
prime-weight distribution in the ten shifted elliptic denominators
`v_(N-k)`, together with the orientations needed to lift each
unoriented contact. A positive proportion of forced inert weight,
or an obstruction to simultaneous low-height oriented lifts,
would add information absent from the divisor-height calculation.
The present argument establishes compatibility of the latter,
not either of those stronger arithmetic assertions.

Independent audits: the root agent and `fresh_algebraic` checked
the geometry, irreducibility, all height constants, the exact
denominator dictionary, and the progression/compactness argument.
The root also reran the exact checker successfully.
