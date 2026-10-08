# Triangle square-root covers: an exact conditional Runge budget

Adjoining square roots of the ten triangle functions creates genuine
ramified covers. In the degree-four pole-orbit case, full specialization
degree and one suitable prime at each triangle already force at least as
many finite places as geometric poles. This is a conditional obstruction
for the covers and function specified below. Degree drops, exceptional
pole orbits, different functions, and unramified Jacobian covers remain
separate questions.

## 1. The fixed curves and their degrees

Use the smooth genus-five cover `C: d^2=Delta(s,t)` and its ten squarefree
triangle quadratics `M_i`. Their roots are pairwise disjoint and avoid
`Delta=0`. Choose a fixed linear form `ell` whose root avoids both sets.
Write

```text
F_i=M_i/ell^2,       psi=ell^20/K,       K=-product_i M_i.
```

These are rational functions on `C`. Let `B_i` be the four geometric
points above `M_i=0`. At each point of `B_i`, `F_i` has valuation one,
and every other `F_j` is a unit. At the two points over `ell=0`, each
`F_i` has valuation minus two. It has no other zero or pole.

For a binary vector space `W <= F_2^10` of dimension `m`, form the smooth
projective normalization `X_W` with function field

```text
Q(C)(sqrt(product_i F_i^(w_i)) : w in W).
```

Set `q=2^m`; call coordinate `i` active if some `w in W` has `w_i=1`,
and let `a` be the number of active coordinates. The odd valuation at
`B_i` detects `w_i`, so no nonzero product is a square even over
`Qbar(C)`. Consequently this cover is geometrically connected of degree
`q`. Its inertia group has order two at each point of every active
`B_i`; it is unramified everywhere else. Even orders at `ell=0` cause
no geometric ramification in characteristic zero.

Riemann--Hurwitz and the fibers of the pole divisor give exactly

```text
genus(X_W) = 1+q(4+a),
number of geometric poles of psi pulled back to X_W = q(40-2a).   (1)
```

Indeed an active base pole has `q/2` preimages, and an inactive one has
`q`. For all ten square roots, the degree is `1024`, the genus is
`14337`, and the pole count is `20480`. Every nontrivial such product
cover is ramified: these particular functions supply no nonconstant
unramified quadratic subcover.

## 2. The field of a specialized lift

Fix a rational point `P` off the above divisors. A lift is defined over

```text
L_0=Q(sqrt(F^w(P)) : w in a basis of W),
[L_0:Q] <= q.                                                  (2)
```

The inequality can be strict. Geometric independence of functions does
not imply independence of their values.

For an exact comparison using geometric poles, fix a finite Galois
extension `k/Q` over which every pole of the pulled-back function is
rational. Put `L=k L_0` and `h=[L:k]`. Then `L/k` is an elementary
abelian extension and `h<=q`. The pole count over `L` remains the full
number in (1), even if `h<q`.

Fix an exception set of rational primes containing two, ramification
of `k/Q`, all bad reduction and coefficient primes for these models,
and primes at which distinct boundary points meet. Enlarge it so the
nonzero unit values of the relevant radicals at the boundary are units
with distinct good specializations. This set is finite and depends on
the fixed frame, `W`, and `k`; it cannot be dropped in a varying-frame
argument. A fixed rescaling of `psi` can also be included here.

Suppose an odd prime `p` outside this set divides exactly `M_i(s,t)`
to depth `e>0`, with primitive `(s,t)`. Then `ell` and the other minors
are units and

```text
v_p(psi(P))=-e.                                                (3)
```

Thus every place of `L` above `p` must occur in the finite denominator
set for `psi(P)`.

Let `v` be a place of `k` at which `P` reduces to a point `Q in B_i`.
If `i` is inactive, every generator is a unit at `Q` whose square root
belongs to `k`: evaluate at any `k`-rational point of `X_W` above `Q`.
The good-reduction exclusions and Hensel's lemma make all their values
at `P` squares in `k_v`. Hence `v` splits into exactly `h` places of `L`.

If `i` is active, choose a basis of `W` with just its first member
having coordinate `i` equal to one. The other `m-1` radicals are local
squares by the same argument. Only one quadratic extension remains:

```text
local degree <=2;
if e is odd: ramification degree 2, residue degree 1,
             and exactly h/2 places of L above v.              (4)
```

The unramifiedness of `p` in `k` preserves the odd valuation. For even
`e`, the residual extension can be trivial or unramified quadratic;
there are at least `h/2` places. No odd-depth assumption is needed for
that lower bound.

Over `Q_p` alone, without the boundary-splitting field, the remaining
unit radicals need not be squares. The local multiquadratic extension
can have degree four: at odd depth it has ramification degree two and
residue degree one or two. Counting `q/2` places directly over `Q`
would therefore be unjustified.

## 3. A conditional equality obstruction

Assume each of the ten partitions has one degree-four rational pole
orbit, equivalently its number `D_i` from the triangle splitting-field
formula has squareclass neither `1` nor `-1`. Suppose that for this
point `P` every partition has a distinct good rational prime as in (3).
Then

```text
number of required finite places of L >= h(40-2a).              (5)
```

To prove this, start with a place `v` of `k` where `P` reduces to one
point of `B_i`. Galois conjugation fixes `P` and carries this reduction
to every point of its degree-four orbit. Outside the exception set
these four reductions are distinct, so they require at least four
distinct places of `k` above the same rational prime. Formula (4)
therefore supplies at least `2h` places of `L` for each active partition,
and `4h` for each inactive partition. The rational primes were distinct,
so these contributions add.

In particular, **if the specialization degree is full, `h=q`, the finite
places alone number at least all the geometric poles**. Standard Runge
requires a strict surplus of pole orbits over the places, including
archimedean places. It fails here even before charging the latter.

A concrete sufficient condition for full degree is that the chosen
prime for each active partition has odd depth. The valuation at such a
prime applied to `F^w(P)` detects exactly `w_i`. A product which becomes
a square in `k` must therefore have every active coordinate zero. The
valuation rows detect every nonzero element of `W`, proving `h=q`.
It suffices more generally that the odd-depth coordinate evaluations
span `W^*`.

This distinguishes two different roles of odd depth: it forces local
ramification, but its crucial global use is to prevent specialization
degree from dropping. If `h<q`, the comparison in (5) scales with `h`
while (1) scales with `q`. A genuine Runge surplus might then occur,
subject to the entire denominator and archimedean count. This audit
does not dispose of such points.

The orbit hypothesis is also substantive. If `D_i` has squareclass
`1` or `-1`, there are two degree-two rational orbits. Conjugating one
reduction only forces the two points in its visited orbit. One cannot
charge all four base points to that prime. The proof of (5) does not
apply unchanged in those exceptional cases.

Finally, logarithmic prime masses in a fair profile imply neither odd
depth nor the existence of a prime outside the fixed exception set for
each partition. The frame exception integers already have height
`U^O(1)` and may absorb relevant prime mass. No assertion of those
additional hypotheses is made here.

## 4. Twists and unramified covers are separate arithmetic problems

One can force a rational lift by twisting each equation by the
squareclass of its value. Since `ell(P)^2` is a square, the individual
twists are squarefree representatives of `M_i(s,t)`. They can carry
primes from the moving triangle values. The elementary available height
bound is `|twist_i| <= U^O(1) H^2`, not a bound in `U` alone. Products
have analogous bounds of degree at most twenty in `H`. The twisted
curve, its boundary fields, and its exceptional primes now vary with
`P`; rational lifts cannot be obtained for free on one fixed cover.

Levin's unramified-cover theorem uses rational Jacobian torsion,
functions with divisible divisors, and lift-degree bounds involving
bad places, units, and ideal classes. Its superelliptic application also
requires affine integrality and a factorization of the defining
polynomial. The odd divisors of the triangle functions do not supply
those inputs, and `K=product M_i` is not the defining `Delta`.
See Sections 3--4, Theorem 6 and Theorem 11 of
[Aaron Levin, *Variations on a theme of Runge*](https://arxiv.org/pdf/0805.1345).

The separate `six_point_branch_two_torsion_audit.md` proves that two
explicit covers have no rational Jacobian two-torsion and hence no
geometrically connected unramified double cover defined over `Q`.
It also records their degree-1024 multiplication-by-two pullback cover:
this is an unramified cover defined over `Q` with nonconstant deck group
scheme. Thus absence of rational two-torsion does **not** exclude larger
unramified covers. An application through a branch-splitting field must
still control its places, lift degree, unit and class groups. None of
those bounds follows from (1)--(5).

The result proved here is the conditional no-surplus comparison (5)
for the fixed product-radical family and the pulled-back triangle pole
function. It is not a closure of all covering versions of Runge, nor
a polynomial endpoint height bound.
