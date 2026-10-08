# Rational two-torsion and unramified double covers of the circle cover

Both existing first-ruling genus-five witnesses have `J(Q)[2]=0`. This has
an exact, short modular certificate, and disproves any claim that the
six-point origin universally supplies a nonconstant unramified double cover
over Q. It does **not** rule out unramified covers of higher degree with a
nonconstant arithmetic deck group, covers after field extension, or special
weights with rational two-torsion. No uniform endpoint bound follows.

## 1. The correct Galois module

Let `C: y^2=f(x)` be a smooth degree-twelve hyperelliptic model over Q,
`Omega` its twelve geometric branch points, and `J=Jac(C)`. Then

```
J[2] = {even-cardinality subsets S of Omega}/(S ~ Omega\S),    (1)
```

with addition symmetric difference and its natural Galois action. In
particular this vector space has dimension ten over F_2. The even-degree
exact sequence underlying (1) is stated in [Poonen, Lectures on rational
points on curves, Exercise 6.2](https://math.mit.edu/~poonen/papers/curves.pdf).
Here is also a direct explanation of the particular formula being used.

Let F_infinity be the degree-two fiber over infinity. The subset S gives

```
D_S=sum_(i in S) P_i - (|S|/2) F_infinity.
```

The divisor of `product_(i in S)(x-alpha_i)` is `2D_S`, and the full subset
has principal divisor class because `div(y)=sum_i P_i-6 F_infinity`.
There are no further identifications: an element of the quadratic field
`Qbar(x,y)` whose square lies in `Qbar(x)` is either in `Qbar(x)` or in
`y Qbar(x)`. The polynomial `product_(i in S)(x-alpha_i)` can be a square
in that field only for the empty or full subset, by its root valuations.
The resulting `2^10` classes exhaust J[2], which has order `2^(2g)`.
The degree-two fiber class is intrinsic, so (1) is Galois equivariant even
if its two points are not rational.

A rational two-torsion class is thus a Galois-invariant **unordered even
partition**, not necessarily a Galois-invariant individual subset. If f
is irreducible, the only possible nonzero invariant partitions have sizes
6+6, with some group elements swapping the two halves. Transitivity alone
therefore does not exclude two-torsion. For example, a transitive 12-cycle
preserves the unordered partition into alternating positions.

A particularly convenient certificate is stronger:

> If the branch Galois group contains a permutation consisting of exactly
> two odd cycles, then `J(Q)[2]=0`.

To prove it, an invariant unordered partition for this permutation either
has each part preserved or has the parts swapped. In the first case its
parts are unions of the two cycles; the only even unions are empty and
full. The second case is impossible because swapping the parts makes every
cycle even. No irreducibility hypothesis is required for this certificate.

## 2. Exact witness certificates

Use the first ruling and the polynomial constructed by
`check_six_point_isotropic_circle_cover.probe`. Clear rational coefficient
denominators and remove the integer content. Its squarefree factorizations
at the following primes have these exact degree patterns:

| Primitive central weights | First prime and degrees | Second prime and degrees |
| --- | --- | --- |
| `(1,-20,250,-1000,1445,-676)` | `59: (1,11)` | `73: (5,7)` |
| `(-45,616,-3850,15895,-27500,14884)` | `43: (3,9)` | `29: (4,4,4)` |

At each prime the degree stays twelve and the polynomial is coprime to its
derivative. Thus the Galois group contains a permutation of the displayed
cycle type, by [Milne, Fields and Galois Theory, Theorem 4.28](https://www.jmilne.org/math/CourseNotes/FT.pdf).
For a nonmonic integral polynomial with leading coefficient a prime to p,
apply the monic theorem to `a^11 f(X/a)`; scaling its roots by a preserves
the factor degrees and squarefreeness modulo p.

The first prime in each row already proves `J(Q)[2]=0` by the two-odd-cycle
argument. Independently, these pairs give an elementary irreducibility
certificate: the degree of a rational factor must be a subset sum of each
modular factor pattern. The common subset sums in either row are only
`{0,12}`. This also verifies the previous irreducibility assertions without
requiring a characteristic-zero factorization oracle or identifying the
full Galois group.

Run:

```
UV_CACHE_DIR=/tmp/uv-cache uv run --offline --with sympy python \
    docs/check_six_point_branch_two_torsion.py
```

The checker rebuilds the actual polynomials, verifies every modular factor
by the Frobenius-power/gcd irreducibility criterion, multiplies the factors
back, checks the subset-sum intersection, and enumerates all 1024 even
partition classes for the indicated Frobenius permutation. All checks pass.
It also verifies that a transitive 12-cycle retains a nonzero class, guarding
against the irreducibility-only shortcut. These computations concern the
specified first-ruling witnesses, not every ruling or every central vector.

## 3. What this excludes, and what survives

If adjoining `sqrt(h)` to Q(C) gives an unramified quadratic extension,
then `div(h)=2D`; the class of D is a rational point of J[2]. If `J(Q)[2]=0`,
D is geometrically principal, so `h=c g^2` with `g in Q(C)` and `c in Q*`
(by descent of a principal rational divisor). Hence the extension is a
constant quadratic field extension or trivial, and is geometrically
disconnected. Thus the witnesses admit no geometrically connected
unramified double cover over Q. This conclusion does not require the
converse descent assertion, although these witnesses also have rational
points and consequently no Picard descent obstruction.

There are always 1023 nontrivial geometric double-cover classes. Over a
splitting field K of f, all are available, with `[K:Q]<=12!`. An explicit
basis of ten independent unramified squareclasses is

```
h_i=(x-alpha_i)/(x-alpha_12),                 1<=i<=10.
```

Each divisor is `2(P_i-P_12)`. Their independence follows from (1), since
no nonempty combination of these ten pairs is the full branch set.
Consequently their compositum has geometric degree `2^10=1024`, unramified
everywhere over C_K. Every elementary abelian geometric double-cover
compositum has degree at most 1024. This bound is independent of coefficient
height, but does not control a Runge place count or its integral models.
It also does not bound towers allowing new covers of the covering curves.

Moreover, when C has a rational base point P0, the pullback of
`[2]:J->J` along `P |-> [P-P0]` is an unramified degree-1024 cover defined
over Q. It is geometrically connected and has geometric deck group J[2],
but its deck group scheme over Q need not be constant. This standard
construction is described in Poonen's same notes, Section 7.2 (pullback
of an isogeny along the Albanese map); connectedness follows because the
Abel-Jacobi map identifies first homology modulo 2. In particular,
`J(Q)[2]=0` does not rule this cover out or give it ten rational quadratic
intermediate covers. Rational lifts, their number fields, and the places
needed for a covering-Runge application still require their own audit.

## 4. Scope of a generic conclusion and of the fair-weight hypothesis

There is an arithmetic generic statement on the rational-configuration
family containing either witness. Parameterize six rational unit-circle
nodes by six indeterminates, choose a rational normalization of the weights
(for example set the first weight to one), and keep one fixed nonzero minor
in constructing the hyperbolic frame. On a normal connected open parameter
space containing the witness, the resulting first-ruling branch cover is
a smooth family. A nonzero rational two-torsion section over its function
field would extend over this open space: the closure of the generic
section in the finite etale group scheme J[2] is finite birational over the
normal base, hence a section. It cannot meet the zero section, since equality
of sections of a finite etale scheme is both open and closed. Specializing
at the witness would contradict `J(Q)[2]=0`. Thus this family's generic
Jacobian has no nonzero two-torsion over its rational function field.

This is arithmetic genericity in that specified parameter family. It is
not a statement over the algebraic closure of its function field, not a
claim about every smooth specialization, and not a proof that the family
dominates every possible weight component. The witnesses also disprove
any universal construction of a nontrivial Q-defined unramified double
cover that is valid on every nonresonant circle configuration.

The full fair-core hypothesis currently proves nonvanishing proper subsums
and hence twelve simple branch points. The checked notes give no Galois
orbit or partition condition on those points for all such weights. These
witnesses satisfy nonvanishing proper subsums, but are not asserted to be
members of an asymptotic full fair family. Therefore neither zero nor
positive rational two-torsion has been established for every full fair
weight, and no uniform covering-Runge conclusion can be based on either
assertion at present.
