# Full numerical relations surject onto every active restriction space

At every active positive depth of the
[all-cut invariant hierarchy](all_cut_invariant_relation_hierarchy.md),
the full space of numerical zero relations has full restriction image at
each individual cut. Its integral image can have finite index, and that
index is exactly a ratio of numerical evaluation gcds. The quotient is
cyclic. These statements hold at any configuration of distinct rational
directions; the integral formulation uses integer binary representatives.

This identifies the distinction between all numerical relations and the
short relations to which the arithmetic minor bounds apply. Full rank
does not force a zero restriction image. A further argument would have
to control how short relations span these full images, or use additional
arithmetic information. No new uniform circle-point bound follows here.

## 1. Integral surjectivity of a single restriction

Let `m=2q`, `1<=k<=floor(m/6)`, and use the hierarchy's space `K_k`
of simultaneous `SL_2` invariants of degree two in each binary row.
Let `V_Z` be the lattice of invariants in `K_k` with integral polynomial
coefficients. Fix an inside set `S` of cardinality `q-k`. Specialize its
rows to `(0,1)` and all outside rows to `(1,z_j)`. The restriction map is

```text
T_S: V_Z -> W_Z,
W_Z={integral squarefree homogeneous polynomials of degree 2k
     in q+k variables, annihilated by sum_j partial_j}.
```

The hierarchy supplies an integral matching basis of `W_Z`. Each basis
element is a product of `2k` differences on `4k` distinct outside rows.
For each of those matching edges choose a different inside row, and
replace the edge by the triangle bracket product on its two endpoints
and that inside row. Pair the remaining `q-3k` inside rows with the
remaining `q-3k` outside rows, using squared brackets. Call the resulting
invariant `Q`.

Every row has degree two. Its connected components are `2k` triangles
and `q-3k` doubled edges, so a set on which the graph restriction does
not vanish contains at most one row from each component. It therefore
has size at most `q-k`. Thus `Q` belongs to `K_k`.

For a triangle with inside row `s` and outside edge `a,b`,

```text
[s,a][a,b][b,s] -> -(z_b-z_a).
```

There are `2k` triangles, so their signs cancel. Every squared crossing
bracket becomes one. Consequently `T_S(Q)` is the desired matching
basis vector. In particular the map is surjective **integrally**:

```text
T_S(V_Z)=W_Z.                                         (1)
```

## 2. The restriction kernel contains a nonvanishing monomial

There are at least two triangles because `k>=1`. In a lifted graph
from Section 1, exchange the inside vertex of the second triangle with
one outside vertex of the first triangle. At the fixed cut `S`, the
first triangle now has two inside vertices, and the second has none.
Every other component still has one inside vertex. The graph remains
a product of `2k` disjoint triangles and `q-3k` squared brackets, hence
still belongs to `K_k`, but its restriction at `S` is zero.

Denote this invariant by `Z`. At any distinct rational directions,
every bracket in `Z` is nonzero. Numerical evaluation `f` consequently
satisfies

```text
Z in ker T_S,       f(Z)!=0.                           (2)
```

Given any basis lift `Q`, the integral invariant

```text
f(Z) Q-f(Q) Z
```

is a numerical zero and restricts to `f(Z) T_S(Q)`. Therefore, writing
`Lambda=ker(f:V_Z->Z)`, the full rational numerical-zero space obeys

```text
T_S(Lambda tensor Q)=W_Z tensor Q.                     (3)
```

This is an individual-cut assertion for every `S`; it makes no claim
that all the restrictions can be prescribed independently at once.

## 3. Exact integral index and cyclic quotient

The following lattice identity applies without invariant-theory
assumptions. Let `V_Z` be a free integer lattice, let `f:V_Z->Z` be
nonzero, and let `T:V_Z->I` be onto a free integer lattice. Suppose
`f` is nonzero on `K=ker T`. Define positive generators

```text
f(V_Z)=b Z,       f(K)=a Z,       b divides a.
```

Evaluation induces a surjective homomorphism

```text
I -> b Z/a Z,       T(v) -> f(v) mod a Z.
```

This is well defined because two lifts differ by an element of `K`.
Its kernel is exactly `T(ker f)`: if `f(v)` is divisible by `a`,
there is `u in K` with `f(u)=f(v)`, and `v-u` is a zero lift.
Thus one has the exact cyclic quotient

```text
I/T(ker f) ~= b Z/a Z ~= Z/(a/b)Z.                    (4)
```

Applied to (1)--(2), this gives

```text
[W_Z:T_S(Lambda)]
   = gcd{f(Q): Q in ker T_S}/gcd{f(Q): Q in V_Z}.       (5)
```

All gcds refer to the entire indicated integral lattices, not just a
chosen collection of monomials. Common numerical contents cancel by
the displayed quotient; their prime powers are retained exactly.

If `t=rank W_Z`, express restrictions of an integral basis of `Lambda`
in the matching basis. The gcd of all `t` by `t` minors is exactly
`a/b`. The Smith invariant factors of this image lattice are
`1,...,1,a/b`. These assertions also cover `a=b`, when the integral
image is the whole target lattice. They do not assert that one selected
minor equals this gcd or that a small generating basis exists.

Thus full-lattice restriction-minor divisibility records precisely the
restricted evaluation-content ratio in (5). The short-relation proper
image bounds in the hierarchy are compatible with (3): the remaining
issue is control of the coefficient heights needed to attain full rank.
The formula does not rule out stronger covolume or collective-minor
arguments.

## Verification

The [exact checker](check_full_invariant_restriction_image_and_index.py)
constructs every matching basis lift at all active depths for
`m=6,8,10,12`, verifies its exact signed restriction and the required
vanishing cuts, and checks the nonzero numerical kernel witnesses.
It separately builds full integer numerical-relation bases by
unimodular operations, then verifies maximal-minor gcds and the cyclic
Smith criterion for primitive and nonprimitive evaluations. For target
lattices of nontrivial ambient index, it additionally checks the scaled
maximal-minor gcd and hence the relative index. The proofs above apply
in all the stated dimensions.

## 4. A realized single-cut spanning height can equal the cyclic index

Cyclicity does not by itself replace the index by its `t`th root in a
spanning-height bound. This already occurs in the actual six-row invariant
evaluation system, with target rank two. Take integer binary directions

```text
(s,1), (0,1), (1,1), (2,1), (3,1), (4,1),
s=5 mod 6, s>=5,
```

labelled `0,...,5`, and the cut `S={0,1}` at depth one. Write `[ab]`
for a binary bracket, and define eight integral bracket monomials:

```text
M4 =[01][02][12][34][35][45],
M6 =[01][02][13][24][35][45],
M7 =[01][02][14][24][35]^2,
M8 =[01][03][13][24][25][45],
M9 =[01][03][14][24][25][35],
M11=[02]^2[13][14][35][45],
M13=[02][03][13][14][25][45],
M14=[02][03][14]^2[25][35].
```

An integral basis of the saturated balanced kernel `K_1` is

```text
E1=M6-M11+M13,         E2=M6-M9-M11+M14,
E3=M4,                E4=M8,                E5=M6-M7+M9.
```

In the matching basis
`((z2-z3)(z4-z5),(z2-z4)(z3-z5))`, the restriction is `[I_2 0]`.
The checker certifies saturation by unimodular column operations on the
full fifteen-dimensional standard-monomial lattice and then an integral
restriction splitting. Direct bracket evaluation gives

```text
f(s)=(14s^2-38s+24, 14s^2-74s+96,
      2s^2-2s,     12s^2-24s,    20s^2-56s).
```

The gcd of the last three entries is exactly `4s`: the gcd of the
first two of those, after dividing by `s`, is
`2 gcd(s-1,6)=4`, and the third is divisible by four. All five entries
are divisible by four, while the first divided by four is congruent
to six modulo `s`. Since `gcd(s,6)=1`, the full evaluation gcd is
exactly four. Formula (5) therefore gives index `s`. More precisely,

```text
T_S(Lambda_s)={(x,y) in Z^2 : x+4y=0 mod s}.             (6)
```

This follows because the first two primitive evaluation entries reduce
to `6,24` modulo `s`. The quotient is cyclic, but any two rationally
independent vectors of (6) include one with sup norm at least `s/5`:
if both coordinates have absolute value less than `s/5`, the congruence
forces `x+4y=0` as an integer equality. All such vectors lie on one line.

The scale is attained by actual zero relations. In the displayed basis,

```text
u=(-4,1,3,3,0),       v_s=(s,0,12-7s,0,0)
```

both evaluate to zero, and their images `(-4,1)` and `(s,0)` form a
basis of (6). Thus the least coefficient sup norm in this fixed basis
needed for relations spanning this cut is `Theta(s)`, with bounds
`s/5` and `7s`. There is
also the independent constant relation `(0,0,8,-3,1)`, whose restriction
vanishes. The full relation-lattice covolume is `Theta(s^2)`, by its
primitive numerical evaluation vector; this is compatible with the
spanning-height assertion.

These are actual distinct rational directions and can be realized by
Gaussian denominator clearing as six points on an integer circle.
They are not a shrinking bounded-endpoint family: the fixed directions
`0,1,2,3,4` keep a positive angular span. No full uniform conductor
profile is asserted. The example excludes only a single-cut shortcut
from cyclic Smith factors or redundant lifts to a root-of-index spanning
bound; it does not exclude a collective argument using other cuts.

The checker additionally verifies the saturated basis and exact
polynomial evaluation, both displayed zero relations, and twenty
integer specializations of the index and image generators.
