# A pair-reflection domain has bounded endpoint capacity

The exact integrality domain of a pair-bisector reflection is the Gaussian
sublattice computed in
[the rational reflection route](rational_reflection_route.md).  This note
adds an endpoint consequence that applies to the domain subset, without
reflecting or clearing denominators on the whole source cluster.

Section 4 of
[the midpoint-scaling obstruction](midpoint_reflection_scaling_obstruction.md)
left open whether short-arc cardinality can force original rows into one
pair's integrality domain.  The bounds below resolve the stronger
positive-fraction version negatively and reduce the literal third-row
question to pair denominators bounded in terms of `C`.  The whole-union
scaling obstruction and its scope are unchanged.

For a fixed normalized arc constant, every individual reflection is integral
on only boundedly many selected rows.  The rows which it carries back into
the selected cluster satisfy a stronger bound.  Thus large cardinality cannot
force one pair reflection to act integrally on a positive fraction of the
cluster.  This does not bound the original cluster and does not control a
number of different reflections growing with its cardinality.

## 1. Exact domain, including prime powers

Let `Z` be a set of Gaussian integers on the circle of radius `R`, so their
common squared norm is `N=R^2`.  Suppose `Z` lies in an arc of angular width
`Delta<pi/2`, and put

```text
C = Delta sqrt(R).
```

Choose distinct `a,b in Z`.  They are nonantipodal under the angular
hypothesis.  Their midpoint reflection is

```text
sigma(z) = rho conjugate(z),       rho=ab/N.
```

Write the reduced norm-one coefficient in the canonical form

```text
rho = eta d/conjugate(d),
eta in {1,i},
gcd_G(d,conjugate(d))=1.                              (1)
```

As proved in `rational_reflection_route.md`, cancellation in (1) gives the
exact statement

```text
sigma(z) in Z[i]   iff   d divides z.                  (2)
```

Put `D=Norm(d)`.  Since the domain contains the anchors, `D` divides `N`.
On the full radius-`R` lattice circle, the largest possible integrality set
is therefore not merely contained in a congruence class; it is exactly

```text
d {w in Z[i] : Norm(w)=N/D}.                           (3)
```

In particular its full-circle cardinality is `r_2(N/D)`.  For several
reflections with reduced factors `d_1,...,d_q`, the set on which all of them
are integral is exactly

```text
lcm_G(d_1,...,d_q) Z[i].                               (4)
```

There is also a direct nested prime-power description.  First make the
equal-norm tuple Gaussian primitive.  At a split rational prime
`p=pi conjugate(pi)`, write its row valuations as

```text
e_z=v_pi(z),       0<=e_z<=W,
v_conjugate(pi)(z)=W-e_z.
```

If the two anchor values are `e_a,e_b` and `s=e_a+e_b`, then (2) is
equivalent at this prime to

```text
max(0,s-W) <= e_z <= min(W,s).                         (5)
```

Indeed the two valuations of `rho` are `s-W` and `W-s`; its positive
valuation is assigned to `d` and its negative valuation to
`conjugate(d)`.  Formula (5) retains arbitrary depths.  In the squarefree
case `W=1`, it says that a domain row must agree with the anchors at every
column where the anchors agree, while a column where they differ imposes no
condition.

## 2. Uniform capacity of the domain inside an endpoint arc

Let

```text
S_ab = {z in Z : sigma(z) is Gaussian integral}
     = Z intersect d Z[i].
```

More precisely,

```text
ceil(|S_ab|/2)-1 <= C^2/(2 sqrt(D)).                  (6a)
```

In particular,

```text
|S_ab| <= 2 + C^2/sqrt(D).                             (6)
```

Since the left side is integral, the right side may of course be replaced
by its floor.

To prove (6), divide only the domain points by `d`:

```text
X=S_ab/d subset Z[i],       |x|=r=R/|d|=R/sqrt(D).
```

Under this division the rational reflection becomes the lattice reflection

```text
tau(x)=eta conjugate(x),                                  (7)
```

because `sigma(dx)/d=eta conjugate(x)`.  The midpoint axis of the anchors is
inside the original arc.  Hence

```text
U = X union tau(X)
```

lies in an arc centered on that axis, of angular width at most `2 Delta`.
It is a finite integral equal-radius set invariant under the coordinate or
diagonal lattice reflection (7).

Projection onto its axis has levels separated by at least one.  For a
coordinate axis these levels are integers.  For a diagonal they are
`(x+y)/sqrt(2)` or `(x-y)/sqrt(2)`.  Every point of `U` has squared norm
`N/D`, and

```text
x+y == x-y == x^2+y^2 == N/D       (mod 2).
```

Thus all diagonal projection numerators have the same parity, so distinct
levels are separated by `2/sqrt(2)=sqrt(2)`.  This argument needs no
primitive normalization.  Each level contains at most two circle points.
If `M=|U|` and its centered containing arc has width `Theta<pi`, radial
projection therefore gives

```text
r (1-cos(Theta/2)) >= ceil(M/2)-1.
```

Using `1-cos t<=t^2/2`,

```text
M <= 2 + (Theta sqrt(r))^2/4.                           (8)
```

Here `Theta<=2 Delta`, so

```text
Theta sqrt(r)
 <= 2 Delta sqrt(R/sqrt(D))
 = 2C/D^(1/4).
```

Equations (8) and this inequality give `|S_ab|<=|U|` and then (6).
Keeping the integer ceiling before weakening gives (6a).  No primitivity of
the source cluster was assumed.

There is also a stronger form containing no angle approximation.  If `t` is
the dot product of the two endpoints of the source arc, then
`cos(Delta)=t/N`.  Keeping the sagitta before using
`1-cos x<=x^2/2` gives

```text
|S_ab|-2 <= 2(N-t)/sqrt(ND).                           (8a)
```

Thus (8a) can be checked using integers alone by squaring:

```text
(|S_ab|-2)^2 N D <= 4(N-t)^2.                          (8b)
```

If the reflection is not globally lattice-preserving, then `D>=5`.  Indeed
the only possible positive Gaussian norms strictly between one and five are
two and four.  Norm-two elements are associates of `1+i`, while norm-four
elements are associates of `2`; neither can be coprime to its conjugate.
Thus every nonunit pair reflection satisfies
the slightly sharper uniform statement

```text
|S_ab| <= 2 + C^2/sqrt(5).                              (9)
```

A useful denominator corollary is

```text
|S_ab|>=3   implies   D<=C^4/4.                         (9a)
```

Indeed the left side of (6a) is at least one.  Consequently every pair with
`D>C^4/4` has exactly its two anchors in the source domain.  Any attempt to
force even one third row can therefore restrict attention to the finitely
many Gaussian denominator ideals of norm at most `C^4/4`.

## 3. Rows actually paired inside the source

The set on which the reflection acts internally is

```text
A_ab = {z in Z : sigma(z) in Z}.
```

It is `sigma`-invariant and is contained in `S_ab`.  After division by `d`,
the set `A_ab/d` is invariant under (7).  Unlike the union used above, it is
already contained in the original rotated arc.  Symmetry about the anchor
axis shows that its centered containing width is at most `Delta`: if the
original lifted interval is `[-L_-,L_+]`, both `t` and `-t` must lie in it,
so `|t|<=min(L_-,L_+)` and `2 min(L_-,L_+)<=Delta`.

Applying the same radial-projection argument first gives

```text
ceil(|A_ab|/2)-1 <= C^2/(8 sqrt(D)),                   (10a)
```

and hence

```text
|A_ab| <= 2 + C^2/(4 sqrt(D)).                          (10)
```

In particular,

```text
|A_ab|>=3   implies   D<=C^4/64.                       (10b)
```

The corresponding exact sagitta bound is

```text
|A_ab|-2
 <= 2 sqrt(N/D)-sqrt(2(N+t)/D).                        (10c)
```

Writing `k=|A_ab|-2` and `B=2(N-t)-k^2D`, (10c) is equivalent to the two
integer inequalities

```text
B>=0,       B^2>=8k^2D(N+t).                           (10d)
```

This counts fixed rows as well as moved rows, so it also bounds the number
on which the reflection acts nontrivially within the source.

Consequently, for every fixed `C` and every endpoint cluster of `m` rows,

```text
max_(a!=b) |S_ab|/m <= (2+C^2)/m,
max_(a!=b) |A_ab|/m <= (2+C^2/4)/m.                    (11)
```

Both ratios tend to zero along any hypothetical sequence with `m` tending
to infinity.  More generally, the union of the domains of `q` chosen pair
reflections contains at most `q(2+C^2)` source rows.  Thus a bounded number
of such reflections cannot touch a positive fraction of a growing cluster.

There is a direct descent interpretation.  Discarding the rows outside
`S_ab` and dividing the retained rows by `d` produces an integral cluster
of radius `R/sqrt(D)` with the same angular width, but (6) says that this
radius decrease retains only boundedly many rows.  Requiring simultaneous
integrality for several reflections gives the intersection (4), which is a
subset of every individual domain and hence cannot improve the retained
cardinality.  Combining differently scaled domain images by a further
common normalization is a different operation and is not covered by this
argument.

In the squarefree binary code, (5) says precisely that a third row lies in
the pair domain when it belongs to the coordinatewise descendant of the two
anchors: it agrees with their common bit at every coordinate where they
agree.  Asking every pair domain to contain only its anchors is the binary
`2`-frameproof condition.  This condition alone gives no row bound here.
For example, normalized mutually orthogonal Hadamard rows are already
`2`-frameproof: for any three rows, on exactly one quarter of all coordinates
the first two agree and the third has the opposite bit. Short-arc
arithmetic, rather than this incidence condition by itself, remains
essential.

The angular restriction is explicit.  For endpoint arcs of fixed `C`, one
has `Delta=C/sqrt(R)<pi/2` once `R>(2C/pi)^2`.  The argument makes no claim
for unconstrained full-circle prime-code profiles, and it does not exclude
an adaptive construction using a number of different reflections that grows
with the row count.  It supplies a negative answer only to the proposed
single-reflection or bounded-many-reflection forcing step; it is not a proof
of the uniform endpoint bound.

## Exact finite verification

Run

```text
python3 docs/check_pair_reflection_domain_capacity.py
```

The checker exhausts cyclic consecutive arcs of angular width less than
`pi/2` on every actual Gaussian circle of squared norm at most 500.  For
every anchor pair it reduces `ab/N` by exact Gaussian Euclidean arithmetic,
checks the exact integrality domain independently from the reflected
quotient, and constructs `S_ab` and `A_ab`.  It then verifies the stronger
integer inequalities (8b) and (10d), which imply (6) and (10).  The cyclic
angular order and the condition `Delta<pi/2` are also decided by exact
half-plane, cross-product, and dot-product comparisons.  The current scan
checks 2,584 arcs and 7,064 anchor pairs, including 4,248 nonunit domains and
2,336 domains containing at least three selected rows.
