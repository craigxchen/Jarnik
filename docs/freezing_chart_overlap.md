# Overlapping reconstruction charts and a one-row compatibility construction

This note tests whether changing the frozen 49-block family can force
incompatibility of the reconstructed candidates. Two new facts narrow
that route. First, 49 is optimal for the strict whole-cut majority
criterion along **any connected graph** on eight labels, not only an
anchor star. Distinct minimal charts therefore lose a necessary
block-scale margin on their overlap. Second, there are actual Gaussian
integers of imaginary part one with arbitrarily many independent,
nearly equal-weight factors. Every majority chart for such a row is
simultaneously feasible. A successful gluing argument must consequently
use how the same factors occur in different rows; consistency among
one row's many charts is not enough.

These are scoped statements about reconstruction. The construction
below is not an eight-row endpoint configuration, and the uniform
lattice-point bound remains unproved.

## 1. The connected-graph optimum is also 49

There are eight labels and 128 unoriented cuts, represented by the
subsets not containing label zero. Let `T` be any spanning tree on
the labels. The map

```text
cut S -> (indicator that edge e crosses S : e in T)
```

is a bijection from these 128 cuts to `{0,1}^7`. Indeed the bit at
each vertex is recovered from the edge bits by taking their parity
along the unique path from zero. Consequently the numbers of cuts
crossing `j` tree edges are exactly `binom(7,j)`, independently of
the shape of the tree.

To obtain the strict majority margin on every tree edge at equal
block weights, at least 33 of its 64 separating cuts must be frozen.
But the maximum possible total number of tree-edge incidences of
48 cuts is

```text
1*7+7*6+21*5+19*4=230 < 7*33.                    (1)
```

Every connected graph contains a spanning tree. Therefore **no
family of 48 whole cuts meets the strict-majority criterion on any
connected graph**. The existing 49-cut anchor design attains this
minimum. This extends the restricted optimality statement in
[fano_freezing_audit.md](fano_freezing_audit.md).

In particular, if `F` and `F'` are two distinct 49-cut charts and
`I=F intersect F'`, the graph of edges with at least 33 separating
cuts in `I` is disconnected. The same conclusion holds for the
intersection of any collection of such charts containing two
distinct ones. The common data do not form a connected
reconstruction chart by this criterion.

The restriction matters. At near-uniform weights, an edge with
exactly 32 frozen cuts might acquire a favorable smaller weight
error. Equation (1) rules out a fixed positive block-scale margin
on all connecting edges; it does not forbid every possible
fine-weight improvement. Nor does it rule out additional arithmetic
information beyond the Gaussian ideal window.

## 2. The exact overlap determinant and the needed extra divisor

Suppose two candidate completions agree on all blocks in `I`.
For a pair of labels write `P_ij` for the common oriented product
of the separating blocks in `I`. The two primitive pair numerators
`h=x+i t` and `h'=x'+i t'` are multiples of `P_ij`. Therefore

```text
N(P_ij) divides x t'-x't,
|x t'-x't| <= 2H T,                              (2)
```

when both candidates have modulus at most `H` and residue modulus
at most `T`. If `N(P_ij)>2HT`, their directions coincide. On a
connected graph these equalities glue the configurations exactly.
This is the determinant form of the sparse reconstruction lemma.

For eight uniform rows, `log H=32w+o(w)` and `log T=o(w)`. If an
overlap edge has 32 frozen separating blocks, then

```text
log N(P_ij)=32w+o(w).                            (3)
```

Thus (2) is precisely at the leading threshold. To propagate
uniqueness across this edge by the same integer-determinant
mechanism, one needs additional information. A concrete sufficient
new input would be a positive integer divisor `J_ij` of the
determinant, coprime to `N(P_ij)`, with

```text
log J_ij >= delta w,
```

for fixed `delta>0`. An equivalent improvement in the available
height window could serve the same purpose. Neither independent
chart uniqueness nor the list of shared blocks supplies this
additional divisor.

The loss of strictness is not just an artifact of the upper bound.
For every positive even integer `n`, put

```text
P=n+i,       q=N(P)=n^2+1,
h_+=q+2n+2i=P(n+2-i),
h_-=q-2n-2i=P(n-2-i).                            (4)
```

Both `h_+` and `h_-` are conjugate-primitive: their real coordinates
are odd and their imaginary coordinates are `2` and `-2`.
Their Gaussian gcd is exactly `P` up to a unit, since the two
cofactors differ by four and have odd norms. They satisfy

```text
|h_+|~|h_-|~q,       |Im h_+|=|Im h_-|=2,
Im(h_- bar(h_+))=-4q.                            (5)
```

Thus a critical shared divisor with `N(P)~H` permits two distinct
primitive directions with bounded residues and asymptotically
equal moduli. Each full divisor `h_+` or `h_-`, by itself, is in
a strict sparse window with cofactor one. This is a two-candidate
window example, not a full 49-block construction.

## 3. What two charts actually recover in one row

For one row suppose `Q=D A` and `Q'=D B` are its two fixed
divisors, with `gcd_G(A,B)=1`. If their candidates glue to the
same primitive numerator `h`, their cofactors satisfy exactly

```text
h=Q U=Q' U',
U=B V,       U'=A V,
h=lcm_G(Q,Q') V.                                (6)
```

There is only one copy of each common factor in (6). If the row
has 64 uniform incident blocks, the charts each freeze 33, and
their incident overlap has `d` blocks, then their union has
`66-d` blocks and

```text
log|V|=(d-2)w/2+o(w)                             (7)
```

in an actual completion, with an upper bound of the same form
when only the common window is specified. At the smallest
possible overlap `d=2`, the union leaves just the small row
correction. This is a useful exact recovery of factors, but it
does not yet force that correction or the residue to be large.
The next construction shows that even all these one-row recovery
conditions can hold simultaneously.

## 4. Balanced independent factors with imaginary part one

For every fixed positive integer `B`, there is a sequence of
conjugate-primitive Gaussian integers

```text
h=K product_(s=1..B) H_s,       Im h=1,
```

whose blocks are pairwise coprime and coprime to all their
conjugates, and which satisfy

```text
w=(1/B)sum_s log N(H_s) -> infinity,
max_s |log N(H_s)/w-1| -> 0,
log|K|=o(w).                                    (8)
```

Here is a construction retaining the full arithmetic denominator
and common-factor costs.

### Gaussian cyclotomic factors

Choose an odd squarefree integer `D`. Over `Q(i)`, the polynomial

```text
x^D+i
```

is a product of monic, pairwise distinct polynomials indexed by
the divisors `d` of `D`, with degrees `phi(d)`. To see this
directly, each root has order `4d` for a unique `d|D`. Among the
primitive `4d`-th roots, exactly `phi(d)` satisfy the required
equation `alpha^D=-i`. They form one orbit over `Q(i)`: the
automorphisms fixing `i` act transitively on either specified
value `alpha^d=i` or `-i`. The orbit polynomial has coefficients
in `Q(i)` which are algebraic integers, hence in `Z[i]`.

Its degrees sum to `sum_(d|D)phi(d)=D`. Since `D` is squarefree,

```text
max_(d|D) phi(d)=phi(D).
```

Take `D` to be a product of successively more odd primes. Then
`phi(D)/D=product_(p|D)(1-1/p)` tends to zero. This last fact
needs only the elementary divergence of the harmonic series:
the reciprocals of the integers supported on those primes sum
to the inverse Euler product, and eventually include every odd
integer below any fixed bound.

### Grouping by degree and removing overlaps

Greedily distribute the polynomial factors among `B` groups,
always adding the next factor to a group of least total degree.
The final group degrees `d_s` have range at most `phi(D)`, so

```text
|d_s-D/B| <= phi(D).                             (9)
```

Once `phi(D)<D/B`, every group is nonempty. Write its monic
Gaussian polynomial as `F_s(x)`. Distinct groups have no common
complex root. For every integer `n`, their evaluated Gaussian
gcds divide their fixed nonzero polynomial resultants. Therefore
there is a finite constant `C_D`, independent of `n`, such that

```text
sum_(s<t) log N(gcd_G(F_s(n),F_t(n))) <= C_D.      (10)
```

Evaluate at a sufficiently large even positive integer `n`.
The complete product is

```text
h=n^D+i,
```

which is conjugate-primitive: its coordinates are coprime, and
its norm is odd. In particular no factor of one `F_s(n)` shares
a Gaussian prime with the conjugate of another.

At each Gaussian prime occurring in more than one group, retain
its largest group exponent in one block and transfer every other
group exponent into `K`. Call the resulting blocks `H_s`.
For nonnegative exponents, the sum other than their maximum is
at most the sum of all pairwise minima. Thus (10) proves

```text
log N(K) <= C_D.                                (11)
```

The `H_s` are pairwise coprime, remain coprime to all conjugates,
and satisfy `h=K product H_s` exactly. Since the polynomials
are monic,

```text
log N(H_s)=2d_s log n+O_D(1),
w=2D log n/B+O_D(1).                            (12)
```

First let `D` tend through the above sequence, and choose each
even `n=n(D)` large enough that `C_D=o(w)`. Equations (9)--(12)
then prove (8). No overlap or resultant factor has been dropped
without being included in `K`.

### Every majority chart is simultaneously feasible

Fix a subset `F` of more than `B/2` incident blocks. Put

```text
Q_F=product_(s in F) H_s,
U_F=K product_(s outside F) H_s.
```

In the common block-weight and correction window from (8),

```text
2|U_F|<|Q_F|
```

eventually, uniformly over all these finitely many subsets. Thus
every such chart has the same primitive candidate `h`, with
residue exactly one. For `B=64`, all 33-block and larger windows
are simultaneously feasible. Their mutual intersections, their
Gaussian lcm identities, and all the corresponding last-convergent
conditions are consequently consistent with negligible correction
and residue heights.

This construction concerns one row and its incident factors. It
does not assign these factors to the overlapping incidences of
seven different rows, does not produce the required cross-anchor
residues, and does not give a counterexample to endpoint uniformity.
It shows why combining that row's reconstruction charts, even
using their exact factor arithmetic, cannot alone force a positive
block-scale residue or correction.

## 5. The remaining gluing input is genuinely across rows

On a connected component determined by common frozen data, all
relative phases are known. The missing rows still permit different
assignments of prime factors among cuts having the same restriction
to that component. Such assignments leave the component's relative
configuration unchanged. Independent charts do not automatically
produce a second completion or prove that a proposed splice remains
inside the small-angle window.

A useful new theorem would have to charge this freedom: for example,
prove the extra determinant divisor described after (3), or give a
positive lower bound for a cross-row residue arising from two
different assignments with the same known component. The one-row
construction (8) shows that the theorem must use shared incidence
across rows. The spanning-tree calculation (1) identifies why the
existing strict-majority windows do not already supply it through
their intersections.

## Verification

Exact finite checks verified the tree-to-cut bijection and the value
230 for a star, a path, and a branching tree; the proof above covers
every tree. They also verified 100 members of (4)--(5), and the
Gaussian cyclotomic factorization and divisor degrees for 28 factors
with odd squarefree `D` up to 165, using monic polynomial division
over `Z[i]`. The arbitrarily fine balance in (8) is supplied by the
proved limiting construction, not by those small finite examples.
