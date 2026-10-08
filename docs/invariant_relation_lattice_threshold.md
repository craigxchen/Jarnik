# The lattice of exact degree-two invariant relations

Complementary balanced cuts strengthen the proved coefficient threshold
from `w` to `2w`. The balanced evaluation map has exactly
`s=binomial(m,m/2)/2` independent rows. Nevertheless, the unconditional
successive-minima estimate does not force an exact relation outside its
kernel below that threshold. Both the full evaluation height and the
height restricted to the kernel can be calculated, and an explicit CRT
model realizes those heights and all the balanced congruences with its
first outside relation at precisely `exp(2w+O_m(1))`.

The new threshold and height statements concern the actual integer
system. The CRT model concerns only its linear arithmetic constraints;
it is not asserted to satisfy the nonlinear Plücker identities or to
come from points on a circle. No uniform bound follows here.

The later simultaneous smaller-cut theorem in
[two_small_invariant_relations.md](two_small_invariant_relations.md)
does exclude this CRT model from the actual relation system: the model
has several constant relations in the balanced kernel, whereas the
actual six-row small-relation space has dimension at most one.
The precise balanced-only scope of the model below is unchanged.

The exact certificate is
[check_invariant_relation_lattice_threshold.py](check_invariant_relation_lattice_threshold.py).

## 1. The invariant space and exact balanced rank

Let `m=2q>=4`. Let `V` be the rational simultaneous `SL_2` invariant
space of degree two in every binary row. Use the standard monomial
lattice: its basis consists of bracket products encoded by semistandard
two-row tableaux, each label appearing twice. Integer straightening
of bracket products follows from the integral Plücker relation. Every
basis element is a product of `m` brackets and has polynomial coefficient
sum at most `2^m`.

The standard monomial and multidegree descriptions are recalled in
[Patrias--Pechenik--Striker, Section 2](https://www.mat.univie.ac.at/~slc/wpapers/FPSAC2022/55.pdf).
The calculations below are independent of the choice of this basis.

Write `r=dim V`. The weights of `Sym^2(Q^2)` are `-1,0,1` after dividing
the usual weights by two. In a semisimple `SL_2` representation, each
nontrivial even highest-weight summand contributes equally to weights
zero and one, and the trivial summands contribute only to weight zero.
Thus

```text
r=[x^m](1+x+x^2)^m-[x^(m-1)](1+x+x^2)^m.                (1)
```

These dimensions are `Theta(3^m/m^(3/2))`.

For each unordered balanced partition `S|S^c`, evaluate at `e_1` on
`S` and `e_2` on its complement. These are integer linear forms `c_S`.
Swapping the two vectors multiplies an invariant by `(-1)^m=1`, so
`c_S=c_(S^c)`. Let `C:V->Q^s` be this evaluation map, where

```text
s=binomial(m,q)/2.
```

**The map has rank exactly `s`.** To prove this, consider the invariants
`product_(ij in M) Delta_ij^2`, indexed by perfect matchings. Their
balanced values are zero or one, according as the matching fails or
succeeds in crossing the cut. The Gram matrix of these incidence rows
has entries

```text
G_(S,T)=a!(q-a)!,       a=|S intersect T|.                (2)
```

Indeed a matching crossing both cuts pairs `S intersect T` with the
complement of their union, and pairs the other two parts with each other.

Here is a direct positive-definiteness proof, avoiding a rank assertion
from numerical examples. Extend a vector on unordered cuts to a
complement-symmetric vector `b` on all `q`-subsets. Use

```text
a!(q-a)!=(q+1)! integral_0^1 t^a(1-t)^(q-a) dt.
```

For `t>1/2`, put `z=t/(1-t)>1`. The matrix
`K_z(S,T)=z^|S intersect T|` satisfies

```text
K_z=sum_U (z-1)^|U| 1_(U subset S)1_(U subset T)
    >=(z-1)^q I.
```

The last terms, with `|U|=q`, supply the diagonal bound. Pairing `t`
with `1-t` on the complement-symmetric subspace and integrating shows

```text
G >= (q!/2) I                                           (3)
```

on unordered cuts. In particular it is positive definite. Since these
squared matching invariants lie in `V`, the full balanced map has rank
`s`, as claimed. Put `K=ker C`, of dimension `k=r-s`.

## 2. Complementary moduli multiply

Use the actual primitive central numerators `P_i=K_i A_i`. The exact
balanced divisibility already proved in
[small_invariant_relations_at_balanced_cuts.md](small_invariant_relations_at_balanced_cuts.md)
says that an integer invariant `Q` with numerical value zero satisfies

```text
M_S | c_S(Q),
M_S=N(H_S)/N(gcd_G(H_S, product_(j outside S)K_j^2)).     (4)
```

The norm blocks are pairwise coprime, so `M_S` and `M_(S^c)` are coprime.
Because their balanced scalars agree, (4) gives the stronger exact fact

```text
M_S M_(S^c) | c_S(Q).                                   (5)
```

In particular, if `Q` is outside `K`, its polynomial coefficient sum
obeys

```text
log C_Q >= min_balanced(log M_S+log M_(S^c))
        >=2(1-eta)w-4 sum_i log|K_i|
        >=2(1-eta)w-4m sigma.                            (6)
```

Here `log|K_i|<=sigma`; the separate letter distinguishes this
correction-height bound from the balanced rank `s`.
No residue-height loss occurs in (6), and arbitrary prime powers are
retained. Every exact relation of coefficient height below this bound
lies in `K`.

## 3. The two primitive numerical heights

Let the actual brackets have the exact factorization

```text
|Delta_ij|=b_ij product_(S containing i,j)n_S,
(1-eta)w<=log n_S<=(1+eta)w,       log b_ij<=beta.
```

For a cut of cardinality `a`, the smallest number of internal edges
among degree-two-each bracket graphs is

```text
g_0(a)=max(0,2a-m).
```

For `m>=6`, the corresponding minimum for invariants in `K` is

```text
g_1(a)=g_0(a)+1_(a=q).                                  (7)
```

The lower bound at `a=q` is exactly the vanishing balanced coefficient.
The other lower bounds follow by invariant weight balance. All minima
are attained by actual bracket monomials: for `g_0`, use a squared
perfect matching with as many crossing pairs as possible. For `g_1`,
use two disjoint triangle products and square the matching edges on the
remaining rows. At a balanced cut, distribute the triangle vertices as
one and two inside, producing exactly one internal edge. Below balance,
use at most one inside vertex per triangle and pair remaining inside
vertices across the cut. The cases with zero or one inside vertex are
handled by putting the other triangle wholly outside. Above balance,
use `g(a)-g(m-a)=2a-m`.

Every such two-triangle monomial belongs to `K` and is nonzero at a
configuration of distinct directions. Hence the numerical evaluation
restricted to `K` is a nonzero functional.

Let `G_0` be the gcd of numerical values of an integral basis of `V`,
and `G_1` the gcd of all values on the saturated integer kernel `K`.
Define

```text
D_j=product_S n_S^g_j(|S|),       B=product_(all edges)b_ij.
```

The exact gcd bounds are

```text
D_j | G_j,       G_j | D_j B^2,       j=0,1.              (8)
```

For the upper bound, at each prime choose one of the monomials attaining
the corresponding minimum. Its extra bracket-residue valuation is at
most twice the sum over all edges; an edge appears at most twice. At a
prime outside the core, the same estimate applies. For the lower bound,
weight balance gives the indicated Gaussian ideal powers for every
invariant in the relevant lattice; its integer value is real, so the
conjugate divisor supplies the full norm power.

Both the whole space and its kernel have raw numerical maximum of
logarithmic size `m 2^(m-2)w+o(w)`. For the kernel use any of the nonzero
two-triangle monomials for the lower bound; expressing it in a fixed
kernel basis costs only a constant depending on `m`. Therefore their
primitive evaluation heights are

```text
h_V=A_m w+o(w),       h_K=(A_m-2s)w+o(w),
A_m=m 2^(m-2)-(m/2)binomial(m,q).                       (9)
```

The binomial identity in (9) is
`sum_a binomial(m,a)max(0,2a-m)=(m/2)binomial(m,q)`.
Equations (8) and the raw maximum bounds give finite estimates if needed;
the errors are bounded explicitly by the core-weight errors, twice
`sum log b_ij`, and fixed basis constants. No cancellation estimate for
a sum of graph monomials is being assumed.

## 4. What successive minima actually force

Fix the integral invariant basis and divide its numerical value vector
by its gcd, obtaining a primitive vector `f in Z^r`. Put

```text
Lambda={a in Z^r : f dot a=0},       n=r-1.
```

Its Euclidean covolume is `||f||_2`. The restricted relation lattice
`Lambda intersect K` has rank `k-1`, since evaluation on `K` is nonzero.
Thus `k` independent exact relations must include one outside `K`.

Let `lambda_1<=...<=lambda_n` be the Euclidean successive minima.
Every nonzero integer vector has norm at least one. Minkowski's second
theorem, with `v_n` the volume of the Euclidean unit ball, gives

```text
product_i lambda_i <= (2^n/v_n)||f||_2,
lambda_k <= ((2^n/v_n)sqrt(r)||f||_infinity)^(1/s).       (10)
```

There are exactly `s` factors at or above index `k`; the denominator
is neither `r-1` nor `s-1`. The corresponding polynomial coefficient
sum costs at most an additional factor `r 2^m`.

For `m>=6`, the guaranteed outside coefficient exponent from (10) is
`A_m/s`. It exceeds the arithmetic threshold `2`:

| `m` | `r` | balanced rank `s` | `A_m` | `A_m/s` |
|---|---:|---:|---:|---:|
| 6 | 15 | 10 | 36 | 3.6 |
| 8 | 91 | 35 | 232 | 6.628571... |
| 10 | 603 | 126 | 1300 | 10.317460... |
| 12 | 4213 | 462 | 6744 | 14.597402... |

Indeed `A_m/s=m(2^(m-1)/binomial(m,q)-1)`, asymptotic to a positive
constant times `m^(3/2)`. For `m=4`, `K=0`; use the first minimum in
rank `r-1=2`, giving `A_4/(r-1)=4/2=2`, exactly the arithmetic threshold.
The expression `A_4/s=4/3` is inapplicable.

One can guarantee many small relations, but the relevant rank deficit
remains. For a coefficient-norm threshold `H>1`, (10) implies at least

```text
n-floor(log((2^n/v_n)||f||_2)/log H)
```

independent relations of norm at most `H`, if this is positive. At
`log H=epsilon w`, its leading lower bound is
`r-1-A_m/epsilon`. Exceeding the allowed kernel-relation rank `r-s-1`
would require `epsilon>A_m/s`, which is incompatible with the proved
small-relation range `epsilon<2`.

The ratio of the two covolumes has leading logarithm `2s w`, by (9)
and fixed-subspace covolume formulas. This does not by itself bound
short lifts from the quotient lattice: the internal lattice can have
long successive minima. In particular one may not replace (10)
unconditionally by a quotient estimate of exponent `2`.

## 5. An exact integer model at the complementary threshold

There is a more concrete warning than a dimension count. It uses the
actual integer balanced matrix `C`, but does not impose nonlinear
identities on its numerical evaluation vector.

Choose a unimodular column basis with
`C=[D 0]`, where `D` is an invertible `s` by `s` integer matrix and
the last `k` coordinates span `K`. Each row of `D` is primitive: each
balanced evaluation takes value one on some squared matching invariant.
Assume `k>=2`, as holds for the degree-two spaces with `m>=6`.

Choose pairwise-coprime integers `M_j>=2`, put `M=product_j M_j`, and
write the CRT idempotent

```text
E_j=(M/M_j)u_j,       1<=u_j<M_j,
E_j=1 mod M_j,       E_j=0 mod M_l for l!=j.
```

For any integer `N>=1`, define an evaluation vector

```text
f=(sum_j E_j D_j,       MN, M, 0,...,0).                 (11)
```

This vector is primitive. Its gcd divides `M` because one coordinate
is `M`, while modulo every prime divisor of `M_j` its first `s`
coordinates are a unit multiple of the primitive row `D_j`.
Its global and restricted primitive heights are

```text
h(f)=log M+log N+O_m(1),       h(f|K)=log N.              (12)
```

Any exact relation `f dot a=0` obeys
`D_j a_first=0 mod M_j`. If
`s||D||_infinity ||a||_infinity<min M_j`, all these integers vanish,
so the relation lies in `K`. Conversely, for any `j`, take

```text
a_first=M_j adj(D)e_j,
a_(s+1)=0,       a_(s+2)=-u_j det(D),
all other coordinates zero.
```

This is an exact relation outside `K`, of norm at most a fixed
`D`-dependent constant times `M_j`. Therefore its first outside-relation
scale is bounded above and below by constant multiples of `min M_j`.

With `log M_j=2w+O_m(1)` and
`log N=(A_m-2s)w+O_m(1)`, (12) reproduces both heights in (9), every
complementary congruence, and the exact `2w` coefficient threshold.
All these data are ordinary integers. They do not assert that the
coordinates in (11) are simultaneous values of the invariant basis
at any binary-row configuration. That missing nonlinear compatibility
is precisely what rank, covolume, and these congruences have not used.

The subsequent [two-relation theorem](two_small_invariant_relations.md)
does use simultaneous first-smaller-cut congruences and excludes this
particular CRT model for six rows. Its restricted functional has three
independent constant-coefficient zero relations, whereas actual rows
admit at most one relation below the new `w/2-o(w)` threshold.
The CRT construction remains a countermodel to the linear package
specified here; it is not a countermodel to those additional constraints.

## 6. Fixed-depth higher jets do not remove the dimension issue

At a balanced cut, put the inside rows at `e_1+t_i e_2` and the
outside rows at `e_2+z_j e_1`. Diagonal weight balance forces every
monomial `t^alpha z^beta` in a degree-two-each invariant to satisfy

```text
|alpha|=|beta|,       0<=alpha_i,beta_j<=2.
```

Let `K_J` impose vanishing of all terms with `|alpha|<J` at every
balanced cut. Put `a_(q,j)=[x^j](1+x+x^2)^q`. Its codimension is at most

```text
d_J <= s sum_(j=0..J-1) a_(q,j)^2.                       (13)
```

For fixed `J`, this is `O_J(2^m m^(2J-5/2))`, whereas `r` is of order
`3^m/m^(3/2)`. The kernel therefore still has almost the full dimension.
For `m>=6J`, a product of `2J` disjoint triangle invariants and squared
matching edges on the remaining rows lies in `K_J` and is nonzero at
distinct directions. Every triangle has a monochromatic edge; at a
balanced cut the internal and external edge counts agree, so there are
at least `J` internal edges.

There is also an explicit lattice estimate after imposing all these jets
and the actual numerical zero equation. Assume `m>=6J`, so the preceding
nonzero triangle product proves that the numerical row is independent of
the jet constraints. Select `d_J` independent jet
rows in the invariant basis. Their entries have absolute value at most
`2^m`; the numerical row has maximum `Z=||f||_infinity`. Their integer
kernel has rank `n_J=r-d_J-1` and covolume at most

```text
r^((d_J+1)/2) 2^(m d_J) Z.                              (14)
```

This follows from the row-norm determinant bound, dividing by the gcd
of maximal minors if necessary. Minkowski's product estimate applied
to (14) gives the same explicit count of many short independent
relations, now already in `K_J`. Constants involving dimension must
be retained when `m` grows; for fixed `m` they are independent of `w`.

No implication from an actual small zero relation to membership in
`K_J` for `J>=2` has been proved. Even granting a fixed-depth jet
restriction would leave the large kernel and the finite lattice
estimates (13)--(14). A uniform argument needs additional arithmetic
compatibility, not merely the existence of many small exact relations.

## Independent verification

The root agent independently checked the complementary modulus, the
Gram rank argument, both core minima, the successive-minima index and
the CRT construction. The exact checker tests the actual tableau
balanced matrices through eight rows and constructs unimodular kernel
coordinates for the six-row matrix before checking the CRT identities.
Degree variations are kept separately in
[invariant_degree_variation.md](invariant_degree_variation.md).
