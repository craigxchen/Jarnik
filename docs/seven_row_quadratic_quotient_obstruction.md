# Seven rows: quadratic quotients and forgetful symmetries

Two independent statements distinguish seven rows from the balanced
elliptic six-row construction. First, a connected full Boolean fiber
product of four rational quadratic maps has genus at least five.
Second, a balanced degree-one elliptic seven-curve has at most three
degree-two single-label forgetful maps. Neither statement proves that
balanced elliptic seven-curves do not exist, or a uniform lattice-arc
bound.

Throughout, curves and function fields are in characteristic zero. For
geometric genus and automorphism arguments we pass to an algebraic
closure. A moduli curve is assumed to meet the open moduli space, and
its source is smooth, projective, and birational onto its image.

## 1. The exact genus of independent quadratic covers

In this section take `k` algebraically closed; the independence below
is geometric independence, not independence contributed by nonsquare
constants over a smaller field. Let

```text
f_i:P^1_(x_i) -> P^1_T,       deg(f_i)=2,       1<=i<=m.
```

Let `B_i` be the two geometric branch values of `f_i`. Write the
corresponding quadratic extension as `k(T)(sqrt(Delta_i))`. Suppose the
`m` squareclasses `Delta_i` are independent in
`k(T)^*/k(T)^{*2}`, so that their compositum `K` has degree `2^m`.
Its smooth projective curve is the connected normalization of the
fiber product. Put

```text
B = union B_i,       r=|B|.
```

Then

```text
2g(K)-2 = 2^(m-1)(r-4),
g(K) = 1 + 2^(m-2)(r-4).                              (1)
```

Indeed, every branch inertia group has order two. Locally at a branch
value, all ramified quadratic classes become the same class of a
uniformizer: the residue field is algebraically closed and units have
square roots. Thus each branch value contributes `2^(m-1)` to the
ramification divisor, and Riemann--Hurwitz gives (1).

There is also an exact combinatorial description. Make a multigraph
whose vertices are `B` and whose `i`th edge joins the two members of
`B_i`. The divisor-parity vector of `Delta_i` is its edge incidence
vector over `F_2`. On a projective line over an algebraically closed
field, an even divisor of a rational function is a square divisor and
the constant factor is also a square. Consequently squareclass
independence is exactly independence of these incidence vectors.
An edge incidence matrix has independent columns exactly when the
multigraph is a forest; a repeated edge is already a dependence.

If this forest has `c` nonempty connected components, then

```text
r=m+c,
g(K)=1+2^(m-2)(m+c-4),
g(K)>=1+2^(m-2)(m-3).                                (2)
```

The elementary alternative proof of `r>=m+1` is that all branch
parity vectors lie in the dimension-`r-1` even-weight subspace of
`F_2^r`. The graph gives the sharper equality and describes its case
of equality: the branch graph is a tree.

For three independent maps, genus one requires exactly four branch
values, hence a three-edge tree. The mixed path in
[the disjoint-support six-curve](disjoint_mixed_elliptic_six_curve.md)
has precisely this form. For four independent maps, `r>=5` and
`g>=5`. For five independent maps, `g>=17`. These lower bounds are
sharp for general systems of rational quadratic covers with prescribed
branch pairs; that does not assert sharpness inside a specified pencil
or a balanced moduli construction.

Independence is essential. Four maps using only four branch values
have geometric squareclass rank at most three. Their total fiber
product still has sixteen inverse-root choices over a general target,
but is disconnected: one connected component realizes at most eight
distinct coordinate tuples. In particular, displaying sixteen choices
on the union of the components does not produce a degree-sixteen
elliptic source. It is enough to require that one source realizes all
sixteen distinct coordinate tuples at one unramified fiber, since its
generic degree is then at least sixteen.

## 2. Balance itself forces independence if a quadratic quotient exists

Let

```text
f:C -> Mbar_(0,n),       n=m+3,
deg(f^*D)=1 for every stable boundary divisor D.        (3)
```

Normalize three labels as `infinity,0,1`, and denote the remaining
coordinate functions by `x_1,...,x_m`. Every `x_i` is a four-point
forgetful map. The pullback of a boundary point of `Mbar_(0,4)` is the
sum of `2^(n-4)=2^(m-1)` boundary divisors: each forgotten label may
be assigned to either side of the retained `2|2` partition. All these
partitions remain stable and are distinct. Therefore

```text
[k(C):k(x_i)]=2^(m-1),       k(C)=k(x_1,...,x_m).        (4)
```

Suppose, as an additional assumption, that there are rational maps of
degree two with

```text
T=f_1(x_1)=...=f_m(x_m).                               (5)
```

The field tower through any coordinate now gives

```text
[k(C):k(T)]=2^(m-1)*2=2^m.                            (6)
```

By (4), this field is the compositum of the `m` quadratic coordinate
extensions. Thus its squareclasses are independent automatically, and
(2) applies. No separate Boolean-fiber hypothesis is needed once
balance and (5) are known.

**Corollary.** A balanced degree-one genus-one curve in `Mbar_(0,n)`
with `n>=7` cannot admit a common quadratic quotient (5). In
particular, extending the six-row construction by a fourth moving
coordinate using another member of the same quadratic pencil cannot
give a balanced degree-one elliptic seven-curve.

The corollary only needs (5), not the stronger assertion that the maps
belong to a linear pencil. Conversely, neither the degrees (4) nor
equal boundary degrees have been shown to imply (5).

## 3. An obstruction from forgetful deck transformations

There is a restriction that does not assume any common quotient.
Continue with (3), now with `C` of genus one and `n>=7`.
For a set `S` of forgotten labels, let `L_S` be the function field of
the normalized image of the corresponding forgetful map, and put
`r_S=[k(C):L_S]`. If at least four labels remain, the map is
nonconstant and

```text
r_S divides 2^|S|.                                    (7)
```

To prove this, choose any boundary divisor downstairs. Its pullback
is the sum of `2^|S|` distinct stable divisors upstairs, all of degree
one on `C`. On the normalization of the image its pullback has
positive integral degree; multiplication by `r_S` gives `2^|S|`.

For a single label `j`, this gives `r_j in {1,2}`. When `r_j=2`, let
`tau_j` be the nontrivial deck involution of `k(C)/L_j`.
These involutions have the following properties.

1. They are distinct. For two chosen labels, normalize three other
   labels as anchors. The involution dropping either label fixes all
   the other normalized coordinates. If both involutions agreed,
   their common value would fix every coordinate, contradicting (4).
2. They commute pairwise. The group generated by `tau_i,tau_j` fixes
   `L_{i,j}`, so it has order at most `r_{i,j}<=4` by (7). Two
   distinct involutions generate at least four elements, so their
   group is exactly the Klein four group. In fact `r_{i,j}=4` and
   `L_{i,j}=k(C)^{<tau_i,tau_j>}`.
3. Any collection of at most `n-3` such involutions is linearly
   independent over `F_2`. Choose three anchors outside its labels.
   If a nonempty product were the identity, solve the relation for
   one `tau_i`. Every other factor fixes `x_i`, so `tau_i` would fix
   `x_i` as well. It already fixes all other coordinates, which is
   impossible.

An elementary abelian two-subgroup of the automorphism group of a
genus-one curve has rank at most three. To see this, choose an origin.
Its translation subgroup lies in `C[2]` and has rank at most two.
Its image among origin-preserving automorphisms has rank at most one:
in characteristic zero that automorphism group is cyclic.

If four labels had `r_j=2`, there would be three available anchors
outside those four labels because `n>=7`. The preceding properties
would give an elementary abelian subgroup of rank four, a
contradiction. We have proved:

**Theorem.** On a balanced degree-one genus-one curve in `Mbar_(0,n)`
with `n>=7`, at most three single-label forgetful maps have degree
two onto the normalization of their image. At least `n-3` of these
maps are birational. For seven labels, at least four six-point images
therefore have genus-one normalization and common boundary degree two.

For a degree-two forgetful map, the image has common boundary degree
one. Its normalization has genus one if its involution is a
translation, and genus zero if it is a reflection. The theorem does
not exclude either kind. The six-row result in
[the forgetful obstruction note](six_row_elliptic_forgetful_obstruction.md)
is stronger at six labels, but uses additional geometry: four chosen
labels there leave only two anchors, so the direct independence
argument for four generators is unavailable.

## 4. Why the six-row quotient argument does not automatically lift

At six labels, three degree-two forgetful involutions imply a common
quadratic quotient. After normalizing the other three labels, each
coordinate has degree four. Its field is therefore exactly the fixed
field of the other two involutions, and the rank-three group gives a
degree-two map from each coordinate line to its common rational
quotient. This is the mechanism used in the six-row note.

At seven labels, four such involutions are impossible by the theorem.
Even three do not provide the same field identities. To make the
remaining gap explicit, normalize as

```text
(infinity,0,1,a,b,c,d)
```

and suppose `tau_a,tau_b,tau_c` exist. They generate a rank-three
group `G`. Since `d` is fixed by all three and has degree eight,

```text
k(d)=k(C)^G.                                          (8)
```

For example, double forgetting gives

```text
k(a,d)=k(C)^{<tau_b,tau_c>},
[k(C):k(a,d)]=4,
[k(a,d):k(a)]=[k(a,d):k(d)]=2.                        (9)
```

Thus `a` and `d` are the two degree-two functions on an intermediate
curve. They do not satisfy `d=f_a(a)` with a rational quadratic map:
indeed (9) says `d` is not in `k(a)`. The intermediate curve has genus
one if the pair consists of translations, and genus zero otherwise.
The same statements hold for `b,d` and `c,d`.

Equations (8)--(9) provide correspondences of bidegree `(2,2)`, not
the rational quadratic maps of (5). Establishing or excluding balanced
seven-curves through these correspondences would require additional
boundary-incidence analysis. Curves with fewer than three degree-two
forgetful maps require still other arguments. No classification of all
balanced degree-one elliptic seven-curves is obtained here.

Finally, if a chosen geometrically connected four-map full Boolean
construction is defined over `Q`, its fixed smooth source has genus at
least five and hence finitely
many rational points by Faltings. This is finiteness for that fixed
curve. The curve, its coefficients, and its height constants may vary
with the proposed configuration; no uniform bound over those varying
families, and no global uniform `sqrt(R)` arc bound, follows.

The [independent audit](seven_row_quadratic_quotient_audit.md) checks
the geometric independence, boundary counts, and fixed-field steps.
