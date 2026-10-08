# Eight-point Pfaffian content and the missing short-arc coupling

There is an exact divisibility involving all twenty-eight primitive chord
numerators. It retains prime powers and has a positive local occupancy
formula. It does not yet prove the balanced-layer bonus in
[finite_local_test_attack.md](finite_local_test_attack.md).

The limitation is concrete: a family of actual integral circle points can
have growing balanced conductor mass while this new common-content
quantity stays bounded. That family is not at endpoint scale. Thus this
route still needs a coupling to the shrinking angular width, beyond cyclic
order and the Pfaffian identity.

## 1. An exact integer identity

Use eight distinct points in the common-unit model of
[ordered_residue_growth.md](ordered_residue_growth.md), on a proper arc
of angular width less than `pi`. For each pair choose the sign of its
conjugate-primitive Gaussian numerator so that

```text
h_ij=x_ij+i t_ij,       x_ij>0,
t_ij/x_ij=tan((theta_j-theta_i)/2),       i<j.
```

Here `t_ij` is signed; its absolute value is the positive integer residue
of that note. In particular `gcd(x_ij,t_ij)=1`. The same algebra below
allows nonzero signed real parts when the points are not in a proper arc.
Let `E` be the twenty-eight pairs, and let `M` range over the 105 perfect
matchings of eight labels, with its Pfaffian sign `epsilon_M`. Then

```text
sum_M epsilon_M product_(ij in M)t_ij
                     product_(ij not in M)x_ij
 = product_(ij in E)t_ij.                                (1)
```

This is Schur's Pfaffian identity after putting `q_i=exp(i theta_i)`.
For completeness, multiply
`Pf((q_j-q_i)/(q_j+q_i))` by `product_(i<j)(q_i+q_j)`.
The result is an alternating polynomial of degree 28, hence is a
constant times the Vandermonde polynomial. With `q_i=t^i` and
`t -> infinity`, every upper-triangular matrix entry tends to one;
the Pfaffian of this matrix is one by its first-row recurrence.
The Vandermonde quotient has the same limit. The constant is therefore
one. Since `(q_j-q_i)/(q_j+q_i)=i t_ij/x_ij`, both the Pfaffian and
the product acquire the same factor (`i^4=i^28=1`). Clearing the
real-part denominators gives (1).

This standard identity is also recorded in the primary research paper
[Okada, *An elliptic generalization of Schur's Pfaffian identity*](https://arxiv.org/abs/math/0412038).
The proof above supplies the specialization and signs used here.

Define positive integers

```text
X=product_E |x_ij|,             T=product_E |t_ij|,
L_X=lcm_M(product_(ij in M)|x_ij|),
B_X=X/L_X.
```

Every matching product divides `X`, so `B_X` is an integer. Moreover,
`B_X` divides each complementary product `product_(ij not in M)x_ij`
in (1). Consequently

```text
B_X divides T.                                           (2)
```

At every prime this says exactly

```text
sum_E v_p(x_ij)-max_M sum_(ij in M)v_p(x_ij)
 <= sum_E v_p(t_ij).                                    (3)
```

Thus substantial shared real-part divisibility outside a single
matching forces primitive imaginary-residue growth. There is no
coefficient denominator or loss from replacing an exponent by its
radical.

## 2. A positive local identity at every prime

Assume the input circle is primitive. At a split conductor prime `p`,
first partition the rows by their allocation `a_i(p)`. In the group
with allocation `a`, write

```text
z_i=pi^a bar(pi)^(e-a) u_i.
```

The `u_i` are units modulo `p`, with the same norm `N/p^e`.
At a prime outside the conductor use one group and put `u_i=z_i`.
This covers all odd primes for a primitive circle, since an inert
conductor factor would divide every point.

If two rows have different allocations, their primitive pair numerator
has positive valuation at exactly one Gaussian prime above `p`.
Its real and imaginary parts are therefore both `p`-units. Within one
allocation group, the exact identities
`u_i+u_j=2 epsilon g_ij x_ij` and
`u_i-u_j=2i epsilon g_ij t_ij`, with `g_ij` a unit modulo `p`, give

```text
v_p(x_ij)>=l  iff  u_i == -u_j modulo p^l,
v_p(t_ij)>=l  iff  u_i ==  u_j modulo p^l.                 (4)
```

For each depth `l` and allocation group, let `n_c` count its rows
in the residue class `c` of `Z[i]/p^l`. Opposite classes `c,-c` are
distinct because the norm is a unit and `p` is odd. Empty classes
are included with count zero. Then the exact identity is

```text
v_p(T)-v_p(B_X)
 =sum_(l, allocation group, unordered {c,-c})
       binom(|n_c-n_-c|,2).                              (5)
```

To prove it, the graph of edges with `v_p(x_ij)>=l` is a disjoint
union of complete bipartite graphs on opposite classes. Its maximum
matching has `sum_{c,-c} min(n_c,n_-c)` edges. These maximum sizes
can be achieved simultaneously for every depth: start at the largest
nonempty depth, match each opposite pair as fully as possible, and
continue downwards. A previously chosen edge removes one vertex from
each side of every ancestor component. The remaining vertices can
therefore be matched until the ancestor reaches its maximum size.
At depth zero extend the resulting matching arbitrarily to a perfect
matching of all eight rows; the additional edges have zero weight.
It follows that

```text
max_M sum_M v_p(x_ij)
 =sum_(l,groups,{c,-c}) min(n_c,n_-c).
```

The total real-part valuation at a depth is `sum n_c n_-c`; the
total residue valuation is `sum [binom(n_c,2)+binom(n_-c,2)]`.
Subtract and use

```text
binom(n,2)+binom(m,2)-nm+min(n,m)=binom(|n-m|,2).
```

This proves (5) at odd primes with arbitrary prime powers.

At the prime two use a single group and, for valuation threshold
`l>=1`, use residue classes of `z_i` modulo `2^(l+1)` instead of
modulo `2^l`. A primitive common-unit circle has odd norm, so each
`g_ij` has odd norm. The exact factors `2 epsilon` and `2i epsilon`
in the sum and difference identities show that (4) holds with this
one-step shift. Opposite classes are distinct: an odd-norm Gaussian
integer cannot satisfy `2z=0 modulo4`. The nested matching proof is
unchanged. Thus (5) holds at two as well, with the stated shifted
classes. The same argument works for any even number of rows.

## 3. The exact sufficient estimate still missing

For a common-unit eight-point endpoint cluster, the exact chord formula
and the identity `sum_(i<j)d_ij=16W-D` imply

```text
D+2log T <= 2W+56log(C/2).                              (6)
```

Indeed `|z_i-z_j|=2|t_ij| exp((W-d_ij)/2)` and each chord is at
most `C exp(W/4)`. Thus a bound of the form

```text
2log B_X >= W_4-B                                       (7)
```

on the actual endpoint cluster would imply the proposed intrinsic
balanced-layer inequality, with an adjusted uniform additive constant.
No estimate (7) has been proved. Formula (5) alone gives no reason why
balanced conductor layers should force large real-part overlap.

## 4. A balanced prime can contribute no local real or imaginary content

There is an exact family illustrating that local gap. Put

```text
pi=3+2i,        p=13,
(X_1,...,X_8)=(26,2,4,10,52,28,30,36),
H_i=X_i+i,     S={1,2,3,4},
h_i=pi^e H_i if i in S, and h_i=H_i otherwise,
G=pi^e product_i H_i,       z_i=bar(G) h_i/bar(h_i).
```

Every `z_i` is a Gaussian integer, of common norm
`13^e product_i N(H_i)`. The rows are distinct. Each `H_i` is
conjugate-primitive and a unit at both primes above 13. The Gaussian
allocation at 13 is exactly `e` on `S` and zero outside `S`, even
after removing any common factor. Thus this prime supplies balanced
layer mass `e log13`.

Within each group, the four residues of `X_i` modulo 13 are
`0,2,4,10`. Their pairwise differences and every `X_i X_j+1` are
13-units, and each `X_i^2+1` is a 13-unit. Hence both primitive
pair coordinates are 13-units within either group. For a pair across
the groups, the primitive numerator has a nonzero allocation at only
one orientation of 13, so both coordinates are again 13-units. Thus

```text
13 does not divide X T B_X,       for every e>=1.         (8)
```

This is not an endpoint family: two rows in the same group have a
fixed positive angular separation, while the radius tends to infinity.
It demonstrates that the required growth must sometimes occur at
other primes, through an additional global arithmetic argument.

## 5. In this family even the global B_X stays bounded

Write `pi^e=A_e+i B_e`, with `gcd(A_e,B_e)=1`. The real parts for
pairs within a group are fixed integers. For a cross pair, after its
fixed Gaussian common norm is removed, its real part is a primitive
integer linear form

```text
L_ij(A_e,B_e)=c_ij A_e+d_ij B_e.                          (9)
```

The common norm is fixed because no `H_i` has a factor above 13.
Primitive coefficient pairs which are proportional differ only by
sign. Edges with the same form up to sign constitute a matching:
if two such cross edges shared a vertex, two distinct fixed `H_i`
would have the same projective complex direction. This is impossible
for the displayed `X_i+i`.

Let `L_1,...,L_q` represent the distinct forms up to sign and write
`Delta_ab=det(coeff(L_a),coeff(L_b))!=0`. For every prime and every
primitive pair `(A_e,B_e)`,

```text
min(v_p(L_a(A_e,B_e)),v_p(L_b(A_e,B_e)))
 <= v_p(|Delta_ab|).                                    (10)
```

This follows by the two linear combinations giving `Delta_ab A_e`
and `Delta_ab B_e`. Let `X_fixed` be the product of the absolute
within-group real parts. At a prime choose a form with largest
valuation and place all of its edges in one matching; this is
possible by the preceding disjointness observation. Extend it to
a perfect matching. Each of the at most sixteen cross edges left
outside has valuation bounded by its minimum with that chosen form.
Consequently (10) gives the uniform bound

```text
B_X <= X_fixed product_(a<b)|Delta_ab|^16.               (11)
```

The right side is independent of `e`. Thus neither integer circle
geometry alone nor the exact Pfaffian common-content test supplies
(7) from balanced conductor mass. The failure is not a counterexample
under the endpoint hypothesis; the fixed within-group span excludes
that hypothesis. A successful application must exploit the shrinking
span in a way absent from (1)--(5).

## 6. Nearly equal pair-sum products: the exact scalar

Let

```text
P_M=product_(ij in M) (z_i+z_j)/(2 epsilon)
   =product_(ij in M) g_ij x_ij.
```

These are Gaussian integers. All have one complex phase, since a
matching uses each point angle exactly once. Let `Gamma_0` be the
Gaussian gcd of the core products `product_M g_ij`. Sorting the eight
allocations at each prime gives

```text
log|Gamma_0|=(1/2)sum_layers |r-4| log p.                 (12)
```

Indeed, matching the four lowest allocations to the four highest
simultaneously minimizes both Gaussian valuations: the first is the
sum of the four lowest allocations, and the conjugate is the sum of
four copies of `e` minus the four highest allocations. At a threshold
layer their sum is `|r-4|`. The two orientation differences are
independent of the matching, so every quotient of a core matching
product by `Gamma_0` is an integer up to sign.

The actual gcd of the `P_M` can therefore be written

```text
Gamma_+=Gamma_0 delta,       delta in Z_(>0).
```

There is a useful exact refinement:

```text
delta divides B_X, and hence delta divides T.           (13)
```

At an ordinary prime choose a perfect matching `F` maximizing the
sum of its real-part valuations. Among core-minimizing matchings,
there is one avoiding every edge of `F`: partition the allocations
into four low and four high rows, and use a perfect matching of
`K_(4,4)` after removing the edges of `F`. At most one edge at each
vertex has been removed, so this graph has a perfect matching (for
example, by Hall's condition). Its core excess is zero and its real-part
valuation is at most `sum_E v_p(x)-sum_F v_p(x)=v_p(B_X)`.
Thus the minimum defining `v_p(delta)` has the required upper bound.
At a nonconductor prime the core excess is zero for every matching;
the same argument applies. This also handles the prime two in the
primitive common-unit normalization.

Choose the unit of `Gamma_+` so that `c_M=P_M/Gamma_+` are positive
integers. Their gcd is one. For an angular width `Theta<=1`, set

```text
U=R^4/(|Gamma_0| delta).
```

The exact cosine formula and `1-cos x<=x^2/2` give

```text
(1-Theta^2/2)U <= c_M <= U.                              (14)
```

They cannot all be equal. For four distinct labels,

```text
(z_1+z_2)(z_3+z_4)-(z_1+z_3)(z_2+z_4)
 =(z_1-z_4)(z_3-z_2),                                  (15)
```

and multiplication by two remaining nonzero pair sums gives a nonzero
difference of matching products. Consequently `U Theta^2/2>=1`.
This near-equality argument supplies the upper bound

```text
delta <= R^4 Theta^2/(2|Gamma_0|).                       (16)
```

It does not give a lower bound for `B_X`. In the critical fair
eight-row profile, `E|r-4|=35/32`, so

```text
log|Gamma_0|=(35/64)W+o(W),
log max_M c_M=(93/64)W-log delta+o(W).                   (17)
```

The near-equality lower height is only `W/2+o(W)`. Its remaining
gap is `61W/64` before subtracting `log delta`. In the small-residue
regime (13) already gives `log delta=o(W)`, much stronger than the
upper bound (16). Thus scalar normalization reverses the tempting
inference that nearly equal integers should force large `B_X`.

## 7. Higher fixed linear cancellation in this space

A fixed integer linear combination of the 105 unnormalized products
is a polynomial of total degree four in the eight `z_i`. If it is
not the zero polynomial, its order `v` after translating all
`z_i=z_0+u_i`, for nonzero `z_0`, satisfies `v<=4`. Homogeneity
shows that the translated coefficients of degree `j` are integer
coefficients times `z_0^(4-j)`, so this statement does not conceal
a dependence on the common phase or radius.

For coefficients with sum of absolute values `C_A`, the corresponding
integer linear combination of the primitive `c_M` has the bound

```text
|sum_M a_M c_M|
 <=16 C_A R^4 Theta^v/(|Gamma_0| delta).                 (18)
```

The fixed factor 16 bounds the binomial expansion after the common
`(2 epsilon)^4` has been removed. On a fair endpoint profile, even
the maximal generic order four leaves the upper exponent

```text
(29/64)W-log delta+log C_A+o(W).                         (19)
```

Thus, with `log C_A+log T=o(W)`, this particular finite linear
space does not give a nonzero-integer contradiction merely from
its Taylor order and common content.

The order-four cancellation has a familiar exact realization:
multiply two identities (15) on disjoint groups of four rows. The
result is a four-term linear combination of plus-matching products
equal to a product of four ordinary chords. It recovers the original
matching-chord coordinates and their `29W/64` primitive height.
This explains the remaining exponent without identifying a new
balanced-layer bonus.

The scope of (18)--(19) is the universal Taylor estimate in this
105-coordinate space. A zero polynomial gives only an identity;
special cancellation at the actual tuple can make a nonzero
polynomial's value smaller or zero. Such special arithmetic effects,
and expressions outside this finite space, are not excluded.

## 8. The exact 35-dimensional polynomial flag

The finite space in Section 7 has an exact connection to the balanced
incidence rank in
[invariant_relation_lattice_threshold.md](invariant_relation_lattice_threshold.md).
The coefficient of a squarefree monomial `product_(i in S)z_i`,
`|S|=4`, in a plus-matching product is one precisely when every
matching edge crosses `S`. Complementary coefficients agree. The
35 incidence rows have Gram entries `a!(4-a)!`, where
`a=|S intersect T|`; the positive-definiteness proof in that note
therefore proves that the plus-product polynomial space has dimension
35.

At a nonzero diagonal point its Taylor filtration has graded dimensions

```text
order 0: 1,       order 2: 20,       order 4: 14.         (20)
```

For a matching product translated by `z_i=z_0+u_i`, the constant
and linear terms are independent of the matching. The coefficient
of `u_i u_j` is `4z_0^2` if `ij` is not a matching edge and zero
otherwise. The span of these coefficient rows has rank 21. Indeed
the matching edge-incidence matrix has rank 21: its annihilator
consists exactly of weights `w_ij=a_i+a_j` with `sum_i a_i=0`.
Two-edge swaps give
`w_ab+w_cd=w_ac+w_bd=w_ad+w_bc`, which force these additive
potentials; their sum condition follows by summing on a matching.
This annihilator has dimension seven among 28 edge coordinates.
Replacing incidence by its complement preserves the rank because
the constant row is the normalized sum of either set of rows.

For a triple the cubic coefficient is
`2z_0(1-number of matching edges inside the triple)`, so it adds
no rank beyond the constant and quadratic terms. The full polynomial
rank is 35, giving the remaining 14 independent directions at order
four and proving (20).

This is a polynomial-space splitting, not a short-lift theorem for
the lattice of exact numerical relations. The evaluation functional
and its large common divisors still depend on the actual input.
The rank calculation supplies neither a small nonzero quotient
coefficient nor a bounded numerical representative of that quotient.

## Verification and scope

[check_eight_point_pfaffian_content.py](check_eight_point_pfaffian_content.py)
checks the 105 signed terms on exact integer examples, the matching
LCM divisibility, the positive occupancy identity at odd primes both
outside and inside the conductor and at two with its shifted depth,
the pair-sum Gaussian gcd, and the explicit balanced-prime family through
exponent 20. These checks supplement the proofs above.

The Pfaffian numerator itself is the Vandermonde polynomial after
clearing its rational denominators. It is therefore not an independent
new projective equation. The useful extra formulation is the integral
denominator overlap and its exact positive local decomposition. The
uniform endpoint bound and the balanced-layer bonus remain unproved.

The later [Pfaffian--Hafnian audit](pfaffian_hafnian_arithmetic.md)
retains a new positive numerator in the nonlinear squared-sum identity.
It also proves that final denominator cancellation can occur to
arbitrary prime-power depth away from every old residue, on actual
nonendpoint circle tuples. Termwise denominator bounds and reduced
denominators must therefore remain distinct in this route.
