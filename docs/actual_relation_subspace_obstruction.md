# Actual short relation subspaces with deficient smaller-cut images

The condition that all relations vanish at the **same distinct rational
configuration** does not turn the first-smaller-cut rank bounds into a
global independence theorem. A fixed triangle factor gives an explicit
large subspace of the actual numerical relation ideal, with every such
restriction image deficient. The construction also persists at the
coefficient-height threshold supplied by the doubled restriction-minor
theorem: already on forty rows, successive minima force many sufficiently
short relations inside this subspace, conditional on the actual uniform
profile hypotheses.

This is an obstruction to using the first layer alone. It is not an
unbounded endpoint construction and does not preclude descent to an
irreducible relation with smaller support. The latter must retain the
amplified block weights after deleting rows.

The standard-library certificate is
[check_actual_relation_subspace_obstruction.py](check_actual_relation_subspace_obstruction.py).
The underlying arithmetic restriction estimates are in
[smaller_cut_relation_rank.md](smaller_cut_relation_rank.md) and
[complement_reflection_invariant_relations.md](complement_reflection_invariant_relations.md).

## 1. A subspace of the actual numerical ideal

Let `m=2q>=8`, let `V_m` be the simultaneous `SL_2` invariant space of
degree two in each row, and fix any distinct rational directions
`p_1,...,p_m`. Put

```text
T=Delta_12 Delta_23 Delta_31,
J={4,...,m},             n=m-3,
U=T * ker(ev_(p_j,j in J): V_n -> Q).                    (1)
```

Here all factors and vector spaces are over `Q`; choose a fixed integral
invariant basis when measuring coefficient heights. Since the three
directions are distinct, `T(p)!=0`. Multiplication by `T` is injective,
and evaluation on `V_n` is nonzero: any Hamilton-cycle bracket monomial
is nonzero on distinct directions. Consequently

```text
dim U=r_n-1,
r_n=[x^n](1+x+x^2)^n-[x^(n-1)](1+x+x^2)^n.              (2)
```

Every member of `U` is an exact numerical zero at the given full
configuration. This is not merely a subspace with the right formal
cut incidences.

Let `K` be the balanced-cut kernel in `V_m`. Then

```text
U subset K,
dim rho_S(U) <= q-2  whenever |S|=q-1.                  (3)
```

The entire first-smaller restriction target has dimension

```text
t=(q+1)(q-2)/2 > q-2.                                  (4)
```

To prove (3), collapse the rows of `S` to `e_1`, and write each other
row as `(z_j,1)`. If two triangle rows are collapsed, the triangle
vanishes. If at most one triangle row is collapsed at a balanced cut,
strictly more than half of the remaining `n` rows are collapsed, so
every residual invariant vanishes by weight balance. This proves
membership in `K`.

For a first-smaller cut, the only potentially surviving case has
exactly one triangle row inside. It then has

```text
a=q-2=(n-1)/2
```

residual rows inside. A degree-two-each invariant on these `n` rows
restricts to a homogeneous polynomial of degree `n-2a=1` in the
`n-a=q-1` outside variables. Invariance under simultaneous translation
makes the sum of its linear coefficients zero. Thus its restriction
space has dimension at most `q-2`. The restricted triangle is one
fixed outside difference, so multiplying by it cannot increase this
dimension. If no triangle row is inside, too many residual rows are
collapsed and the restriction is zero. This proves (3) for every cut,
without a generic-position hypothesis.

In particular, all maximal minors of size `t` vanish identically on
every tuple from `U`. Every modulus dividing those minors, including
the doubled complementary modulus, is therefore respected.

## 2. An eight-row exact example

Take the actual primitive Gaussian rows

```text
P_j=4j+i,       j=1,...,8.                              (5)
```

They have odd norms and distinct rational directions. The usual Gaussian
lcm construction realizes their ratios on an integer circle; no endpoint
or full-profile assertion is made for this fixed example.

For each of the six triples

```text
456, 457, 458, 467, 468, 478,
```

let `B_j` be the triangle product on that triple times the square of
the bracket on its complementary pair in `{4,...,8}`. These six
polynomials form a basis of `V_5`: direct exact expansion gives rank
six, agreeing with (2). All six numerical values `v_j=B_j(P)` are
nonzero. The five polynomials

```text
T * (v_1 B_j-v_j B_1),       j=2,...,6,                 (6)
```

are independent exact numerical zeros. Divide their polynomial
coefficients by their ordinary gcds. Their coefficient sums are

```text
1776, 912, 4032, 4752, 984.                             (7)
```

The checker verifies directly:

* all five numerical values are exactly zero and their polynomial rank
  is five;
* all seventy balanced restrictions are zero;
* among all fifty-six first-smaller restrictions, twenty-six have image
  rank zero and thirty have image rank two.

Thus even at eight rows there is a `t=5` dimensional subspace of an
actual evaluation hyperplane with all first-smaller images proper.
This strengthens the earlier purely formal common-factor example.

## 3. The residual primitive height, including odd row counts

Now assume the actual full `m`-row uniform central-profile hypotheses,
with inherited block weight `w`, all individual block logarithms
`w+o(w)`, and corrections and primitive pair-residue heights `o(w)`.
Throughout this section `m` is fixed and `w` tends to infinity.

Delete the first three rows but retain the inherited normalization.
For each nonempty `A subset J`, combine the eight old blocks whose
intersection with `J` is `A`. The resulting residual block has norm logarithm

```text
w'=8w+o(w).                                             (8)
```

The full-support block is retained in the raw arithmetic and removed
only through the exact common-content calculation. The empty-intersection
product never divides a retained row and needs no weight assertion; it
has only seven factors if the original convention omits `H_empty`.
We do not renormalize the inherited radius or silently discard the
full-support common factor. In (9)--(11), the empty-set exponent is zero,
so that block's presence and weight are immaterial.

For any `n>=3`, including odd `n`, define

```text
g_0(a)=max(0,2a-n),
A_n=n 2^(n-2)-sum_(a=0)^n binomial(n,a)g_0(a).           (9)
```

The primitive numerical evaluation vector on `V_n` has height

```text
h_(V_n)=A_n w'+o(w)=8 A_n w+o(w).                       (10)
```

Here height is the logarithm of the largest absolute coordinate after
dividing an integral evaluation vector by its gcd. To check the odd
case, use the exact residual factorization

```text
|Delta_ij|=b_ij product_(A containing i,j) n_A,
log b_ij=o(w),       log n_A=w'+o(w).
```

Every degree-two-each bracket graph has at least `g_0(a)` internal
edges on an `a`-vertex subset. A Hamilton cycle attains this minimum
for each specified subset: when `a<=n/2`, place all its vertices
nonadjacently around the cycle; the complementary case follows from
degree balance. This works for odd `n` as well as even `n`.

Let `D=product_A n_A^g_0(|A|)`, let `G` be the gcd of the values of a
fixed integral basis of `V_n`, and put `B=product_(i<j)b_ij`. Integral
bracket straightening and the preceding attaining monomials give

```text
D | G,                 G | D B^2.                     (11)
```

The lower bound applies to every monomial. For the upper bound, at
each prime choose an attaining cycle. Its correction valuation is
bounded by twice the sum of all edge correction valuations. The same
argument applies to primes outside the core. No cancellation assertion
for a sum of monomials is needed. Raw nonzero graph values have
logarithm `n 2^(n-2)w'+o(w)`, so (11) gives (10).

For odd `n`, the binomial sum simplifies to

```text
A_n=n 2^(n-2)-n binomial(n-1,(n-1)/2).                  (12)
```

The dimension formula (2) and the primitive-height convention agree
with [invariant_relation_lattice_threshold.md](invariant_relation_lattice_threshold.md).
Equation (11) supplies the required odd-row extension rather than
assuming an even-row formula unchanged.

## 4. Short vectors inside the deficient-image subspace

The exact residual numerical relation lattice has rank `r_n-1` and
Euclidean covolume equal to the Euclidean norm of its primitive
evaluation vector. For its successive minima, Minkowski's second
theorem and the fixed lower bound on nonzero integral vectors imply

```text
log lambda_j <= 8 A_n w/(r_n-j)+o(w),
                       1<=j<=r_n-1.                   (13)
```

The denominator is `r_n-j`: these are precisely the minima from
indices `j` through `r_n-1`. Changes between fixed integral-coordinate
and polynomial-coefficient norms cost `O_m(1)` in logarithmic height.

Multiplication by the triangle changes coefficient sums by at most a
factor eight. More generally, its inverse on its image is a fixed
rational linear map, so comparison in either direction costs only
`O_m(1)`. Thus (13) produces independent **actual** numerical relations
inside `U` with the same leading coefficient-height bounds.

The doubled smaller-cut theorem forces a proper image for any `t`
relations when their individual coefficient heights are less than
`2w/t-o(w)`; its fixed finite correction losses are absorbed in this
notation. The common-triangle subspace itself supplies such a tuple
whenever

```text
8 t A_(m-3)/(r_(m-3)-t) < 2.                            (14)
```

The exact arithmetic comparison is:

| m | n=m-3 | r_n | A_n | t | Ratio in (14) |
|---:|---:|---:|---:|---:|---:|
| 30 | 27 | 18,985,057,351 | 625,153,464 | 104 | 27.396688 |
| 36 | 33 | 10,331,450,919,456 | 51,031,307,514 | 152 | 6.006327 |
| 38 | 35 | 85,317,692,667,643 | 218,971,493,020 | 170 | 3.490498 |
| 40 | 37 | 707,854,577,312,178 | 935,530,313,516 | 189 | 1.998323 |
| 42 | 39 | 5,897,493,615,536,452 | 3,981,653,897,208 | 209 | 1.128840 |
| 44 | 41 | 49,320,944,483,427,000 | 16,888,280,687,788 | 230 | 0.630045 |

Forty is the first even row count satisfying this comparison. The
strict margin at forty is small but positive; forty-two gives a more
comfortable asymptotic margin. At forty rows the local target has
dimension `189`, while every image of `U` has dimension at most `18`.

In fact (13) gives more than `t` such vectors. For any fixed integer

```text
j < r_n-4t A_n,                                        (15)
```

all first `j` minima lie strictly below `2w/t` for sufficiently accurate
profiles and sufficiently large `w`. At forty rows the right side is
`593,660,294,082`. This is a statement conditional on the actual profile
being present; it does not construct such a profile.

## 5. What a further argument must use

This rules out a proposed algebraic implication of the form

```text
U is contained in the actual evaluation hyperplane,
dim U>=t, and all first-smaller images are proper
    => contradiction.
```

It also shows why simply imposing short coefficient heights does not
repair that implication in large fixed dimension. The same actual
relation lattice necessarily contains a large common-factor subspace
at the relevant scale, if the full-profile configuration exists.

The triangle is nonzero at the configuration, so each numerical zero
in (1) is really a zero on the remaining rows. Factoring is therefore
the appropriate next step. But deleting three rows amplifies the block
weight by eight, as in (8); the residual coefficient height is preserved
up to a fixed constant, rather than rescaled with that weight. A viable
ideal-descent or irreducible-support theorem must exploit this precise
height-to-block comparison, or information from deeper restriction
layers. No independence between the separate cut conditions has been
assumed anywhere in this example.
