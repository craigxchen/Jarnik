# Maximal minors of smaller-cut relation matrices

At the first smaller cuts of an arbitrary even-row actual central
configuration, a fixed large ordinary integer divides every maximal
minor of the coefficient matrix of exact invariant relations.
Consequently a sufficiently small relation space has a proper
restriction image at every such cut. Applying the exact complement
reflection doubles the leading logarithmic modulus to `2w`.
The six-row case supplies the
two-relation theorem; for more rows this rank condition by itself
does not force proportionality or a common factor.

The exact algebra and an explicit common-factor example are checked by
[check_smaller_cut_relation_rank.py](check_smaller_cut_relation_rank.py).
This note asserts no uniform count.

## 1. An integral restriction basis with no conversion loss

Let `m=2q>=6`, and let `K` be the balanced-cut kernel in the invariant
space of degree two in each row. For an inside set `S` of size `q-1`,
there are `n=q+1` outside rows. Relabel them temporarily `1,...,n`.
The restriction `F_(S,Q)` of `Q in K` is a homogeneous quadratic in
their `z` variables. It is squarefree and translation invariant.
Indeed its square coefficient at `z_j^2` is the vanishing balanced
coefficient at `S union {j}`; translation invariance follows from
the unipotent subgroup of `SL_2` fixing the inside vector `e_1`.

Let `W_n` be this quadratic space. If its coefficients are `a_ij`
on `z_i z_j`, translation invariance says

```text
sum_(j!=i) a_ij=0       for every i.
```

The coefficients with both indices below `n` therefore have sum
zero, and those incident to `n` are determined by them. This gives
the following integral basis, of dimension `t=n(n-3)/2`:

```text
q_ij=(z_i-z_n)(z_j-z_n)-(z_1-z_n)(z_2-z_n),
1<=i<j<n,       (i,j)!=(1,2).                            (1)
```

The coordinates in this basis are exactly the coefficients `a_ij`
for the indicated pairs. In particular, for an integer invariant,

```text
F_(S,Q)=sum c_(ij,Q) q_ij,
c_(ij,Q) in Z,       sum |c_(ij,Q)|<=C_Q.                (2)
```

Here `C_Q` is its coefficient sum in the original binary variables.
Direct substitution into the restriction does not inflate the
coefficients, so (2) has no basis-dependent loss.

Each basis element is a sum of at most two matching products on four
distinct outside rows. For example, when `3<=i<j<n`,

```text
q_ij=(z_i-z_1)(z_j-z_n)+(z_1-z_n)(z_j-z_2).
```

The distinguished element is a single matching:

```text
q_13=(z_1-z_n)(z_3-z_2).                                (3)
```

## 2. A common norm modulus for every maximal minor

Use the same actual primitive Gaussian data as
[near_balanced_invariant_congruences.md](near_balanced_invariant_congruences.md),
with `log|K_i|<=sigma` and `0<|t_ij|<=T`. Fix `H=H_S`.
Set `z_j=Y_j/P_j` on the outside rows, and define

```text
v_ij=(product_outside P_j) q_ij(z),
kappa=product_outside K_j.                              (4)
```

The matching decompositions in Section 1 show that all `v_ij` are
Gaussian integers. Also, up to an irrelevant sign,

```text
v_13=D_1n D_32 product_(j outside {1,2,3,n}) P_j !=0.    (5)
```

The nonvanishing uses distinct actual directions and nonzero rows.
If `Q` is a numerical zero, its specialized restriction, coprimality
of the inside `Y_i` with `H`, and coprimality of the outside core
factors with `H` give

```text
H | kappa sum c_(ij,Q) v_ij.                            (6)
```

Take any `t` such numerical zero relations in `K`, and let `delta`
be the determinant of their `t` by `t` coefficient matrix in (2).
Multiplying its equations (6) by the integer adjugate gives

```text
H | kappa delta v_ij       for every basis index ij.
```

This argument does not invert the coefficient matrix modulo `H`.
Since `delta` is an ordinary integer and `H` contains one Gaussian
orientation over each of its odd split primes, we obtain

```text
M_S^* | delta,
M_S^*=N(H)/N(gcd_G(H,kappa v_13)).                       (7)
```

The modulus is fixed for the actual configuration and this inside
set; it is independent of the chosen relations. It therefore divides
every maximal minor of any rectangular restriction matrix of
numerical zeros in `K`.

For completeness the finite loss is

```text
log M_S^* >= log N(H_S)-(4n-4)sigma-2log T.               (8)
```

The factors in `kappa` cost at most `2n sigma` in logarithmic norm.
The remaining `n-4` outside row factors in (5) cost at most
`2(n-4)sigma`; their core parts are units at `H`. The exact
outside-bracket gcd bound from the preceding congruence note bounds
the two determinants' contribution by `4sigma+2log T`. Adding these
costs gives (8), and retains arbitrary prime powers.

Hadamard's inequality and (2) give

```text
|delta|<=product_(j=1..t) C_(Q_j).                       (9)
```

Thus, if a family of actual numerical zeros in `K` has
`log C_Q<=h`, and

```text
t h+(4n-4)sigma+2log T < log N(H_S),                    (10)
```

its restriction image at `S` has rank at most `t-1` over `Q`.
This is an actual arithmetic assertion. The scalar norm in (7) is
available even though, for `n>4`, the individual `v_ij` need not be
real. For six rows, `t=2`, and the joint injectivity proved in
[two_small_invariant_relations.md](two_small_invariant_relations.md)
then forces proportionality globally.

The common reflection in
[complement_reflection_invariant_relations.md](complement_reflection_invariant_relations.md)
preserves the actual zero relations and the same labelled coefficient
matrix. Its modulus is supported on the complementary block and is
coprime to the first modulus. Consequently every maximal minor is
divisible by a common ordinary integer whose logarithm is at least

```text
log N(H_S)+log N(H_(S^c))-2m(4n-3)sigma-4log T.          (10a)
```

Thus the stronger sufficient condition for a proper restriction image is

```text
t h+2m(4n-3)sigma+4log T < 2(1-eta)w.                  (10b)
```

The reflection proof charges all removed row contents and core powers.
It does not assume that reflected imaginary coordinates remain small.

## 3. A concrete shared-factor obstruction to pencil rigidity

Already on ten rows, independent invariant polynomials can have
proportional first-smaller restrictions at every cut. Let

```text
P=T_123 T_456,
Q=P Delta_78^2 Delta_(9,10)^2,
R=P Delta_79^2 Delta_(8,10)^2.                           (11)
```

They have degree two in each row and belong to the balanced kernel.
They are independent polynomials: the last four rows give two
nonproportional squared matchings. For an inside set of size four,
if either triangle contains two inside rows, both restrictions vanish.
Otherwise at most two of the first six rows are inside. At least
two of the last four are then inside. If more than two are inside,
both last-four invariant restrictions vanish by weight balance.
If exactly two are inside, each last-four restriction is a scalar,
so the two complete restrictions are scalar multiples of the same
restriction of `P`.

The checker verifies this at all `binomial(10,4)=210` cuts by exact
polynomial substitution. These are not asserted to be numerical
zero relations in an actual central configuration. At distinct
directions `P` is nonzero, so a numerical zero of either product
would descend to a numerical zero of its last-four factor. This is
why a prospective rigidity theorem must either allow and control
common factors, or use more than proportionality of the first
smaller restrictions. The example does not classify all possible
degeneracies.

## 4. Exact comparison with the number of short relations

Write

```text
r=dim V,
A_m=m 2^(m-2)-(m/2)binomial(m,m/2),
s=binomial(m,m/2)/2,
t=(m/2+1)(m/2-2)/2.
```

The actual full numerical relation lattice has rank `r-1` and
covolume of logarithmic size `A_m w+o(w)`. Minkowski's bound on
its `t`th successive minimum, in a fixed coefficient norm, is

```text
log lambda_t <= A_m w/(r-t)+o(w).                       (12)
```

The exact comparison for the first-smaller-cut determinant threshold
is therefore `t A_m/(r-t)<2`, in addition to membership in `K`.
The latter follows from the balanced `2w` bound whenever the former
holds. Some values are:

| m | r | A_m | t | t A_m/(r-t) |
|---:|---:|---:|---:|---:|
| 18 | 1730787 | 742068 | 35 | 15.006413 |
| 24 | 834086421 | 68213424 | 65 | 5.315844 |
| 30 | 439742222071 | 5726300880 | 104 | 1.354283 |
| 32 | 3602118427251 | 24742452128 | 119 | 0.817395 |
| 36 | 245613376802185 | 455122855224 | 152 | 0.281657 |

For example, at `m=30` the bounds do guarantee at least `104`
independent small numerical zeros whose images are proper at every
first-smaller cut. There is no contradiction in this statement:
the domain of each restriction map is much larger than its image.

Even a dimension count using the stronger requirement that all
first-smaller restrictions vanish does not exhaust the elementary
minima allowance. There are `N=binomial(m,m/2-1)` such cuts, so
their common kernel inside `V` has codimension at most `s+Nt`.
At coefficient exponent approaching `2/t` from below, the product
bound for the minima alone guarantees roughly
`r-1-A_m t/2` small independent numerical relations. At `m=30`,

```text
s+Nt < A_m t/2,
```

so this guaranteed number is smaller even than the lower bound
`r-s-Nt-1` for the dimension of numerical zeros in the common
restriction kernel. This comparison is about what dimension and
covolume arguments force; it does not assert that those specific
kernel directions are short.

To continue, one would need an additional global statement relating
the proper restriction images for different `S`, or a controlled
descent using the fact that actual relations form an ideal. It is
not valid to charge one independent codimension for every cut
without proving independence, nor to replace a proper image by
a zero image. No such global subspace estimate is established in
this note.
