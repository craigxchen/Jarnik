# The exact Specht filtration of the smaller-cut kernels

The algebraic kernels in the [all-cut hierarchy](all_cut_invariant_relation_hierarchy.md)
have an exact, multiplicity-free description. This determines their
dimensions and their last nonzero depth. It does not put short numerical
relations into successive kernels: the restriction criterion at any
depth proves only that an image is proper, not that it is zero.

Let `m=2q`, and let `V_m` be the rational space of simultaneous `SL_2`
invariants of degree two in each of the `m` labelled binary rows. For
`k>=0`, let `K_k` consist of polynomials whose restrictions vanish
whenever at least `q-k+1` rows are collapsed to one common direction.
Thus `K_0=V_m`; for `k>q` the condition includes an empty inside set
and `K_k=0`. Write `[2a,2b,2c]` for the irreducible `S_m` module with
that partition, allowing `c=0` and `b=0`.

**Filtration theorem.** Over `Q`, with `a>=b>=c>=0` and `a+b+c=q`,

```text
V_m = direct_sum_(a>=b>=c>=0, a+b+c=q) [2a,2b,2c],
K_k = direct_sum_(a>=b>=c>=k, a+b+c=q) [2a,2b,2c].       (1)
```

In particular, each constituent has a nonzero restriction with exactly
`q-c=a+b` collapsed rows, and vanishes with any larger number. Its
occurrence in `K_k/K_(k+1)` is therefore characterized by `c=k`.

## 1. The spin-one and Gram-minor model

Set `E=Sym^2(C^2)`, a three-dimensional vector space. Its invariant
nondegenerate symmetric form makes the image of `SL_2` equal to `SO(E)`;
the nonzero null vectors are the squares of binary linear forms. The
space `V_m(C)` is the multilinear part of the degree-`m` invariant
polynomials on `E^m`. Since `m` is even, `SO(E)` and `O(E)` have the same
invariants in this degree: an element of `O(E)` outside `SO(E)` is
`-I` times an element of `SO(E)`, and `-I` acts by `(-1)^m=1`.

Let `v_1,...,v_m` be generic vectors of `E`. The degree-`m` polynomial
space on them has its usual `GL_m` action by mixing vector slots. Cauchy
decomposition followed by `O(E)` invariants gives one `GL_m` irreducible
of highest weight `(2a,2b,2c)` for each partition `a+b+c=q`, and no
others. One can also see the nonzero highest-weight vectors explicitly.
Let `G=(<v_i,v_j>)`, and let `Delta_j` be the leading principal `j` by
`j` Gram determinant, for `j=1,2,3`. Then

```text
f_(a,b,c) = Delta_1^(a-b) Delta_2^(b-c) Delta_3^c       (2)
```

is a nonzero highest-weight vector of weight `(2a,2b,2c)`. The
multilinear weight `(1,...,1)` of that `GL_m` irreducible is the Specht
module `[2a,2b,2c]`. This also follows from [Howard--Millson--Snowden--Vakil,
Proposition 6.5 and its proof](https://ems.press/content/serial-article-files/31806),
which state that the degree-two invariant space is multiplicity-free
with precisely the even partitions of at most three parts and identify
the spin-one representation with `SO_3`.

## 2. Exact maximum number of collapsed null rows

Fix a nonzero isotropic vector `u` and substitute

```text
v_i=w_i+t_i u.
```

Put `p_i=<u,w_i>`. The highest total `t` degrees in the three Gram
minors are `1,2,2`, respectively. More explicitly, their top parts are

```text
Delta_1 : 2t_1 p_1,
Delta_2 : -(t_1 p_2-t_2 p_1)^2,
Delta_3 : det(form) (sum_(i=1)^3 t_i det(w_1,...,u at i,...,w_3))^2.
```

Each displayed polynomial is nonzero for generic `w_i`; therefore
the top part of (2) is nonzero and its degree is exactly

```text
(a-b)+2(b-c)+2c = a+b = q-c.                         (3)
```

This degree bound holds for the entire `GL_m` constituent, not only
for its highest-weight vector. Indeed, after any invertible linear
mixing of the vector variables, the added multiples of `u` are still
linear forms in the `t_i`; taking linear combinations and weight
components cannot raise their total degree.

For clarity, the top degree survives passage to the multilinear
weight space. Split the slots into groups of sizes `2a,2b,2c`, replace
`v_i` in (2) by the sum of the independent variables in group `i`, and
take the coefficient multilinear in all `m` variables. This is a
weight projection of a `GL_m` translate of (2): the sums for the
nonempty groups are independent linear forms and can be completed to
an invertible `m` by `m` mixing matrix. Thus it belongs to the same
irreducible constituent, including when `b=0` or `c=0`. If a top
monomial of (2) has exponents
`e_1,e_2,e_3` in the `t_i`, set `e_i` distinct variables in group `i`
equal to `u` and the others equal to `w_i`. The resulting polarized
evaluation is a nonzero factorial multiple of that top coefficient.
It has exactly `e_1+e_2+e_3=a+b` collapsed slots. The remaining slots
may also be chosen on the null cone: the multilinear polynomial is
nonzero for some arbitrary remaining vectors, and null vectors span
`E` (the squares of three distinct binary linear forms already do).
Nonvanishing is open, so those outside null directions can be taken
distinct from one another and from the common direction if desired.

Conversely, any multilinear polynomial in this constituent vanishes
when more than `a+b` variables equal `u`, since that evaluation would
require a coefficient of total `t` degree greater than (3). These
facts hold for every subset of labels, by the permutation action.
The kernel definition now gives (1): vanishing at `q-k+1` collapsed
rows is equivalent to `a+b<q-k+1`, or `c>=k`. Nonvanishing at the
maximum also gives nonvanishing when any smaller prescribed number of
rows is collapsed, by leaving the extra rows free and choosing them
equal to `u`.

## 3. Dimensions and the exact layers

Let `r_k=dim K_k`. The hook-length formula gives

```text
f^(2a,2b,2c)
 = m! (2a-2b+1)(2a-2c+2)(2b-2c+1)
   / ((2a+2)!(2b+1)!(2c)!),
r_k = sum_(a>=b>=c>=k, a+b+c=q) f^(2a,2b,2c),
D_k = r_k-r_(k+1)
    = sum_(a>=b>=k, a+b+k=q) f^(2a,2b,2k).          (4)
```

There is a constituent with `c=k` exactly when `3k<=q`. Hence the
last nonzero kernel is `K_floor(q/3)`, and
`K_(floor(q/3)+1)=0`, sharpening the dimension statement behind the
termination argument in the all-cut hierarchy. For example,
`dim K_1(6)=5`, `dim K_1(8)=56`, and
`dim K_1(10)=225+252=477`, while `K_2(10)=0`.

## 4. Conditional use in height comparisons

The dimensions in (4) can be used to locate which successive minimum
would force a numerical zero relation outside the next kernel. Suppose
evaluation at the actual distinct directions is nonzero on each
nonzero `K_k`, as witnessed by the triangle/doubled-edge polynomial in
the all-cut hierarchy. Then the numerical-zero lattice in `K_k` has
rank `r_k-1`. If `r_(k+1)>0`, its sublattice in `K_(k+1)` has rank
`r_(k+1)-1`. Thus `r_(k+1)` independent numerical-zero relations in
`K_k` force one outside `K_(k+1)`.

The [higher-kernel evaluation-height calculation](higher_kernel_evaluation_heights.md)
supplies such an aggregate bound under its full core-profile hypotheses.
Its coefficient is

```text
A_0=m*2^(m-2)-q*binomial(m,q),
A_k=A_0-sum_(s=0)^m binomial(m,s) max(0,k-|s-q|).
```

The primitive evaluation vector has logarithmic norm
`A_k w+O_m(epsilon w)+O(log B)+O_m(1)` there. For fixed `m`, the
successive-minima theorem transfers this to the sum of the logarithms
of the `r_k-1` ordered coefficient-lattice minima `lambda_j`, up to
`O_m(1)`. In a regime where the displayed errors are `o(w)`, this
gives `sum_j log lambda_j <= A_k w+o(w)`, and hence

```text
log lambda_(r_(k+1)) <= (A_k/D_k)w+o(w),
       D_k=r_k-r_(k+1),                 if r_(k+1)>0.  (5)
```

This uses the `D_k` minima from index `r_(k+1)` through `r_k-1`, all
at least `lambda_(r_(k+1))`. At the terminal step
`r_(k+1)=0`, one needs the *first* minimum instead, and the
denominator is `r_k-1`, provided `r_k>1`. If `r_k<=1`, evaluation
leaves no nonzero numerical-zero relation in that kernel.

Equation (5) is conditional on the stated aggregate height bound; the
representation-theoretic result alone supplies no short relations.
Even if its exponent is small, a proper restriction-image criterion at
depth `k` does not imply zero image or membership in `K_(k+1)`.
There is no automatic short lift from a quotient kernel.

The finite checks of the top degrees and dimensions are in
[check_kernel_specht_filtration.py](check_kernel_specht_filtration.py).
