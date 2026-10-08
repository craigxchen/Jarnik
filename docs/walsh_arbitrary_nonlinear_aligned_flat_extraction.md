# Every nonzero binary vector field has a growing aligned affine flat

Let `a:F_2^t -> F_2^t minus {0}` be **any** function, and put
`M=2^t`. No degree, affinity, balance, or injectivity is assumed.
There are a subspace `V`, a coset `P=x0+V`, and a nonzero
functional `ell in V*` such that

```text
a(x)|V = ell             for every x in P.                 (1)
```

For `t>=4096`, these can be chosen with

```text
h=|V| >= t / [16 (log_2(t+2))^2].                          (2)
```

This removes the nonlinear-assignment obstruction to a **growing aligned
flat**; it does not produce an affine approximation on almost all rows.
The [balanced 32-row example](walsh_nonlinear_coset_height_obstruction.md)
with no aligned plane is compatible with the large threshold here.

For repeated-Walsh endpoint profiles with arbitrary nonlinear assignments,
the existing aligned-flat phase certificate consequently gives

```text
log R^2 >= c_C M (log M)^2 / (loglog M)^2,                  (3)
M = O_C(1 + log R (logloglog R)^2 / (loglog R)^2),          (4)
```

above a fixed absolute logarithm cutoff. This arithmetic consequence
still assumes the **actual near-equal physical prime logs** in (15)
below. Independent row units and all common Gaussian content are allowed.
Arbitrary physical prime weights remain outside this consequence, and
the uniform endpoint count remains unproved.

## A rank-to-density fact

Let `K` be a symmetric `0/1` matrix of order `n`, with diagonal one
and rank at most `r` over `F_2`. Make a simple graph from its
off-diagonal ones. Any independent vertex set gives an identity principal
submatrix, so the independence number `alpha` is at most `r`.
For degrees `d_i`, a random ordering and selection of vertices earlier
than all their neighbors gives an independent set of expected size
`sum_i 1/(d_i+1)`. Thus Cauchy--Schwarz gives

```text
r >= alpha >= sum_i 1/(d_i+1) >= n^2/sum_i(d_i+1),
number of ones in K >= n^2/r.                             (5)
```

The last count includes the diagonal. Rank is over the binary field;
the entry count is an ordinary integer count. No real positive
semidefiniteness is asserted.

## Starting with many aligned nonzero edges

For `x,y in F_2^t`, set

```text
Q_xy = a(x) dot (x+y)       in F_2,
A_xy = 1 + Q_xy + Q_yx      in F_2.
```

The matrix `A` is symmetric with diagonal one. Expansion of `Q`
as `a(x) dot x + a(x) dot y` gives

```text
rank_F2 A <= r0 := 2t+3.                                  (6)
```

Because every `a(x)` is nonzero, each row of `Q` has `M/2` ones.
Among all ordered pairs, let `n_ij` count cases
`(Q_xy,Q_yx)=(i,j)`. Comparing their total count with the total
entries in `Q` and its transpose yields `n_00=n_11`.
The ones of `A` are exactly these two cases, so (5)--(6) imply

```text
n_11 >= M^2/(2r0).                                       (7)
```

There are no diagonal pairs of type `11`. Averaging by `w=x+y`
over the `M-1` nonzero directions gives a `w` for which

```text
B1={x:a(x) dot w=a(x+w) dot w=1}
```

has size at least `M/(2r0)`. It is invariant under translation by
`w`. Its cosets of `V1=<w>` have the same nonzero restriction
`ell1(w)=1` and union row density

```text
rho1 >= 1/(4t+6).                                        (8)
```

The nonzero restriction begins here and survives all later extensions.
Constructing only zero-restriction flats would not give the required
zero-sum character certificate.

## Extending an aligned family in one common direction

Suppose `V` has size `h=2^i` and `N` distinct `V`-cosets are
aligned with the same nonzero restriction `ell`. Write `rho=Nh/M`
for their union row density. Choose a linear complement `L` to `V`
and the unique representative `x in L` of each selected coset.
For representatives `x,y`, define over `F_2`

```text
K_xy = product_(v in V) [1+a(x+v) dot (x+y)]
                       [1+a(y+v) dot (x+y)].              (9)
```

The matrix is symmetric with diagonal one. For fixed `v`, the first
factor expands as

```text
[1+a(x+v) dot x] + sum_(j=1)^t a(x+v)_j y_j,
```

whose matrix has rank at most `t+1`, hence at most `t+2`.
The second factor is its transpose. An entrywise product has rank
at most the product of the ranks, by multiplying rank-one outer-product
expansions. Consequently

```text
rank_F2 K <= r_i := (t+2)^(2h).                           (10)
```

If `N>=2r_i`, (5) gives at least `N^2/(2r_i)` **ordered
off-diagonal** compatible pairs. For such a pair set `w=x+y`.
The linear complement ensures `w in L minus {0}`, hence `w` is
outside `V`. Every actual assigned label on **both** cosets annihilates
`w` by (9), and restricts to `ell` on `V`. Their union is exactly
the coset `x+V'`, where `V'=V+<w>`, aligned with

```text
ell'|V=ell,       ell'(w)=0.
```

Thus `ell'` is still nonzero. Each factor in (9) uses the assigned
label at that selected row itself, not an opposite-coset label.

There are `M/h-1` possible nonzero directions in `L`. A popular one
has at least `N^2 h/(2r_i M)` ordered compatible pairs. Once this
direction is fixed, each representative has the unique partner `x+w`.
Each new coset is counted exactly twice, and different pairs give
disjoint new cosets. Their number is at least `N^2 h/(4r_i M)`, so

```text
rho' >= rho^2/(2r_i) = rho^2/[2(t+2)^(2h)].               (11)
```

One may choose a new complement at the next step. This changes
representatives but does not change the actual restrictions on a coset.

## Explicit iteration and density retained

Put `L_t=log_2(t+2)` and `A_i=-log_2 rho_i`. For `t>=2`,
(8) gives `A_1<=2L_t`. The recurrence (11) gives

```text
A_(i+1) <= 2A_i+1+2^(i+1)L_t,
A_i <= 2^i i L_t+2^(i-1)-1.                              (12)
```

Divide by `2^i` and sum to obtain the second inequality. Since
the number of selected cosets is `N_i=rho_i 2^(t-i)`, a sufficient
condition for `N_i>=2r_i` is

```text
t >= i+2^i(i+2)L_t+2^(i-1).                              (13)
```

A common aligned family of dimension `d>=1` thus exists whenever
(13) holds for `1<=i<d`. It suffices to impose it at `i=d`,
because the right side increases with `i`.
For `t>=4096`, choose

```text
d=floor log_2[t/(8L_t^2)],       h=2^d.                   (14)
```

Then `d>=1`, `d<=L_t`, and `t/(16L_t^2)<=h<=t/(8L_t^2)`.
The right side of (13) at `i=d` is at most

```text
L_t+t[1/8+1/(4L_t)+1/(16L_t^2)] < t.
```

Indeed `L_t<=t/4`, `L_t>=1`, and the bracket is at most `7/16`.
Also `t>=16L_t^2` at `4096` and persists afterward, since
`t/(log_2(t+2))^2` is increasing there. This proves `d>=1` and
(2). For smaller `t>=1` the start still gives an aligned edge;
no aligned-plane claim is made.

The family, not just one flat, is retained: (12) supplies
`rho_d>=2^(-h d L_t-h/2+1)`. In the range (14), it in particular
gives `rho_d>=M^(-1/4)`, since `h d L_t+h/2<=t/4`.
Thus at least `M^(3/4)` rows lie in disjoint aligned cosets for the
same `V,ell`. This may be useful for a future weighted extraction.

## Actual near-equal prime weights

Take `b>=5` physical copies of each nonzero Walsh label, flip one
assigned entry at every row, and require all assigned physical columns
distinct. The label map can be completely nonlinear and can repeat
labels within physical copy capacity. Set `r=b(M-1)`. Use arbitrary
orientations of pairwise distinct split Gaussian primes with rational
norms `p_j`, and assume the actual logs obey

```text
u <= log p_j < u+1/r,       u>=log 5.                      (15)
```

Set `W0=sum_j log p_j`, `D=log Norm(g)>=0` for any common Gaussian
factor `g`, and `W=log R^2=W0+D`. Independent row units are allowed.
Assume the actual rows lie on an arc of length at most `C sqrt(R)`.
For (1), use `c_(x0+v)=(-1)^ell(v)` and zero elsewhere. These
coefficients have sum zero and absolute sum `h`.

The baseline physical half-coefficient at label `a` is
`(h/2)(-1)^(a dot x0)` when `a|V=ell`, and zero otherwise.
Exactly `M/h` labels qualify, none zero. Each selected row's actual
assigned label qualifies, and its distinct physical flip lowers the
absolute coefficient by one. Hence exactly

```text
sum_j |v_j| = bM/2-h.                                    (16)
```

The associated Gaussian integer `Beta` is conjugate-primitive and a
nonunit: physical rational primes are distinct, and `bM/2-h>0`.
The prime-log error is less than `(bM/2-h)/r<1`. Therefore

```text
log Norm(Beta)-W/2 < -(h-b/2)u+1-D/2.                    (17)
```

The exact signed source product is `unit * Beta/bar(Beta)`. Lifting
arguments on the arc, the zero row sum cancels the central argument;
the argument of `Beta` is within `h C exp(-W/4)/4` of `(pi/4)Z`.
A conjugate-primitive nonunit Gaussian integer has a nonzero integer
distance numerator from every axis or diagonal. The elementary gap gives

```text
log Norm(Beta) >= W/2+log 8-2log(hC).                     (18)
```

Thus, with the common content retained in the source radius,

```text
(h-b/2)u+D/2 < 1+2log(hC)-log 8.                          (19)
```

For `h>=H_C=max(16,ceil(8log^+ C/log 5))`, this excludes `b<=h`:
the left side would be at least `(h/2)log 5`, while
`2log h<=(h/4)log 5`, `2log^+ C<=(h/4)log 5`, and
`1-log 8<0`. Hence `b>h` for sufficiently large `t`, depending
only on `C`. Distinctness of the `r` physical prime norms then gives

```text
W >= W0 >= log((r+1)!) >= (r/2)log(r/2)
                         >= (hM/4)log M.                 (20)
```

The last inequality holds in the large-`M`, `h>=16` regime here.
Insert (2) and `t=log_2 M` to prove (3). Inverting
`M(log M)^2/(loglog M)^2` proves (4). This improves the growth
scale `M=O(W/log W)` within the stated family for all nonlinear maps.

## The remaining arbitrary-weight gap

For general positive physical weights the same certificate satisfies
exactly

```text
log Norm(Beta)=(h/2) W_(V,ell)-F_P,                       (21)
```

where `W_(V,ell)` totals weights of labels restricting to `ell`,
and `F_P` totals flipped-prime weights on the flat. The individual
label-mass bound from [the weighted affine theorem](walsh_three_family_weighted_affine_growth.md)
only gives roughly `(h/2)W_(V,ell)<=W0`, not `W0/2`.
The extraction gives a dense family for one nonzero restriction, not
balanced families over all restrictions. It therefore does not close
this factor-of-two gap, and (3)--(4) are not proved for arbitrary weights.

There is nevertheless a weighted necessary condition. Distinct assigned
primes on the at least `M^(3/4)` retained rows have total weight at
least `c M^(3/4)log M`. Averaging over their disjoint aligned flats
gives a flat with `F_P>=c h log M` (one may use `c=1/4` in this
large regime). Together with (18) and (21),

```text
W_(V,ell) >= W/h + (1/2)log M - O_C(log h/h).             (22)
```

This keeps arbitrary actual prime weights and common content, but the
excess is in a single restriction fibre and is not a growth bound.

The [checker](check_walsh_arbitrary_nonlinear_aligned_flat_extraction.py)
checks binary ranks and counts, quotient representatives and merges,
direction and ordered-pair losses, the recurrence cutoff, and exact
physical-column heights. Fixtures are not asserted to be endpoint arcs.
