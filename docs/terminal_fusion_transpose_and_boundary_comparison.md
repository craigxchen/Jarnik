# Transpose symmetry and the two boundary families

The higher-rank determinant numerator has the same coalition orders
after interchanging the rank `n` and the rectangular width `d`.
This follows from a tableau trace identity, without assuming
rank-level duality. It proves the strengthened height comparison for
every `n=2` or `d=2`. The comparison for all interior pairs
`n,d>=3` remains open.

The later [tensor-product floor](higher_rank_tensor_multiplicity_floor.md)
strengthens the union threshold `U` to an occupancy multiplicity `M`.
It gives `M(n,2)=B(n,2)` in every rank, so the entire width-two
family becomes equality for that stronger comparison. The transpose
identities below and the stated inequalities against `U` remain valid.

Use the notation and conformal-block identification of
[the higher-rank note](higher_rank_terminal_fusion_grid.md):

```text
m=n*d,                 kappa=n+d-1,
h(n,d)=level-(d-1) fundamental conformal-block rank,
c_n(omega_1)=(n^2-1)/n,
U(n,d)=2^(nd)-2*(2^d-1)^n+(2^d-2)^n.
```

Here `B(n,d)` is the determinant-covolume exponent from that note,
and `B(n,d)/h(n,d)` is its first-minimum upper exponent. The
arithmetic exclusion exponent is `U(n,d)`, under the stated
primitive Gaussian profile and all-coordinate residual hypotheses.
The present note is a comparison of these exponents; it does not
assert a radius-independent arc bound.

## 1. Fusion paths are tableaux with an ordered pair of corners

Lift a fundamental fusion path to partitions by adding one box at
each step. Nonnegative Dynkin coordinates say that the row lengths
remain weakly decreasing. A path of length `m=nd` ending at zero
has final shape the rectangle `d^n`. Thus the unrestricted paths
are standard Young tableaux of this rectangle.

During such a path all row lengths lie in `[0,d]`. The level
condition is

```text
lambda_1-lambda_n <= d-1.
```

It fails precisely when `lambda_1=d` and `lambda_n=0`, that is,
after the top-right box has appeared but before the bottom-left
box appears. Consequently

```text
h(n,d) = number of tableaux of shape d^n in which
         the bottom-left label precedes the top-right label.   (1)
```

Transposition exchanges the two distinguished corners, which have
different labels. If `F(n,d)` is the number of all tableaux, then

```text
h(n,d)+h(d,n)=F(n,d),
F(n,d)=(nd)! product_(i=0)^(n-1) i!
              / product_(i=0)^(n-1) (d+i)!.                    (2)
```

Only the alcove fundamental-step rule from the higher-rank note is
used here. There is no generic-position qualification in (1).

## 2. Casimir deficits are symmetric under transposition

For a shape `lambda` of size `s`, write its content as
`ct(lambda)=sum_(i,j in lambda)(j-i)`. Its `sl_n` Casimir is

```text
c_n(lambda)=n*s+2*ct(lambda)-s^2/n.
```

This expression is unchanged by adding full columns of height `n`,
as required when passing from partitions to `sl_n` weights.
Since transposition reverses content,

```text
c_n(lambda)+c_d(lambda^t)
       =(n+d)*s*(m-s)/m.                                      (3)
```

Let `C_s(n,d)` be the mean of `c_n(lambda_s)` over the tableaux
counted by (1). The corresponding mean over *all* rectangular
tableaux is

```text
mu_s(n,d)=s*(m-s)*c_n(omega_1)/(m-1).                           (4)
```

One proof of (4) is a trace calculation on the full invariant
space in `(C^n)^(tensor m)`. Each individual factor has Casimir
`c_n(omega_1)`. The total Casimir is zero, and permutations of the
tensor factors make the traces of all pair Casimir operators equal.
The partial Casimir on any `s` factors therefore has mean
`s*c_n(omega_1)-s*(s-1)*c_n(omega_1)/(m-1)`.
Branching successively by tensor factors identifies its eigenvalues
with the Casimirs of the intermediate tableau shapes, giving (4).

Define the unnormalized deficit

```text
J_s(n,d)=h(n,d)*(mu_s(n,d)-C_s(n,d)).                            (5)
```

Then

```text
J_s(n,d)=J_s(d,n)                  for every 0<=s<=m.           (6)
```

To see this, partition all tableaux into the admissible set in (1)
and its complement. The latter set transposes to the admissible
set for `(d,n)`. Its Casimir sum is
`h(d,n)*[(n+d)s(m-s)/m-C_s(d,n)]`, by (3).
The total sum is `F(n,d)*mu_s(n,d)`, by (4). Finally,
`mu_s(n,d)+mu_s(d,n)=(n+d)s(m-s)/m`. Rearranging proves (6).
No sign or positivity of the individual `J_s` is needed.

## 3. Every coalition order, and the whole height numerator, agree

The fusion residue can now be written as

```text
tau_s(n,d)=[-h(n,d)*s*(s-1)*c_n(omega_1)/(m-1)-J_s(n,d)]
              /(2*kappa).
```

Pair subtraction cancels its first term exactly. Hence, on the
coalitions where the residue computation applies directly,

```text
nu_s(n,d)=[binom(s,2)*J_2(n,d)-J_s(n,d)]/(2*kappa).             (7)
```

The Euler-degree formula likewise becomes

```text
delta(n,d)=(m-1)*J_2(n,d)/(2*kappa).                            (8)
```

Equations (6)--(8) prove equality of `delta` and of all orders with
`s<=m/2`. Complementary binary-coordinate symmetry then gives
equality of every remaining order. In fact (7) holds for those
orders as well: fusion duality gives `C_(m-s)=C_s`, and thus
`J_(m-s)=J_s`; the difference of the two right sides of (7) is
exactly `delta*(s-m/2)`.

We obtain the all-parameter identities

```text
delta(n,d)=delta(d,n),
nu_s(n,d)=nu_s(d,n)                 for every s,
B(n,d)=B(d,n).                                             (9)
```

There is also a useful cancellation formula:

```text
B(n,d)=sum_s binom(m,s)*J_s(n,d)/(2*kappa),
B(n,d)/h(n,d)
  =[m*c_n(omega_1)*2^(m-2)-sum_s binom(m,s)*C_s(n,d)]
       /(2*kappa).                                         (10)
```

Indeed `sum_s binom(m,s)*binom(s,2)=m*(m-1)*2^(m-3)`,
so the `J_2` term cancels the raw determinant degree exactly.

## 4. Explicit formulas when one parameter is two

At `d=2` the level is one and there is a unique admissible fusion
path. Its intermediate weight is `omega_(s mod n)`, with the
zero weight at multiples of `n`. Its Casimir is
`j*(n-j)*(n+1)/n`, where `j=s mod n`. Substitution gives

```text
h(n,2)=1,        delta(n,2)=1,
nu_s(n,2)=max(0,s-n).
```

The elementary binomial identity
`sum_(s>n) binom(2n,s)*(s-n)=n*binom(2n,n)/2` then yields

```text
B(n,2)=n*4^(n-1)-n*binom(2n,n)/2.                            (11)
```

By (2), the Catalan count for a two-row rectangle, and (9),

```text
h(2,d)=binom(2d,d)/(d+1)-1,
B(2,d)=d*4^(d-1)-d*binom(2d,d)/2.                            (12)
```

Thus the same numerator is divided by very different ranks in
the two boundary families.

## 5. All-size comparison with the combined arithmetic threshold

For every `n,d>=2` with `min(n,d)=2`,

```text
B(n,d)/h(n,d) >= U(n,d).                                     (13)
```

Equality holds exactly at `(n,d)=(2,2)` and `(3,2)`.

First consider `n=2`. Here `U(2,d)=2`. For `d>=3`, the ratio
`binom(2d,d)/4^d` is decreasing and is already `5/16<1/3`
at `d=3`. Write `C=binom(2d,d)`. Equations (12) give

```text
B(2,d)-2*h(2,d)
 =d*4^(d-1)+2-C*(d/2+2/(d+1))
 >4^d*(d/12-2/(3*(d+1)))+2 >0,
```

since `d*(d+1)>8`. The case `d=2` is equality.

Now consider `d=2`. In this case

```text
U(n,2)=4^n-2*3^n+2^n,
B(n,2)-U(n,2)
 =4^n*(n/4-1)-n*binom(2n,n)/2+2*3^n-2^n.                    (14)
```

For all integers `j>=1`,

```text
binom(2j,j)/4^j <= 1/sqrt(3j+1).                             (15)
```

The case `j=1` is equality. Induction follows from

```text
((2j+1)/(2j+2))^2 <= (3j+1)/(3j+4),
4*(j+1)^2*(3j+1)-(2j+1)^2*(3j+4)=j >=0.
```

For `n>=8`, (15) bounds the central binomial coefficient by
`4^n/5`. The right side of (14) is consequently at least

```text
4^n*(3n/20-1)+2*3^n-2^n >0.
```

The remaining exact differences for `n=2,...,7` are
`0,0,6,80,670,4522`. This proves (13) with its stated equality
cases.

The [exact checker](check_terminal_fusion_transpose_and_boundary.py)
builds the two tableau path spaces directly, verifies the Casimir
trace and transpose identities on a finite grid, and checks the
boundary formulas and integer inequalities through size 100.
The arguments above, rather than that finite range, prove the
two all-size boundary statements. No claim of (13) for all
`n,d>=3` is made here.
