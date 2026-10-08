# Positive sum kernels and their exact denominator cost

The symmetric sum kernel does exploit the shrinking endpoint angle:
its normalized entries approach one quadratically. The denominator cost
can also be determined exactly. In the primitive common-unit model,

```text
LCD of all normalized secant permutation products = L_X^2.       (1)
```

Here `L_X` is the matching least common multiple of the primitive real
parts from [eight_point_pfaffian_content.md](eight_point_pfaffian_content.md).
That note proves `X/L_X | T`, where `X=product x_ij` and
`T=product |t_ij|`. Thus (1) has logarithmic height `14W+o(W)` in a
critical fair eight-row profile. It is not merely an upper estimate by
a convenient product denominator.

As a consequence, every fixed-degree positive Taylor remainder, when
cleared term by term in the natural way, produces a large positive
integer rather than an integer smaller than one. A further cancellation
in the denominator of the final rational expression is not excluded.

## 1. The real normalization of actual sums

Take eight distinct points in a primitive common Gaussian-unit class on
an arc of angular width `Theta<pi`. Choose their primitive half-angle
numerators with positive real part:

```text
h_ij=x_ij+i t_ij,       x_ij>0,
q_ij=N(h_ij)=x_ij^2+t_ij^2,
gcd(x_ij,t_ij)=1,
q_ij=product_p p^|a_i(p)-a_j(p)|.
```

Set diagonal values `x_ii=q_ii=1`. With consistent square roots along
the arc, put

```text
S_ij=2sqrt(z_i z_j)/(z_i+z_j)
    =sec((theta_i-theta_j)/2)=sqrt(q_ij)/x_ij,
S_ii=1.
```

Every permutation product is positive and rational. Indeed

```text
Y_sigma=product_i S_(i,sigma(i))=A_sigma/X_sigma,
A_sigma=product_p p^[(1/2)sum_i |a_i-a_sigma(i)|] in Z_(>0),
X_sigma=product_i x_(i,sigma(i)) in Z_(>0).                       (2)
```

The exponents in `A_sigma` are integers because the sum of positive
allocation increments equals the sum of negative increments around
the permutation cycles. This retains all prime powers.

Let `D_sum` be the least positive integer clearing each of the `8!`
rational numbers `Y_sigma` separately. It is not defined to be the
reduced denominator of their sum.

## 2. Proof of the exact denominator formula

Fix a rational prime `p` and give edge `ij` weight `v_p(x_ij)`.
The graph of its positive-weight edges is bipartite, component by
component. At an odd conductor prime, edges between different allocation
groups have weight zero. Within one allocation group, divide out its
common local Gaussian factor: positive edges join opposite nonzero
classes modulo `p`. At an odd nonconductor prime the same statement
holds without an allocation partition.

The prime two requires the shifted modulus. The primitive circle has
odd norm, its Gaussian gcd factors are units modulo every power of two,
and

```text
z_i+z_j=2 epsilon g_ij x_ij.
```

Thus `2|x_ij` means `z_i==-z_j modulo4`. For an odd-norm class `c`,
the classes `c` and `-c` are distinct modulo four. Hence this positive
edge graph is bipartite too. Higher thresholds use `2^(l+1)` for
`v_2(x_ij)>=l`; opposition modulo two alone would give the wrong graph.

Let `nu_p` be the maximum total edge weight of a perfect matching on
the eight labels. A permutation cycle cover has positive-weight edges
forming paths and even cycles, with double edges allowed. Its edges can
be split into two matchings, so its total weight is at most `2nu_p`.
The denominator exponent of (2) is at most that weight, since the
numerator `A_sigma` can only cancel factors. Therefore

```text
v_p(D_sum)<=2nu_p.
```

For the opposite inequality, choose a maximum-weight matching, transpose
the endpoints of its positive-weight edges, and fix every remaining
vertex. On such an edge, `p|x_ij` implies `p` does not divide `q_ij`,
because `q_ij=x_ij^2+t_ij^2` and `gcd(x_ij,t_ij)=1`. The transposition
contributes precisely `q_ij/x_ij^2`, with no `p` cancellation. Fixed
points contribute one. The reduced denominator exponent is `2nu_p`.
Using zero-weight matching edges instead of fixed points could introduce
unwanted norm factors; the construction avoids that issue explicitly.

Since `v_p(L_X)=nu_p`, this proves at every prime that

```text
D_sum=L_X^2=(X/B_X)^2,       B_X=X/L_X | T.                      (3)
```

In particular `X^2/T^2<=D_sum<=X^2` as real numbers. This exact
description retains the real-part contents and all odd, inert, and
ramified residue primes.

## 3. The positive quadratic deviation

For `Theta<=1`, each term in (2) satisfies

```text
1<=Y_sigma<=exp(2Theta^2).
```

For example, `log sec u<=u^2` for `|u|<=1/2`, and each permutation
uses eight differences of magnitude at most `Theta/2`. At least one
transposition has value strictly larger than one. Therefore

```text
0<per(S)-8! <=7*8!*Theta^2,
U=D_sum(per(S)-8!) in Z_(>0).                                  (4)
```

The near-equality is real and arithmetic: the scaled remainder is a
nonzero positive integer. It only gives
`D_sum>=1/(7*8!*Theta^2)`, however.

The symmetric Cauchy identity also gives, for eight rows,

```text
det(S)=product_(i<j)tan^2((theta_i-theta_j)/2)=(T/X)^2.           (5)
```

Consequently its reduced denominator is
`(X/gcd(X,T))^2`. In a fair endpoint profile this has the same leading
height as (3). The standard squared-kernel identity
`det(S entrywise squared)=det(S) per(S)` then introduces the same
permanent as (4); it does not remove the real-part denominator.

## 4. Every fixed positive Taylor degree has an explicit budget

Put

```text
u_ij=t_ij^2/q_ij=sin^2((theta_i-theta_j)/2).
```

Each permutation product is
`product_edges(1-u_ij)^(-m_ij/2)`, where `m_ij` is zero, one, or two.
All coefficients of its convergent multivariate Taylor series are
nonnegative rational numbers. Let `P_d(u)` be the total-degree-`d`
Taylor polynomial of `per(S)`, and put

```text
R_d=per(S)-P_d(u)>0.                                           (6)
```

For each fixed `d`, the coefficient of `u_ij^(d+1)` is positive:
permutations transposing `i,j` already contribute the geometric series
`(1-u_ij)^(-1)`. Thus when every pair has
`log u_ij=-W/2+o(W)`,

```text
log R_d=-(d+1)W/2+o(W).                                       (7)
```

Now define `H_d` to be the least common denominator of every permutation
product in (2) and every individual evaluated monomial of `P_d`. This
specifies the clearing procedure exactly. Since the circle is primitive,

```text
Q=lcm_(i<j)q_ij=R_primitive^2.                                (8)
```

For a fixed degree, the coefficient denominators and numerators cost
only a fixed factor. The monomials `c_d u_ij^d` have a nonzero coefficient
`c_d` independent of the pair, while `gcd(t_ij,q_ij)=1`. Hence their
common denominator is at least `Q^d` divided by a fixed factor depending
on `d`. Conversely every degree-at-most-`d` monomial denominator divides
a fixed multiple of `Q^d`.

Suppose the inherited conductor distribution is fair on eight rows, with
pair residues of logarithmic size `o(W)`, as in the critical reduction.
Here `W` denotes the inherited weight, before removal of the two constant
patterns. Each of those patterns has weight `W/256+o(W)`. Common Gaussian
division preserves every primitive pair numerator, so
`log Q=(127/128)W+o(W)`, while `log q_ij=W/2+o(W)` for each pair.
The endpoint angle gives
`log x_ij=W/4+o(W)`. Therefore

```text
log X=7W+o(W),       log T=o(W),
log D_sum=14W+o(W).
```

Combining (3), (7), and (8), for every fixed integer `d>=0`, gives

```text
log H_d >= max(14,127d/128)W+o(W),
H_d R_d in Z_(>0),
log(H_d R_d)
 >=[max(14,127d/128)-(d+1)/2]W+o(W)
 >=(13/2)W+o(W).                                             (9)
```

For integer `d`, the minimum coefficient still occurs at `d=14`.
For `d<=14` the first term gives the bound; for `d>=15`, the coefficient
is at least `63d/128-1/2>13/2`. Thus increasing the Taylor
order does not produce a small nonzero integer under this natural
termwise denominator clearing. The failure persists at every fixed
degree, even though the remainder is positive and exploits the actual
shrinking angle to its full order.

This is not a theorem about the reduced denominator of the final rational
number `R_d`. Arithmetic cancellation in that sum might make its reduced
denominator much smaller than `H_d`. Proving and exploiting such a
reduction would be additional information, beyond the identities,
positive coefficients, and denominator bounds proved here. Neither (9)
nor the preceding Cauchy formula rules out other rational identities or
the desired uniform endpoint bound.

Two subsequent tests clarify this remaining scope. The
[local cancellation construction](reduced_sum_kernel_cancellation.md)
shows that arbitrarily deep final cancellation can occur at a prime
which divides no old pair residue, on actual nonendpoint circles.
The [shrinking-arc Pell construction](shrinking_arcs_with_bounded_pair_residues.md)
shows that the final reduced denominator itself can grow while all
pair residues stay bounded; its three-point specialization is already
at endpoint scale. These results do not exclude a special eight-point
endpoint estimate, but require such an estimate to retain the relevant
cardinality and conductor hypotheses.

## Exact checks and audit

[check_sum_kernel_denominator_budget.py](check_sum_kernel_denominator_budget.py)
checks all 40,320 permutation products on each of three actual equal-norm
Gaussian eight-tuples. It verifies the exact common denominator (3),
including nontrivial two-adic factors, the normalized symmetric
Cauchy--Borchardt relation, and positivity and integrality of the first
Taylor remainder. The test tuples are not asserted to be endpoint families.

The root agent independently checked the fixed-point completion in the
denominator-attaining permutation, the shifted modulus at two, and the
scope of the Taylor clearing argument. The inherited-versus-primitive
normalization in (8)--(9) includes the two constant fair-profile patterns;
replacing `127W/128` by `W` there would silently lose that distinction.
