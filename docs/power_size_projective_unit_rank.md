# Full projective-unit rank on a power-sized subset

The uniform cut extraction previously gave full multiplicative rank for
`k=O(log M)` selected rows. This note obtains full rank on
`k` of order `M^(1/4)`, without requiring every cut to occur or every
individual block to have a controlled weight. The proof keeps the
nearest-endpoint rounding errors as rational height errors, and uses
only a selected Gram matrix of quadratic sign functions.

This is a new structural consequence of an endpoint cluster. It does
not prove a uniform point count or improve the current general growth
bound. In particular, the resulting rank is of order `sqrt(M)`;
it cannot simply be converted into a stronger count estimate without
another arithmetic input.

## 1. Statement and inherited normalization

Take `M` distinct Gaussian lattice points in an arc of length at most
`C sqrt(R)`, where `0<C<=sqrt(2)`. Divide out the common Gaussian gcd
of the **whole** tuple. Its endpoint constant can only decrease.
Write the resulting primitive radius and split-prime allocations as

```text
N=R^2>1,       W=log N,
z_i=epsilon_i product_p pi_p^(a_ip) bar(pi_p)^(e_p-a_ip),
E_i=sum_p min(a_ip,e_p-a_ip) log p.
```

Let `k>=4`, `K=binom(k,2)`, and suppose

```text
M>=48,       K<=sqrt(M)/64.                    (1)
```

There is a same-unit subset of `k` rows for which the following holds.
For **every** real vector `c=(c_ij)` of zero row degrees, put

```text
sum_(j!=i)c_ij=0,
||c||_1=sum_(i<j)|c_ij|,
||c||_2=(sum_(i<j)c_ij^2)^(1/2).
```

For every integer such vector, the rational projective unit

```text
Lambda_c=product_(i<j)(z_i-z_j)^c_ij
```

satisfies the simultaneous lower bound

```text
h(Lambda_c) >= W ||c||_2/(32 sqrt(K)),          (2)
h(n/q)=log max(|n|,q)  for gcd(n,q)=1, q>0.
```

There is no cutoff on the integer coefficients. Consequently the
evaluated cross-ratios on these `k` rows generate a subgroup of
`Q^*/{1,-1}` of rank exactly `k(k-3)/2`.

For example, for sufficiently large `M`, one may take
`k=floor(M^(1/4)/8)`. All heights in (2) remain inherited from the
original primitive tuple; no new subset gcd is removed.

## 2. Good rows and two simultaneous quadratic estimates

The interior-allocation estimate proved in
[endpoint_allocation_rounding.md](endpoint_allocation_rounding.md)
and used in central truncation is

```text
sum_i E_i <= W sqrt(M)/2.
```

At least `M/2` rows therefore satisfy `E_i<=W/sqrt(M)`.
Take the largest literal Gaussian-unit class among them, of size
`L>=M/8`. The rows in this class have threshold signs

```text
f_i(p,t)=2 indicator(a_ip>=t)-1,
weight(p,t)=log p/W,
mu_ij=E_layers(f_i f_j)=1-2d_ij/W<=0.
```

All are unit vectors in the weighted real Hilbert space. For the
nonzero signed primitive residues set `u_ij=log|t_ij|>=0`.
The exact chord bound gives

```text
0<=u_ij<=-W mu_ij/4.                          (3)
```

As before, pairwise obtuseness gives the Bessel estimate
`sum_i <f_i,h>^2<=2||h||^2`. Conditioning on all but one selected
row consequently proves, for a uniformly sampled set of `r`
distinct row positions,

```text
E_rows |E_layers product_(r positions) f_i|^2
 <=2/(L-r+1).                                (4)
```

Sample an ordered `k`-tuple of distinct rows of this class. Index
the `K` edge products `f_i f_j` by unordered pairs, and let `A`
be their Gram matrix. Its diagonal is one. The symmetric difference
of two distinct edges has two or four vertices. Thus (4) gives

```text
E_rows ||A-I||_F^2 <= 2K(K-1)/(L-3).           (5)
```

There is a sharper estimate for the residues than bounding them all
by a common maximum. Since `-1<=mu_ij<=0`,
`mu_ij^2<=-mu_ij`. Also

```text
0<=||sum_(i=1..L) f_i||^2
  =L+sum_(i!=j)mu_ij,
mean_pair(-mu_ij)<=1/(L-1).
```

Together with (3), this gives

```text
E_rows sum_edges u_ij^2 <= K W^2/[16(L-1)].    (6)
```

Divide the two nonnegative quantities in (5)--(6) by their displayed
positive upper bounds and add. The expectation is at most two, so
some tuple has sum at most two. Fix that tuple. In particular both
estimates hold with twice their right sides:

```text
||A-I||_F <=2 sqrt(K(K-1)/(L-3)),
sum_edges u_ij^2 <=K W^2/[8(L-1)].             (7)
```

This is one common choice of rows for all subsequent coefficient
vectors, not a new selection depending on `c`.

## 3. Rounded sign products retain a positive Gram matrix

At each prime, round each allocation to its nearest endpoint:

```text
b_ip in {0,e_p},
|a_ip-b_ip|=min(a_ip,e_p-a_ip),
g_i(p)=2 b_ip/e_p-1.
```

Either choice is allowed at a tie. Regard `g_i` as constant over
the threshold layers of its prime. Then

```text
E_layers |f_i-g_i|=2E_i/W<=2/sqrt(M).
```

For a product of at most four signs, telescoping the product changes
its expectation by at most `8/sqrt(M)`. Let `B` be the Gram matrix
of the rounded edge products `g_i g_j`, with prime weights
`e_p log p/W`. Its diagonal still equals one. From (7),

```text
||B-I||_op
 <=2K/sqrt(L-3)+8K/sqrt(M)
 <=16K/sqrt(M)<=1/4.                          (8)
```

Here `M>=48` and `L>=M/8` imply `L-3>=M/16`.
In particular, for every real edge vector `c`,

```text
E_p (sum_edges c_ij g_i(p)g_j(p))^2
 >=(3/4)||c||_2^2.                            (9)
```

Only this Gram bound is asserted. The selected prime cuts need not
have full support, equal weights, or small error in every pattern.

## 4. The conductor height and its exact rounding error

For the zero-degree edge lattice, the exact rational factorization is

```text
Lambda_c=product_p p^beta_p(c) product_edges t_ij^c_ij,
beta_p(c)=sum_edges c_ij min(a_ip,a_jp).         (10)
```

It is proved, with all prime powers and residue signs, in
[coupled_crossratio_height_norm.md](coupled_crossratio_height_norm.md).
At a rounded prime let `S_p={i:b_ip=e_p}` and put

```text
C_p(c)=sum_(i<j, i,j in S_p)c_ij,
U_c=product_p p^(e_p C_p(c)).                   (11)
```

For integer `c`, this is an ordinary rational number. Zero row degrees
give the exact identity

```text
C_p(c)=(1/4) sum_edges c_ij g_i(p)g_j(p).       (12)
```

For a rational product of prime powers, its height is the maximum
of the positive and negative logarithmic valuation sums. It is
therefore at least half their total absolute sum. Using (12), then
`|sum c_ij g_i g_j|<=||c||_1` and (9), gives, for nonzero `c`,

```text
h(U_c) >= (W/8) E_p |sum_edges c_ij g_i g_j|
 >= (3W/32) ||c||_2^2/||c||_1
 >= (3W/(32 sqrt(K)))||c||_2.                  (13)
```

The crucial rounding estimate uses each row's own interior allocation.
The elementary Lipschitz inequality for the minimum gives

```text
|min(a_ip,a_jp)-min(b_ip,b_jp)|
 <=|a_ip-b_ip|+|a_jp-b_jp|.
```

Thus the height of the rational correction from the first product in
(10) to `U_c` is at most

```text
sum_edges |c_ij|(E_i+E_j)
 <=2W||c||_1/sqrt(M)
 <=2W sqrt(K)||c||_2/sqrt(M).                  (14)
```

No central truncation and no integral row multiplier is claimed here.
The error is explicitly a rational height error. In particular,
replacing (14) by a bound involving the sum of all `k` interior
weights would unnecessarily lose another factor of `k`.

The residue factor in (10) has height at most

```text
sum_edges |c_ij|u_ij
 <=||c||_2 (sum_edges u_ij^2)^(1/2)
 <=sqrt(2)W sqrt(K)||c||_2/sqrt(M),             (15)
```

where (7) and `L-1>=M/16` were used. This is the only residue
estimate needed; it holds uniformly for unbounded coefficient
vectors by Cauchy--Schwarz.

Equations (14)--(15) bound the combined correction height by
`4W sqrt(K)||c||_2/sqrt(M)`. The rational height triangle inequality
and (13) now prove

```text
h(Lambda_c)
 >= [3/32-4K/sqrt(M)] W||c||_2/sqrt(K)
 >= W||c||_2/(32 sqrt(K)),                     (16)
```

which is (2).

## 5. Consequence and the remaining growth-rate problem

The integer zero-degree edge lattice has rank `k(k-3)/2` and is
generated integrally by quadruple cross-ratio vectors. Since (16)
is strictly positive for every nonzero integer vector, evaluation
into `Q^*/{1,-1}` is injective. This proves the stated full rank.

The result needs only the original endpoint assumptions, the interior
allocation bound, and quadratic/fourth sign moments. It extends the
size of the subset with proved full multiplicative rank from
logarithmic to a fixed power of the original point count. It does
**not** extend the full independent-block normal form to that size:
no such exponential collection of nearly equal blocks has been
extracted here.

For `k` of order `M^(1/4)`, the rank has order `sqrt(M)`. Neither
that rank nor (16) is a contradiction, and the rational cross-ratios
still satisfy their universal additive exchanges. A further argument
must use those actual arithmetic values to improve the count. No
bound uniform in `R` is claimed.

## Verification

Two independent proof audits checked the common sample selection,
the all-unit normalization `L>=M/8`, the quadratic Gram perturbation,
and the two rational correction budgets. They verified the final
constant `1/32` in (16), including the `sqrt(2)` residue-error factor.

An exact Gaussian-arithmetic check on 120 configurations with four
through seven rows verified the rational identity (10), the zero-degree
cut identity (12), and the prime-by-prime minimum-rounding inequalities.
It used split primes `5,13,17,29,37`, exponents one through three,
and integer combinations of four-cycle coefficient vectors, with
seed 27891. These configurations test the unconditional algebraic
steps; they are not claimed to satisfy the endpoint selection bounds.
The latter are supplied by the audited proof, not by enumeration.
No Lean formalization is claimed.
