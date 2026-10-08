# A positive Hankel determinant at the full Boolean endpoint scale

The weighted moment matrix below detects repeated cotangent residues in a
way the positive even/odd Bezoutian does not. Its integer minors have
exact, simultaneous divisibility by every Boolean core block. Those
guaranteed divisors nevertheless fall short of the archimedean size for
every proper minor. The top minor is the ordinary squared Vandermonde and
is exactly critical. This is a budget audit, not an exclusion of extra
arithmetic divisibility or a proof of the uniform circle bound.

## 1. Positive integer minors from the actual Gaussian rows

Fix `m>=2` distinct actual conjugate-primitive anchor numerators

```text
P_i=x_i+i y_i,                 1<=i<=m,
P_i=K_i product_(T containing i) H_T,
n_T=Norm(H_T),
N(P_i)=x_i^2+y_i^2,
Delta_ij=x_i y_j-x_j y_i=G_ij t_ij !=0,
G_ij=product_(T containing i,j) n_T.
```

Assume every `y_i!=0`. The `n_T` are pairwise coprime core norms indexed
by nonempty `T subset [m]`, and `n_T|N(P_i)` whenever `i in T`.
Correcting factors may overlap core primes. Choose any positive integer
`L` divisible by every `|y_i|` and every `|t_ij|`, and put

```text
X_i=L x_i/y_i in Z,       C_i=X_i^2+L^2=(L/y_i)^2 N(P_i) in Z_(>0).
```

The `X_i` are distinct because `Delta_ij!=0`. Their pair cotangents
are also integers: writing
`z_ij=Re(conjugate(P_i)P_j)/G_ij in Z`, one has

```text
(X_i X_j+L^2)/(X_i-X_j)=L z_ij/t_ij in Z.           (1)
```

The signs of the `X_i` are irrelevant for the matrix argument. For
`i in T`, `n_T|C_i`. For `i,j in T`, the exact difference formula is

```text
X_i-X_j=L G_ij t_ij/(y_i y_j).                        (2)
```

At any prime `p|n_T`, conjugate primitivity gives `p` coprime to
`y_i y_j`: an integer-coordinate primitive row has
`gcd(y_i,N(P_i))=1`. Thus (2) proves `n_T|(X_i-X_j)`, including
arbitrary prime-power depth and primes dividing `L`. No assumption
that corrections avoid the core is used.

For `1<=r<=m`, form the real symmetric Hankel matrix

```text
H_r[a,b]=sum_(i=1)^m X_i^(a+b)/C_i,     0<=a,b<r.
```

It is positive definite, since it is the Gram matrix of the
independent monomial-evaluation columns with positive weights `1/C_i`.
Cauchy--Binet gives the positive integer

```text
D_r=(product_i C_i) det H_r
   =sum_(S subset [m], |S|=r)
      [product_(i<j in S)(X_i-X_j)^2]
      [product_(i notin S) C_i].                       (3)
```

## 2. Exact core divisor of each minor

For a core block `T` and a summand indexed by `S`, write
`s=|T|`, `t=|T intersect S|`. Equation (2) and `n_T|C_i` show that
the summand is divisible by `n_T` to the exponent

```text
f_T(S)=(s-t)+t(t-1).                                  (4)
```

The first term counts norm factors outside `S`; the second counts
twice the pairs inside `T intersect S`. Therefore, without relying
on cancellation between summands,

```text
product_(nonempty T) n_T^c_T(r) | D_r,
c_T(r)=min_(|S|=r) f_T(S).                             (5)
```

Because the core norms have disjoint rational-prime support, their
divisibilities multiply. Formula (5) is exact at each block prime
power as a *guaranteed* divisor; an individual `D_r` may have extra
divisibility.

In a balanced limiting profile `log n_T=w+o(w)`, write

```text
c_m(r)=sum_(nonempty T) c_T(r).
```

Then the guaranteed divisor in (5) has logarithm `c_m(r)w+o(w)`.
For `m=4` and `r=1,2,3,4`, these exponents are respectively
`17,18,25,48`.

## 3. The critical min-versus-sum calculation

For any fixed `S` of size `r`, summing (4) over all nonempty `T`
gives

```text
sum_T f_T(S)
 = (m-r)2^(m-1)+r(r-1)2^(m-2)
 =: a_m(r).                                            (6)
```

Indeed, each index outside `S` belongs to `2^(m-1)` blocks; each
ordered pair inside `S` belongs to `2^(m-2)` blocks. Consequently

```text
c_m(r)=sum_T min_(|S|=r) f_T(S)
       <= min_(|S|=r) sum_T f_T(S)=a_m(r).             (7)
```

The inequality is **strict for `r<m`**. For every singleton block
`T={i}`, its minimum in (5) is zero, attained by selecting `i in S`.
No set `S` of size `r<m` selects every singleton. Hence for every
such `S`, at least one term `f_{ {i} }(S)` exceeds its individual
minimum, and `c_m(r)<a_m(r)`. At `r=m`, there is only one possible
`S`, so equality holds. In that case (3) is simply the squared
Vandermonde product.

Now impose the genuine balanced small-imaginary profile, with

```text
N(P_i)=Norm(K_i) product_(T containing i) n_T,
log n_T=w+o(w),
log max(1,|y_i|,|t_ij|)=o(w),
log Norm(K_i)=o(w).
```

Taking `L=lcm(|y_i|,|t_ij|)` gives `log L=o(w)`. The actual row
and pair formulas imply

```text
log C_i=2^(m-1)w+o(w),
log|X_i-X_j|=2^(m-2)w+o(w).                                (8)
```

The first follows from `C_i=(L/y_i)^2 N(P_i)`; the second from (2),
since `G_ij` contains `2^(m-2)` balanced blocks and its remaining
factors have subpower height. Every positive summand in (3), and
therefore `D_r`, has logarithm

```text
log D_r=a_m(r)w+o(w).                                  (9)
```

For `m=4` the size exponents `a_4(r)` are `24,24,32,48`. The
guaranteed core divisor exponents `17,18,25,48` leave positive
linear slack in each proper minor and meet the size only at the top
Vandermonde minor. Thus this positive matrix records one-orientation
conductor clusters, but the prescribed Boolean divisibilities alone
do not produce a fixed positive residual-height gap. Any improvement
would have to prove additional divisibility or a smaller archimedean
size using the actual residues and joint norm identities.

The [exact checker](check_positive_hankel_full_profile_budget.py)
verifies the min-versus-sum calculation for `2<=m<=9` and all four
positive minors on a primitive fifteen-block Gaussian fixture. In
that fixture each correction factor overlaps its incident singleton
core block. These finite checks supplement the general identities.

## 4. Leading residue at a singleton core prime

There is also a sharp local test for possible *extra* divisibility.
Let `p` divide only the singleton core block `n_{ {i} }`, suppose
`p` is odd, `p` does not divide `L`, and suppose `p` divides none of
the other `C_j`. Write `alpha=X_i mod p`; then
`alpha^2=-L^2 mod p` and all `X_j+alpha` for `j!=i` are units.
Set

```text
rho_j=(X_j-alpha)/(X_j+alpha) mod p,     j!=i,
V_U=product_(j<k in U)(X_j-X_k).
```

In (3), every summand with `i notin S` contains `C_i` and vanishes
modulo `p`. For `S={i} union U`, the cross factors simplify by

```text
(X_i-X_j)^2/C_j = rho_j mod p.
```

Consequently the *exact* leading-residue formula is

```text
D_r = (product_(j!=i) C_j)
      sum_(U subset [m]\{i}, |U|=r-1)
          V_U^2 product_(j in U) rho_j       (mod p).   (10)
```

For `r=1` the sum is one, so `D_1` is necessarily a `p`-unit under
these hypotheses. For `r=2` the sum is `sum_(j!=i) rho_j`.
It need not vanish: in the checker's actual primitive fifteen-block
Gaussian fixture, the singleton primes `13,29,61` meet all the
hypotheses and give respective `D_2 mod p` values `12,6,1`.
Thus the Gaussian norm equations and Boolean divisibilities do not
automatically contribute an extra singleton-core factor to these
minors. An endpoint argument could still exploit small residues and
simultaneous constraints across different core primes; this local
calculation does not exclude that possibility.
