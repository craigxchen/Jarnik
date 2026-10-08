# Shrinking arcs with bounded primitive pair residues

Small primitive chord residues and a shrinking angle do not by themselves
bound a reduced rational denominator. The following actual Gaussian
integer family makes that limitation explicit. It has arbitrarily many
points if the arc exponent is allowed to exceed one half. Its three-point
specialization lies in arbitrarily small endpoint arcs.

This is not a counterexample to the uniform endpoint conjecture. For a
fixed number `m>4` of rows its normalized endpoint constant tends to
infinity. The conductor growth is part of the proof, not suppressed.

## 1. A unimodular Gaussian pencil

Let positive integers `r,b` run through

```text
r+b sqrt(2)=(3+2 sqrt(2))^k,       k>=1.
```

Thus `r^2-2b^2=1`, `r` is odd, and `b` is even. Put

```text
a=r+b,       G=a-i b,
d=N(G),     n=a^2-b^2+2ab,
B_j=G+j(1-i)bar(G),       j=0,...,m-1.
```

Direct multiplication gives

```text
(1+i)G^2=n+i,       n^2+1=2d^2,
(1+i)G B_j=n+2jd+i,
Im(B_j bar(B_l))=l-j.                                  (1)
```

In particular `gcd_G(B_0,B_1)=1`. Each `B_j` is coprime to its
conjugate and has odd norm. To see this, the real part of
`n+2jd+i` is odd, so its Gaussian gcd with its conjugate is exactly
`1+i`, up to a unit. Removing the factor `1+i` in (1) leaves
`G B_j` conjugate-coprime, and therefore also `B_j`.

If `g_jl=gcd_G(B_j,B_l)`, (1) implies

```text
N(g_jl) divides l-j,
t_jl=(l-j)/N(g_jl) in {1,...,l-j}.                       (2)
```

Here `t_jl` is exactly the absolute imaginary coordinate of the
conjugate-primitive numerator of the phase ratio between rows `j,l`.
All their primitive pair residues, including their product, are bounded
in terms of `m` alone:

```text
1<=T=product_(j<l)t_jl<=product_(j<l)(l-j).              (3)
```

## 2. Actual circles and their full radius

Let `L` be a Gaussian least common multiple of the `B_j`, and set

```text
z_j=bar(L) B_j/bar(B_j),       R=|L|.                    (4)
```

Every `z_j` is a Gaussian integer and has modulus `R`. They are distinct
by (1). Their common Gaussian gcd is a unit. Indeed, at an oriented
Gaussian prime let the orders of `B_j` and `bar(B_j)` be `u_j,v_j`.
At most one is positive for each row, and the gcd of all rows is a unit.
The orders of `z_j` are
`max_l v_l+u_j-v_j`; their minimum is zero. Apply the same argument at
the conjugate prime. There are no ramified or inert factors in the
conjugate-primitive `B_j`.

Primewise comparison of products, pairwise gcds, and the least common
multiple gives the finite bounds

```text
product_j |B_j| / sqrt(product_(j<l)(l-j))
 <= R <= product_j |B_j|.                               (5)
```

For example, after sorting prime orders, the sum of all but the largest
is at most the sum of all pairwise minima. This proves the lower bound
using (2). It does not presume pairwise coprimality.

As `k` tends to infinity, `n/d` tends to `sqrt(2)`. Equation (1) yields

```text
|B_j|/sqrt(d) -> 1+j sqrt(2),
R asymp_m d^(m/2).                                     (6)
```

The exact angular width, for the continuous arguments supplied by (1), is

```text
Theta=2[arctan(1/n)-arctan(1/(n+2(m-1)d))]
      asymp_m 1/d,       m>=2.                         (7)
```

Consequently

```text
Theta sqrt(R) asymp_m d^(m/4-1).                        (8)
```

For `m=3` this tends to zero; for `m=4` it stays between positive
constants; for every fixed `m>4` it tends to infinity. Thus this family
retains the exact radius distinction required by the endpoint problem.
It does not possess the fair conductor profile extracted from a failure
of uniformity.

The same construction can be kept in one fixed Gaussian-unit class.
Retain `B_0,B_2,...,B_(2m-2)` instead. Since
`2(1-i)bar(G)` is divisible by `(1+i)^3`, these rows have the same
unit in the primary Gaussian-prime convention. Their phase quotients
therefore also have the same unit. The determinants in (1) become
`2(l-j)`, and the bound in (3) is multiplied by
`2^(m(m-1)/2)`. The first two retained rows still have unit gcd:
its norm divides two and is odd. All arguments for primitivity and
the radius carry over. In (6) the constants become `1+2j sqrt(2)`,
and (7) uses `n+4(m-1)d`; the exponents in (8) are unchanged.

## 3. Even the reduced Taylor denominator can grow with bounded T

For these `m` points use the normalized sum kernel

```text
S_jl=sec((theta_j-theta_l)/2),       S_jj=1,
u_jl=sin^2((theta_j-theta_l)/2).
```

Every permutation product of `S` is rational: writing
`S_jl=sqrt(N(B_j)N(B_l))/Re(B_j bar(B_l))`, the square roots
cancel around the permutation. For any fixed integer `h>=0`, subtract
the total-degree-`h` Taylor polynomial in the variables `u_jl` from
`per(S)`, and call the resulting rational number `E_h`.

All coefficients of the series are nonnegative. A transposition
contributes `(1-u_jl)^(-1)`, so the first omitted coefficient is
strictly positive. Every nonzero pair angle in (7) is of order `1/d`.
It follows, with constants depending only on `m,h`, that

```text
E_h>0,       E_h asymp_(m,h) d^(-2(h+1)).                (9)
```

If `E_h=A_h/D_h` in lowest terms, positivity and `A_h>=1` imply

```text
D_h >= 1/E_h >>_(m,h) d^(2(h+1)).                       (10)
```

This concerns the reduced denominator of the final sum, not a termwise
common denominator. Together with (3), it excludes a bound for this
denominator solely in terms of `m,h,T` on arbitrary shrinking arcs.
The same exclusion already holds for three-point endpoint arcs of any
fixed positive constant, including within one Gaussian-unit class by
the preceding even-index selection. It does not exclude a special eight-point
endpoint estimate, nor an estimate using the full conductor distribution.

The construction and integer identities are checked in
[check_shrinking_bounded_residues.py](check_shrinking_bounded_residues.py).
The asymptotic assertions follow from the displayed Pell identities;
finite computations are not being substituted for those proofs.
The denominator-audit agent independently checked the construction,
primitivity, complete radius, and reduced-denominator argument.
