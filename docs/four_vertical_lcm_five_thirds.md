# Four vertical cofactors: a five-thirds bound and the remaining obstruction

For any positive integer `d` and distinct positive integers
`t_1<t_2<t_3<t_4`, put `w_j=-d+i t_j` and let `L` be their Gaussian
lcm. Then

```text
|L| >= d^(-1) t_1^(2/3) t_2^(1/2) t_3^(1/3) t_4^(1/6)
     >= d^(-1) t_1^(5/3).                                  (1)
```

No comparability assumption is needed. This improves the direct
product-over-all-pair-gcd estimate, but does **not** prove the proposed
quadratic bound. It does not use ordinary primitivity of complements,
and it does not address varying contents or the uniform endpoint count.

## 1. A better way to combine the pair gcds

For four arbitrary nonzero elements in a unique factorization domain,
write `g_ij=gcd(w_i,w_j)`. Prime by prime, sort their exponents as
`e_1>=e_2>=e_3>=e_4`. The exponent in the product of the six pair gcds is

```text
e_2+2e_3+3e_4.
```

Consequently

```text
3e_1 + (e_2+2e_3+3e_4) - 2(e_1+e_2+e_3+e_4)
 = e_1-e_2+e_4 >= 0.
```

Thus there is the divisibility, and hence modulus inequality,

```text
(product_j w_j)^2 | L^3 product_(i<j) g_ij,
|L|^3 >= (product_j |w_j|)^2 / product_(i<j) |g_ij|.      (2)
```

Using the exact vertical-line gcd bound from
[the endpoint note](endpoint_vertical_line_gcd.md), (2) implies

```text
|L|^3 >= (product_j (d^2+t_j^2))
          / (d^3 sqrt(product_(i<j) (t_j-t_i))).          (3)
```

Since all heights are positive and ordered,

```text
product_(i<j) (t_j-t_i) <= t_2 t_3^2 t_4^3.
```

Replace `d^2+t_j^2` by `t_j^2` in the numerator of (3) and take cube
roots to obtain (1).

## 2. Exact identification of the high-order gcd correction

Define

```text
H = lcm_G(gcd_G(w_i,w_j,w_k) : i<j<k),
G = gcd_G(w_1,w_2,w_3,w_4),
S = product_j [w_j / gcd_G(w_j,lcm_G(w_k : k != j))].
```

At each Gaussian prime the exponents of `H`, `G`, and `S` are
respectively `e_3`, `e_4`, and `e_1-e_2`. Therefore the following is an
exact identity up to a unit:

```text
L^2 H G = (product_j w_j) S.                             (4)
```

In particular,

```text
|L| >= sqrt(product_j |w_j| / (|H| |G|)).                (5)
```

If no Gaussian prime divides three of the cofactors, then `H=G=1`,
and this immediately gives `|L|>=t_1^2`. More generally, a bound for
`|H||G|` depending only on `d` would imply the desired quadratic bound.
Such a bound has not been established, and (4) shows why discarding
the numerator `S` can also lose information. Large higher gcds alone
do not establish a counterexample.

For reference the complete inclusion-exclusion identity is

```text
L = (product_j w_j)
      (product_(i<j<k) gcd_G(w_i,w_j,w_k))
      / ((product_(i<j) g_ij) G),                        (6)
```

again up to a unit. The four triple gcds must not be omitted when
trying to understand simultaneous cancellation.

## 3. Why incidence and pairwise size bounds alone stop at five-thirds

Here is an exact obstruction in the *relaxed valuation model*, not a
construction of Gaussian integers on the required vertical line.
Give each two-element subset of `{1,2,3,4}` a separate formal prime
block, and likewise each three-element subset. There are ten blocks.
Let `w_j` be the product of the blocks whose subsets contain `j`.
Assign each block logarithmic size `1/6`, in units of `log T`.

Then every cofactor contains three pair blocks and three triple blocks:

```text
log_T |w_j| = 1.
```

Every pair gcd contains its one pair block and the two triple blocks
containing that pair, so

```text
log_T |g_ij| = 1/2.
```

Every triple gcd consists of its one triple block; the fourfold gcd
and all the unique-factor terms in `S` are units. Finally,

```text
log_T |L| = 10/6 = 5/3,
log_T |H| = 4/6 = 2/3.
```

All exact lcm/gcd incidence identities hold, including (2), (4), and
(6), with equality in the exponent bound (2). The pair-gcd exponents
also saturate `|g_ij|=O_d(T^(1/2))` for comparable, separated heights.
Thus no argument using only these incidence identities and the
coarse pairwise size bounds can force a quadratic exponent. It must
also use the additive compatibility of the six differences, the
common real part, or the primitive complements. The formal blocks
have **not** been realized with those constraints; this is not a
counterfamily to the quadratic conjecture.

The existing four-cofactor Pell construction attains exponent two;
it does not settle whether a smaller exponent is realizable. See
[the Pell-family note](vertical_lcm_pell_counterfamily.md).

## 4. The same argument for arbitrary cardinality

Let `k>=2`, `r=floor(k/2)`, `B=r(r+1)/2`, and `K=k(k-1)/2`.
For any integer multiplicity `s`,

```text
r s - s(s-1)/2 <= B,
B - r s + s(s-1)/2 = (s-r)(s-r-1)/2 >= 0.
```

Apply this at every prime-power level, where `s` counts the cofactors
divisible by that power. Summing gives

```text
(product_j w_j)^r | L^B product_(i<j) gcd_G(w_i,w_j).
```

For `0<t_1<...<t_k`, the same vertical-line bound and
`product_(i<j)(t_j-t_i)<=product_j t_j^(j-1)` now give

```text
|L| >= d^(-K/(r(r+1)))
       product_(j=1)^k t_j^((2r-j+1)/(r(r+1)))
     >= d^(-K/(r(r+1))) T^(2-1/(r+1)),       T=t_1.      (7)
```

All individual height exponents are nonnegative, even for odd `k`,
so **comparability is unnecessary** here too. Because `d>=1`, the
weaker but simpler consequence is

```text
|L| >= (T/d)^(2-1/(floor(k/2)+1)).                      (8)
```

For every fixed `k` this exponent is strictly smaller than two. Inserting
the endpoint condition `T/d>=sqrt(R)/C` and `|L|<=2R` therefore gives
no contradiction at large `R`, and proves no uniform endpoint bound.
Letting `k` increase does not turn (8) into a fixed-cardinality endpoint
theorem.

Run `python3 docs/check_four_vertical_lcm_five_thirds.py` for exact
finite checks of (2)--(6), the height bound, and the formal incidence
obstruction. The proofs above establish the infinite statements.

## 5. The same exponent holds with varying chord contents

The exponent in (8) does not actually require equal real parts. Let

```text
w_j=-d_j+i t_j,     d_j,t_j positive integers,
x_j=t_j/d_j,       0<x_1<...<x_k,
L=lcm_G(w_1,...,w_k),       r=floor(k/2).
```

For any two Gaussian integers `u,v` not on the same real line through zero,
write `u=g u'`, `v=g v'`, where `g=gcd_G(u,v)`. Then

```text
det(u,v)=N(g) det(u',v'),
```

and the determinant on the right is a nonzero integer. Hence

```text
N(gcd_G(w_i,w_j)) <= |det(w_i,w_j)|
                         =d_i d_j (x_j-x_i).           (9)
```

Apply the same primewise divisibility from Section 4 and use
`|w_j|>=d_j x_j`. The product of the six, or more generally `k(k-1)/2`,
pair bounds in (9) contains each `d_j` to the power `(k-1)/2` after
taking Gaussian moduli. The resulting inequality is

```text
|L| >= (product_j d_j)^((2r-k+1)/(r(r+1)))
       product_(j=1)^k x_j^((2r-j+1)/(r(r+1)))
     >= X^(2-1/(r+1)),                  X=min_j x_j.   (10)
```

All powers used for the last inequality are nonnegative, and `d_j>=1`.
For even `k=2r` the content exponent is `1/(r(r+1))`; for odd `k=2r+1`
it is zero. Thus varying chord contents do not obstruct this particular
subquadratic bound. For four factors the more precise version is

```text
|L| >= (d_1 d_2 d_3 d_4)^(1/6)
       x_1^(2/3) x_2^(1/2) x_3^(1/3) x_4^(1/6).      (11)
```

Inserting a common `d` into (11) recovers (1) exactly.

For genuine anchored endpoint factors, `L|2P` and `|P|=R`. On a one-sided
arc of length at most `C sqrt(R)`, the earlier elementary endpoint estimate
gives `x_j>=sqrt(R)/C` for sufficiently large `R`. With `k=2r`, (10)
therefore implies the conditional content restriction

```text
product_(j=1)^(2r) d_j
 <= 2^(r(r+1)) C^(r(2r+1)) R^(r/2).                 (12)
```

Its geometric-mean content scale is at most a constant times `R^(1/4)`.
Small contents are possible in actual endpoint configurations, so (12)
does not give a uniform number of points. The remaining exponent deficit
is present even after removing the equal-content assumption.
