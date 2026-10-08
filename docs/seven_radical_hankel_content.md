# Simultaneous radical forms, exact content, and the middle determinant

The seven actual Gale radicals give a natural family of simultaneous
small forms. Their determinant contents can be normalized exactly.
The middle determinant then factors into the 21 pairwise quadratic
factors: its full formal sign product is an explicit integer power, not an
additional independent arithmetic constraint. This explains the exact
zero exponent in the refined middle-determinant norm budget.

All notation `q,Q_i,D_i,tau,c_i,sigma_i,X` is as in
[the positive radical note](seven_positive_radical_height.md).

## 1. Simultaneous forms with their denominators retained

For `k=1,...,7`, put `r=k-1` and use the ordinary binary Veronese column

```text
u_r(x,y)=(x^r,x^(r-1)y,...,y^r)^T.
H_k=sum_i sigma_i sqrt(c_i) u_r(v_i)u_r(v_i)^T/Q_i^r,
h_k=det H_k.                                       (1)
```

The entries have shared coefficients from the actual marked vectors;
they are not arbitrary radical linear forms. With `Q(t)=q(1,t)`, entry
`a,b`, indexed from zero, is exactly

```text
(H_k)_(a,b)=1/(X sqrt(tau))
  [t_1,...,t_7](t^(a+b) Q(t)^(7/2-k)).              (2)
```

Thus all entries are `O(exp(-77w+o(w)))` on the stated full fair profile.
Their rational coefficients have denominators involving the `Q_i`;
small archimedean coefficient size does not make the entries integral.

To obtain a nonzero determinant asymptotic, normalize the real metric
and rotate the cluster to parameter zero. This is used only for the
archimedean estimate; its matrix and inverse have subpower size and it
changes determinants by subpower factors. The divided differences in
(2) converge to the matrix

```text
(C_k)_(a,b)=[t^(6-a-b)](1+t^2)^(7/2-k),             (3)
```

with negative coefficient indices interpreted as zero. Each `C_k` has
nonzero rational determinant. Consequently

```text
log|h_k|=-77k w+o(w).                               (4)
```

This gives several simultaneously small determinants, with exact shared
coefficients and no numerical-coefficient relaxation.

There are also genuinely independent linear forms here. For fixed
`r=k-1`, the distinct entries are the `2r+1` forms
`L_j=sum_i sigma_i sqrt(c_i)t_i^j/Q(t_i)^r`. Their coefficient matrix
has rank `min(2r+1,7)`, by the distinct-node Vandermonde. The families
for `r=0,1,2` lie in the five-dimensional family with denominator
`Q(t)^2`; for `r=3` the seven forms span all seven source radicals.

Their simultaneous denominator cost is substantial. For primitive
marked vectors, the rational gcd of the squares of **all individual
terms in all these forms** is exactly

```text
e_r=gcd_Q(c_i/Q_i^(2r): i=1,...,7).                 (3a)
```

To see this, at every prime at least one coordinate of each primitive
`v_i` is a unit, and the family includes both pure-coordinate monomials.
Dividing every form by `sqrt(e_r)` makes it an algebraic integer. For
`r=0,1,2,3`, the full fair profile gives respectively

```text
log sqrt(e_r)/w = 0, -17.5, -63, -122.5,
log|L_j/sqrt(e_r)|/w <= -77, -59.5, -14, 45.5.      (3b)
```

All errors are `o(w)`. In particular the first seven independent forms
are not a family of seven small algebraic integers after this exact
termwise normalization.

The archimedean conditioning also balances exactly. For `r=3`, let
`ell_i(t)=P(t)/((t-t_i)P'(t_i))=sum_j ell_(i,j)t^j` be the Lagrange
basis. The exact inverse is

```text
sigma_i sqrt(c_i)=Q(t_i)^3 sum_(j=0)^6 ell_(i,j)L_j. (3c)
```

Its last coefficient is `Q(t_i)^3/P'(t_i)`, of size
`exp(96w+o(w))` because the six parameter gaps have exponent `-16w`.
Multiplying this by the `exp(-77w)` form scale recovers the original
`exp(19w)` radical scale. Treating the seven nearly singular moment
forms as a uniformly conditioned change of basis would lose exactly
this necessary cost.

## 2. Exact algebraic-integer normalization

For a `k`-subset `I`, let `V_I=product_(i<j in I)d_ij`, where
`d_ij=det(v_i,v_j)`. Cauchy--Binet gives

```text
h_k=sum_(|I|=k) w_I,
w_I=(product_(i in I)sigma_i sqrt(c_i))
       V_I^2 / product_(i in I)Q_i^(k-1).
g_k=gcd_Q(w_I^2: |I|=k),
Z_k=h_k/sqrt(g_k).                                 (5)
```

Every `w_I^2/g_k` is a positive integer. Thus `Z_k` is an algebraic
integer, as is every formal choice of signs of its seven source
radicals. No claim that an arbitrary entry of `H_k` was already integral
is needed.

At a clean minority cut of size `s<=3`, a term indexed by `I` with
`a=|I intersect S|` has valuation

```text
v_p(w_I)=a(a+7/2-s-k)e.                             (6)
```

The allowed range is `max(0,k-7+s)<=a<=min(k,s)`. Minimizing twice (6)
is exactly the valuation of `g_k`. Fixed finite minima preserve the
full-profile error budget, including unselected primes. The resulting
table uses `L_k=log sqrt(g_k)/w`:

| k | L_k | log absolute value of Z_k / w | Uniform other-branch bound / w |
|---:|---:|---:|---:|
| 1 | 0 | -77 | 19 |
| 2 | -17.5 | -136.5 | 23.5 |
| 3 | -63 | -168 | 24 |
| 4 | -140 | -168 | 24 |
| 5 | -248.5 | -136.5 | 23.5 |
| 6 | -385 | -77 | 19 |
| 7 | -539 | 0 | 0 |

Every entry in this asymptotic table allows `o(w)` error. The uniform
branch bound before subtracting `L_k` is `35k-16k^2`: each term has
`k` radicals of size `exp(19w)` and a squared normalized Vandermonde
with angular exponent `-16k(k-1)w`.

The stronger rank-of-sign-flips accounting is given in
[the conjugate-budget audit](seven_radical_determinant_conjugate_budget.md).
For `k=3` it gives normalized exponents `-168,-72,-8,24` for zero,
one, two, three flips. Their multiplicities `1,7,21,35` sum to an
exact zero total exponent; a strict saving is not obtained by these
simultaneous forms.

## 3. The apparent symmetry is exact complement duality

Set `A=product_i Q_i`. Multiplying the factors of `D_i` inside `I`
cancels its internal squared Vandermonde. Since `k(7-k)` is even, the
two oriented cross products for `I` and its complement agree. This gives

```text
w_I(k)=-(tau A)^((7-2k)/2) w_(I^c)(7-k).            (7)
```

Taking rational gcds of squares and summing proves

```text
Z_k=-Z_(7-k) for k=1,...,6,       Z_7=-1.            (8)
```

Thus these determinant combinations reduce to three quantities rather
than six independent small algebraic integers.

## 4. The middle determinant is explicitly a conic determinant

First take `Q(t)=1+t^2`, choose either branch `y_i^2=1+t_i^2` at each
node, and let

```text
M_(a,b)=sum_i y_i t_i^(a+b)/P'(t_i),   0<=a,b<=2.
```

Let `E` be the seven-by-seven evaluation matrix with columns
`1,t,t^2,t^3,y,ty,t^2y`, and let `V(t)=product_(i<j)(t_j-t_i)`.
Laplace expansion and Cauchy--Binet give `det M=-det E/V(t)`.

Put `u_i=y_i+t_i`. Then
`t_i=(u_i-u_i^(-1))/2` and `y_i=(u_i+u_i^(-1))/2`.
The seven functions span Laurent powers `-3,...,3`; their exact
coefficient determinant is `1/512`. Therefore

```text
det M = -2^12 product_i u_i^3 /
                   product_(i<j)(1+u_i u_j),
(det M)^2 = 8 / product_(i<j)(1+t_i t_j+y_i y_j).    (9)
```

These determinants are nonzero for every branch choice when the nodes
are distinct: `u_i` is nonzero, and `1+u_i u_j=0` would imply
`t_i=t_j`. No generic nonvanishing assumption remains in the middle
case.

Flipping a branch sends `u_i` to `-1/u_i`. For each pair, the product
of its four denominator factors is `-4(t_i-t_j)^2`. Hence the product
of `det M` over all 128 formal branches is exactly

```text
2^192/V(t)^64.                                     (10)
```

For `Q(t)=A+2Bt+Ct^2`, the corresponding formula is
`2^192 det(q)^96/V(t)^64`. This follows by completing the square and
an affine change of parameter, retaining both the moment-matrix and
Vandermonde scaling factors.

## 5. Exact norm base and pairwise factorization

Write `r=sqrt(det(q))>0`, `D=product_(i<j)|d_ij|`, and

```text
B_0=8r^3/(tau^3 D g_3).                            (11)
```

Since `H_3=(X sqrt(tau))^(-1)M` and `D=X^6|V(t)|`, (10) gives the
exact universal identity

```text
product_(64 sign classes modulo simultaneous reversal) Z_3^2
  = B_0^64.                                        (12)
```

The left side uses formal branches; if square classes are dependent it
must not automatically be identified with the minimal-field norm.
Every factor is nonetheless an algebraic integer. The rational positive
number `B_0` therefore has integral 64th power and is itself an integer.

There are two useful expressions for it. With `I` ranging over triples,

```text
G=gcd_Q((product_(i in I)Q_i) V_I^2 V_(I^c)^2),
g_3=G/(tau^3 D^2),
B_0=8r^3 D/G.                                      (13)
```

In a saturated Gale basis with primitive minors `p_ij`, normalized metric
values `c_i`, and square-star content `s`, the equivalent formulas are

```text
g_3=gcd_Q((p_ij p_ik p_jk)^4/(c_i c_j c_k)^3),
B_0=8r^3 s^5/(product_i c_i * g_3).                 (14)
```

Here `r` in (14) is the determinant square root of that basis's
normalized metric; it must not be substituted from another frame.

Finally define the positive pairwise quadratic factors

```text
U_ij=(q(v_i,v_j)+sqrt(Q_i Q_j))/(r |d_ij|).
```

The bilinear quadratic form satisfies
`q(v_i,v_j)^2-Q_iQ_j=-r^2 d_ij^2`. Thus the formal quadratic conjugate
of `U_ij` is `-U_ij^(-1)`. The covariant version of (9) gives the exact
identity

```text
h_3^2=8r^24/(tau^3 product_(i<j)(q(v_i,v_j)+sqrt(Q_iQ_j))).
```

The power `r^24` is essential for metric scaling. Substituting the
displayed definitions of `B_0` and `U_ij` gives the factorization

```text
Z_3^2 = B_0 / product_(i<j) U_ij.                   (15)
```

The `U_ij` are norm-minus-one factors when their quadratic field is
nontrivial. They need not be algebraic integers: their denominators
remain in the displayed `r|d_ij|`. For example `(1+sqrt(10))/3` has
norm minus one but is not integral. Equation (15) does not silently
promote them to global units.

On the full fair profile each `U_ij` has logarithm `16w+o(w)`, while
`log B_0=o(w)` by (11) and the content table. This recovers
`log Z_3^2=-336w+o(w)` and its quadratic-character conjugate profile.
Thus the middle determinant's norm is already an exact product of
pairwise factors; simultaneous smallness has not supplied an extra
independent height inequality.

## 6. The integer base is not uniformly bounded

Take `q=x^2+y^2` and `v_i=(1,m i)` for `i=0,...,6`. If a prime
`p>7` divides `m` to depth `e`, then every `Q_i` is a `p`-unit,
every bracket has valuation `e`, and

```text
v_p(D)=21e,      v_p(G)=18e,      v_p(B_0)=3e.       (16)
```

Thus `B_0` is unbounded; arbitrary sufficiently large split or inert
primes can occur. This family has uncontrolled collision mass and is
not a fair endpoint family. It shows why replacing the subpower error
bound on `B_0` by an absolute bound would be incorrect.

[check_seven_radical_hankel_content.py](check_seven_radical_hankel_content.py)
verifies the content table, all coalescent determinants, Cauchy--Binet,
the simultaneous linear-form content and exact Lagrange inverse,
integral normalization and complement duality, all 128 branches of a
rational conic example, both metric formulas for the integer base, the
pairwise factorization, and the unbounded-base valuation example.
