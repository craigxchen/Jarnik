# Uniform reciprocal-square tail tightness on critical real-rooted supports

For every `epsilon>0` there is an integer `K(epsilon)`, independent of
the degree, with the following property. Let `m=2s+1`, `s>=1`, and let
`z_1,...,z_m` be nonzero real numbers with distinct absolute values and

```
sum_j z_j^(2r+1)=0,       0<=r<s.
```

Order `a_j=|z_j|` increasingly. Then

```
sum_(j>K(epsilon)) a_j^(-2) <= epsilon a_1^(-2).       (1)
```

The statement is uniform over all such polynomials and degrees. It is
stronger than the bound on the total reciprocal-square mass. The proof
here is elementary compactness for finite root products; no classification
theorem for entire functions is needed. This compactness proof supplies no
explicit rate. The subsequent [harmonic argument](critical_moment_effective_reciprocal_tail.md)
proves the stronger explicit `O(1/K)` tail and linear magnitude growth
under these exact critical-support hypotheses. The finite approximate
coefficient version in Section 5 below remains non-effective.
There is no exclusion of the full eight-row system or transfer to a
general lattice-circle bound.

## 1. Normalized root products

Scale the roots so that `a_1=1`, and put `y_j=1/z_j`, ordered by decreasing
absolute value. Extend this list by zeros after its finite length. The
root polynomial is odd plus a nonzero constant, so its normalized form

```
E(t)=product_j(1-y_j t)
```

satisfies the exact identity

```
E(t)+E(-t)=2.                                        (2)
```

Write `p_l=sum_j y_j^l`. The quadratic coefficient in (2) gives
`p_1^2=p_2`. The [critical paired sign order](critical_odd_moment_root_order.md)
and decreasing reciprocal magnitudes give

```
|p_1|<2,       p_2<4,       |y_j|<=2/sqrt(j).          (3)
```

The total-mass argument is recorded in the
[reciprocal concentration note](critical_moment_reciprocal_concentration.md).
The weak versions `|p_1|<=2`, `p_2<=4` suffice below.

## 2. Every possible escaped mass has an explicit Gaussian factor

Take any sequence of these normalized products `E_n`. By successive
subsequence selection, assume

```
y_(j,n) -> y_j for every fixed j,
p_(1,n) -> ell,       p_(2,n) -> P.
```

The bounds in (3) give `sum_j y_j^2<=P<=4`. Define

```
a=(P-sum_j y_j^2)/2 >=0,
E_1(u)=(1-u)exp(u).
```

Then the products converge uniformly on every bounded complex disk to

```
E_infinity(t)=exp(-ell t-a t^2) product_(j>=1) E_1(y_j t).  (4)
```

Here are the estimates justifying this assertion, including the tail.
For `|u|<=1/2`, the power-series logarithm satisfies

```
log E_1(u)=-u^2/2+R(u),       |R(u)|<=(2/3)|u|^3.     (5)
```

On `|t|<=L`, choose `J` so that `2L/sqrt(J+1)<=1/2`. Uniformly in `n`,

```
sum_(j>J) log E_1(y_(j,n)t)
 =-(t^2/2)[p_(2,n)-sum_(j<=J)y_(j,n)^2]+error,
|error| <= (16/3)L^3/sqrt(J+1).                       (6)
```

This follows from (3), since the sum of the tail squares is at most
four. The same estimate applies to the limiting `y_j`. The finite first
`J` factors converge, and their tail logarithms are uniformly bounded
on the disk. First passing `n` to infinity and then `J` to infinity
proves (4). In particular the infinite product in (4) converges locally
uniformly. This argument remains valid if limiting roots coalesce or
some limiting reciprocals vanish.

Taking limits in (2) gives

```
E_infinity(t)+E_infinity(-t)=2.                       (7)
```

## 3. The constant even part forbids any escaped mass

For every real `u`,

```
|E_1(u)|<=exp(u^2).                                  (8)
```

For `u<1`, this follows from `log(1-u)<=-u`. For `u>1`, use
`log(u-1)<=u-2`, followed by `2u-2<=u^2`. The case `u=1` is immediate.

Put `G(x)=exp(-ell x) product_j E_1(y_j x)` for real `x`. Given any
`eta>0`, choose `J` with `sum_(j>J)y_j^2<eta`. Equation (8) bounds the
tail product by `exp(eta x^2)`. Each of the finitely many initial factors
has modulus at most `(1+|y_j x|)exp(|y_j x|)`. Consequently there is
a finite constant `C_J` such that

```
|G(x)| <= exp(eta x^2+C_J|x|)(1+|x|)^J.              (9)
```

If `a>0`, choose `eta<a`, for example `eta=a/2`. Equations (4) and
(9) then imply `E_infinity(x)->0` as `x->+infinity` and as
`x->-infinity`. This contradicts (7). Therefore

```
a=0,       P=sum_j y_j^2.                            (10)
```

Thus reciprocal-square mass cannot disappear in a normalized limiting
sequence. This proof passes the whole polynomial identity (2) to the
limit. Section 5 shows that for each fixed tail tolerance, a suitable
finite collection of approximate coefficient identities already suffices;
the collection and its error threshold may depend on that tolerance.

## 4. Why absence of mass loss gives a uniform index tail

If (1) failed, there would be an `epsilon>0`, integers `K_n->infinity`,
and normalized products with
`sum_(j>K_n)y_(j,n)^2>epsilon`. Extract the subsequence of Section 2.
For each fixed `J`, eventually `K_n>=J`, and therefore

```
sum_(j<=J)y_j^2 <= P-epsilon.
```

Letting `J` increase gives `sum_j y_j^2<=P-epsilon`, in contradiction
with (10). This proves (1), and scaling restores arbitrary `a_1`.

For each real `r>=2`, the same result also gives
`sum_(j>K) a_j^(-r)<=epsilon a_1^(-r)`, because `a_j>=a_1`.

There is also uniform improvement on the elementary square-root growth
of magnitudes. If `N(T)=#{j:a_j<=T a_1}`, (1) gives

```
N(T)<=K(epsilon)+epsilon T^2,
sup_(critical supports) N(T)/T^2 ->0 as T->infinity.  (11)
```

Equivalently, among supports with at least `j` roots,
`a_j/(a_1 sqrt(j))->infinity` uniformly as `j->infinity`.
Indeed the tail after `floor(j/2)` contains at least `j/2` entries
whose reciprocal squares are at least `a_j^(-2)`, so

```
(j/2)a_j^(-2) <= sum_(h>floor(j/2)) a_h^(-2).
```

This concerns the index growth of magnitudes within a critical support;
it is not an estimate in the original circle radius or a point-count
bound for lattice arcs.

## 5. Finite approximate coefficient conditions suffice at each tolerance

For every `epsilon>0` there exist positive integers `K,J` and a real
`delta>0` with the following property. Let `y_1,...,y_m` be any finite
list of real numbers, sorted by decreasing absolute value, such that

```
max_j |y_j|<=1,       |sum_j y_j|<=2,       sum_j y_j^2<=4.
```

Write `e_l(y)` for its elementary symmetric functions, taking `e_l=0`
when `l>m`. If

```
|e_(2r)(y)|<=delta       for 1<=r<=J,
```

then

```
sum_(j>K) y_j^2<=epsilon.                            (12)
```

This version assumes the displayed bounds. It requires neither the
paired root sign order nor exact forward odd moments; the entries may
repeat or vanish. The parameters `K,J,delta` are uniform over the finite
list length, but the proof supplies no effective values for them.

To prove it, suppose some `epsilon>0` admits no such triple of parameters.
For each integer `n>=1`, choose a violating list with `K=J=n` and
`delta=1/n`. Its tail after `n` has square mass greater than `epsilon`.
Apply the subsequence and canonical-product construction of Section 2;
that construction uses only the displayed bounds, which give
`|y_(j,n)|<=2/sqrt(j)` as before.

The finite product `E_n(t)=product_j(1-y_(j,n)t)` has coefficient
`e_(2r)(y_n)` at `t^(2r)`. Local uniform convergence on complex disks
implies convergence of each Taylor coefficient, for example by Cauchy's
coefficient formula on a fixed circle. For each fixed `r>=1`, the bound
`|e_(2r)(y_n)|<=1/n` eventually applies, so the limiting coefficient is
zero. The constant term is one. Consequently the limiting entire
function again satisfies `E_infinity(t)+E_infinity(-t)=2`.

Section 3 forbids positive Gaussian mass. On the other hand the violated
tail condition, by the fixed-prefix argument of Section 4, forces
`a>=epsilon/2>0`. This contradiction proves (12).

The dependence on the number of coefficient conditions is essential.
For every fixed `J`, [exact Gaussian-cloud constructions](critical_moment_finite_coefficient_gaussian_cloud.md)
match all coefficients through degree `2J+1` while retaining positive
reciprocal-square mass beyond every fixed index. Even distinct absolute
values do not prevent this phenomenon. Those examples do not satisfy
the full critical moment hypotheses.

## 6. Consequence for the nonconstant columns of an eight-row system

Suppose every cross pair between two groups of four rows satisfies the
full critical moments on its differing support. Discard columns constant
on all eight rows; every remaining column belongs to at least one of
the sixteen cross supports. Let `a_*` be the smallest remaining magnitude.

In fact each such column belongs to at least four cross supports. If its
positive-sign counts in the two groups are `r,t`, the number is
`r(4-t)+(4-r)t`, whose minimum over nonconstant columns is four.
Summing (3) over the cross supports consequently gives the useful bound

```
sum_(active h) a_h^(-2) < 16 a_*^(-2).                (13)
```

On each support retain its first `K(epsilon/16)` magnitudes. Their union
has at most `16K(epsilon/16)` columns, and the total reciprocal-square
mass outside that union is at most `epsilon a_*^(-2)`: sum the sixteen
tail estimates and use that each support's minimum is at least `a_*`.
Retaining instead the globally smallest `16K(epsilon/16)` magnitudes
can only decrease the omitted mass. Hence the nonconstant columns of
the full system have a uniform reciprocal-square tail as well.

Constant columns are explicitly excluded because cross-pair differences
do not see them. The [effective follow-up](critical_moment_effective_reciprocal_tail.md)
also gives a uniform `O(1/K)` tail for these active columns. Neither
result is a finite decision procedure: no robust exclusion of all
finite prefixes has been proved.

The degree itself cannot be bounded under the one-support hypotheses.
For every odd `m>=3`, the polynomial `T_m(X)-1/2`, where
`T_m(cos(theta))=cos(m theta)`, has `m` simple real roots in `(-1,1)`.
It is odd plus a nonzero constant, so Newton's identities give the
required odd moments through degree `m-2`. No root is zero, and no two
roots are opposite because `T_m` is odd. Its roots therefore have
distinct absolute values. These unbounded-degree examples explain why
tail tightness controls mass, rather than the total number of roots.

## Verification scope

The [exact checker](check_critical_moment_uniform_reciprocal_tail.py)
tests the finite polynomial identity and reciprocal bounds on rational
critical fixtures, and the rational constants in the logarithmic tail
estimate. These checks do not certify compactness by finite enumeration;
Sections 2--4 give that argument directly.
