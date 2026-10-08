# Uniform counts when distinct affine slopes have distinct 2-adic valuations

The common-slope theorem extends to arbitrarily many slope values if their
odd-harmonic frequency sets are disjoint. This permits arbitrary odd parts
of the slopes: only their 2-adic valuations must be distinct. The count is
independent of all slopes, intercepts, layer indices and multiplicities.
It is still a theorem about fixed growing Fibonacci templates, with no
classification or uniform entry threshold for arbitrary lattice circles.

## 1. Statement and hypotheses

Use the exact cyclotomic layers and phase polynomials of
[the mixed-class reduction](mixed_affine_cyclotomic_rank_scope.md). Classes
have rates

```text
d_c(n)=A_c n+B_c,     A_c>0 even, B_c odd,
```

on a fixed parity subsequence. All data are fixed as n grows. Merge
identical classes and retain only active coordinates, with integer widths
`W_(c,e)>=1` and odd indices e. Assume that distinct slope values have
distinct 2-adic valuations:

```text
A_c != A_c'  implies  v_2(A_c) != v_2(A_c').           (1)
```

Several classes may share one slope and have different intercepts.
Condition (1) includes a common slope and also all slopes of the form
`A_0 2^j`, with no bound on the number of j's. It is more general than
the latter condition because the odd parts need not agree.

Assume fixed equal-norm Gaussian prefactors, positive effective degree

```text
L=sum_(c,e) A_c phi(e) W_(c,e)>0,
```

and distinct eventual rows with aligned limiting phases. Alignment forces
common total base parity. After division by the entire Gaussian gcd,
assume the rows lie on an arc of length at most `C sqrt(R_n)` along an
unbounded index subsequence, for fixed `C>0`. A family with at most one
row is trivial. Otherwise let tau be the minimum first nonzero pair
frequency in the normalized phase expansion. Then

```text
L<=4tau.                                             (2)
```

With `D_0=2^32`, the number of distinct eventual rows is at most

```text
2(D_0+1125)+1 = 8,589,936,843.                         (3)
```

No optimization of this large absolute constant is intended. For one
slope, the [common-slope theorem](common_slope_cyclotomic_uniform_count.md)
gives the slightly better bound `2D_0+1`, and at most 2251 rows when
its unweighted degree is at least D_0. The significance of (3) is that
the number and sizes of the slopes in (1) are unrestricted.

## 2. Separating the frequency sets

At frequency t, a class with slope A contributes only if

```text
t=A k,     k positive odd.
```

Consequently every such t has `v_2(t)=v_2(A)`. Under (1), no two distinct
slope values contribute at the same frequency. The low-frequency
equations therefore separate by slope before any arithmetic argument.

For a slope value A, set

```text
h_A = least positive odd integer h with A h>=tau,
D_A = sum_(c:A_c=A) sum_e phi(e) W_(c,e),
Q_A = ceil(L/A).
```

For every actual pair difference, the slope group's coupled moments
vanish at every odd k<h_A. Notice that this is exactly the condition
`A k<tau`, including the case when tau/A is itself an odd integer.
We have

```text
D_A<=L/A<=Q_A<=4h_A.                                 (4)
```

The last inequality follows from (2) and the integrality of `4h_A`.
It does not require tau/A or L/A to be integral.

## 3. Large relative budgets force classwise vanishing

Call a slope group large-budget when `L/A>=D_0`. Apply the arithmetic
decoupling lemma in the common-slope note with the integer upper budget
`Q_A`, rather than requiring the actual degree to equal that budget.
By (4) its hypotheses hold. For every actual pair difference,

```text
S_(c,k)=sum_e b_(c,e)c_e(k)=0
             for every c in the group and every odd k<h_A.  (5)
```

Here c_e(k) denotes the Ramanujan sum. The lemma uses the following
two facts. An integer Laurent polynomial with coefficient mass at most
Q_A cannot vanish at `q^k`, `q=-phi^(-2)`, if `phi^(2k)>Q_A`, unless
it is the zero polynomial: conjugation gives a root of modulus greater
than its coefficient mass. The remaining logarithmically many moments
are recovered by three distinct primes, using
`c_e(dp)=c_e(d) mod p` and the bound `|S_(c,d)|<=Q_A`.

This assertion applies to the bounded integer differences of actual
rows. It does not bound the full rational coupled kernel, which can
have arbitrarily large dimension even in the common-slope case.

Ramanujan inversion turns (5) into the separate truncated Mobius
equations in each class:

```text
B_(c,d)=sum_(e:d|e) mu(e/d)b_(c,e)=0,   d odd, d<h_A.  (6)
```

## 4. Weighted colored compression bounds the sum of the nullities

Let E_c be the active support of class c in the large-budget groups.
Apply [prime-power compression](cyclotomic_uniform_rank.md#1-a-prime-power-compression-preserving-nullity)
separately to its matrix in (6), obtaining a divisor downset F_c.
The proof preserves or increases nullity and does not increase the
ordinary totient cost `sum_e phi(e)`. There is no extra base charge
here: only total base parity is assumed, not parity within each class.
Multiplying the cost in class c by its fixed slope A_c preserves the
inequality. Thus

```text
sum_(large-budget c) A_c sum_(e in F_c) phi(e) <= L.    (7)
```

For a divisor downset its nullity equals the number of high indices,
now specified by `e>=h_(A_c)`, or equivalently `A_c e>=tau`.
We show that the total number of these colored high indices is at most
1125, irrespective of their slopes.

Give class c a separate color. For each root of unity whose exact
order belongs to F_c, include A_c labeled copies. The total cardinality
of this disjoint colored universe is the left side of (7). For any
n in F_c, the complete nth-root group, with these copies, has size
`A_c n`. Two such groups of the same color intersect in
`A_c gcd(n,m)` points; different colors have empty intersection.

Build a graph on the colored high indices. Join two distinct vertices
only when their copied root groups intersect in at least `2tau/15`
points. Neighbors have the same color. Writing
`n=a g, m=b g, g=gcd(n,m)`, their odd positive ratios satisfy

```text
a,b <= L/(2tau/15) <=30.                              (8)
```

There are at most 225 such ratio pairs, including `(1,1)` for the
vertex itself. Each pair determines m from n, so the degree is at most
224. If there were at least 1126 high vertices, greedy selection would
give six independent vertices. Their copied root groups would obey

```text
size(union of the six groups)
 >= sum group sizes - sum pair intersection sizes
 > 6tau - 15(2tau/15) =4tau >= L,                     (9)
```

contradicting (7). Thus the total classwise kernel dimension on all
large-budget groups is at most 1125.

## 5. The other groups cost only a bounded number of coordinates

Every coordinate in a remaining group has `A_c>L/D_0`, and contributes
at least A_c to the total degree, since `phi(e)W_(c,e)>=1`. If q is
the number of all such coordinates, the conclusion `q<D_0` is immediate
when q=0. When q>0,

```text
q L/D_0 < sum_(remaining c,e) A_c phi(e)W_(c,e) <=L,
q<D_0.                                                (10)
```

No separate moment restriction on these coordinates is needed. Combining
(6), the nullity bound, and (10), every actual pair difference lies in
a fixed rational linear space of dimension at most `1125+D_0`.

The global coefficientwise reciprocity from the mixed-class note uses
the automorphism fixing i and sending phi to `-1/phi`. Total base
parity gives a common reciprocal scalar for the two canceled pair
polynomials. It implies

```text
sum_(c,e) A_c phi(e)|b_(c,e)| >=2tau.                  (11)
```

Map rows to the weighted box vectors

```text
f_i(c,e)=sqrt(A_c phi(e)W_(c,e)) (2a_(i,c,e)/W_(c,e)-1).
```

Their pairwise inner products are at most `L-4tau<=0`. Their affine
dimension is at most `1125+D_0` by the preceding linear containment.
The affine-projection/Gram-matrix argument for an obtuse set then gives
at most `2(1125+D_0)+1` rows, proving (3).

One must not replace (11) by an assertion that every individual class
has canceled degree at least twice its contact. A class with odd base
difference can violate that assertion. The global parity and global
reciprocity are essential to the proof used here.

## 6. Scope and the next obstruction

All constants in (3) are independent of the template and the fixed arc
constant C. The threshold at which the fixed template reaches its
asymptotic regime may depend on it. The theorem does not bound arbitrary
isolated evaluations of varying templates, and no reduction of arbitrary
centered-circle configurations to these templates is known here.

If two distinct slopes have the same 2-adic valuation, their odd-harmonic
frequency sets overlap. Then a single frequency combines different
harmonic orders in different classes; the polynomial in q used in
Section 3 no longer has one common power q^k. Multiplying a frequency
by a prime can introduce new contributing slopes as well. Neither
classwise decoupling nor the constant bound has been proved in that
case. The general radius-independent circle bound remains unproved.
