# Independent audit of the quartet matching gcd and balanced-cut budget

This note checks the local global-gcd formula and the balanced-cut exponent in
`quartet_matching_gcd_cut_budget.md`.  It also separates three quantities that
must not be conflated: the Gaussian gcd of all quartet contents, the integer
gcd of their norms, and the minimum norm of one quartet content.

## 1. Local valuation of the global Gaussian gcd

Fix a Gaussian prime `pi`, let

```text
r=min_(i<j) v_pi(z_i-z_j),
```

and partition the indices by the residues of
`(z_i-z_1)/pi^r` modulo `pi`.  A product of two disjoint differences has
valuation at least `2r`.

If the largest residue class has size at most `m-2`, the complete multipartite
graph between distinct classes has a two-edge matching.  Both selected
differences have valuation exactly `r`, so the lower bound is attained:

```text
v_pi(G)=2r.                                           (1)
```

If one class `B` has size `m-1`, it is unique.  Put

```text
r_B=min_(i<j in B)v_pi(z_i-z_j)>r.
```

Every disjoint two-edge matching either contains an outlier-to-`B` edge and a
within-`B` edge, or consists of two within-`B` edges.  Its valuation is at
least `r+r_B`.  A pair in `B` attaining `r_B`, together with an edge from the
outlier to a third member of `B`, attains the bound.  The third member exists
because `m>=4`.  Hence

```text
v_pi(G)=r+r_B.                                        (2)
```

This verifies both cases and their exact attainment.  The quantity `r_B` is
a minimum within the unique majority class; it is not a global second-smallest
pair valuation unless that description happens to agree in a special case.

## 2. Balanced-cut exponent

Let `m=2n>=6` and consider one balanced cut of size `n`.  For a quartet meeting
one side in `j` labels, its norm contribution is `|j-2|`.  Summing over all
quartets gives

```text
J_m=4 binom(n,4)+2n binom(n,3)
   =n(n-1)^2(n-2)/2.                                  (3)
```

Since

```text
binom(2n,4)=n(n-1)(2n-1)(2n-3)/6,
```

the normalized contribution is

```text
beta_m=J_m/binom(2n,4)
      =3(n-1)(n-2)/((2n-1)(2n-3)).                   (4)
```

It tends to `3/4`; the checks `beta_6=2/5`, `beta_8=18/35`, and
`beta_10=4/7` agree.

If all `binom(m,n)` balanced cuts have equal logarithmic norm weight, symmetry
makes the total contribution the same for every quartet.  Dividing that
quartet contribution by the total source norm gives exactly (4), so

```text
Gamma_Q >= N^(beta_m-o(1))                            (5)
```

simultaneously for every quartet.

## 3. Why the common gcd can nevertheless lose every source prime

Fix one balanced cut.  Because `n>=3`, there is a quartet with two labels on
each side.  Pairing across the cut twice gives valuation zero at `pi`; the
same is true at its conjugate.  Therefore that source prime divides neither
the global Gaussian gcd `G` nor the integer gcd of all `Gamma_Q`.

This argument may use a different witnessing quartet for each cut.  It does
not produce one quartet small at every source prime.  Formula (5) shows that
the fair profile can instead make every individual quartet large by retaining
different subsets of the source support.

For arbitrary Gaussian tuples the unconditional logical relations are only

```text
Norm(G) divides gcd_Q Gamma_Q <= min_Q Gamma_Q.       (6)
```

The first divisibility need not be equality in general: a rational prime can
divide every norm `Gamma_Q` while the contributing Gaussian orientation
alternates between `pi` and `bar(pi)` from one quartet to another.  Likewise,
a gcd across quartets is not the minimum of their norms.

For same-circle Gaussian tuples the first divisibility is equality; Section 6
proves this using the common-norm valuation levels.

For the balanced source-cut construction, both the Gaussian global gcd and the
integer norm gcd have zero valuation at every deliberately introduced source
prime.  Accidental congruence primes outside the source norm are outside this
claim.

## 4. Audit of the explicit realization

For

```text
kappa_j=2jMq+i,       n_j=Norm(kappa_j),
```

with `M` divisible by every nonzero `j^2-k^2`, the integers `n_j` are pairwise
coprime.  A common prime divisor of `n_j,n_k` would divide their difference
`4M^2q^2(j^2-k^2)` but cannot divide `2Mq`; divisibility of `M` then gives the
contradiction.  Also `kappa_j` and `bar(kappa_j)` are Gaussian-coprime because
their common divisor would divide `2i`, whereas `n_j` is odd.

Allocating `kappa_j` or its conjugate according to the balanced cut therefore
creates equal-norm rows with disjoint source supports and Gaussian gcd one.
Every pair of labels is separated by some balanced cut, so unique prime
support makes the rows distinct.  The construction realizes (5) exactly at
the valuation-profile level. Its angular width is comparable to `1/q`:
the root audit in the main note uses the unique largest power-of-two
denominator to show the leading angle coefficients are nonzero, and their
sum is zero. Since `R` has order `q^K`, `K=binom(m,n)`, its endpoint
constant grows as `q^(K/2-1)`. Thus it gives no endpoint counterexample.

## 5. Relation to triangle content

The exact identity

```text
Gamma_Q=(product_(I subset Q, |I|=3)|D_I|)/s_Q^2
```

shows that quartet content redistributes triangle content and the pencil
factor `s_Q`.  At a residue-free two-level cut, a nontrivial cut has zero
valuation in the gcd of all triangle determinants, while (3) assigns positive
average quartet weight. This compares product weights, not global gcds.
It neither proves nor refutes a bound for the full `Norm(G)` in terms of
the common triangle content; no such comparison is established here.

After removing the common difference factor, (1)--(2) show precisely what an
additional bound would have to control: primes for which all but one reduced
point occupy one residue class.  Source cuts having at least two labels on
each side contribute nothing to `G`; any proposed comparison with triangle
content must separately account for common-difference primes, these unique
`(m-1)+1` concentration primes, and accidental primes outside the source
conductor.  The present audit makes no stronger triangle-content comparison.

## 6. Circle tuples synchronize the conjugate gcd minima

For any same-circle Gaussian tuple,

```text
Norm(G)=gcd_Q Gamma_Q.                                (7)
```

A common Gaussian factor `H` multiplies every matching product, and hence
every `gamma_Q`, by `H^2`.  Both sides of (7) scale by `Norm(H)^2`, so reduce
to a Gaussian-primitive tuple of common norm `N`.

Fix a split rational prime `p=pi bar(pi)`.  If `p` does not divide `N`, equal
norms imply that the `pi`- and `bar(pi)`-valuations of every difference agree.
Thus any quartet attaining one global minimum attains the conjugate minimum.

Suppose `p^e || N`. Put `t_i=v_pi(z_i)` and order the levels as

```text
0=t_1 <= a=t_2 <= ... <= b=t_(m-1) <= t_m=e.         (8)
```

The extrema follow from Gaussian primitivity at `pi` and `bar(pi)`. For two
unequal levels `u<v`, their difference has valuation pair

```text
(v_pi,v_bar(pi))=(u,e-v).                             (9)
```

First assume `a<e` and `b>0`. Every two-edge matching has `pi`-valuation at
least the sum of the two smallest endpoint levels, namely `a`; at the
conjugate prime the analogous lower bound is `e-b`. Choose the quartet with
levels `0,a,b,e`. The matching

```text
(0,b), (a,e)
```

has, by (9), valuation `(a,e-b)` simultaneously. Hence this one quartet
attains both conjugate global minima. Equalities such as `a=b` cause no issue,
because both displayed edges join unequal levels.

If `a=e`, one point is the unique level-zero outlier and the other `m-1`
points have level `e`. Choose a within-majority pair attaining the minimum
residual valuation `r`, a third majority point, and the outlier. Pair the
minimum pair together and the remaining majority point with the outlier.
Equal norms make the within-majority valuation pair `(e+r,r)`, so the quartet
attains both global exponents simultaneously. The exceptional case `b=0` is
the conjugate situation and gives `(r,e+r)`.

These cases exhaust (8). Inert and ramified primes have only one
conjugation-stable orientation, so joint attainment is automatic. Therefore
at every rational prime some quartet realizes the norm of the Gaussian global
gcd, proving (7).

The witnessing quartet can depend on `p`. Equality (7) does not imply that a
single quartet has small norm, and it does not weaken the fair-cut obstruction.

For a concrete distinction, the primitive five-point circle tuple

```text
(-11,-2), (-10,-5), (-5,-10), (-2,-11), (2,-11)
```

has common norm 125. Its five quartet norms, in lexicographic subset
order, are `100,500,40,100,200`. Their common gcd is 20 and their minimum
is 40. Thus even on one actual circle the common gcd need not be attained
as the norm of a single quartet. This finite example has no asserted
endpoint asymptotic significance.
