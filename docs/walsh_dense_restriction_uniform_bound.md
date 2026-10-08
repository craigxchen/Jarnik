# Arbitrary arc constants for positive-density Walsh profiles

For fixed `C>0`, a fixed number `k` of additional cut classes, and fixed
`delta>0`, there is an effective bound on the ambient Walsh dimension of
the actual profiles defined below. The bound is independent of the circle
radius and allows arbitrary Gaussian row units. In particular, the Walsh
family with at most a fixed number of deleted rows and extra cuts has a
uniform point bound for every fixed endpoint constant `C`.

The proof uses an elementary dense affine-cube argument and exact
threshold-layer restriction. It does not require a Ramsey theorem or a
uniform version of Roth. The effective arithmetic input is the previously
proved [Hadamard family bound](deleted_hadamard_uniform_family_bound.md).

## 1. The exact family

Let `G=F_2^t`, `N=2^t`, and let `A subset G` index distinct nonzero
Gaussian integers `z_x` of one modulus `R`. Assume `|A|>=delta N`.
For each `v in G`, write `chi_v(x)=(-1)^(v dot x)`. Let `W(A)` be the
**set of distinct nonconstant unoriented cut classes** obtained by
restricting all the `chi_v` to `A`.

The hypothesis is that the exact threshold-layer whole-pattern support
of the source points contains every class in `W(A)` and at most `k`
other classes. Every present class is assigned its actual
conjugate-primitive nonunit Gaussian block. Repeated restrictions of
ambient Walsh characters are one class, not separately charged blocks.
Constant restrictions are not present classes.

This refers to the actual split-prime allocations of the source, after
removing their complete common Gaussian divisor. An arbitrary formal
product of Gaussian blocks with overlapping prime support is not a
substitute for that exact threshold factorization.

Assume all points lie in a circular arc of length at most `C sqrt(R)`.
Set

```text
J=ceil(C/sqrt(2)).                                      (1)
```

Subdivide the containing arc into `J` pieces, each of length at most
`sqrt(2) sqrt(R)`. Assign each subdivision endpoint to just one piece.
This colors `A` with at most `J` colors. One color class has density at
least `delta/J` in the ambient group. The source radius is unchanged.

## 2. An elementary dense affine-flat lemma

If `B subset F_2^t` has density at least `rho>0`, then it contains an
affine subspace of dimension `s>=1` whenever

```text
N >= 2^(s-1) (2/rho)^(2^(s-1)).                        (2)
```

Here is a constructive proof. Start with `L_0={0}`, `B_0=B`.
At step `i`, maintain an `i`-dimensional subspace `L_i` and the
translation-invariant set

```text
B_i=intersection_(u in L_i) (B+u).
```

Write `b_i=|B_i|` and `alpha_i=b_i/N`. The exact identity

```text
sum_(h outside L_i) |B_i intersect (B_i+h)|
   = b_i^2-2^i b_i                                    (3)
```

holds because the sum over all `h` is `b_i^2`, while each of the `2^i`
translations in `L_i` fixes `B_i`. If `b_i>=2^(i+1)`, the right side
is at least `b_i^2/2`. Averaging over the at most `N` admissible
translations gives an `h_i outside L_i` with

```text
|B_i intersect (B_i+h_i)| >= b_i^2/(2N).
```

Take `L_(i+1)=span(L_i,h_i)` and `B_(i+1)=B_i intersect (B_i+h_i)`.
Then

```text
alpha_(i+1)>=alpha_i^2/2,
alpha_i>=2(rho/2)^(2^i).                              (4)
```

For `0<=i<s`, condition (2) implies
`N alpha_i>=2^(i+1)`: the sufficient expressions
`2^i(2/rho)^(2^i)` increase with `i`, and (2) is the last one.
Thus all `s` steps are valid and `B_s` is nonempty. Choose `x in B_s`;
then `x+L_s subset B`. This proves (2), including its weak endpoint.

## 3. Restriction preserves every surviving actual block

Fix a split rational prime `p`, a chosen Gaussian prime `pi_p`, and
the original residual allocations `a_x(p)` on `A`. Its raw threshold
sets

```text
T_(p,l)={x in A:a_x(p)>=l}
```

form a chain. Restrict the rows to any nonempty subset `A' subset A`.
The sets `T_(p,l) intersect A'` remain nested. Two such sets cannot be
nonconstant complements: a nonempty set cannot be contained in its
disjoint complement. Therefore, within one nonconstant restricted
unoriented cut class, all layers at the same rational prime use the
same Gaussian orientation.

It follows that merging original blocks whose restricted cuts agree
up to sign produces a conjugate-primitive block, with no rational
prime cancellation. If even one old layer survives in this class,
the merged block is a nonunit. This remains true when different
original classes share a prime through nested threshold layers.

The constant restricted layers are exactly additional common Gaussian
content. More explicitly, if the old residual width at `p` is `e_p`,
and the restricted allocations have minimum `m_p` and maximum `h_p`,
the additional common factor is

```text
pi_p^m_p bar(pi_p)^(e_p-h_p).
```

The remaining allocations have width `h_p-m_p` and complete common
gcd one. Thus restriction gives the literal factorization

```text
z_x=g' epsilon_x product_B Gamma'_B^((1+s_B(x))/2)
                              bar(Gamma'_B)^((1-s_B(x))/2),
log R^2=log Norm(g')+sum_B log Norm(Gamma'_B).            (5)
```

Orienting a restricted cut differently just conjugates its block.
The literal row units `epsilon_x` and the source radius `R` remain
unchanged. We do not divide the source by `g'` or renormalize its arc.

This argument is specific to actual threshold layers. For contrast,
formal columns `(+,+,-)` and `(+,-,-)` with blocks `pi` and `bar(pi)`
merge, upon retaining the first and third rows, into the rational
factor `p`. These are not the threshold columns of their actual row
allocations. Formal overlapping blocks alone do not have the
noncancellation property proved above.

## 4. An affine flat recovers a complete smaller Walsh core

Let `A'=x_0+U subset A` be an affine subspace of dimension `s`.
Under the identification of its rows with `U`,

```text
chi_v(x_0+y)=chi_v(x_0) chi_(v|U)(y).
```

The restrictions `v|U` run over all linear functionals on `U`.
For every nonzero functional, any extending ambient character is
nonconstant already on `A`, so its class belongs to `W(A)` and is
present. Section 3 proves that its surviving merged block is nonunit.
Consequently the exact support on `A'` consists of a complete
`2^s`-row Walsh core plus at most `k` additional classes. Extras may
become constant or merge with the core or with one another; they
cannot increase the number beyond `k`.

Choose an integer `s>=3` such that

```text
2^s >= max_(0<=r<=k) U(r,0),                           (6)
```

where the explicit integer `U` is defined in equation (16) of
[the effective Hadamard family theorem](deleted_hadamard_uniform_family_bound.md).
A monochromatic affine `s`-flat would now be an actual endpoint cluster
of at least `U(r,0)` rows, for some `r<=k`, in an arc of length at most
`sqrt(2) sqrt(R)`. That theorem excludes it, with arbitrary row units,
all moving prime powers, and the complete common factor in (5).

Combining (1), (2), and (6), no profile in Section 1 exists when

```text
2^t >= 2^(s-1) (2J/delta)^(2^(s-1)).                   (7)
```

All constants are effective and depend only on `C,k,delta`. In
particular, if at most `d` ambient rows are missing, then `N>=2d`
gives `delta=1/2`. Substitution in (7) yields a bound depending only
on `C,k,d`, and hence a uniform bound on the number of actual points.
This is an arbitrary-`C` extension for **Walsh** profiles, not a
claim about every Hadamard matrix.

## 5. Why this does not cover the extracted full-fair profile

Any `M`-row full-fair sign matrix can be written as the Walsh
restrictions to

```text
A={0,e_1,...,e_(M-1)} subset F_2^(M-1).
```

Indeed, `chi_v(0)=1` and `chi_v(e_j)=(-1)^(v_j)`, so the nonzero
characters supply exactly all `2^(M-1)-1` nonconstant unoriented cuts.
But these row labels are affinely independent: four distinct labels
cannot sum to zero, so they contain no affine two-dimensional flat.
Their density is

```text
|A|/|G|=M/2^(M-1) -> 0.
```

Thus an ambient Walsh representation by itself provides no extraction
of a large affine flat. The fixed positive-density hypothesis in (7)
is essential to the proved family theorem. No uniform bound or new
growth estimate for general endpoint arcs follows here.

## Verification

[check_walsh_dense_restriction_uniform_bound.py](check_walsh_dense_restriction_uniform_bound.py)
checks the affine-flat construction and its exact averaging recurrence,
the restricted character support, and literal Gaussian threshold
factorizations with nested layers, orientation changes, arbitrary units,
and enlarged common content. It also checks the sparse full-fair simplex
example. These Gaussian fixtures have uncontrolled arguments; their role
is to verify the exact restriction dictionary, not to exhibit short arcs.
