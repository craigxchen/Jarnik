# Audit of the CM denominator phase lift

This note records what is actually justified in
`cm_denominator_phase_lift.md`, with particular attention to vertical
divisors and moving correction phases.  It does not prove the uniform
circle short-arc bound.

## 1. The exact local correction

Put

```text
A=(ell-i)P,  Abar=(ell+i)P,  I=iP,  R=NP,
d_v(T)=max(0,-v(x(T))/2).
```

At every odd Gaussian prime the displayed Weierstrass equation has good
reduction.  If `f_ell` is normalized by `f_ell(O)=1`, intersection with
its divisor gives

```text
v(f_ell(R))
 = d_v(R-A)+d_v(R-I)-d_v(R-Abar)-d_v(R+I)
   -d_v(-A)-d_v(-I)+d_v(-Abar)+d_v(I).               (A1)
```

The last line is the base-section correction.  The four moving terms are
the denominator depths of, respectively,

```text
Q_(N-ell),  bar(Q_N),  bar(Q_(N-ell)),  Q_N.
```

The two terms involving `I` vanish, and
`d_v(-A)-d_v(-Abar)=v(gamma_(-ell))`.  Consequently (A1) is exactly

```text
v(gamma_(N-ell)/gamma_N)
 = v(gamma_(-ell) f_ell(NP))                         (A2)
```

at every odd Gaussian prime.  This explains why the fixed factor
`gamma_(-ell)` is necessary; it is not an arbitrary normalization.
Collisions of horizontal sections cause no loss because their complete
intersection multiplicities are the depths in (A1).

Let the quotient of the two sides of (A2) be `q_(N,ell)`.  Every factor
has norm one, so `q_(N,ell) bar(q_(N,ell))=1`.  The curve has bad
reduction only above two, and `1+i` is fixed by conjugation.  Hence

```text
2 v_(1+i)(q_(N,ell))=0.
```

Together with (A2), this says that `q_(N,ell)` has zero valuation at
every finite prime.  Therefore it is one of the four Gaussian units and

```text
gamma_(N-ell)/gamma_N
 = epsilon_(N,ell) gamma_(-ell) f_ell(NP),
epsilon_(N,ell) in {1,-1,i,-i}.                       (A3)
```

Thus no unbounded bad-prime factor is missing in this example.  More
generally, on a fixed proper regular model, the vertical part of the
divisor of a fixed function has finite support and finitely many fixed
component coefficients.  A varying section can change which component
it meets, but cannot make those coefficients unbounded.  For a fixed
finite list of shifts, an infinite pigeonhole subsequence freezes all
such component choices and all unit factors simultaneously.  No
subsequence is needed for (A3).

The exceptional input `N=-ell`, at which the two unreduced line formulas
both vanish at their common third intersection, must be evaluated after
cancelling that common intersection.  It is one isolated value for each
fixed shift and is irrelevant to the asymptotic use.  The checker skips
this value, so it tests (A3) extensively but does not prove the
intersection argument or this cancelled evaluation.

## 2. Bounded trimming and Boolean words

Replacing a denominator generator by the explicit Gaussian gcd changes
it by a unit times `(1+i)^e` with `e` in a fixed finite interval.  On
taking `b/bar(b)`, this contributes only a Gaussian unit.  Dividing by
common factors drawn from fixed denominator ideals likewise gives only
finitely many phase corrections.  Hence bounded trimming preserves the
finite-target conclusion.

The divisor formula for a zero-sum Boolean word is also correct when a
shift is zero.  There is no nonvertical `f_0`; instead the fixed
`[iP]-[-iP]` terms from the nonzero shifts contribute

```text
-c_0([iP]-[-iP])=c_0([-iP]-[iP]),
```

which is precisely the `k=0` term in the displayed divisor in the main
note.  Distinctness of all `(k-i)P` and `(k+i)P` then makes every nonzero
word nonconstant.

Complement reflection does not remove arbitrary row corrections.  If
`kappa_i=K_i/bar(K_i)`, any monomial in row phases has correction factor
`product_i kappa_i^(r_i)`.  Independence from all `kappa_i` forces every
`r_i=0`, in which case the corresponding row-phase monomial is already
the trivial identity.  Products around pairwise cycles cancel both the
corrections and the elliptic functions tautologically.  The phase lift
therefore supplies no new correction-free symmetry obstruction.

## 3. What David's theorem gives for a moving target

For fixed targets, the use of elliptic logarithms in the main note has
the right shape.  The primary source is Sinnou David, *Minorations de
formes lineaires de logarithmes elliptiques*, Memoire SMF 62 (1995),
Theorem 2.1, DOI `10.24033/msmf.376`.  Taking logarithms of `P`, a fixed
zero `Q` of `F-a`, and two periods gives four fixed elliptic logarithms;
the coefficient heights are `O(log N)`.  The theorem therefore yields
the stated

```text
exp(-C log N (log log N)^5)
```

lower bound when the target belongs to a fixed finite set.

The explicit dependencies also show the limit of the same method for a
moving correction.  Let

```text
s_N=max_i log max(2,|K_(N,i)|).
```

The phase target `a_N` has height `O(s_N)`.  Every point in the fixed
finite fiber `F^(-1)(a_N)` has uniformly bounded degree and canonical
height `O(s_N+1)`.  In David's Theorem 2.1 this varying point contributes
one factor `O(s_N+1)` through its `log V` parameter, while the number
field degrees remain bounded.  There is a second dependence which must
not be hidden: condition (1) of that theorem requires `D log B` to
dominate the largest `log V`.  Thus `log B` must be at least a constant
multiple of `s_N`, as well as dominating the coefficient height
`O(log N)`.  The argument gives, away from exact equality, a uniform
bound of the form

```text
|F(NP)-a_N|
 >= exp(-C (s_N+1) max(s_N,log(N+2))
              (log(s_N+log(N+3)))^5).                 (A4)
```

This is useful under the stronger hypothesis

```text
(s_N+1) max(s_N,log N)
    (log(s_N+log N))^5=o(N^2),                         (A5)
```

and then the residue conclusion of the main note still has its full
quadratic leading term.  In the usual range `s_N >> log N`, this is
roughly the requirement `s_N^2(log s_N)^5=o(N^2)`.  The available
hypothesis `s_N=o(N^2)` does not imply (A5), even on a subsequence: for
example `s_N=N^(3/2)` is little-oh of `N^2` but violates (A5) along
every subsequence.  Thus David's verified dependencies do not close the
moving-correction gap.  Exact equality is eventually impossible under
`s_N=o(N^2)`, since functorial height gives
`h(F(NP))=deg(F) N^2 hhat(P)+O(N)`, whereas equality with the correction
phase would bound this by `O(s_N)`; the quantitative near-equality is
unresolved by the single-function estimate.  The collective argument in
the next section supplies the missing information for this fixed family.

The exact checker was rerun successfully.  It verifies 181 finite
instances of (A3), 40 generator comparisons, and the three-row Boolean
identities.  These computations support the formulas but do not replace
the divisor-model proof or provide a moving-target estimate.

## 4. The collective phase map removes the moving targets

First identify the effect of complement symmetry.  Index the blocks by
`b in {0,1}^m`, use row `m` as anchor, and put

```text
c_b=(b_1-b_m,...,b_(m-1)-b_m),
A_b=(k_b-i)P,  B_b=(k_b+i)P.
```

Let `F_1,...,F_(m-1)` be the row-to-anchor phase functions.  Their
simultaneous valuation vector is

```text
ord_(A_b)(F_1,...,F_(m-1))=c_b,
ord_(B_b)(F_1,...,F_(m-1))=-c_b.                      (A6)
```

The all-zero and all-one labels have zero vector.  For every other
label, `c_b` is primitive, and each nonzero vector in (A6) occurs at
exactly two points:

```text
A_b and B_(bar b).
```

All these points are distinct because a nonzero Gaussian endomorphism
cannot annihilate the nontorsion point `P`.

Let `C` be the smooth curve with function field
`Q(i)(F_1,...,F_(m-1))`, and let `pi:E->C` have degree `d`.  Pullback
multiplies a valuation vector by the ramification index.  The vectors in
(A6) are primitive, so the displayed points are unramified, and their
fibers contain at most the two points just listed.  Hence `d<=2`.

If `d=2` and `C` has genus one, the two points in a fiber differ by a
fixed torsion point.  But

```text
B_(bar b)-A_b=(k_(bar b)-k_b+2i)P
```

is nontorsion, a contradiction.  If `C` has genus zero, the deck
involution has the form `Q -> T-Q`.  Pairing `A_b` with `B_(bar b)`
then forces

```text
k_b+k_(bar b)=constant                               (A7)
```

for every nontrivial complementary pair.  Conversely, (A7) makes all
divisors in (A6) invariant under that involution.  Their descended
degree-zero divisors on `P^1` are principal, so all the `F_j` factor
through the degree-two quotient.  Therefore, for `m>=3`,

```text
(F_1,...,F_(m-1)) is birational to E
iff the nontrivial complement sums are not all equal. (A8)
```

Both function-field cases contain a rational degree-two coordinate.  In
the birational case take `u=x` and set `T=O` in the height formula below.
In the complement-symmetric case let
the common sum in (A7) be `C` and put `T=CP`.  If `T!=O`, set

```text
u(Q)=(y(Q)+y(T))/(x(Q)-x(T)).                           (A9)
```

This is the slope of the line through `Q` and `-T`.  The points
`Q,T-Q,-T` are collinear, so `u(Q)=u(T-Q)`.  Its poles are `O` and `T`,
and it is a rational degree-two parameter for the quotient.  If `T=O`,
take `u=x`.  Thus in either case `u` is defined over `Q`, has degree two,
and belongs to `Q(i)(F_1,...,F_(m-1))`.  Along the real orbit,

```text
h(u(NP))=hhat(NP)+hhat(NP-T)+O(1)
        =2 hhat(P) N^2+O(N).                          (A10)
```

Now let the row corrections have `s_N=o(N^2)` and assume the resulting
directions are distinct.  Put

```text
L_N=log max_(j<m)|t_(j,m)|.
```

The exact primitive-pair phase identity gives norm-one targets

```text
a_N=(a_(N,1),...,a_(N,m-1)) in Q(i)^(m-1),
h(a_N)=O(s_N),
max_j |F_j(NP)-a_(N,j)|
 <= exp(-2^(m-1) hhat(P) N^2+L_N+o(N^2)).             (A11)
```

The coefficient is exact at leading order: a pair numerator contains
`2^(m-1)` separated blocks, and each has
`log|H_(N,b)|=hhat(P)N^2+O(N)`.  Finite units and trimming factors are
absorbed in `a_N`.

Suppose, on a subsequence and for some fixed `epsilon>0`, that

```text
L_N <= (2^(m-1)-2)hhat(P)N^2-epsilon N^2.             (A12)
```

Then (A11) is at most
`exp(-(2hhat(P)+epsilon)N^2+o(N^2))`.

Take fixed generators `G` for the ideal of the image curve in the
algebraic torus `(G_m)^(m-1)`.  The norm-one complex locus is compact,
which gives uniform Lipschitz constants for these fixed polynomials.
Equation (A11), under (A12), makes `G(a_N)` exponentially small with a fixed positive
quadratic exponent.  If it were nonzero, the product formula and
`h(a_N)=o(N^2)` would give the incompatible lower bound
`|G(a_N)|>=exp(-o(N^2))`.  Thus every defining relation vanishes at
`a_N` for all large `N`: the moving target lies exactly on the fixed
image curve.

Write `u=R(F_1,...,F_(m-1))=H(F)/D(F)` using fixed Laurent polynomials
`H,D`, with `D` nonzero in the image function field, and put
`z_N=R(a_N) in Q(i)`.  Exclude the finite set of zeros of `D` on the
image, including removable, base, and singular values.  Away from this
set, height functoriality and the product formula give

```text
h(z_N)=o(N^2),
|u(NP)-z_N|
 <= exp(-(2hhat(P)+epsilon)N^2+o(N^2)).                (A13)
```

Indeed `|D(a_N)|>=exp(-o(N^2))`.  The error in (A11) is much smaller,
so also `|D(F(NP))|>=|D(a_N)|/2`; applying `R` loses only an
`exp(o(N^2))` factor.  Repeated exceptional values are fixed targets and
are excluded by applying the fixed-target elliptic-logarithm bound to a
nonconstant coordinate function.

Write the reduced rational number `u(NP)=a/b` with `a,b in Z`, and write
`z_N=c/d` with coprime `c,d in Z[i]`.  If the values are unequal, then
`ad-bc` is a nonzero Gaussian integer and

```text
|u(NP)-z_N| >= 1/(b|d|)
 >= exp(-2 hhat(P)N^2-o(N^2)).                         (A14)
```

Equality is also impossible for large `N`, because (A10) gives
`h(u(NP))=2hhat(P)N^2+O(N)`, whereas `h(z_N)=o(N^2)`.  Equations
(A13)--(A14) contradict each other.  Any violation of the following
bound would supply such an `epsilon` and subsequence, so

```text
liminf_(N->infinity) log max_(i<j)|t_(i,j)| / N^2
 >= (2^(m-1)-2) hhat(P).                              (A15)
```

The maximum over all pairs is at least the anchor maximum used above.
For `m>=3`, (A15) is a positive quadratic residue-growth theorem for
this fixed family.

The checker uses the symmetric binary assignment
`k_b=0,...,2^m-1`, for which `k_b+k_(bar b)=2^m-1`.  Its phase map is
genus zero, but the quotient coordinate (A9) shows that it is not an
exception to the height argument.  Therefore the whole fixed CM
denominator family is excluded for `m>=3`, even with arbitrary moving
corrections of height `o(N^2)`: the extracted endpoint system requires
`L_N=o(N^2)`, contrary to (A15).  This remains an obstruction to this
family, not a proof of the uniform circle short-arc theorem.

## 5. Actual arc growth for this fixed family

There is also a direct conclusion for centered Gaussian lattice circles,
without assuming in advance that all primitive residues have subquadratic
logarithmic height. Write `B=2^m` and `q=hhat(P)>0`. Keep the fixed shifts,
bounded trimming, and conjugate-primitive rows of Section 4, with
`sum_i log|K_(N,i)|=o(N^2)`. Let `delta_N` be the angular width, in radians,
of an arc containing their `m` unit-circle phases
`P_(N,i)/bar(P_(N,i))`.

The same prime-count calculation as in (A11), for every pair, gives

```text
log|h_(N,ij)|=(B/2)qN^2+o(N^2),
|P_(N,i)/bar(P_(N,i))-P_(N,j)/bar(P_(N,j))|
 =2|t_(N,ij)|/|h_(N,ij)|.
```

The width of a containing arc is at least the distance between any two
of its unit-circle points. Applying (A15) to a pair attaining the largest
residue therefore gives

```text
delta_N >= exp(-2qN^2-o(N^2)).                         (A16)
```

To compute the least possible lattice radius, put

```text
L_N=lcm_G(P_(N,1),...,P_(N,m)),
G_N=gcd_G(P_(N,1),...,P_(N,m)).
```

For clarity, the elementary radius formula is

```text
R_min=|L_N|/|G_N|.                                    (A17)
```

Indeed a common rational Gaussian multiplier `beta` makes all
`beta P_i/bar(P_i)` Gaussian integers exactly when, at each Gaussian
prime `pi`,

```text
v_pi(beta) >= max_i(v_pi(bar(P_i))-v_pi(P_i))
           =v_pi(bar(L_N))-v_pi(G_N).
```

The equality uses conjugate-primitivity: if any row has positive
`pi`-valuation, that row has zero `bar(pi)`-valuation, and conversely.
Thus the permitted multipliers form the fractional ideal
`(bar(L_N)/G_N) Z[i]`, whose least nonzero modulus is (A17).
Any lattice realization has such a rational multiplier, since its
ratio to one of the rational Gaussian phases is rational Gaussian.
This is also the exact formula proved in
[the minimal-radius note, Section 1](fresh_algebraic_parallel.md).

Initially omit the corrections `K_i`. Independence of the block supports
means that the lcm contains every nonempty cut block once, and the gcd
contains exactly the all-rows block. Consequently

```text
log R_min^(0)=sum_(empty != b != full) log|H_(N,b)|
            =(B-2)qN^2+O(N).
```

The cost of restoring corrections can be proved prime by prime. For
a fixed Gaussian prime let `e_i` be its uncorrected row valuations and
`f_i=v_pi(K_i)>=0`. Its contribution to `log(|L|/|G|)` is
`(max_i e_i-min_i e_i) log Norm(pi)/2`. Since

```text
|range(e+f)-range(e)| <= range(f) <= sum_i f_i,
```

summing over all Gaussian primes gives the full-depth bound

```text
|log R_min-log R_min^(0)| <= sum_i log|K_i|.
```

No radical replacement or uncharged opposite-orientation correction is
being used. Hence, for the original `m` points,

```text
log R_min=(B-2)qN^2+o(N^2),
delta_N sqrt(R_min)
 >= exp(((B-6)/2)qN^2-o(N^2)).                        (A18)
```

For `m>=3`, `B>=8`, and the right side tends to infinity. In particular
these fixed CM families cannot lie on arcs of length at most
`C sqrt(R)` for any fixed `C`: such an arc would require bounded
`delta_N sqrt(R)`, whereas every lattice realization has `R>=R_min`.

If one additionally includes the relative anchor phase `1` as an actual
point of the configuration, its primitive representative is `1`, so
`G_N=1` and the least radius is `|L_N|`. Explicitly, if its Gaussian
coordinate is `z_0`, the other coordinates are `z_0 P_i/bar(P_i)`;
conjugate-primitivity makes integrality equivalent to `bar(L_N)|z_0`.
The choice `z_0=bar(L_N)` realizes the least radius. The relative phase
`1` does not require `z_0` to lie on the real axis.

Here the uncorrected lcm has all `B-1` nonempty blocks. Corrections
alter its logarithmic modulus by at most `sum_i log|K_i|`, giving

```text
log R_min,anchored=(B-1)qN^2+o(N^2),
delta_N,anchored sqrt(R_min,anchored)
 >= exp(((B-5)/2)qN^2-o(N^2)).                        (A19)
```

Gaussian units change none of these integrality ideals or radius
formulas; the actual phase choices are already retained in Section 4.
Equations (A18)--(A19) concern this fixed denominator construction only.
They give no uniform improvement for arbitrary lattice-circle families.

The exact [radius checker](check_least_radius_formula.py) verifies the
Gaussian lcm/gcd normalization, the added-anchor distinction, and the
primewise correction inequality. The
[phase checker](check_cm_denominator_phase_lift.py) also verifies 102
rational secant-quotient reflection identities. These finite checks
supplement the algebraic and height arguments above.

The [slowly varying shift theorem](cm_slow_shift_phase_exclusion.md)
extends (A15), (A18), and (A19) to
`max|k_b(N)|<=N^(1/2-epsilon)` for fixed `m` and `epsilon>0`.
Its separate proof retains the moving coefficient, exceptional-point,
and translated height constants; the fixed-family statements here
alone do not justify that extension.
