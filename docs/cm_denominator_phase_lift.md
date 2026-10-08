# Actual Gaussian phases in the fixed CM denominator family

The CM norm construction has a rigid phase law. Equal-cardinality
products of its Gaussian denominator generators have relative phases
given by fixed elliptic functions of `NP`, up to finitely many fixed
Gaussian constants. Consequently this fixed family cannot give small
primitive residues when its row corrections have bounded height.

For arbitrary moving corrections, `log|K_i|=o(N^2)` is much weaker
than a fixed bound. The single-function argument below does not handle
them. The subsequent [collective phase audit](cm_denominator_phase_lift_audit.md)
does exclude them for every fixed shift assignment with at least three
rows. The collective image contains a rational degree-two coordinate,
including in the constant-complement-sum case, and denominator
separation gives a positive residue-growth exponent. This result
concerns the fixed CM family and does not prove the uniform lattice-arc
bound.

The [slowly varying shift extension](cm_slow_shift_phase_exclusion.md)
also permits `max|k_b(N)|<=N^(1/2-epsilon)` for fixed row count and
fixed `epsilon>0`, with all parameter dependence charged explicitly.

## 1. An explicit generator rather than only its norm

Use

```text
E:y^2=x^3-2x,       P=(2,2),       I=[i]P=(-2,2i),
Q_n=nP+I,
nP=(A_n/B_n^2,C_n/B_n^3),       B_n>0,
gcd(A_n,B_n)=gcd(C_n,B_n)=1.
```

For `n!=0`, set

```text
D_n=A_n+2B_n^2>0,
G_n=gcd_G(D_n,C_n+2i B_n^3).                           (1)
```

Let `b_n` be the Gaussian denominator ideal of `Q_n`, as in
[elliptic_split_inert_norm_route.md](elliptic_split_inert_norm_route.md).
At every odd Gaussian prime, `G_n` generates exactly `b_n`.
In particular its odd norm equals the odd part of `D_n`.

Indeed

```text
(C_n+2iB_n^3)(C_n-2iB_n^3)
 =D_n (A_n^2-2A_n B_n^2+2B_n^4).                     (2)
```

At an odd prime dividing `D_n`, the two factors on the left are
coprime because their difference is `4iB_n^3` and
`gcd(D_n,B_n)=1`. Thus the factor `C_n+2iB_n^3` takes one entire
orientation of every prime power of `D_n`. At that orientation the
addition slope `(C_n-2iB_n^3)/(B_n D_n)` has valuation
`-v(D_n)`, so `x(Q_n)` has valuation `-2v(D_n)` and the denominator
ideal has exponent `v(D_n)`. At the opposite orientation the
addition slope is integral and no pole occurs. If a prime divides
`B_n`, then `D_n` is a unit there and translation by the integral
point `I` sends the reduction of `nP=O` to `I`, so again no pole
occurs. This proves the exact odd-prime claim, including arbitrary
prime powers and the prime five.

At two, `v_2(D_n)<=2`. If `B_n` is even, `D_n` is odd. If `B_n`
is odd and `8|D_n`, then `A_n=6 mod8`, and the curve equation
would give `C_n^2=12 mod16`, impossible. The earlier formal-group
argument gives `v_(1+i)(b_n)<=2`. Thus (1) differs from an actual
denominator generator only by a Gaussian unit and one of finitely
many powers of `1+i`; their exponent lies between `-2` and `4`.
All these factors have bounded height. Put `G_0=1` when needed.

This proves the suggested gcd formula with its phase and ramified
costs retained. An arbitrary choice of a representation of `D_n`
as a sum of two squares would not retain the same orientation.

## 2. An exact elliptic-function law for relative phases

Choose generators `b_n=(beta_n)` and put

```text
gamma_n=beta_n/bar(beta_n) in Q(i)^*.
```

Changing a Gaussian associate changes `gamma_n` by a sign. For a
fixed nonzero integer shift `l`, let `ell_l` be the line through
`(l-i)P` and `iP`, normalized to have coefficient one on `y`.
It is nonvertical: the two points sum to `lP!=O`. Its conjugate
is the line through `(l+i)P` and `-iP`. Both lines have the same
third intersection `-lP`. Therefore the rational function

```text
f_l=ell_l/bar(ell_l),       f_l(O)=1,
div(f_l)=[(l-i)P]+[iP]-[(l+i)P]-[-iP]                 (3)
```

is nonconstant and has norm one on real points. The apparent ratio
`0/0` at the common third intersection is interpreted after
cancellation as a rational function. None of its zeros or poles is
rational over `Q`.

There is an exact finite-unit relation

```text
gamma_(N-l)/gamma_N
    =epsilon_(N,l) gamma_(-l) f_l(NP),
epsilon_(N,l) in {1,-1,i,-i}.                         (4)
```

Here is the arithmetic justification, including the exceptional
primes. At an odd prime the curve has good reduction. The divisor
of `f_l` on the smooth model is the four horizontal sections in
(3), plus a constant multiple of the special fiber. Intersecting
with the section `NP` gives the differences of the four contact
depths, plus a constant independent of `N`. Translation on the
smooth group identifies those depths with the denominator ideal
exponents of `Q_(N-l)`, `Q_N`, and their conjugates. Comparing with
the section `O` removes that vertical constant and gives the same
odd-prime valuations on both sides of (4). This argument includes
colliding reductions of the four sections: it uses their full
intersection multiplicities, not only distinct residue classes.

The quotient of the two sides has complex norm one. At the unique
prime over two, any norm-one element of `Q(i)` has valuation zero.
Consequently this quotient is a unit at every Gaussian prime and
is a Gaussian unit. This proves (4). The needed contact/formal
depth identification is the one in
[Streng, Section 2](https://pub.math.leidenuniv.nl/~strengtc/cmeds.pdf);
the chord divisor calculation and the comparison with `O` give
the exact normalization here.

For `l=1`, the line is particularly simple:

```text
L=2y+x+2+i(x-2),       f_1=L/bar(L).                  (5)
```

The use of `G_n` from (1) in place of `beta_n` only changes the
finite unit in (4). The bounded common-factor trimming in the CM
norm note enlarges the finite set of constants by bounded Gaussian
factors. It does not introduce a factor whose height grows with `N`.

## 3. Every equal-cardinality Boolean phase word is algebraic

The construction and its bounded trimming work for any fixed finite
list of distinct integer shifts, not just the original ten. Assign
distinct shifts `k_S` to all `B=2^m` subsets of `{1,...,m}` and
write the resulting independent, odd, conjugate-primitive blocks
as `H_(N,S)`. Then

```text
log N(H_(N,S))=2(N-k_S)^2 hhat(P)+O(1),
```

where the harmless `+2hhat(P)` from `Q_(N-k)` is absorbed in
`O(1)`. The constants may depend on this entire fixed list.

For integer coefficients `c_S` with `sum_S c_S=0`, (4) gives

```text
product_S (H_(N,S)/bar(H_(N,S)))^c_S
        =a_N F_c(NP),                                 (6)
```

where `a_N` belongs to a fixed finite subset of the norm-one group
of `Q(i)`, and `F_c` is a fixed rational function over `Q(i)` with

```text
div(F_c)=sum_S c_S ([(k_S-i)P]-[(k_S+i)P]).            (7)
```

The zero-sum condition cancels the common `gamma_N` and is exactly
the Abel-sum condition for this divisor to be principal. All the
points in (7) are distinct: equality would give a nonzero Gaussian
endomorphism annihilating the nontorsion point `P`. Thus every
nonzero word `c` gives a nonconstant function `F_c`.

In particular, for two distinct Boolean rows `i,j`, take

```text
c_S=1_(i in S)-1_(j in S).
```

Both rows contain `B/2` blocks, so this is a nonzero zero-sum word.
Their relative phases are constrained by a nonconstant fixed
elliptic function. Unit choices or bounded trimming cannot make
this function constant.

## 4. Quantitative residue growth for bounded corrections

The following fixed-data consequence uses a proved elliptic-logarithm
bound. If `F` is a fixed nonconstant rational function on `E` and
`A` is a fixed finite set of nonzero algebraic constants, then,
apart from finitely many exact equalities,

```text
min_(a in A) |F(NP)-a|
 >= exp(-C log(N+2) (log log(N+3))^5).                (8)
```

To justify the precise scope, the zeros of the finitely many
functions `F-a` are fixed algebraic points. Near such a point,
their values are bounded below by a fixed power of the analytic
distance. Write a lift of that distance as
`N u_P-u_Q-r omega_1-s omega_2`, with `r,s=O(N)`. The fixed-target
linear-form theorem
[David's theorem as stated in Streng, Theorem 3.13](https://pub.math.leidenuniv.nl/~strengtc/cmeds.pdf)
gives (8) with four fixed elliptic logarithms. Away from the finite
zero set, compactness supplies a constant lower bound. Exact
equalities occur only finitely often because `P` is nontorsion.
The constant `C` depends on all the fixed data and targets.

Now form actual conjugate-primitive rows

```text
P_(N,i)=K_(N,i) product_(S containing i)H_(N,S),
|K_(N,i)|<=K_0,
```

with a fixed `K_0`. Their correction phases form a finite set.
Combine (6) with (8). Every pair of distinct Boolean rows has
relative phase distance from one bounded below by the right side
of (8), after changing `C`. The exact primitive pair numerator

```text
h_ij=P_(N,i) bar(P_(N,j))/N(gcd_G(P_(N,i),P_(N,j)))
```

has

```text
log|h_ij|=(B/2) hhat(P) N^2+O(N).                    (9)
```

There are exactly `B/2` separated blocks; their factors survive
with the appropriate orientation. Bounded corrections can change
only bounded norm factors. Since

```text
|h_ij/bar(h_ij)-1|=2|Im(h_ij)|/|h_ij|,
```

equations (8)--(9) give

```text
log|t_ij| >= (B/2) hhat(P) N^2
              -O(N)-O(log N (log log N)^5),
t_ij=Im(h_ij)!=0,                                   (10)
```

for all sufficiently large `N`. In particular the primitive pair
residues have the same leading logarithmic size as the pair
numerators. They are not subpower in the common block norm.
Thus bounded corrections to this fixed CM family cannot yield
the endpoint system, regardless of the Gaussian-unit choices.

## 5. The remaining moving-coefficient gap

The endpoint reduction permits corrections whose logarithmic
heights are `o(N^2)` but are unbounded. In that case (6) still
holds for the uncorrected blocks, but the desired small-angle
condition compares `F_c(NP)` with a moving norm-one rational
constant built from `K_(N,i)` and `K_(N,j)`.

The fixed-target constant in (8) is not uniform over these targets.
Allowing it to depend on `N` and then calling the error `o(N^2)`
would be invalid. The norm construction, the exact gcd generator,
and the elliptic-function phase law do not by themselves prove the
needed moving-target bound. A new theorem controlling that
dependence, or an exact relation eliminating the moving correction
phases, is still required for this single-function argument. The
collective argument in the linked audit supplies an alternative for
every fixed shift assignment with at least three rows, including the
symmetric case.

The phase law is also not a construction of small residues: it
identifies the nonconstant function which such a construction
would have to approximate with its correcting coefficients.

Exact generator, chord, phase, and Boolean-word identities are
checked in
[check_cm_denominator_phase_lift.py](check_cm_denominator_phase_lift.py).
