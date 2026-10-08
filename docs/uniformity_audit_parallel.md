# Independent audit of the fixed-profile uniformity reduction

This is a prose audit of `uniform_profile_extraction.md`. Its reduction
is valid, but its missing arithmetic assertion remains unproved. The
main additional result here is a quantitative obstruction: the existing
complex-factor relaxation reaches **every fixed neighborhood** of the
extracted profile while satisfying all the product and separation tests
already established. Thus those tests do not acquire a new consequence
merely by being combined with extraction.

## 1. Audit of the extraction and its finite constants

The following steps in the original note check out.

* With one Gaussian-unit class and the common Gaussian gcd removed,
  the primitive-chord formula gives
  `d_ij >= W/2 + 2 log(2/C)`. Hence the threshold sign vectors are
  pairwise obtuse when `C<=2`. No lower bound on the growth rate of
  the number of rows is used.
* The Bessel constant two follows by separating positive and negative
  inner products. Conditioning on `r-1` distinct row choices leaves
  at least `M-r+1` possible last rows, so the moment estimate
  `E mu_S^2 <= 2/(M-r+1)` has the correct denominator. The weighted
  product of the previously chosen sign functions has norm one.
* Fourier inversion and Cauchy--Schwarz give the stated simultaneous
  relative error `(2^k-1) sqrt(2/(M-k+1))`.
* Passing to a tuple's own Gaussian gcd removes exactly its initial
  all-positive and final all-negative threshold layers. The surviving
  layer signs and weights are unchanged. This remains true for prime
  powers; squarefreeness is not required.
* For eight rows, inherited error `eta` gives intrinsic unoriented-cut
  relative error at most `128 eta/(127-eta)`. Taking `eta=eta_0/2`
  with `0<eta_0<=1` is sufficient, and yields the advertised threshold
  `M >= 7+520200/eta_0^2`.
* At `C=1/2`, the strict obtuse gap gives
  `W >= 4(M-1) log 4`. The inherited constant mass is at most `1/64`
  when `eta<=1`, so `W' >= (63/64)W`. Therefore the height condition
  `M >= 1+16 max(W_0,0)/(63 log 4)` is sufficient as claimed.
* The eight-point uniform intrinsic excess is `21W'/127`.
  The neighborhood estimate `217W'/2540` at relative error `1/20`
  follows from coefficient sum `21` and absolute coefficient sum `203`.
* For even `k`, the positive-bonus criterion is exactly
  `delta binom(k,k/2) > k(k-1)/2`. In particular the stated
  `k=32`, `delta=10^(-6)` comparison is correct.

No gap was found in these reductions or in the stated verification
boundary. The full extraction is still a prose proof; the Hilbert-space
ingredient in `GaussianChain/ObtuseBessel.lean` is the formalized part.

One useful scale must not be forgotten: if the removed squared norm is
`exp(A)`, the new radius and endpoint constant are

```text
R' = R exp(-A/2),       C' = C exp(-A/4).
```

Discarding this improvement is legitimate for the stated reduction.
Keeping it makes the normalized tuple a still shorter arc, but does
not contradict the product-rigidity results: the relaxation below has
no multiplicative relations of any length.

## 2. Complete half-cuts approach every prescribed intrinsic neighborhood

Fix `k>=2` and `0<eta<=1`. Choose an even integer

```text
M >= 32 k^2/eta.
```

Give equal weight to all `M/2`-subsets of `M` row labels, as in
`integer_combination_separation.md`. Restrict to any `k` distinct rows.
For one specified oriented sign pattern with `r` positive entries,
its inherited weight fraction is

```text
P_M(r) = binom(M-k,M/2-r)/binom(M,M/2)
       = (M/2)_r (M/2)_(k-r)/(M)_k,
```

where the subscript denotes a falling factorial. The formula holds for
each individual pattern, not the total weight of its Hamming level.
Every pattern has positive weight because `M>=2k`.

For `0<=x<=1/2`, use `|log(1-x)|<=2x`. Expanding the falling
factorials gives

```text
|log(2^k P_M(r))|
 <= [2r(r-1)+2(k-r)(k-r-1)+k(k-1)]/M
 <= 3k(k-1)/M
 <= eta/8.
```

Thus, uniformly in the pattern,

```text
|2^k P_M(r)-1| <= exp(eta/8)-1 <= eta/4.
```

Remove the two constant patterns and normalize the remaining weight.
Put `a_0=2^(1-k)` and `epsilon=eta/4`. Write the inherited weight
of each pattern as `2^(-k)(1+delta_sigma)`, and the combined constant
weight as `a_0(1+delta_0)`, with all these errors at most `epsilon`.
For each nonconstant pattern the relative intrinsic error is at most

```text
epsilon/[1-a_0-a_0 epsilon] <= 2eta/3 < eta,
```

because `a_0<=1/2` and `epsilon<=1/4`.

Consequently **every** intrinsic `k`-row profile in this one model is
within the requested relative tolerance of the uniform distribution on
the `2^k-2` nonconstant patterns. Spatial selection of a different
tuple cannot avoid this conclusion: the formula is identical for every
tuple of the given size.

## 3. The simultaneous product, phase, and determinant tests survive

For this same even `M`, choose `0<C<=1/2` and then choose `W` as
large as necessary. Section 3 of `integer_combination_separation.md`
constructs compatible short phases satisfying, for every nonzero
integer vector `c` with coordinate sum zero,

```text
exp(H(c)/2) |sin((sum_i c_i phi_i)/2)| >= 1.
```

This is a simultaneous assertion over all such `c`, not a finite list.
The factor-phase realization there has the complete-half-cut incidence
matrix used in Section 2. Therefore the following coexist:

1. arbitrarily large fixed `M` and arbitrarily large `W`;
2. an actual real-circle arc of length at most `C sqrt(R)`;
3. compatible complex-factor products with the prescribed valuations;
4. every integer-combination imaginary-part separation lower bound;
5. every exact multiplicative product-rigidity conclusion, of every
   length, since a nonzero zero-sum combination cannot have phase zero
   modulo `2 pi`;
6. every chord-retaining subset determinant inequality, by the pair
   separation calculation in Section 7 of `finite_local_test_attack.md`;
7. the intrinsic profile conclusion of Section 2 for every `k`-tuple.

This also addresses the shrinking constant `C'` after normalization.
The subset chord identities remain true after division by the formal
common factor, and a common factor cancels from every zero-sum product
comparison. Absence of relations of every length is stronger than any
larger finite `B_K` conclusion obtained from the smaller `C'`.

The qualification is essential: these are **complex factors**, and
their imaginary parts are not asserted to be integers. This construction
is not a configuration of Gaussian lattice points. Indeed the canonical
factor realization has a proved failure of Gaussian integrality.

### Consequence for every proposed fixed-size positive bonus

For even `k`, let

```text
lambda = [delta binom(k,k/2)-k(k-1)/2]/(2^k-2).
```

If the proposed bonus satisfies the sufficient comparison, then
`lambda>0`. The intrinsic excess per unit weight at the uniform
profile is exactly `lambda`. For the coefficient on a pattern of
size `r`, write

```text
q(r) = (r-k/2)^2 + delta 1_(r=k/2) - k/4.
```

On a relative-`eta` profile neighborhood, the excess per unit weight
is at least `lambda-eta A`, where

```text
A = (1/(2^k-2)) sum_(r=1..k-1) binom(k,r) |q(r)|.
```

Choose `eta>0` sufficiently small that `eta A<=lambda/2`, and then
apply Section 2. Letting `W` grow makes the intrinsic excess at least
`lambda W'/2`, and hence larger than every fixed additive allowance.
All the simultaneous tests above still hold.

Thus moving to a larger fixed tuple or an arbitrarily tiny positive
bonus does not circumvent this combined relaxation. A proof of the
bonus must retain additional Gaussian arithmetic.

## 4. What integral reflections add at the incidence level

For squarefree split factors, a pair-bisector reflection sends a third
row to a Gaussian integer exactly when that third row agrees with the
anchors at every factor where the anchors agree. At any nonconstant
full-support profile, this criterion fails for every triple: some
factor has the anchors equal and the third opposite. Thus the necessary
reflection-avoidance condition is already supported by the limiting
profile, rather than conflicting with it.

There is a useful quantitative version. In a primitive `k`-row
squarefree model with weights within relative error `eta<1` of the
uniform nonconstant profile, put `W=log N`. For two anchors, the norm
`D` of their reflection-domain generator obeys

```text
(1-eta) (2^(k-1)-2)/(2^k-2) W
 <= log D
 <= (1+eta) (2^(k-1)-2)/(2^k-2) W.
```

The count `2^(k-1)-2` enumerates nonconstant patterns where the
anchors agree. For a distinct third row (`k>=3`), let `Q` be the norm
of the reduced Gaussian denominator of its reflected image. Then

```text
(1-eta) 2^(k-2)/(2^k-2) W
 <= log Q
 <= (1+eta) 2^(k-2)/(2^k-2) W.
```

Here the contributing patterns are precisely `001` and `110` on the
three specified rows, with the other entries free. They number
`2 * 2^(k-3)`. Each such factor occurs once in the reduced denominator,
and no other factor occurs there. For large fixed `k` the two logarithmic
fractions are close to `1/2` and `1/4`, respectively.

These are exact statements about genuine Gaussian factors when such a
squarefree profile is given. They impose no claim about its arc length.
In the complex relaxation they describe the corresponding formal prime
valuations only. The independent symmetry audit establishes a stronger
whole-orbit obstruction for complete half-cuts.

### Prime-power caution

The last two formulas must not be inferred from threshold-layer
statistics in the general prime-power setting. For one split prime with
norm exponent `e`, anchor allocations `a,b` and third allocation `c`,
the reflected allocation is `a+b-c`, so the reduced denominator's
exponent is

```text
max(0,c-a-b,a+b-c-e).
```

Its reflection-domain exponent is `|a+b-e|`. These quantities depend
on how threshold layers are grouped into the same prime. A linear
count of threshold patterns can lose cancellation inside one prime.
The extraction theorem is valid for arbitrary exponents, but its
profile alone does not determine these reflection quantities.

## 5. Result of the audit

The fixed-profile exclusion would indeed prove a bound independent of
`R`, including for arbitrarily slowly growing clusters. No error in
that reduction was found. However, almost-uniform profile extraction,
all established product-rigidity conclusions, all integer-combination
separation inequalities, and the incidence restrictions of rational
reflections do not produce the missing exclusion. The simultaneous
relaxation reaches arbitrarily small fixed profile neighborhoods.

The arithmetic that remains absent is simultaneous integrality of the
actual Gaussian factors and residues, beyond their separate magnitude
lower bounds. None of the statements here proves uniformity or improves
the growth rate.

## 6. Independent audit of a new fixed-projective-template theorem

The subsequent projective investigation supplies a genuinely different
restricted result. Fix `M>=5` distinct rational projective directions,
represented by primitive integer vectors `v_1,...,v_M`. Let `A` range
over `GL_2(Q)`, write its images as complex numbers `U_j=A v_j`, and
allow any `beta in Q(i)` for which all

```text
z_j = beta U_j/conjugate(U_j)
```

are Gaussian integers. For each fixed endpoint constant `C`, the radii
of all such configurations on arcs of length at most `C sqrt(R)`
are bounded in terms of `C` and the fixed template. This does not
bound configurations whose projective shape varies with the radius.

The independent checks below support the theorem developed in the
parallel projective investigation.

### Normalizing the matrix is legitimate

First clear rational denominators of `A`; real scalar multiplication
does not change `U_j/conjugate(U_j)`. View its integer columns as
Gaussian integers and divide both by their Gaussian gcd `g`. The
result again has integer columns, now Gaussian-coprime, and nonzero
integer real determinant `D`. Its change to every quotient is the
common modulus-one factor `conjugate(g)/g`, which is absorbed in
`beta`. Thus the radius `R=|beta|` is unchanged.

Put

```text
Delta_ij = det(v_i,v_j),
S = product_(i<j) |Delta_ij|,
d_* = min_(i<j) |Delta_ij|,
H = max_j |U_j|,       m = min_j |U_j|.
```

All `Delta_ij` are nonzero integers and `m>0`. Gaussian coprimality
of the two normalized columns implies

```text
gcd_G(U_i,U_j) divides Delta_ij.
```

Indeed the gcd divides `Delta_ij` times each column; a Gaussian
Bezout identity for the columns removes the latter factors.

### Rational contents have only one determinant's worth of common growth

Let `r_j=gcd(|Re U_j|,|Im U_j|)`. The integer adjugate identity
`adj(A) U_j=D v_j`, together with primitivity of `v_j`, gives `r_j|D`.
For every rational prime `p`, at least one entry of the normalized
matrix is a `p`-unit. If `p^t` divides both contents `r_i,r_j`, the
adjugate identity for the two template columns shows that `p^t`
divides `Delta_ij` times every entry of `A`. Hence

```text
min(v_p(r_i),v_p(r_j)) <= v_p(Delta_ij).
```

The elementary inequality
`sum_j b_j <= max_j b_j + sum_(i<j) min(b_i,b_j)` for nonnegative
integers now proves the divisibility, and consequently the bound,

```text
product_j r_j divides |D| S,
product_j r_j <= |D| S.
```

Using only the separate bounds `r_j|D` would lose the theorem's
exponent. The product bound above is therefore a necessary check.

### Fractional beta is covered by an exact ideal calculation

For each `j`, divide `U_j` by `r_j` and, if its two primitive
coordinates are odd, additionally divide by `1+i`. Call the result
`d_j`. It is Gaussian-coprime to its conjugate, and

```text
U_j/conjugate(U_j) = epsilon_j d_j/conjugate(d_j),
|d_j| >= |U_j|/(sqrt(2) r_j),
```

where `epsilon_j` is a Gaussian unit. Let `L=lcm_G(d_j)` and
`G=gcd_G(d_j)`. The permissible common multipliers form exactly

```text
(conjugate(L)/G) Z[i].
```

To check this even when `beta` is fractional, fix a Gaussian prime
and put `x_j=v_pi(d_j)`, `y_j=v_pi(conjugate(d_j))`. For every `j`,
`x_j y_j=0`. Therefore

```text
max_j(y_j-x_j) = max_j y_j - min_j x_j,
```

which is precisely the required lower bound on `v_pi(beta)`.
Negative lower bounds are allowed, so this calculation has not assumed
that `beta` itself is Gaussian. It gives `R>=|L|/|G|`.

Pairwise Gaussian gcds of the `d_j` divide the fixed `Delta_ij`.
Consequently `|G|<=d_*` and
`|L|>=product_j |d_j|/S`. Combining with the content-product bound,

```text
R >= 2^(-M/2) product_j |U_j| / (|D| S^2 d_*).
```

### All but the shortest image are comparable to the longest

Choose a fixed `K>=1` bounding the absolute values of the two
coefficients that express any `v_l` in the basis of any two other
distinct template directions. Then

```text
H <= K (|U_i|+|U_j|)       for every distinct i,j.
```

Apply this to the two shortest images. Every image except possibly
the shortest has magnitude at least `H/(2K)`. Thus, for a positive
template constant `c_0`,

```text
R >= c_0 m H^(M-1)/|D|,
c_0 = 2^(-M/2) (2K)^(-(M-1))/(S^2 d_*).
```

### The angular inequality and the final exponent

Suppose the containing circle arc has angular width `Delta<pi`.
The projective directions of the `U_j` then lie in a half-angle
interval of width `Delta/2`. Take a shortest image and any distinct
one. Since
`|det(U_i,U_j)|=|D Delta_ij|`, their sine separation gives

```text
Delta >= 2|D| d_* /(m H).
```

Consequently, writing `c_1=2 d_* sqrt(c_0)>0`,

```text
Delta sqrt(R)
 >= c_1 sqrt(|D|) H^((M-3)/2)/sqrt(m)
 >= c_1 sqrt(|D|) H^((M-4)/2)
 >= c_1 H^((M-4)/2).
```

The last step uses the integer condition `|D|>=1`. This confirms
the candidate's exponent, including the use of `m<=H`.

When the arc constant is at most `C`, it follows that

```text
H <= (C/c_1)^(2/(M-4)).
```

There is also a direct radius bound. The angular inequality implies
`Delta>=2d_*/H^2`, so

```text
R <= C^2 H^4/(4d_*^2)
  <= C^2 (C/c_1)^(8/(M-4))/(4d_*^2).
```

If the given containing arc has width at least `pi`, its endpoint
length bound already gives `R<=(C/pi)^2`. Taking the maximum of
these two bounds covers all cases. Since bounded radii contain only
finitely many centered Gaussian lattice points, this proves finiteness
of the endpoint configurations from the fixed template.

The cutoff `M=5` is substantive. At `M=4`, the sequence of columns
giving `U_j=t+j+i` for four fixed distinct integers `j` has constant
nonzero determinant, image magnitudes of order `t`, and half-angle
span of order `t^(-2)`. The common multiplier
`beta=product_j conjugate(U_j)` gives Gaussian points of radius of
order `t^4`, so their endpoint constants stay bounded. This matches
the zero power of `H` at `M=4`.

Exact checks over 2,112 invertible integer matrices with entries from
`-3` to `3`, for the five directions `(1,j)`, `j=0,...,4`, verified
Gaussian column normalization, the content-product divisibility,
integrality of the exact minimum multiplier, and the displayed radius
lower bound. These checks supplement the proofs above; they do not
replace them. No normalization, denominator, or exponent gap was found.

## 7. Cross-ratio heights forced by the extracted profile

The fixed-projective-template theorem is compatible with the extraction
for a precise arithmetic reason: the extracted rational invariants
necessarily acquire height approximately `R^(1/4)`. This statement
holds for arbitrary prime exponents, not only squarefree conductors.

Take four distinct rows `a,b,c,d` in one common Gaussian-unit class
and form the real rational cross-ratio

```text
x = (z_a-z_c)(z_b-z_d)/((z_a-z_d)(z_b-z_c)).
```

Let `t_ij` denote the signed nonzero integer primitive chord residue.
At one split prime `p`, with row allocations denoted `a_p,b_p,c_p,d_p`,
define

```text
b(p) = min(a_p,c_p)+min(b_p,d_p)
       -min(a_p,d_p)-min(b_p,c_p).
```

The exact common Gaussian factors in the four differences cancel to
the same exponent `b(p)` at `pi_p` and its conjugate. Hence

```text
x = (t_ac t_bd)/(t_ad t_bc) product_p p^(b(p)).       (17)
```

Signs are included in the signed residues. This identity follows
directly from `z_i-z_j=2i epsilon t_ij G_ij`, with

```text
G_ij = product_p pi_p^min(a_i,a_j)
                     conjugate(pi_p)^(e_p-max(a_i,a_j)).
```

### The conductor numerator and denominator cannot cancel

Write `A=product_(b(p)>0) p^(b(p))` and
`B=product_(b(p)<0) p^(-b(p))`. These integers are coprime.
If `b(p)>0`, then the denominator pairs `(a,d)` and `(b,c)` both
have different allocations at `p`. A primitive chord residue for
two different allocations is a `p`-unit: its reduced Gaussian factor
is divisible by exactly one of `pi_p,conjugate(pi_p)`, so its
difference from its conjugate is divisible by neither. Therefore
`p` divides neither `t_ad` nor `t_bc`.

Similarly, if `b(p)<0`, then `p` divides neither `t_ac` nor `t_bd`.
It follows from (17) that, for the reduced rational expression `x=n/q`
with `q>0`,

```text
A divides n,       B divides q.                    (18)
```

Thus residue cancellation cannot erase either of these conductor
products. In particular

```text
max(A,B) <= height(x)
 <= max(A |t_ac t_bd|, B |t_ad t_bc|),              (19)
```

where `height(n/q)=max(|n|,q)`.

To verify the asserted signs with arbitrary exponents, use

```text
b(p) = sum_(t=1..e_p)
  (1_(a_p>=t)-1_(b_p>=t)) (1_(c_p>=t)-1_(d_p>=t)).
```

Each of the two differences has one fixed sign on its interval of
nonzero thresholds. Their product consequently has one fixed sign
where nonzero; there is no cancellation between positive and negative
threshold contributions inside a prime. A positive product means
that the two allocation intervals overlap with the same orientation,
which forces `a_p!=d_p` and `b_p!=c_p`. A negative product similarly
forces `a_p!=c_p` and `b_p!=d_p`. This proves the `p`-unit assertions
used in (18).

### Quantitative consequence for an inherited nearly uniform profile

Suppose the inherited distribution on a fixed set of `k>=4` rows has
every oriented sign-pattern weight within relative error `eta` of
`2^(-k) W`, where `W=log(R^2)`. On the four specified rows, the
positive threshold contribution in (17) occurs in the patterns `1010`
and `0101`; the negative contribution occurs in `1001` and `0110`.
Each pair accounts for exactly one eighth of the uniform weight.
The absence of within-prime sign cancellation gives

```text
(1-eta)W/8 <= log A, log B <= (1+eta)W/8.           (20)
```

Every pair cut distance is at most `(1+eta)W/2`. The endpoint chord
bound therefore gives

```text
log |t_ij| <= log(C/2)+eta W/4.
```

Combining with (19)--(20),

```text
R^((1-eta)/4) <= height(x)
 <= (C^2/4) R^((1+5eta)/4).                        (21)
```

The statement holds simultaneously for every ordered quadruple in
the selected tuple. In a hypothetical unbounded sequence, the
extraction gives `eta -> 0` and `R -> infinity`, so (21) yields

```text
log height(x)/log R -> 1/4.
```

Dividing the selected tuple by its own common Gaussian gcd leaves
its cross-ratios unchanged. If `R'` is its intrinsic radius, then
`log R'/log R -> 1-2^(1-k)`, so the corresponding intrinsic exponent
is `1/[4(1-2^(1-k))]`.

For comparison, in the exact complete-half-cut model on `M` ambient
rows, each of the two conductor logarithms has fraction

```text
M(M-2)/(8(M-1)(M-3))
```

of the ambient `W`. This tends to `1/8`, consistently with (20).

### Implication for the new projective route

Normalizing three selected rational directions to `infinity,0,1`
expresses the remaining projective coordinates as such cross-ratios
or their permutations. All their reduced rational heights therefore
grow like `R^(1/4+o(1))` in the extracted configuration. A bound on
their real values would not change this height conclusion.

Thus the fixed-shape projective theorem rules out bounded arithmetic
shape, while the profile extraction predicts a specific unbounded
arithmetic shape. For instance a coarse radius bound of the form
`R <= constant(C) B^96` is fully compatible with `B=R^(1/4+o(1))`.
The exponent gap is substantial. No averaging of the mere fixed-height
statement bridges it; a stronger estimate exploiting relations among
the several changing invariants would be needed.

## 8. Independent audit of growing-tuple central truncation

The reduction in
[endpoint_central_truncation.md](endpoint_central_truncation.md)
has been checked independently, including arbitrary prime exponents,
ties in endpoint rounding, the common-unit selection, and the distinction
between the old and core radii. Its conclusion is valid: for every fixed
`0<c<1/2`, a hypothetical sequence with `M -> infinity` supplies
`k=floor(c log_2 M)` rows and `B=2^(k-1)` independent Gaussian core
blocks with relative weight error tending to zero. The exact primitive
anchor numerators are small Gaussian multiples of the core row products,
and their nonzero imaginary parts are small in logarithmic height relative
to every individual core block.

The following points are material to the conclusion.

1. The retained good rows have size at least `M/2`; taking a common
   Gaussian-unit class subsequently gives `L>=M/8`. This constant loss
   does not change the permissible range of `c`. Keeping the original
   conductor weight `W` throughout is essential.

2. The simultaneous averaging is valid: if `E A<=a` and `E P<=b`
   for the pattern and pair-moment square sums, averaging `A/a+P/b`
   gives one tuple with `A<=2a` and `P<=2b`. It yields
   `eta=2(2^k-1)/sqrt(L-k+1)` and
   `zeta=2sqrt(binom(k,2)/(L-1))` simultaneously. A bound derived
   only from the full pattern error would lose an unnecessary additional
   factor `2^k` in the comparison with individual block sizes.

3. For a prime exponent `e`, put `u_i=min(a_i,e-a_i)`,
   `r=max_i u_i`, and `e'=e-2r`. The retained threshold layers
   `r<t<=e-r` all have the rounded endpoint pattern. Thus they are
   literally retained original layers. If `q_S` is the original cut
   weight and `v_S` its removed weight, then
   `w'_S=q_S-v_S`, `v_S>=0`, and `sum v_S=2s`, where
   `s=sum_p r_p log p<=kW/sqrt(M)`. This proves directly

   ```text
   eta_core <= [eta+2(B-1)s/W]/[1-2s/W].
   ```

   The removed threshold masses `v_S` need not be the removed rounded
   block-divisor masses. The distinction is correctly retained in the
   central-truncation note.

4. For the original primitive anchor numerator
   `h_i=product pi^((a_i-a_0)_+) bar(pi)^((a_0-a_i)_+)`, the
   core row product `A'_i` divides `h_i` in the correct orientation.
   Opposite rounded labels leave exponent `2r-u_i-u_0`; equal labels
   leave `|u_i-u_0|`. Both lie in `[0,2r]`. Consequently

   ```text
   h_i=K_i A'_i,        log|K_i|<=s,
   W/4-s <= log|A'_i| <= (1+zeta)W/4.
   ```

   This factorization preserves conjugate coprimality. It handles ties
   because a tie forces `e'=0` at that prime.

5. The original arc, rather than an asserted core arc, gives the exact
   small-value bound

   ```text
   0<|Im h_i| <= (C/2)exp(zeta W/4).
   ```

   No correcting-factor height must be added to this exponent: `h_i`
   is the original primitive numerator. The alternative correction
   formula `B_i/D=C_i/bar(C_i)` is also exact, with
   `log N(C_i)=s-E_i`; using `B_i bar(B_0)` as a conjugate-quotient
   numerator directly would instead square the desired ratio.

For `C<=sqrt(2)`, hence `C/2<1`, the logarithmic coefficient and
nonzero-value bounds have no positive additive constant. Dividing by
`log N(H'_S)=(1+o(1))W/B` gives uniformly in `i,S`

```text
log max(1,|Re K_i|,|Im K_i|,|Im h_i|)/log N(H'_S)
  = O(k 2^k/sqrt(M)) = o(1).
```

For a fixed `C<sqrt(2)`, the strict pair-distance gap additionally
gives `W>=4(M-1)log(sqrt(2)/C)`. In particular every block's
logarithmic norm tends to infinity. Such a strict fixed `C` is available
by subdividing endpoint arcs before beginning the reduction.

As an exact finite check of the primewise arithmetic, all 138,464
allocation tuples with `2<=k<=5` and `1<=e<=8` were enumerated.
The check verified the central-layer cut, nonnegative correcting
exponents, the identity `|c_i-r|=r-u_i`, oriented divisibility
`A'_i|h_i`, and both formulas for the remaining exponent. This is
verification of the formulas, not a replacement for the argument above.

The audit finds no gap in this reduction. It strengthens the arithmetic
system to be excluded; an incompatibility theorem for that growing
system remains necessary for the requested uniform bound.

## 9. Independent audit of the converse Gaussian reconstruction

Suppose independent, conjugate-coprime Gaussian blocks `H_S` are given
for all subsets of `[m]`. Put `A_i=product_(S contains i) H_S`, and
suppose `h_i=K_i A_i` is conjugate-coprime with nonzero imaginary part
`t_i`. Set

```text
P=product_(S nonempty)|H_S|,
E=sum_i log|K_i|,
a=min_i log|A_i|,
T=max_i |t_i|,
sigma=a-(1/2)log P,
L=lcm_G(h_1,...,h_m).
```

The converse construction proposed in
[independent_block_reconstruction.md](independent_block_reconstruction.md)
is exact:

```text
z_0=bar(L),
z_i=bar(L) h_i/bar(h_i),
R=|L|,
P<=R<=P exp(E).                                   (26)
```

Integrality follows from `bar(h_i)|bar(L)`. For the radius bounds,
`product_(S nonempty)H_S` divides `L`, since every nonempty block
occurs in some row. Conversely
`L` divides `product_(S nonempty)H_S product_i K_i` by checking each
Gaussian prime valuation.

This is also the **least possible radius** for these ratios with the
anchor included. If an integral anchor `b` makes every
`b h_i/bar(h_i)` integral, conjugate coprimality gives
`bar(h_i)|b` for every `i`. Thus `bar(L)|b`, so `|b|>=|L|`.
The constructed array is itself primitive: at a Gaussian prime
dividing `bar(L)`, choose a row attaining the maximal corresponding
conjugate valuation in `L`. That row's constructed point has zero
valuation at the prime. No Gaussian prime divides all the points.

For distinctness, equality of two ratios forces `h_i=+h_j` or
`h_i=-h_j`. Indeed the equality and conjugate coprimality first force
the numerators to be associates, and equality of the conjugate
quotients then restricts the unit to `+1` or `-1`. A block with
`i in S`, `j notin S` would therefore divide `K_j`. In particular

```text
max_i log|K_i| < min_(S nonempty) log|H_S|          (27)
```

ensures pairwise distinct nonanchor points. Every point differs from
the anchor because `t_i!=0`.

If `beta=T exp(-a)<=1/2`, then `|t_i|/|h_i|<=beta`. Every ratio
`h_i/bar(h_i)` consequently has its argument in
`[-2 arcsin(beta),2 arcsin(beta)]`, regardless of the sign of
`Re h_i`. Thus one common arc containing the anchor has angular
width at most `4 arcsin(beta)<=8 beta`. Its length divided by
`sqrt(R)` is at most

```text
8 T exp(E/2-sigma).                               (28)
```

No assertion that each `h_i` itself lies near the positive real axis
is needed.

The exact balanced-weight calculation, formally putting
`log N(H_S)=w` for all blocks, is `sigma=w/4`. Actual independent
nonunit blocks cannot have literally identical integer norms, so the
asymptotic conclusion should use actual marginal balance
`sigma=w/4+o(w)`. Relative uniformity of all block weights alone does
not supply this stronger additive accuracy when `m` grows.

The central-truncation output does supply it through its independent
pair-moment estimate. If `s=log D` is the truncation cost and
`w_empty=log N(H_empty)`, the audited row bounds give

```text
w_empty/4-s/2 <= sigma
 <= w_empty/4+s/2+zeta W/4,
E<=(k-1)s.                                        (29)
```

Here `s=O(kW/sqrt(M))`, `zeta=O(k/sqrt(M))`, and
`w_empty=(1+o(1))W_core/2^(k-1)`. Hence for every fixed `c<1/2`,
`k=floor(c log_2 M)` yields
`sigma=w/4+o(w)` and `E+log T=o(w)`, with
`w=W_core/2^(k-1)->infinity`. Formula (28) then tends to zero.
Conversely any abstract sequence satisfying these smallness and
balance hypotheses reconstructs actual endpoint clusters. This
establishes a converse reduction; it does not assert that such an
abstract sequence exists.

An exact Gaussian-integer check covered 1,000 constructions with two
or three nonreference rows, independent blocks of varying prime-power
exponents, and correcting factors including conjugate primes where
primitivity allowed them. It verified both divisibility bounds in
(26), common radius, primitive reconstructed arrays, and criterion
(27). The angle estimate and minimal-radius assertion follow from
the exact arguments above.
